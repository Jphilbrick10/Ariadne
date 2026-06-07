"""Transit search + Equation-of-One vetting of TESS light curves.

Finding a periodic dip is the easy, commodity part (Box Least Squares). The hard
part -- the actual bottleneck in every transit survey -- is deciding whether a dip
is a real planet or one of the impostors:

  * an ECLIPSING BINARY (two stars): far too deep, often a secondary eclipse, and
    V-shaped (grazing) rather than the flat-bottomed U of a planet crossing a disk;
  * a BLENDED binary at half the period: alternating (odd vs even) transit depths;
  * stellar variability / systematics: low signal-to-noise, no coherent shape.

A real planet's transit is a tightly constrained physical object: depth below a few
percent, no secondary, equal odd/even depths, a flat U-shaped floor, good SNR. That
is a coherence BASIN. We score each candidate's incoherence energy against it (the
S_One E_c term) -- no training set needed, which is the point: this works on the
candidates a trained classifier never saw and a standard pipeline DISCARDED.

Honest scope: this surfaces and ranks candidates; it does not confirm planets
(confirmation needs radial-velocity or imaging follow-up we do not have). It reuses
the same coherence engine that beat hard rules for asteroid-tracklet vetting.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

# A real transiting planet, as physical penalties (the S_One E_total = sum of
# energy sectors). Each impostor signature ADDS incoherence energy; we sum rather
# than average so a single unambiguous disqualifier (e.g. a 9% "transit" can only
# be a star) is not diluted by the benign axes -- the failure mode of a mean.
# `tol` is the 1-sigma scale; `w` the sector weight.
PLANET_PENALTIES = {
    # depth is ONE-SIDED: planets are rarely deeper than ~3% (Jupiter on a Sun is
    # ~1%); deeper means a stellar companion. Shallow is fine, so no floor penalty.
    "depth_excess": {"knee": 0.03, "tol": 0.02, "w": 1.6},
    "secondary_ratio": {"knee": 0.0, "tol": 0.15, "w": 1.7},  # any secondary eclipse
    "odd_even": {"knee": 0.0, "tol": 0.18, "w": 1.2},  # alternating depths
    "v_shape": {"knee": 0.0, "tol": 0.30, "w": 1.4},  # grazing V floor
    "low_snr": {"knee": 7.0, "tol": 4.0, "w": 0.5},  # ONE-SIDED below 7
}
PLANET_COH_TAU = 0.42  # coherence = exp(-E/2) above this -> planet-candidate


@dataclass
class TransitCandidate:
    period: float
    t0: float
    duration: float
    depth: float
    snr: float
    features: dict = field(default_factory=dict)
    energy: float = float("nan")
    coherence: float = float("nan")
    verdict: str = "?"


# --------------------------------------------------------------------------- #
#  transit search (commodity: astropy Box Least Squares)
# --------------------------------------------------------------------------- #
def flatten_lightcurve(time, flux, *, window_days=0.4):
    """High-pass detrend: divide out a running median on a window LONGER than a
    transit, removing slow stellar variability (spots, M-dwarf activity) while
    preserving the brief transit dip. Recovers planets whose folded shape was
    contaminated by stellar variability (e.g. GJ 1132 b)."""
    from scipy.ndimage import median_filter

    t = np.asarray(time, float)
    f = np.asarray(flux, float)
    good = np.isfinite(t) & np.isfinite(f)
    t, f = t[good], f[good]
    if t.size < 100:
        return t, f
    cad = np.median(np.diff(np.sort(t)))
    win = max(11, int(window_days / max(cad, 1e-6)) | 1)  # odd window in points
    trend = median_filter(f, size=min(win, f.size // 2 * 2 - 1), mode="nearest")
    trend = np.where(np.abs(trend) > 1e-9, trend, np.nanmedian(f))
    return t, f / trend


def search_transits(
    time, flux, *, min_period=0.5, max_period=20.0, n_durations=8
) -> TransitCandidate | None:
    """Box Least Squares periodogram. Returns the strongest candidate (period, t0,
    duration, depth, SNR) or None. Flux is normalized to ~1 internally."""
    from astropy.timeseries import BoxLeastSquares

    t = np.asarray(time, float)
    f = np.asarray(flux, float)
    good = np.isfinite(t) & np.isfinite(f)
    t, f = t[good], f[good]
    if t.size < 100:
        return None
    f = f / np.nanmedian(f)
    durations = np.linspace(0.05, 0.4, n_durations)
    bls = BoxLeastSquares(t, f)
    span = t.max() - t.min()
    pg = bls.autopower(
        durations,
        minimum_period=min_period,
        maximum_period=min(max_period, span / 2.0),
        objective="snr",
    )
    i = int(np.argmax(pg.power))
    period = float(pg.period[i])
    t0 = float(pg.transit_time[i])
    duration = float(pg.duration[i])
    depth = float(pg.depth[i])
    stats = bls.compute_stats(period, duration, t0)
    snr = (
        float(pg.depth_snr[i])
        if hasattr(pg, "depth_snr")
        else (depth / (np.nanstd(f) / math.sqrt(max(stats["transit_times"].size, 1))))
    )
    return TransitCandidate(
        period=period, t0=t0, duration=duration, depth=max(depth, 0.0), snr=abs(snr)
    )


# --------------------------------------------------------------------------- #
#  vetting features (the physics that separates planets from impostors)
# --------------------------------------------------------------------------- #
def vet_features(time, flux, cand: TransitCandidate) -> dict:
    """Physical features of a folded candidate that distinguish a planet from an
    eclipsing binary / systematic. All derived from (time, flux) + the ephemeris."""
    t = np.asarray(time, float)
    f = np.asarray(flux, float)
    good = np.isfinite(t) & np.isfinite(f)
    t, f = t[good], f[good] / np.nanmedian(f[good])
    P, t0, dur = cand.period, cand.t0, max(cand.duration, 1e-3)
    phase = ((t - t0) / P + 0.5) % 1.0 - 0.5  # transit at phase 0
    half = 0.5 * dur / P
    intransit = np.abs(phase) < half
    oot = np.abs(phase) > 3 * half  # out of transit
    if intransit.sum() < 5 or oot.sum() < 20:
        return {}
    base = np.nanmedian(f[oot])
    scatter = np.nanstd(f[oot]) or 1e-6
    depth = base - np.nanmedian(f[intransit])
    depth = max(depth, 1e-6)
    # secondary eclipse near phase 0.5 -- significance-gated so photon noise on a
    # shallow real planet does not masquerade as a secondary (the low-SNR bug).
    sec_mask = np.abs(np.abs(phase) - 0.5) < half
    if sec_mask.sum() >= 5:
        sec_depth = base - np.nanmedian(f[sec_mask])
        sec_err = scatter / math.sqrt(sec_mask.sum())
        sec_depth = sec_depth if sec_depth > 3 * sec_err else 0.0  # 3-sigma gate
    else:
        sec_depth = 0.0
    secondary_ratio = max(sec_depth, 0.0) / depth
    # odd vs even transit depth (blended-binary-at-half-period tell), also gated.
    epoch = np.round((t - t0) / P)
    odd = intransit & (epoch.astype(int) % 2 == 1)
    even = intransit & (epoch.astype(int) % 2 == 0)
    if odd.sum() >= 3 and even.sum() >= 3:
        d_odd = base - np.nanmedian(f[odd])
        d_even = base - np.nanmedian(f[even])
        diff = abs(d_odd - d_even)
        oe_err = scatter * math.sqrt(1.0 / odd.sum() + 1.0 / even.sum())
        odd_even = diff / (abs(d_odd) + abs(d_even) + 1e-6) if diff > 3 * oe_err else 0.0
    else:
        odd_even = 0.0
    # V-shape vs U-shape: inner third should be ~as deep as the whole if flat (U).
    inner = np.abs(phase) < 0.33 * half
    if inner.sum() >= 3:
        d_inner = base - np.nanmedian(f[inner])
        v_shape = float(np.clip((d_inner / depth) - 1.0, 0.0, 2.0))  # V -> inner deeper
    else:
        v_shape = 0.0
    snr = depth / (scatter / math.sqrt(max(intransit.sum(), 1)))
    return {
        "depth": float(depth),
        "log_depth": math.log10(depth),
        "secondary_ratio": float(np.clip(secondary_ratio, 0, 2)),
        "odd_even": float(np.clip(odd_even, 0, 1)),
        "v_shape": float(v_shape),
        "snr": float(snr),
        "inv_snr": float(1.0 / max(snr, 1e-3)),
    }


# --------------------------------------------------------------------------- #
#  coherence vetting (Equation-of-One)
# --------------------------------------------------------------------------- #
def transit_incoherence(feats: dict) -> tuple[float, dict]:
    """Total incoherence energy of a candidate vs a real planet -- the S_One
    E_total as a SUM of physical penalty sectors (not a mean, so one unambiguous
    disqualifier dominates). Returns (E, per-sector contributions)."""
    p = PLANET_PENALTIES
    sectors = {}
    snr = feats.get("snr", 99.0)
    # SHAPE-FEATURE CONFIDENCE: secondary eclipse, odd-even, and V-shape are only
    # meaningful when the transit is well-measured. On a shallow / low-SNR signal the
    # folded shape is dominated by noise, so trusting it rejects genuine small planets
    # (the completeness failure on TOI-700 d / LHS 3844 b / GJ 1132 b). Ramp these
    # sectors from 0 below SNR~5 to full weight by SNR~20; depth (one-sided) and the
    # low-SNR penalty are always in force, so deep eclipsing binaries are still caught.
    f_shape = min(max((snr - 5.0) / 15.0, 0.0), 1.0)
    # depth excess: one-sided above the planet ceiling (always trusted -- robust)
    de = max(feats.get("depth", 0.0) - p["depth_excess"]["knee"], 0.0)
    sectors["depth_excess"] = p["depth_excess"]["w"] * (de / p["depth_excess"]["tol"]) ** 2
    for key in ("secondary_ratio", "odd_even", "v_shape"):
        sectors[key] = f_shape * p[key]["w"] * (feats.get(key, 0.0) / p[key]["tol"]) ** 2
    # low SNR: one-sided below the knee
    deficit = max(p["low_snr"]["knee"] - snr, 0.0)
    sectors["low_snr"] = p["low_snr"]["w"] * (deficit / p["low_snr"]["tol"]) ** 2
    return sum(sectors.values()), sectors


def vet_candidate(
    time, flux, cand: TransitCandidate, *, coherence_tau: float = PLANET_COH_TAU
) -> TransitCandidate:
    """Score a candidate's coherence with a real planet and label it. coherence =
    exp(-E/2); above `coherence_tau` -> planet-candidate, else likely-FP."""
    feats = vet_features(time, flux, cand)
    cand.features = feats
    if not feats:
        cand.energy, cand.coherence, cand.verdict = math.inf, 0.0, "insufficient-data"
        return cand
    E, sectors = transit_incoherence(feats)
    cand.features["energy_sectors"] = {k: round(v, 3) for k, v in sectors.items()}
    cand.energy = E
    cand.coherence = math.exp(-0.5 * E)
    cand.verdict = (
        "planet-candidate" if cand.coherence >= coherence_tau else "likely-false-positive"
    )
    return cand


def search_and_vet(time, flux, **kw) -> TransitCandidate | None:
    """Full pipeline: BLS search then coherence vetting. The one call a caller
    needs per light curve. BLS frequently locks onto 2x or 1/2 the true period
    (a real planet then folds like an eclipsing binary), so we also vet the
    half- and double-period aliases and keep the most planet-coherent solution --
    a completeness fix that recovers aliased small planets (e.g. LHS 3844 b)."""
    tau = kw.get("coherence_tau", 0.42)
    if kw.get("flatten", True):
        time, flux = flatten_lightcurve(time, flux)  # remove stellar variability
    cand = search_transits(
        time,
        flux,
        **{k: v for k, v in kw.items() if k in ("min_period", "max_period", "n_durations")},
    )
    if cand is None:
        return None
    best = vet_candidate(time, flux, cand, coherence_tau=tau)
    for factor in (0.5, 2.0):
        alt = TransitCandidate(
            period=cand.period * factor,
            t0=cand.t0,
            duration=cand.duration * factor,
            depth=cand.depth,
            snr=cand.snr,
        )
        if alt.period < 0.2:
            continue
        avc = vet_candidate(time, flux, alt, coherence_tau=tau)
        if avc.coherence > best.coherence:
            best = avc
    return best


# --------------------------------------------------------------------------- #
#  data access (free TESS light curves; lightkurve if present, else astroquery)
# --------------------------------------------------------------------------- #
def fetch_tess_lightcurve(target, *, sector=None, author="SPOC"):
    """(time, flux) for a TESS target (TIC id, name). Tries lightkurve, then a
    minimal astroquery+astropy fallback. Returns (np.array, np.array) or (None,None)."""
    try:
        import lightkurve as lk

        sr = lk.search_lightcurve(str(target), mission="TESS", author=author, sector=sector)
        if len(sr) == 0:
            return None, None
        lc = sr[0].download().remove_nans().normalize()
        return np.asarray(lc.time.value, float), np.asarray(lc.flux.value, float)
    except Exception:
        pass
    try:  # fallback: astroquery + FITS
        from astropy.io import fits
        from astroquery.mast import Observations

        obs = Observations.query_criteria(
            target_name=str(target), obs_collection="TESS", dataproduct_type="timeseries"
        )
        if len(obs) == 0:
            return None, None
        prod = Observations.get_product_list(obs[:1])
        lcp = prod[[str(u).endswith("lc.fits") for u in prod["productFilename"]]]
        if len(lcp) == 0:
            return None, None
        dl = Observations.download_products(lcp[:1])
        with fits.open(dl["Local Path"][0]) as h:
            d = h[1].data
            return np.asarray(d["TIME"], float), np.asarray(d["PDCSAP_FLUX"], float)
    except Exception:
        return None, None
