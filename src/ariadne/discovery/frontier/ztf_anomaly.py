"""Physics-coherence novelty scoring for transient / variable light curves -- the
triage layer for an alert firehose (ZTF now, LSST next).

LSST will issue ~10 million alerts a night. No human and no single classifier can
read that. The open problem is TRIAGE: which handful of objects are genuinely new
or weird enough to deserve a follow-up telescope. Trained classifiers are built to
recognize KNOWN classes; by construction they push a never-before-seen phenomenon
into whichever known bin fits least badly. The thing you most want to find is the
thing they are worst at surfacing.

A coherence energy inverts that. We define each well-understood class as a physical
basin (period, amplitude, periodicity strength, light-curve asymmetry -- set by the
physics of RR Lyrae, eclipsing binaries, Cepheids, long-period variables,
supernovae, AGN). An object's NOVELTY is its incoherence with the BEST-matching
known basin: low if it looks like something we understand, HIGH if it coheres with
nothing. No training set -- the basins come from physics -- which is exactly why it
can flag the never-seen.

Honest scope: this is a triage / ranking score, not a classifier and not a
discovery. A high score means "worth a human's attention," nothing more. The basins
are deliberately coarse physical priors; the point is to float the anomalous up, not
to label the ordinary.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from ..imaging.coherence_field import incoherence_energy

# Calibrated basins (fit to a real labeled ALeRCE sample) live here, committed
# alongside the code so the detector works out of the box; regenerate with
# scripts/calibrate_ztf.py. If absent we fall back to the physics priors below.
# Physics sets the STRUCTURE (which axes, the class topology); calibration sets the
# SCALES to the actual survey.
CALIBRATED_BASINS_PATH = Path(__file__).resolve().parent / "ztf_class_basins.json"

# Coarse physical priors for the common classes. mu = typical, sig = spread.
# Periodic classes constrain log_period; the stochastic/transient ones do NOT
# (they leave it unconstrained, so incoherence_energy skips it gracefully).
KNOWN_CLASS_BASINS = {
    "rr_lyrae": {
        "mu": {"log_period": math.log10(0.55), "amplitude": 0.8, "ls_power": 0.75, "abs_skew": 0.6},
        "sig": {"log_period": 0.25, "amplitude": 0.4, "ls_power": 0.15, "abs_skew": 0.5},
    },
    "eclipsing": {
        "mu": {"log_period": math.log10(1.5), "amplitude": 0.45, "ls_power": 0.6, "abs_skew": 1.2},
        "sig": {"log_period": 0.6, "amplitude": 0.35, "ls_power": 0.2, "abs_skew": 0.8},
    },
    "cepheid": {
        "mu": {"log_period": math.log10(7.0), "amplitude": 0.8, "ls_power": 0.75, "abs_skew": 0.5},
        "sig": {"log_period": 0.35, "amplitude": 0.4, "ls_power": 0.15, "abs_skew": 0.5},
    },
    "long_period": {
        "mu": {"log_period": math.log10(250.0), "amplitude": 1.6, "ls_power": 0.5, "abs_skew": 0.6},
        "sig": {"log_period": 0.5, "amplitude": 0.8, "ls_power": 0.25, "abs_skew": 0.6},
    },
    "supernova": {
        "mu": {"amplitude": 1.8, "ls_power": 0.12, "abs_skew": 1.4},
        "sig": {"amplitude": 1.0, "ls_power": 0.12, "abs_skew": 0.8},
    },
    "agn": {
        "mu": {"amplitude": 0.25, "ls_power": 0.15, "abs_skew": 0.4},
        "sig": {"amplitude": 0.2, "ls_power": 0.12, "abs_skew": 0.5},
    },
}
CLASS_W = {
    "log_period": 1.0,
    "amplitude": 1.0,
    "ls_power": 1.2,
    "abs_skew": 0.8,
    "eta": 1.0,
    "frac_bright": 1.0,
}
ANOMALY_TAU = 3.0  # min-basin incoherence above this -> flag as anomalous


# --------------------------------------------------------------------------- #
#  light-curve features (robust, physically meaningful)
# --------------------------------------------------------------------------- #
def lightcurve_features(time, mag, err=None, *, min_period=0.05, max_period=500.0) -> dict:
    """Period (Lomb-Scargle), amplitude, periodicity strength, and asymmetry of a
    light curve. Magnitudes (brighter = smaller). Robust to outliers and gaps."""
    from astropy.timeseries import LombScargle

    t = np.asarray(time, float)
    m = np.asarray(mag, float)
    good = np.isfinite(t) & np.isfinite(m)
    t, m = t[good], m[good]
    if t.size < 20:
        return {}
    span = t.max() - t.min()
    # robust amplitude: 5th-95th percentile span (mag)
    amplitude = float(np.percentile(m, 95) - np.percentile(m, 5))
    # Lomb-Scargle periodogram. Irregular survey sampling resolves sub-day periods
    # despite a slow median cadence, so the search frequency must NOT be capped by
    # the cadence -- bound it by the requested period range (down to ~0.05 d).
    try:
        ls = LombScargle(t, m, err[good] if err is not None else None)
        freq, power = ls.autopower(
            minimum_frequency=1.0 / min(max_period, span),
            maximum_frequency=1.0 / min_period,
            samples_per_peak=5,
        )
        i = int(np.argmax(power))
        period = float(1.0 / freq[i])
        ls_power = float(power[i])
    except Exception:
        period, ls_power = float("nan"), 0.0
    # asymmetry: |skew| of the magnitude distribution (transients/EBs are skewed)
    mc = m - np.mean(m)
    sd = np.std(m) or 1e-6
    skew = float(np.mean((mc / sd) ** 3))
    # von Neumann eta: mean squared successive difference / variance. Low (<<1) for
    # smooth/correlated curves (periodic variables, slow stochastic), ~2 for white
    # noise -- separates structured variability from noise-like.
    order = np.argsort(t)
    ms = m[order]
    eta = float(np.mean(np.diff(ms) ** 2) / (np.var(m) or 1e-9))
    # outburst fraction: points >1 mag BRIGHTER than the median (mag: brighter=smaller).
    # High for cataclysmic/nova outbursts and supernova/transient peaks, ~0 otherwise.
    med = float(np.median(m))
    frac_bright = float(np.mean(m < med - 1.0))
    feats = {
        "amplitude": amplitude,
        "ls_power": ls_power,
        "skew": skew,
        "abs_skew": abs(skew),
        "eta": eta,
        "frac_bright": frac_bright,
        "n": int(t.size),
        "baseline": float(span),
    }
    if math.isfinite(period) and ls_power > 0.05:
        feats["period"] = period
        feats["log_period"] = math.log10(period)
    return feats


# --------------------------------------------------------------------------- #
#  basin loading (calibrated if available, else physics priors)
# --------------------------------------------------------------------------- #
_BASINS_CACHE = None


def load_basins(*, force: bool = False) -> dict:
    """The class basins in use: calibrated (fit to real ALeRCE data) if the cache
    file exists, otherwise the physics priors. Loaded once and memoized."""
    global _BASINS_CACHE
    if _BASINS_CACHE is not None and not force:
        return _BASINS_CACHE
    if CALIBRATED_BASINS_PATH.exists():
        try:
            _BASINS_CACHE = json.loads(CALIBRATED_BASINS_PATH.read_text())["basins"]
            return _BASINS_CACHE
        except Exception:
            pass
    _BASINS_CACHE = KNOWN_CLASS_BASINS
    return _BASINS_CACHE


# --------------------------------------------------------------------------- #
#  novelty / anomaly score (Equation-of-One: incoherence with the best class)
# --------------------------------------------------------------------------- #
def anomaly_score(feats: dict, basins: dict | None = None) -> dict:
    """Novelty = incoherence with the BEST-matching known class. Returns the score,
    the closest class, and the per-class energies. High score = coheres with
    nothing known = worth a look."""
    if not feats:
        return {
            "score": float("inf"),
            "best_class": None,
            "verdict": "insufficient-data",
            "per_class": {},
        }
    basins = basins if basins is not None else load_basins()
    per = {name: incoherence_energy(feats, b, CLASS_W) for name, b in basins.items()}
    best_class = min(per, key=per.get)
    score = per[best_class]
    return {
        "score": float(score),
        "best_class": best_class,
        "verdict": "anomalous" if score >= ANOMALY_TAU else f"known:{best_class}",
        "per_class": {k: round(v, 2) for k, v in per.items()},
    }


def score_lightcurve(time, mag, err=None, *, basins=None, **kw) -> dict:
    """Features + novelty score for one light curve -- the single call a triage
    loop needs per object. `basins` defaults to the loaded (calibrated) set."""
    feats = lightcurve_features(
        time, mag, err, **{k: v for k, v in kw.items() if k in ("min_period", "max_period")}
    )
    out = anomaly_score(feats, basins)
    out["features"] = feats
    return out


# --------------------------------------------------------------------------- #
#  data access (free public ZTF light curves via the ALeRCE broker)
# --------------------------------------------------------------------------- #
def fetch_ztf_lightcurve(oid: str, *, timeout: float = 60.0):
    """(time, mag, err) for a ZTF object id from the ALeRCE broker (public, no
    auth). Uses the g band by default. Returns (None,None,None) on failure."""
    import requests

    try:
        r = requests.get(
            f"https://api.alerce.online/ztf/v1/objects/{oid}/lightcurve", timeout=timeout
        )
        r.raise_for_status()
        det = r.json().get("detections", [])
        g = [(d["mjd"], d["magpsf"], d.get("sigmapsf", 0.1)) for d in det if d.get("fid") == 1]
        if len(g) < 20:
            g = [(d["mjd"], d["magpsf"], d.get("sigmapsf", 0.1)) for d in det]
        if len(g) < 10:
            return None, None, None
        a = np.array(g, float)
        return a[:, 0], a[:, 1], a[:, 2]
    except Exception:
        return None, None, None


def alerce_oids(class_name: str, *, n: int = 20, timeout: float = 60.0) -> list[str]:
    """Object ids of the n highest-probability ALeRCE objects of a given class
    (e.g. RRL, CEP, E, LPV, QSO, SNIa). Free, no auth. [] on failure."""
    import requests

    try:
        r = requests.get(
            "https://api.alerce.online/ztf/v1/objects",
            params={
                "classifier": "lc_classifier",
                "class": class_name,
                "page_size": n,
                "order_by": "probability",
                "order_mode": "DESC",
            },
            timeout=timeout,
        )
        r.raise_for_status()
        return [o["oid"] for o in r.json().get("items", [])]
    except Exception:
        return []


def _robust_basin(feature_dicts: list[dict], axes) -> dict:
    """Median + robust (MAD) spread per axis over a labeled sample -> a basin.
    An axis is included only if enough sample members actually have it."""
    mu, sig = {}, {}
    for ax in axes:
        vals = np.array([f[ax] for f in feature_dicts if ax in f and np.isfinite(f[ax])], float)
        if vals.size >= 5:
            med = float(np.median(vals))
            mad = float(np.median(np.abs(vals - med)))
            mu[ax] = med
            sig[ax] = max(1.4826 * mad, 0.15 * abs(med) + 0.05)  # floor avoids 0
    return {"mu": mu, "sig": sig, "n": len(feature_dicts)}


def calibrate_from_alerce(
    classes=("RRL", "CEP", "E", "LPV", "QSO", "SNIa"), *, per_class: int = 25, save: bool = True
) -> dict:
    """Fit class basins to a REAL labeled ALeRCE sample and (optionally) save them
    to CALIBRATED_BASINS_PATH. Periodic classes get a log_period axis; all classes
    get amplitude / ls_power / abs_skew. This sets the basin scales to the actual
    survey while the physics still sets which axes matter."""
    basins = {}
    for cls in classes:
        feats = []
        for oid in alerce_oids(cls, n=per_class):
            t, m, e = fetch_ztf_lightcurve(oid)
            if t is not None:
                f = lightcurve_features(t, m, e)
                if f:
                    feats.append(f)
        if len(feats) < 5:
            continue
        axes = ["amplitude", "ls_power", "abs_skew", "eta", "frac_bright"]
        if np.mean([("log_period" in f) for f in feats]) > 0.5:
            axes = ["log_period"] + axes
        basins[cls] = _robust_basin(feats, axes)
    if save and basins:
        CALIBRATED_BASINS_PATH.parent.mkdir(parents=True, exist_ok=True)
        CALIBRATED_BASINS_PATH.write_text(json.dumps({"basins": basins}, indent=2))
        load_basins(force=True)
    return basins
