"""Everything computable about Planet Nine from the real extreme-TNO data.

This is the honest, exhaustive version: not just "the cluster mean + 180", but
(1) rigorous significance of the clustering (Rayleigh AND Kuiper, on longitude of
perihelion, node, and the orbital pole), (2) the orbital orientation those imply,
(3) a genuine DYNAMICAL test -- secular-integrate the real TNOs under the giant
planets with and without a trial Planet Nine and measure whether P9 actually
preserves the clustering the giants alone disperse, (4) the implied sky position
(time-weighted along the orbit -> near aphelion) and apparent magnitude, hence why
it is unfound.

The governing honesty: the answer is only as strong as the clustering, and on the
current large sample the longitude-of-perihelion signal is NOT statistically
significant. We report that plainly. What is robust (the node clustering at
a>=150, the shared orbital plane) we use; what is not we flag. Mass and semimajor
axis are NOT determined by orientation -- only the dynamical test and the
literature N-body work constrain them, and we keep that line bright.
"""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass

import numpy as np

from .tno_clustering import TNO, circular_stats, extreme_population, fetch_tnos

OBLIQUITY = math.radians(23.439)
M_EARTH_SOLAR = 5.972e24 / 1.989e30  # Earth mass in solar masses


@dataclass
class P9Orbit:
    a_au: float
    e: float
    i_deg: float
    Omega_deg: float
    omega_deg: float
    mass_earths: float = 6.0

    @property
    def varpi_deg(self):
        return (self.Omega_deg + self.omega_deg) % 360.0

    @property
    def q_au(self):
        return self.a_au * (1 - self.e)

    @property
    def Q_au(self):
        return self.a_au * (1 + self.e)


# --------------------------------------------------------------------------- #
#  (1) rigorous clustering significance
# --------------------------------------------------------------------------- #
def rigorous_clustering(ext: list[TNO]) -> dict:
    """Rayleigh AND Kuiper uniformity tests for longitude of perihelion (varpi),
    node (Omega), argument of perihelion (omega), plus the orbital-pole clustering
    (the shared-plane signal). The honest significance of every directional claim."""
    from astropy.stats import kuiper

    n = len(ext)
    out = {"n": n}
    for name, vals in (
        ("varpi", [t.varpi for t in ext]),
        ("Omega", [t.Omega % 360 for t in ext]),
        ("omega", [t.omega % 360 for t in ext]),
    ):
        a = np.asarray(vals, float)
        st = circular_stats(a)
        try:
            _, kp = kuiper((a % 360) / 360.0)
        except Exception:
            kp = float("nan")
        out[name] = {
            "mean_deg": round(st["mean_deg"], 1),
            "R": round(st["R"], 3),
            "rayleigh_p": round(st["rayleigh_p"], 4),
            "kuiper_p": round(float(kp), 4),
            "significant": bool(st["rayleigh_p"] < 0.05 and kp < 0.05),
        }
    # orbital pole clustering -> the common plane (a real, robust signal)
    inc = np.radians([t.inc for t in ext])
    Om = np.radians([t.Omega % 360 for t in ext])
    poles = np.array([np.sin(inc) * np.sin(Om), -np.sin(inc) * np.cos(Om), np.cos(inc)]).T
    mp = poles.mean(axis=0)
    R = float(np.linalg.norm(mp))
    ph = mp / R
    out["orbit_pole"] = {
        "R": round(R, 3),
        "plane_tilt_deg": round(math.degrees(math.acos(min(1, ph[2]))), 1),
        "plane_node_deg": round(math.degrees(math.atan2(ph[0], -ph[1])) % 360, 1),
    }
    return out


# --------------------------------------------------------------------------- #
#  (2) implied P9 orientation
# --------------------------------------------------------------------------- #
def derive_orientation(ext: list[TNO]) -> dict:
    """P9 orbital orientation implied by the data: longitude of perihelion from
    apsidal anti-alignment (varpi_P9 = varpi_TNO + 180), inclined plane from the
    TNO common-plane tilt. Direction only -- never mass or distance."""
    sv = circular_stats([t.varpi for t in ext])
    rc = rigorous_clustering(ext)["orbit_pole"]
    varpi_p9 = (sv["mean_deg"] + 180.0) % 360.0
    return {
        "varpi_p9_deg": round(varpi_p9, 1),
        "anti_aligned_from_tno_varpi_deg": round(sv["mean_deg"], 1),
        "plane_tilt_deg": rc["plane_tilt_deg"],
        "plane_node_deg": rc["plane_node_deg"],
        "varpi_signal_significant": bool(sv["rayleigh_p"] < 0.05),
    }


# --------------------------------------------------------------------------- #
#  (3) orbit geometry -> sky position
# --------------------------------------------------------------------------- #
def orbit_xyz(o: P9Orbit, nu_deg: float) -> np.ndarray:
    """Heliocentric ECLIPTIC xyz (AU) of a point at true anomaly nu on orbit o."""
    nu = math.radians(nu_deg)
    i = math.radians(o.i_deg)
    Om = math.radians(o.Omega_deg)
    om = math.radians(o.omega_deg)
    r = o.a_au * (1 - o.e**2) / (1 + o.e * math.cos(nu))
    xp, yp = r * math.cos(nu), r * math.sin(nu)  # perifocal
    cO, sO, ci, si = math.cos(Om), math.sin(Om), math.cos(i), math.sin(i)
    co, so = math.cos(om), math.sin(om)
    x = (cO * co - sO * so * ci) * xp + (-cO * so - sO * co * ci) * yp
    y = (sO * co + cO * so * ci) * xp + (-sO * so + cO * co * ci) * yp
    z = (so * si) * xp + (co * si) * yp
    return np.array([x, y, z])


def ecliptic_to_radec(xyz: np.ndarray) -> tuple[float, float, float]:
    """Ecliptic xyz (AU) -> (RA deg, Dec deg, distance AU), rotating by obliquity."""
    x, y, z = xyz
    xe = x
    ye = y * math.cos(OBLIQUITY) - z * math.sin(OBLIQUITY)
    ze = y * math.sin(OBLIQUITY) + z * math.cos(OBLIQUITY)
    d = float(np.linalg.norm(xyz))
    ra = math.degrees(math.atan2(ye, xe)) % 360.0
    dec = math.degrees(math.asin(ze / d)) if d else 0.0
    return ra, dec, d


def sky_position(o: P9Orbit, n: int = 720) -> dict:
    """Time-weighted current-position estimate. With no constraint on where P9 is
    on its orbit, Kepler's second law says it spends the most time near aphelion --
    so the time-weighted mean position IS the best a-priori sky location, and the
    orbit traces a locus. Returns the aphelion point + the time-weighted centroid."""
    nus = np.linspace(0, 360, n, endpoint=False)
    rows = []
    for nu in nus:
        xyz = orbit_xyz(o, nu)
        ra, dec, d = ecliptic_to_radec(xyz)
        # time spent ~ r^2 / sqrt(...) ; dt/dnu proportional to r^2 (vis-viva/areal)
        w = d * d
        rows.append((ra, dec, d, w))
    rows = np.array(rows)
    w = rows[:, 3] / rows[:, 3].sum()
    # circular-weighted mean RA, weighted Dec
    ra_r = np.radians(rows[:, 0])
    ra_mean = math.degrees(math.atan2((np.sin(ra_r) * w).sum(), (np.cos(ra_r) * w).sum())) % 360
    dec_mean = float((rows[:, 1] * w).sum())
    # aphelion (nu=180): the single most-likely point
    apo = orbit_xyz(o, 180.0)
    ra_a, dec_a, d_a = ecliptic_to_radec(apo)
    return {
        "most_likely_RA_hours": round(ra_mean / 15.0, 2),
        "most_likely_Dec_deg": round(dec_mean, 1),
        "aphelion_RA_hours": round(ra_a / 15.0, 2),
        "aphelion_Dec_deg": round(dec_a, 1),
        "aphelion_distance_AU": round(d_a, 1),
        "constellation_hint": _constellation(ra_a, dec_a),
    }


def _constellation(ra_deg: float, dec_deg: float) -> str:
    """Very rough RA/Dec -> constellation region (for orientation, not astrometry)."""
    h = ra_deg / 15.0
    if 1.5 <= h < 3.5:
        return "Cetus/Aries region"
    if 3.5 <= h < 5.5:
        return "Taurus/Orion region (near galactic plane -- crowded, hard)"
    if 5.5 <= h < 7.5:
        return "Gemini/Orion region"
    if 0 <= h < 1.5 or h >= 22.5:
        return "Pisces region"
    return f"RA {h:.1f}h"


# --------------------------------------------------------------------------- #
#  (4) apparent magnitude / detectability
# --------------------------------------------------------------------------- #
def apparent_magnitude(dist_au: float, *, radius_earths: float = 2.5, albedo: float = 0.4) -> float:
    """Reflected-light V magnitude of P9 at heliocentric = geocentric ~ dist_au
    (it is far, so Earth-Sun baseline is negligible). Standard H + 5log10(r*d)."""
    R_km = radius_earths * 6371.0
    # absolute magnitude from radius + albedo (H = 5 log10(1329/(2R_km/1000)/sqrt(p)))
    diam_km = 2 * R_km
    H = 5 * math.log10(1329.0 / (diam_km * math.sqrt(albedo))) if albedo > 0 else 99
    return H + 5 * math.log10(max(dist_au, 1e-3) * max(dist_au - 1.0, 1e-3))


# --------------------------------------------------------------------------- #
#  (5) the dynamical test (secular integration of the REAL TNOs +/- P9)
# --------------------------------------------------------------------------- #
def secular_dispersal_test(
    ext: list[TNO],
    p9: P9Orbit | None,
    *,
    span_gyr: float = 1.5,
    dt_yr: float = 1.5e5,
    snapshots: int = 25,
    n_ring: int = 48,
    n_tp: int = 48,
    epoch: str = "2026-01-01T00:00:00",
) -> dict:
    """Secular-integrate the real extreme TNOs for span_gyr under the giant planets
    (p9=None) or giants+P9, recording the clustering resultant R of longitude of
    perihelion and node vs time. Giants alone differentially precess the orbits and
    DISPERSE any clustering; a real shepherding P9 should PRESERVE it. da/dt is
    reported as the conservation check (doubly-averaged -> ~0)."""
    from ariadne.data.constants import GM_SUN
    from ariadne.dynamics.secular_avg import giant_rings, integrate_secular

    elems = [
        dict(a_au=t.a, e=t.e, i_deg=t.inc, Omega_deg=t.Omega % 360, omega_deg=t.omega % 360)
        for t in ext
    ]
    extra = None
    if p9 is not None:
        extra = [
            dict(
                name="P9",
                gm=p9.mass_earths * M_EARTH_SOLAR * GM_SUN,
                a_au=p9.a_au,
                e=p9.e,
                i=p9.i_deg,
                Omega=p9.Omega_deg,
                omega=p9.omega_deg,
            )
        ]
    rings, _ = giant_rings(epoch, n=n_ring, extra=extra)
    n_steps = int(span_gyr * 1e9 / dt_yr)
    rec_every = max(1, n_steps // snapshots)
    out = integrate_secular(
        [dict(e) for e in elems], rings, dt_yr, n_steps, record_every=rec_every, n_tp=n_tp
    )

    def R(state, key):
        ang = np.radians(
            [
                s.get(key, (s["Omega_deg"] + s["omega_deg"]) % 360 if key == "varpi_deg" else 0)
                for s in state
            ]
        )
        return float(np.hypot(np.cos(ang).mean(), np.sin(ang).mean()))

    times = out.get("times_yr", np.array([0.0]))
    hist = out.get("history", [out["elements"]])
    return {
        "label": "giants+P9" if p9 else "giants-only",
        "times_gyr": [round(t / 1e9, 3) for t in times],
        "R_varpi": [round(R(h, "varpi_deg"), 3) for h in hist],
        "R_Omega": [round(R(h, "Omega_deg"), 3) for h in hist],
        "da_au_per_yr_max": float(out["da_au_per_yr_max"]),
        "span_gyr": span_gyr,
        "n_tno": len(ext),
    }


# --------------------------------------------------------------------------- #
#  (6) open-minded multi-population anomaly scan (no P9 model assumed)
# --------------------------------------------------------------------------- #
def _kuiper_p(angles_deg) -> float:
    from astropy.stats import kuiper

    try:
        return float(kuiper((np.asarray(angles_deg, float) % 360) / 360.0)[1])
    except Exception:
        return float("nan")


def _bin_stats(arr: np.ndarray) -> dict:
    """Clustering of varpi and Omega for an (N, >=7) array [a,e,i,Om,om,q,varpi]."""
    if len(arr) < 8:
        return {"n": len(arr)}
    sv = circular_stats(arr[:, 6])
    so = circular_stats(arr[:, 3] % 360)
    return {
        "n": len(arr),
        "varpi": {
            "mean": round(sv["mean_deg"], 0),
            "R": round(sv["R"], 2),
            "kuiper_p": round(_kuiper_p(arr[:, 6]), 3),
        },
        "Omega": {
            "mean": round(so["mean_deg"], 0),
            "R": round(so["R"], 2),
            "kuiper_p": round(_kuiper_p(arr[:, 3]), 3),
        },
    }


def outer_system_anomaly_scan(
    tnos: list[TNO] | None = None, *, use_cache: bool = True, include_comets: bool = True
) -> dict:
    """Pool every distant population and ask what sticks out -- WITHOUT assuming the
    textbook Planet Nine. Bins the whole TNO catalogue by semimajor axis; flags the
    MODEL-INDEPENDENT anomaly (detached high-perihelion objects Neptune cannot
    emplace); checks whether independent bins agree on a node direction; counts the
    retrograde/polar population; and tests the Oort-cloud comets for any orientation
    imprint. The honest map a perturber hypothesis (any kind) must explain."""
    tn = tnos if tnos is not None else fetch_tnos(use_cache=use_cache)
    A = np.array([[t.a, t.e, t.inc, t.Omega % 360, t.omega % 360, t.q, t.varpi] for t in tn])
    a, e, i, Om, w, q, vp = A.T
    report = {"n_tno": len(tn)}
    # clustering vs semimajor axis (q>30)
    bins = {}
    for lo, hi in [(50, 100), (100, 150), (150, 250), (250, 500), (500, 1e12)]:
        key = f"a{lo}-{'inf' if hi > 1e9 else hi}"
        bins[key] = _bin_stats(A[(a >= lo) & (a < hi) & (q > 30)])
    report["a_bins"] = bins
    # node-direction consistency across independent significant bins
    nodes = [
        b["Omega"]["mean"]
        for b in bins.values()
        if b.get("n", 0) >= 8 and b["Omega"]["kuiper_p"] < 0.05
    ]
    if len(nodes) >= 2:
        nr = circular_stats(nodes)
        report["node_consistency"] = {
            "significant_bins": len(nodes),
            "mean_node_deg": round(nr["mean_deg"], 0),
            "spread_R": round(nr["R"], 2),
        }
    # the robust, model-independent anomaly: detached objects Neptune can't make
    det = A[(q > 40) & (a > 150)]
    report["detached_anomaly"] = {
        "n_decoupled_q_gt_40_a_gt_150": len(det),
        "max_perihelion_AU": round(float(det[:, 5].max()), 0) if len(det) else 0,
        "max_a_AU": round(float(det[:, 0].max()), 0) if len(det) else 0,
        "note": (
            "Neptune at 30 AU cannot lift a perihelion past ~38 AU; these "
            "REQUIRE an external agent (planet, past stellar flyby, birth "
            "cluster, or primordial-disk self-gravity) -- independent of any "
            "orbit clustering."
        ),
        "stats": _bin_stats(det),
    }
    # weird dynamical classes
    report["retrograde_polar"] = {
        "n_i_gt_90": int(((i > 90) & (a > 30)).sum()),
        "n_high_i_40_90": int(((i > 40) & (i <= 90) & (a > 50)).sum()),
    }
    # comets / Oort cloud: any orientation imprint?
    if include_comets:
        try:
            import requests

            d = requests.get(
                "https://ssd-api.jpl.nasa.gov/sbdb_query.api"
                "?fields=full_name,a,e,i,om,w,q&sb-kind=c",
                timeout=90,
            ).json()
            C = np.array(
                [[float(r[1]), float(r[4]), float(r[5])] for r in d["data"] if r[1] not in (None,)]
            )  # a, Om, om
            cvp = (C[:, 1] + C[:, 2]) % 360
            oort = cvp[C[:, 0] > 1000]
            report["comets"] = {
                "n_total": len(C),
                "n_oort_a_gt_1000": len(oort),
                "oort_varpi_R": round(float(circular_stats(oort)["R"]), 2)
                if len(oort) >= 8
                else None,
                "note": "isotropic Oort cloud -> no strong perturber imprint",
            }
        except Exception:
            report["comets"] = {"error": "comet query failed"}
    return report


# --------------------------------------------------------------------------- #
#  (7) discrimination tests -- try to KILL the signal before believing it
# --------------------------------------------------------------------------- #
INVARIABLE_PLANE = (1.578, 107.58)  # (inclination deg, node deg) vs J2000 ecliptic


def _pole(i_deg, Om_deg):
    i, O = math.radians(i_deg), math.radians(Om_deg)
    return np.array([math.sin(i) * math.sin(O), -math.sin(i) * math.cos(O), math.cos(i)])


def plane_warp(ext: list[TNO]) -> dict:
    """Angle between the distant-TNO mean orbital plane and the solar system's
    INVARIABLE plane. Distant objects should settle toward the invariable/Laplace
    plane; a large deviation is a genuine warp (the real, if modest, anomaly) --
    not just objects sitting in the plane they are expected to."""
    inc = np.radians([t.inc for t in ext])
    Om = np.radians([t.Omega % 360 for t in ext])
    P = np.array([np.sin(inc) * np.sin(Om), -np.sin(inc) * np.cos(Om), np.cos(inc)]).T
    mp = P.mean(axis=0)
    ph = mp / np.linalg.norm(mp)
    inv = _pole(*INVARIABLE_PLANE)
    return {
        "n": len(ext),
        "mean_plane_tilt_deg": round(math.degrees(math.acos(min(1, ph[2]))), 1),
        "mean_plane_node_deg": round(math.degrees(math.atan2(ph[0], -ph[1])) % 360, 1),
        "angle_from_invariable_deg": round(
            math.degrees(math.acos(np.clip(np.dot(ph, inv), -1, 1))), 1
        ),
    }


def planetary_perturbation_constraint(
    p9_mass_earths: float,
    p9_dist_au: float,
    *,
    mission_years: float = 13.0,
    ranging_precision_m: float = 100.0,
) -> dict:
    """Can the KNOWN planets' tracking detect or constrain a distant perturber?

    A Planet Nine at hundreds of AU never approaches the planets (its perihelion is
    >> Neptune), so there are no close encounters -- but its TIDAL pull (the
    Sun-vs-planet differential acceleration, ~2 G M9 a_planet / d^3) slowly drifts
    each planet's orbit. Over a precision-tracking mission that drift accumulates as
    ~1/2 a t^2. Saturn, ranged by Cassini to ~100 m, is the gold standard: a drift
    above that precision means planetary ephemerides CONSTRAIN the perturber (the
    method of Fienga et al. 2020). Honest caveat: much of a steady drift is absorbed
    by re-fitting the planet's orbit; the real, weaker constraint is on the
    direction/time-dependent part, so treat 'detectable' as an upper bound on
    sensitivity, not a guaranteed detection."""
    M9 = p9_mass_earths * M_EARTH_SOLAR * 1.989e30  # kg
    G = 6.674e-11
    au = 1.495978707e11
    d = p9_dist_au * au
    t = mission_years * 3.156e7
    planets = {"Jupiter": 5.20, "Saturn": 9.54, "Uranus": 19.2, "Neptune": 30.1}
    out = {}
    for name, a_au in planets.items():
        tidal = 2 * G * M9 * (a_au * au) / d**3
        drift_m = 0.5 * tidal * t**2
        out[name] = {
            "tidal_accel_ms2": tidal,
            "drift_m": round(drift_m, 0),
            "detectable": drift_m > ranging_precision_m,
        }
    return {
        "p9_mass_earths": p9_mass_earths,
        "p9_dist_au": p9_dist_au,
        "p9_perihelion_note": "never approaches the planets; constraint is tidal, not close-encounter",
        "ranging_precision_m": ranging_precision_m,
        "per_planet": out,
        "saturn_constrains": out["Saturn"]["detectable"],
    }


def selection_bias_null_test(
    ext: list[TNO], *, n_trials: int = 200, b_cut: float = 15.0, pool: int = 600, seed: int = 0
) -> dict:
    """Can galactic-plane selection ALONE fake the observed clustering? Build
    populations with ISOTROPIC orientations (random Omega, omega) but the real a/e/i
    distribution, keep only those detectable away from the galactic plane (|b|>b_cut
    at perihelion, where they are brightest), and compare the induced clustering R to
    the observed. P(null>=obs) > 0.05 means the observed clustering is CONSISTENT
    with selection bias -- not, by itself, evidence of a perturber."""
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    rng = np.random.default_rng(seed)
    n = len(ext)
    obs = {
        "R_Omega": circular_stats([t.Omega % 360 for t in ext])["R"],
        "R_varpi": circular_stats([t.varpi for t in ext])["R"],
    }
    aa = np.array([t.a for t in ext])
    ee = np.array([t.e for t in ext])
    ii = np.array([t.inc for t in ext])
    Ro_null, Rv_null = [], []
    for _ in range(n_trials):
        idx = rng.integers(0, n, pool)  # resample real a/e/i marginals
        a, e, i = aa[idx], ee[idx], ii[idx]
        Om = rng.uniform(0, 360, pool)
        om = rng.uniform(0, 360, pool)  # isotropic orientation
        ra = np.empty(pool)
        dec = np.empty(pool)
        for k in range(pool):
            ra[k], dec[k], _ = ecliptic_to_radec(
                orbit_xyz(P9Orbit(a[k], e[k], i[k], Om[k], om[k]), 0.0)
            )
        b = SkyCoord(ra * u.deg, dec * u.deg).galactic.b.deg
        keep = np.where(np.abs(b) > b_cut)[0]
        if len(keep) < n:
            continue
        sel = rng.choice(keep, n, replace=False)
        Ro_null.append(circular_stats(Om[sel] % 360)["R"])
        Rv_null.append(circular_stats((Om[sel] + om[sel]) % 360)["R"])
    Ro_null, Rv_null = np.array(Ro_null), np.array(Rv_null)
    return {
        "n": n,
        "n_trials": len(Ro_null),
        "Omega": {
            "observed_R": round(obs["R_Omega"], 3),
            "null_median_R": round(float(np.median(Ro_null)), 3),
            "p_value": round(float(np.mean(Ro_null >= obs["R_Omega"])), 3),
            "survives_bias": bool(np.mean(Ro_null >= obs["R_Omega"]) < 0.05),
        },
        "varpi": {
            "observed_R": round(obs["R_varpi"], 3),
            "null_median_R": round(float(np.median(Rv_null)), 3),
            "p_value": round(float(np.mean(Rv_null >= obs["R_varpi"])), 3),
            "survives_bias": bool(np.mean(Rv_null >= obs["R_varpi"]) < 0.05),
        },
    }
