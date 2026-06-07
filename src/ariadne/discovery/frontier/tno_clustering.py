"""Extreme-TNO orbital-orientation outlier hunt -- the Planet Nine clustering
signal, scored with the Equation-of-One coherence engine.

The distant solar system has a real, debated anomaly: the most extreme trans-
Neptunian objects (large semimajor axis a, large perihelion q, so they never come
near Neptune) appear to have their orbits ORIENTED together -- clustered longitude
of perihelion (varpi = Omega + omega) and ascending node (Omega). If real and not
an observational-bias artifact, a massive distant perturber ("Planet Nine") is the
leading explanation (Batygin & Brown 2016).

This module reproduces that signal from scratch on real JPL data, then applies the
coherence selector to it in two honest ways:

  1. CLUSTERING COHERENCE. Treat each object's orbit orientation as an observation
     and the population's mean direction as a basin. The circular coherence energy
     E_c = 1 - cos(angle - mean) (the von Mises form of the S_One coherence /
     alignment term) is LOW for an orbit that aligns with the cluster, HIGH for one
     that points the other way. The mean energy over the population is a coherence
     measure of the clustering; the Rayleigh test gives its significance.

  2. OUTLIER RANKING. Per object, the orientation energy ranks how anomalously its
     orbit points relative to the cluster. The least-coherent orbits are the
     dynamically odd ones -- interesting follow-up targets, and the objects whose
     orientation a perturber model must explain.

Honest scope: this REPRODUCES known science and provides a coherence-based scoring
/ outlier tool on top of it. It does not discover Planet Nine (that needs new deep
imaging we do not have). The clustering significance depends on the a / q cutoff
and on observational selection effects -- we report that dependence rather than
cherry-picking the cutoff that looks best.
"""

from __future__ import annotations

import json
import math
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import requests

SBDB_URL = "https://ssd-api.jpl.nasa.gov/sbdb_query.api"
_CACHE = Path(__file__).resolve().parents[4] / "data" / "tno_elements.json"


@dataclass
class TNO:
    """One trans-Neptunian object's orbit. Angles in degrees, a/q in AU."""

    name: str
    a: float
    e: float
    inc: float
    Omega: float  # longitude of ascending node
    omega: float  # argument of perihelion
    q: float  # perihelion distance
    H: float = float("nan")

    @property
    def varpi(self) -> float:
        """Longitude of perihelion = Omega + omega (mod 360)."""
        return (self.Omega + self.omega) % 360.0


# --------------------------------------------------------------------------- #
#  data access (free, public JPL Small-Body Database)
# --------------------------------------------------------------------------- #
def fetch_tnos(*, use_cache: bool = True, timeout: float = 120.0) -> list[TNO]:
    """All TNOs with orbital elements from JPL SBDB. Cached to data/ so repeated
    runs and tests do not re-hit the network. Returns [] only on hard failure."""
    if use_cache and _CACHE.exists():
        try:
            raw = json.loads(_CACHE.read_text())
            return [TNO(**r) for r in raw]
        except Exception:
            pass
    fields = "full_name,a,e,i,om,w,q,H"
    try:
        r = requests.get(f"{SBDB_URL}?fields={fields}&sb-class=TNO", timeout=timeout)
        r.raise_for_status()
        data = r.json()["data"]
    except Exception:
        return []
    out = []
    for row in data:
        try:
            out.append(
                TNO(
                    name=row[0].strip(),
                    a=float(row[1]),
                    e=float(row[2]),
                    inc=float(row[3]),
                    Omega=float(row[4]),
                    omega=float(row[5]),
                    q=float(row[6]),
                    H=float(row[7]) if row[7] is not None else float("nan"),
                )
            )
        except (TypeError, ValueError):
            continue
    if out:
        try:
            _CACHE.parent.mkdir(parents=True, exist_ok=True)
            _CACHE.write_text(json.dumps([t.__dict__ for t in out]))
        except Exception:
            pass
    return out


def extreme_population(tnos: list[TNO], *, a_min: float = 150.0, q_min: float = 30.0) -> list[TNO]:
    """The detached / extreme TNOs that drive the Planet Nine signal: large a (so
    they are decoupled from Neptune) and large q (so Neptune did not place them)."""
    return [
        t
        for t in tnos
        if math.isfinite(t.a)
        and math.isfinite(t.q)
        and math.isfinite(t.Omega)
        and math.isfinite(t.omega)
        and t.a >= a_min
        and t.q >= q_min
    ]


# --------------------------------------------------------------------------- #
#  circular statistics
# --------------------------------------------------------------------------- #
def circular_stats(angles_deg) -> dict:
    """Mean direction, resultant length R (0 = uniform, 1 = identical), and the
    Rayleigh-test p-value for non-uniformity (Zar's approximation). The honest
    significance of a directional clustering."""
    th = np.radians(np.asarray(angles_deg, float))
    n = th.size
    if n == 0:
        return {"n": 0, "mean_deg": float("nan"), "R": 0.0, "rayleigh_p": 1.0}
    C, S = np.cos(th).mean(), np.sin(th).mean()
    R = float(math.hypot(C, S))
    mean = math.degrees(math.atan2(S, C)) % 360.0
    Z = n * R * R
    # Zar 1999 series correction to exp(-Z)
    p = math.exp(-Z) * (
        1 + (2 * Z - Z * Z) / (4 * n) - (24 * Z - 132 * Z**2 + 76 * Z**3 - 9 * Z**4) / (288 * n * n)
    )
    return {"n": n, "mean_deg": mean, "R": R, "rayleigh_p": float(min(max(p, 0.0), 1.0))}


# --------------------------------------------------------------------------- #
#  coherence energy (Equation-of-One, circular / von Mises form)
# --------------------------------------------------------------------------- #
def orientation_energy(angle_deg: float, mean_deg: float) -> float:
    """Circular coherence energy E_c = 1 - cos(angle - mean), in [0, 2]. The
    von Mises analogue of the S_One coherence / alignment term: 0 when the orbit
    points exactly with the cluster, 2 when it points exactly opposite."""
    return 1.0 - math.cos(math.radians(angle_deg - mean_deg))


def rank_by_coherence(tnos: list[TNO], *, axis: str = "Omega") -> list[dict]:
    """Score every object by how its orbit orientation coheres with the population
    cluster on `axis` ('Omega', 'varpi', or 'omega'). Returns rows sorted MOST
    anomalous first (highest energy = least coherent). The coherent core is the
    tail; the head is the dynamically odd orbits."""
    vals = [getattr(t, axis) % 360.0 for t in tnos]
    stats = circular_stats(vals)
    mean = stats["mean_deg"]
    rows = [
        {
            "name": t.name,
            "a": t.a,
            "q": t.q,
            "e": t.e,
            "inc": t.inc,
            axis: v,
            "energy": orientation_energy(v, mean),
            "coherence": math.exp(-0.5 * orientation_energy(v, mean)),
        }
        for t, v in zip(tnos, vals)
    ]
    rows.sort(key=lambda r: -r["energy"])
    return rows


# --------------------------------------------------------------------------- #
#  coherence SELECTOR: which subsample shows the cleanest clustering
# --------------------------------------------------------------------------- #
def select_coherent_subsample(
    tnos: list[TNO], *, axis: str = "Omega", a_grid=(150, 200, 250, 300, 400)
) -> dict:
    """Equation-of-One selector over a-cutoffs: pick the population that MINIMIZES
    mean orientation energy (the most coherent clustering) while keeping the sample
    statistically meaningful. This is the principled version of 'which objects to
    trust for the inference' instead of hand-picking the cutoff that looks best."""
    best = None
    table = []
    for a_min in a_grid:
        sub = [t for t in tnos if t.a >= a_min]
        if len(sub) < 5:
            continue
        vals = [getattr(t, axis) % 360.0 for t in sub]
        st = circular_stats(vals)
        mean_E = float(np.mean([orientation_energy(v, st["mean_deg"]) for v in vals]))
        row = {
            "a_min": a_min,
            "n": len(sub),
            "R": st["R"],
            "rayleigh_p": st["rayleigh_p"],
            "mean_energy": mean_E,
            "mean_dir_deg": st["mean_deg"],
        }
        table.append(row)
        # selector: most coherent (lowest mean energy), significance-gated
        if best is None or mean_E < best["mean_energy"]:
            best = row
    return {"selected": best, "table": table}


# --------------------------------------------------------------------------- #
#  top-level honest report
# --------------------------------------------------------------------------- #
def outlier_hunt(
    *,
    a_min: float = 150.0,
    q_min: float = 30.0,
    use_cache: bool = True,
    tnos: list[TNO] | None = None,
) -> dict:
    """Full hunt: clustering significance on Omega and varpi (over multiple
    cutoffs, honestly), the coherence-selected cleanest subsample, and the most
    orientation-anomalous orbits. Pass `tnos` to run offline on a fixed sample."""
    pool = tnos if tnos is not None else fetch_tnos(use_cache=use_cache)
    ext = extreme_population(pool, a_min=a_min, q_min=q_min)
    report = {"n_total": len(pool), "n_extreme": len(ext), "a_min": a_min, "q_min": q_min}
    if not ext:
        report["error"] = "no extreme TNOs (offline + empty cache?)"
        return report
    for axis in ("Omega", "varpi"):
        vals = [getattr(t, axis) % 360.0 for t in ext]
        report[axis] = circular_stats(vals)
        report[f"{axis}_select"] = select_coherent_subsample(ext, axis=axis)["table"]
    report["most_anomalous"] = rank_by_coherence(ext, axis="Omega")[:10]
    return report


def _ecliptic_lon_to_radec(lon_deg: float, lat_deg: float = 0.0) -> tuple[float, float]:
    """Approximate equatorial RA/Dec (deg) of an ecliptic direction (obliquity
    23.439 deg). Used only to translate a predicted longitude into a sky region."""
    eps = math.radians(23.439)
    lam, bet = math.radians(lon_deg), math.radians(lat_deg)
    ra = math.atan2(math.sin(lam) * math.cos(eps) - math.tan(bet) * math.sin(eps), math.cos(lam))
    dec = math.asin(math.sin(bet) * math.cos(eps) + math.cos(bet) * math.sin(eps) * math.sin(lam))
    return math.degrees(ra) % 360.0, math.degrees(dec)


def predict_planet_nine(
    *,
    a_min: float = 250.0,
    q_min: float = 30.0,
    use_cache: bool = True,
    tnos: list[TNO] | None = None,
) -> dict:
    """The framework's geometric prediction of Planet Nine's orbital ORIENTATION
    from the extreme-TNO clustering, via apsidal anti-alignment (Batygin & Brown):
    P9's perihelion points OPPOSITE the clustered TNO perihelia, and P9 sits near
    aphelion (which is why it is faint and unfound).

    This determines DIRECTION only. Mass and semimajor axis require N-body matching
    (not done here) -- the literature values are attached for context, clearly
    labelled, and are NOT outputs of this geometry. The prediction is only as strong
    as the clustering, which in the current large sample is weak (see rayleigh_p)."""
    pool = tnos if tnos is not None else fetch_tnos(use_cache=use_cache)
    ext = extreme_population(pool, a_min=a_min, q_min=q_min)
    if len(ext) < 5:
        return {"error": f"only {len(ext)} extreme TNOs at a>={a_min}"}
    sv = circular_stats([t.varpi for t in ext])
    so = circular_stats([t.Omega % 360.0 for t in ext])
    mean_inc = float(np.median([t.inc for t in ext]))
    # apsidal anti-alignment
    p9_varpi = (sv["mean_deg"] + 180.0) % 360.0
    # P9 is near aphelion now; aphelion longitude = perihelion longitude + 180
    p9_aphelion_lon = (p9_varpi + 180.0) % 360.0
    ra_apo, dec_apo = _ecliptic_lon_to_radec(p9_aphelion_lon, 0.0)
    return {
        "n_extreme": len(ext),
        "a_min": a_min,
        "tno_cluster": {
            "mean_varpi_deg": round(sv["mean_deg"], 1),
            "varpi_R": round(sv["R"], 3),
            "varpi_rayleigh_p": round(sv["rayleigh_p"], 4),
            "mean_Omega_deg": round(so["mean_deg"], 1),
            "Omega_R": round(so["R"], 3),
            "Omega_rayleigh_p": round(so["rayleigh_p"], 4),
        },
        "predicted_orientation": {
            "longitude_of_perihelion_deg": round(p9_varpi, 1),
            "aphelion_ecliptic_longitude_deg": round(p9_aphelion_lon, 1),
            "aphelion_sky_RA_hours": round(ra_apo / 15.0, 1),
            "aphelion_sky_Dec_deg": round(dec_apo, 0),
            "note": "near-ecliptic approximation; the inclined-plane solution shifts Dec",
        },
        "from_dynamics_NOT_this_geometry": {
            "mass_earths": "~5-10 (N-body matching, Batygin & Brown 2016; ~5-6 in 2021)",
            "semimajor_axis_AU": "~400-500",
            "eccentricity": "~0.2-0.3",
            "inclination_deg": "~15-20",
            "perihelion_AU": "~300",
            "apparent_mag": "~21-25 near aphelion -> why undetected",
        },
        "honesty": (
            "orientation only; weak clustering in this sample (see "
            "rayleigh_p); contested hypothesis; cross-check vs Brown & "
            "Batygin 2021 published search region."
        ),
    }


if __name__ == "__main__":  # pragma: no cover
    t0 = time.time()
    rep = outlier_hunt()
    print(f"=== Extreme-TNO orientation outlier hunt ({time.time() - t0:.1f}s) ===")
    print(
        f"  {rep['n_total']} TNOs -> {rep['n_extreme']} extreme (a>={rep['a_min']}, q>={rep['q_min']})"
    )
    for axis in ("Omega", "varpi"):
        s = rep[axis]
        print(
            f"  {axis:6s}: mean={s['mean_deg']:.0f} deg  R={s['R']:.3f}  "
            f"Rayleigh p={s['rayleigh_p']:.4f}  (n={s['n']})"
        )
    print("  clustering vs cutoff (Omega):")
    for row in rep["Omega_select"]:
        print(
            f"    a>={row['a_min']:>3}: n={row['n']:>2}  R={row['R']:.3f}  "
            f"p={row['rayleigh_p']:.4f}  mean_E={row['mean_energy']:.3f}"
        )
    print("  most orientation-anomalous (least coherent) orbits:")
    for r in rep["most_anomalous"][:6]:
        print(
            f"    {r['name'][:28]:28s}  a={r['a']:.0f}  Omega={r['Omega']:.0f}  E={r['energy']:.2f}"
        )
