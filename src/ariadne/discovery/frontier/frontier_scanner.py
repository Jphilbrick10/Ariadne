"""Automated anomaly scanner over every tractable solar-system population.

Pulls each dynamical class from JPL SBDB (cached) and runs the whole anomaly suite
on it: orientation clustering (with an isotropic-shuffle significance), the plane
warp vs the invariable plane (for distant populations), the Southworth-Hawkins
D-criterion family search (with a chance-excess null), and exotic-origin scoring
(retrograde / polar / co-orbital / unbound = captured-object candidates). Every
flagged anomaly carries a significance so the ranked report does not cry wolf.

Honest scope: osculating elements (noisier than proper elements); the family search
is O(N^2) so it is capped and skipped/subsampled for the largest classes; the main
asteroid belt (~1.3M, families already well-mapped) is intentionally out of scope.
This is a triage layer -- it points, it does not confirm.
"""

from __future__ import annotations

import json
import math
import time
from pathlib import Path

import numpy as np
import requests

from .exotic_orbit_hunt import Orb, _count_below, dsh_matrix, exotic_score, find_pairs
from .planet_nine import INVARIABLE_PLANE, _pole
from .tno_clustering import circular_stats

SBDB = "https://ssd-api.jpl.nasa.gov/sbdb_query.api"
_CACHE = Path(__file__).resolve().parents[4] / "data" / "scan_cache"

# population catalogue: name -> JPL query (sb-class / sb-cdata). Tractable classes
# across the whole solar system; the ~1.3M main belt is deliberately excluded.
POPULATIONS = {
    "NEO Atira (interior)": "sb-class=IEO",
    "NEO Aten": "sb-class=ATE",
    "NEO Apollo": "sb-class=APO",
    "NEO Amor": "sb-class=AMO",
    "Mars-crossers": "sb-class=MCA",
    "Jupiter Trojans": "sb-class=TJN",
    "Centaurs": "sb-class=CEN",
    "Trans-Neptunian": "sb-class=TNO",
    "Hyperbolic asteroids": "sb-class=HYA",
    "Parabolic asteroids": "sb-class=PAA",
    "Comets (all)": "sb-kind=c",
    # derived distance/inclination cuts via constraint data
    "Scattered disk (a 50-200)": 'sb-class=TNO&sb-cdata={"AND":["a|LT|200","a|GT|50"]}',
    "Detached (q>40, a>150)": 'sb-cdata={"AND":["q|GT|40","a|GT|150"]}',
    "Retrograde TNO/Centaur (i>90)": 'sb-cdata={"AND":["i|GT|90","a|GT|5"]}',
}
FAMILY_CAP = 4000  # O(N^2) D-criterion only below this; subsample above


def fetch_population(query: str, *, use_cache: bool = True, timeout: float = 150.0) -> list[Orb]:
    """Pull a population's orbital elements from JPL SBDB (cached to disk)."""
    key = _CACHE / (str(abs(hash(query))) + ".json")
    if use_cache and key.exists():
        try:
            return [Orb(*r) for r in json.loads(key.read_text())]
        except Exception:
            pass
    try:
        url = f"{SBDB}?fields=full_name,a,e,i,om,w,q&{query}"
        data = requests.get(url, timeout=timeout).json().get("data", [])
    except Exception:
        return []
    objs = []
    for r in data:
        try:
            objs.append(
                Orb(
                    r[0].strip(),
                    float(r[1]),
                    float(r[2]),
                    float(r[3]),
                    float(r[4]),
                    float(r[5]),
                    float(r[6]),
                )
            )
        except (TypeError, ValueError):
            continue
    if objs:
        try:
            _CACHE.mkdir(parents=True, exist_ok=True)
            key.write_text(
                json.dumps([[o.name, o.a, o.e, o.inc, o.Omega, o.omega, o.q] for o in objs])
            )
        except Exception:
            pass
    return objs


# --------------------------------------------------------------------------- #
#  per-population analyses (each returns a value + a significance)
# --------------------------------------------------------------------------- #
def _orientation_sig(angles_deg, n_null=500, seed=0) -> dict:
    """Clustering resultant R + isotropic-shuffle significance for an angle set."""
    a = np.asarray(angles_deg, float) % 360
    n = a.size
    R = circular_stats(a)["R"]
    rng = np.random.default_rng(seed)
    null = np.array([circular_stats(rng.uniform(0, 360, n))["R"] for _ in range(n_null)])
    p = float(np.mean(null >= R))
    return {
        "R": round(float(R), 3),
        "mean_deg": round(circular_stats(a)["mean_deg"], 0),
        "p_value": round(p, 4),
        "significant": bool(p < 0.01),
    }


def _family_sig(objs: list[Orb], thr=0.05, n_null=150, seed=0) -> dict:
    """D-criterion close-pair count vs an orientation-shuffled null. The observed
    count EXCLUDES same-parent fragments (split comets etc.) -- a known artifact."""
    obs = len(find_pairs(objs, thr, exclude_same_parent=True))
    rng = np.random.default_rng(seed)
    Om = np.array([o.Omega for o in objs])
    om = np.array([o.omega for o in objs])
    base = [Orb(o.name, o.a, o.e, o.inc, 0, 0, o.q) for o in objs]
    null = np.empty(n_null)
    for k in range(n_null):
        po, pw = rng.permutation(Om), rng.permutation(om)
        for j, o in enumerate(base):
            o.Omega = po[j]
            o.omega = pw[j]
        null[k] = _count_below(dsh_matrix(base), thr)
    mean, std = float(null.mean()), float(null.std() or 1e-9)
    return {
        "observed_pairs": int(obs),
        "null_mean": round(mean, 2),
        "excess_sigma": round((obs - mean) / std, 2),
        "p_value": round(float(np.mean(null >= obs)), 4),
    }


def scan_population(name: str, objs: list[Orb], *, seed: int = 0) -> dict:
    """Run the full anomaly suite on one population."""
    res = {"name": name, "n": len(objs)}
    if len(objs) < 8:
        res["note"] = "too few objects"
        return res
    a = np.array([o.a for o in objs])
    e = np.array([o.e for o in objs])
    res["a_median"] = round(float(np.median(a)), 1)
    # exotic-origin scoring
    exo = [exotic_score(o) for o in objs if exotic_score(o)["flags"]]
    res["n_exotic"] = len(exo)
    res["exotic_frac"] = round(len(exo) / len(objs), 3)
    res["top_exotic"] = sorted(exo, key=lambda x: -x["score"])[:5]
    res["n_unbound"] = int((e >= 1.0).sum())
    # orientation clustering (node + perihelion), with significance
    res["node"] = _orientation_sig([o.Omega for o in objs], seed=seed)
    res["apsidal"] = _orientation_sig([(o.Omega + o.omega) for o in objs], seed=seed + 1)
    # plane warp vs invariable plane (distant populations only)
    if res["a_median"] > 20:
        inc = np.radians([o.inc for o in objs])
        Om = np.radians([o.Omega % 360 for o in objs])
        P = np.array([np.sin(inc) * np.sin(Om), -np.sin(inc) * np.cos(Om), np.cos(inc)]).T.mean(0)
        ph = P / np.linalg.norm(P)
        res["plane_warp_deg"] = round(
            math.degrees(math.acos(np.clip(np.dot(ph, _pole(*INVARIABLE_PLANE)), -1, 1))), 1
        )
    # D-criterion family search (tractable sizes only)
    if 8 <= len(objs) <= FAMILY_CAP:
        res["family"] = _family_sig(objs, seed=seed)
        res["top_pairs"] = [
            {"D": p["D"], "a": p["a"], "b": p["b"]} for p in find_pairs(objs, 0.05)[:5]
        ]
    else:
        res["family"] = {"note": f"skipped O(N^2) at n={len(objs)} (>cap {FAMILY_CAP})"}
    # flags. Orientation clustering is an ANOMALY only for distant populations,
    # where differential precession should isotropize orientations absent a
    # perturber. For inner populations non-uniformity is EXPECTED (selection +
    # secular resonances + Lagrange clouds), so we record it as structure, not a flag.
    distant = res["a_median"] > 30
    fam = res.get("family", {})
    res["flags"] = []
    res["expected_structure"] = []
    for chan, label in (("node", "node"), ("apsidal", "apsidal")):
        if res[chan]["significant"]:
            (res["flags"] if distant else res["expected_structure"]).append(
                f"{label} clustering p={res[chan]['p_value']}"
                + ("" if distant else " (resonant/selection -- expected)")
            )
    if isinstance(fam, dict) and fam.get("excess_sigma", 0) > 3:
        res["flags"].append(
            f"family excess {fam['excess_sigma']} sigma "
            f"({fam['observed_pairs']} pairs, fragments excluded)"
        )
    # exotic objects are reported individually (top_exotic), not as a population flag
    return res


def scan_all(
    populations: dict | None = None, *, use_cache: bool = True, log=lambda m: None
) -> dict:
    """Pull and scan every population; return the per-population results."""
    populations = populations or POPULATIONS
    results = {}
    for name, query in populations.items():
        objs = fetch_population(query, use_cache=use_cache)
        log(f"{name}: {len(objs)} objects")
        results[name] = scan_population(name, objs)
    return results
