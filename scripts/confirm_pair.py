"""Confirm (or refute) a candidate asteroid pair by backward N-body integration.

Two objects on nearly identical orbits but different mean anomalies are the
signature of an asteroid PAIR (a fragment that split). The test: integrate both
backward under the full planetary system; a genuine pair's orbits CONVERGE (the
Southworth-Hawkins D and the node/perihelion differences shrink toward zero) at the
recent epoch they separated. Coincidental look-alikes do not converge.

  python scripts/confirm_pair.py "2021 PH27" "2025 GN1"

Uses the terrestrial planets too (Atiras are dominated by them). Honest limits: no
Yarkovsky or orbit-uncertainty clones here, so this is a first-pass convergence
check, not a publication-grade pair age.
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np
import requests

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from ariadne.discovery.frontier.exotic_orbit_hunt import Orb, dsh_matrix  # noqa: E402
from ariadne.dynamics.secular import (
    YEAR_S,
    add_test_particles,
    build_system,
    elements_to_state,
)
from ariadne.dynamics.secular_fast import integrate_fast_elements  # noqa: E402

ALL_PLANETS = (
    "MERCURY BARYCENTER",
    "VENUS BARYCENTER",
    "EARTH BARYCENTER",
    "MARS BARYCENTER",
    "JUPITER BARYCENTER",
    "SATURN BARYCENTER",
    "URANUS BARYCENTER",
    "NEPTUNE BARYCENTER",
)


def sbdb_elements(des):
    d = requests.get(
        "https://ssd-api.jpl.nasa.gov/sbdb.api",
        params={"sstr": des, "full-prec": "true"},
        timeout=40,
    ).json()
    o = d["orbit"]
    el = {e["name"]: float(e["value"]) for e in o["elements"]}
    return {
        "epoch_jd": float(o["epoch"]),
        "a": el["a"],
        "e": el["e"],
        "i": el["i"],
        "om": el["om"],
        "w": el["w"],
        "ma": el["ma"],
    }


def M_to_nu(M_deg, e):
    M = math.radians(M_deg % 360)
    E = M
    for _ in range(80):
        E -= (E - e * math.sin(E) - M) / (1 - e * math.cos(E))
    return math.degrees(
        2 * math.atan2(math.sqrt(1 + e) * math.sin(E / 2), math.sqrt(1 - e) * math.cos(E / 2))
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("objects", nargs="*", default=["2021 PH27", "2025 GN1"])
    ap.add_argument("--span-kyr", type=float, default=3.0)
    ap.add_argument("--dt-days", type=float, default=2.0)
    args = ap.parse_args()
    objs = args.objects or ["2021 PH27", "2025 GN1"]
    els = [sbdb_elements(o) for o in objs]
    # propagate every object's mean anomaly (2-body) to a COMMON epoch so the
    # initial conditions are simultaneous -- required if SBDB returns them at
    # different epochs (else the pair geometry is wrong from the start).
    target = els[0]["epoch_jd"]
    for e in els:
        if abs(e["epoch_jd"] - target) > 0.5:
            n_deg_day = 360.0 / (e["a"] ** 1.5 * 365.25)
            e["ma"] = (e["ma"] + n_deg_day * (target - e["epoch_jd"])) % 360.0
            e["epoch_jd"] = target
    from astropy.time import Time

    epoch = Time(target, format="jd", scale="tdb").utc.isot
    print(f"=== Pair-convergence test: {objs[0]} <-> {objs[1]} ===")
    print(
        f"  epoch {epoch} | full planetary system (8 planets) | {args.span_kyr} kyr backward",
        flush=True,
    )

    sys_ = build_system(epoch, massive=ALL_PLANETS)
    tps = []
    for e in els:
        nu = M_to_nu(e["ma"], e["e"])
        tps.append(elements_to_state(e["a"], e["e"], e["i"], e["om"], e["w"], nu))
    sys_ = add_test_particles(sys_, tps)

    dt = -args.dt_days * 86400.0
    n_steps = int(args.span_kyr * 1000 * YEAR_S / abs(dt))
    out = integrate_fast_elements(sys_, dt, n_steps, n_snap=120)

    def D_at(snap):
        a = [
            Orb(
                "A",
                s["a_au"],
                s["e"],
                s["i_deg"],
                s["Omega_deg"],
                s["omega_deg"],
                s["a_au"] * (1 - s["e"]),
            )
            for s in snap
        ]
        return float(dsh_matrix(a)[0, 1]), a

    best = (1e9, 0.0, None)
    rows = []
    for t, snap in zip(out["times_yr"], out["elements"]):
        D, orbs = D_at(snap)
        dOm = abs(((orbs[1].Omega - orbs[0].Omega + 180) % 360) - 180)
        dvp = abs(
            (((orbs[1].Omega + orbs[1].omega) - (orbs[0].Omega + orbs[0].omega) + 180) % 360) - 180
        )
        rows.append((t, D, dOm, dvp))
        if best[0] > D:
            best = (D, t, (dOm, dvp))
    print(
        f"  D_SH now: {rows[0][1]:.4f}  (dOmega={rows[0][2]:.2f} deg, dvarpi={rows[0][3]:.2f} deg)"
    )
    for t, D, dOm, dvp in rows[:: max(1, len(rows) // 10)]:
        print(f"    t={t:+8.0f} yr:  D_SH={D:.4f}  dOmega={dOm:5.2f}  dvarpi={dvp:5.2f}")
    print(
        f"  MINIMUM D_SH = {best[0]:.4f} at t={best[1]:+.0f} yr "
        f"(dOmega={best[2][0]:.2f}, dvarpi={best[2][1]:.2f})"
    )
    verdict = (
        "orbits CONVERGE -> consistent with a genuine pair (recent common origin)"
        if best[0] < rows[0][1] * 0.5 and best[0] < 0.02
        else "no clean convergence -> likely coincidental look-alike (or older/perturbed)"
    )
    print(f"  VERDICT: {verdict}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
