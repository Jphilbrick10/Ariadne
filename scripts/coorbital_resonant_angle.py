"""Co-orbital confirmation by RESONANT-ANGLE LIBRATION -- the correct diagnostic that
replaces the retracted semimajor-axis-band proxy.

A 1:1 co-orbital is defined by libration (bounded oscillation, not circulation) of the
resonant angle phi = lambda_obj - lambda_planet (mean longitudes). The a-band proxy
failed because a co-orbital's a can swing while phi stays bounded. Here we integrate
with REBOUND (MERCURIUS) from accurate JPL-Horizons initial conditions, track phi, and
classify: if phi avoids part of the circle (bounded) it LIBRATES (co-orbital); if it
sweeps the full circle it CIRCULATES (not co-orbital).

Retrograde subtlety: for a retrograde orbit (i>90) the prograde phi circulates even when
the body IS resonant; Morais & Namouni show the librating argument uses the retrograde
mean longitude lambda* = M + Omega - omega. We compute BOTH and report whichever
librates -- the resonance is whichever argument is bounded.

  python scripts/coorbital_resonant_angle.py "514107" Jupiter --span-kyr 100

Validate on a known case (Jupiter Trojan = prograde libration ~+/-60; Ka'epaoka'awela =
retrograde libration) BEFORE trusting any new verdict.
"""

from __future__ import annotations

import argparse
import math
import sys

import numpy as np


def libration_diagnosis(phi_deg: np.ndarray, nbins: int = 36) -> dict:
    """Does phi librate (bounded, avoids part of the circle) or circulate (covers it)?"""
    occ = np.zeros(nbins, dtype=bool)
    idx = (np.asarray(phi_deg) % 360 / (360 / nbins)).astype(int) % nbins
    occ[idx] = True
    frac = occ.mean()
    # libration center/amplitude (circular)
    th = np.radians(phi_deg % 360)
    C, S = np.cos(th).mean(), np.sin(th).mean()
    center = math.degrees(math.atan2(S, C)) % 360
    R = math.hypot(C, S)
    return {
        "circle_fraction_visited": round(float(frac), 2),
        "librates": bool(frac < 0.85),
        "center_deg": round(center, 0),
        "concentration_R": round(R, 2),
    }


def main():
    import rebound

    ap = argparse.ArgumentParser()
    ap.add_argument("object")
    ap.add_argument("planet")
    ap.add_argument("--span-kyr", type=float, default=100.0)
    args = ap.parse_args()
    print(f"=== Resonant-angle co-orbital test: {args.object} vs {args.planet} ===", flush=True)
    sim = rebound.Simulation()
    sim.units = ("yr", "AU", "Msun")
    bodies = ["Sun", "Jupiter", "Saturn", "Uranus", "Neptune", args.object]
    for b in bodies:
        sim.add(b)
    pidx = bodies.index(args.planet)
    oidx = len(bodies) - 1
    nom = sim.particles[oidx].orbit(primary=sim.particles[0])
    retro = math.degrees(nom.inc) > 90
    print(
        f"  {args.object}: a={nom.a:.3f} e={nom.e:.2f} i={math.degrees(nom.inc):.0f}deg "
        f"({'RETROGRADE' if retro else 'prograde'})",
        flush=True,
    )
    sim.integrator = "mercurius"
    sim.dt = 0.5
    sim.move_to_com()

    nsnap = 600
    phi_pro, phi_retro, a_obj = [], [], []
    for t in -np.linspace(0, args.span_kyr * 1000, nsnap):
        sim.integrate(t)
        sun = sim.particles[0]
        op = sim.particles[pidx].orbit(primary=sun)
        oo = sim.particles[oidx].orbit(primary=sun)
        lam_p = math.degrees(op.l)  # planet mean longitude
        lam_pro = math.degrees(oo.l)  # prograde mean longitude
        # retrograde mean longitude: lambda* = M + Omega - omega
        lam_ret = math.degrees(oo.M + oo.Omega - oo.omega)
        phi_pro.append((lam_pro - lam_p) % 360)
        phi_retro.append((lam_ret - lam_p) % 360)
        a_obj.append(oo.a)
    dp = libration_diagnosis(np.array(phi_pro))
    dr = libration_diagnosis(np.array(phi_retro))
    print(
        f"  prograde   phi=lam-lam_p : visits {dp['circle_fraction_visited'] * 100:.0f}% of circle "
        f"-> {'LIBRATES' if dp['librates'] else 'circulates'} (center {dp['center_deg']:.0f})",
        flush=True,
    )
    print(
        f"  retrograde phi=lam*-lam_p: visits {dr['circle_fraction_visited'] * 100:.0f}% of circle "
        f"-> {'LIBRATES' if dr['librates'] else 'circulates'} (center {dr['center_deg']:.0f})",
        flush=True,
    )
    a = np.array(a_obj)
    print(
        f"  (a ranged {a.min():.2f}-{a.max():.2f} AU over {args.span_kyr} kyr -- note a alone is "
        "not the diagnostic)",
        flush=True,
    )
    co = dr["librates"] if retro else dp["librates"]
    print(
        f"  VERDICT: {'CO-ORBITAL -- resonant angle LIBRATES' if co else 'NOT co-orbital -- resonant angle circulates'}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
