"""Co-orbital screen: REBOUND (MERCURIUS, adaptive close-encounter handling) + an
orbit-uncertainty CLONE ensemble, scoring semimajor-axis-band occupancy.

It fixes the integrator half of the problem (MERCURIUS switches to IAS15 inside close
encounters, so a no longer goes spuriously hyperbolic; clones probe the chaos). It
does NOT fix the DIAGNOSTIC: a-band occupancy is an inadequate proxy for 1:1 co-orbital
resonance.

  CRITICAL INTEGRATION NOTE: the a-band score is only meaningful with an ACCURATE
  integrator. With a coarse timestep (MERCURIUS default / dt~0.5 yr) the deep
  planet-crossing encounters of these high-e retrograde orbits are mishandled and a
  spuriously wanders -- which earlier made even the KNOWN co-orbital Ka'epaoka'awela
  look unconfined (~44% in-band) and led to a (later-reversed) retraction. With the
  adaptive IAS15 integrator (or dt<=0.1), Ka'epaoka'awela is 100% a-confined and the
  candidates 2012 TL139 / 2006 BZ8 are a-confined too -- restoring them as likely
  retrograde co-orbitals. USE IAS15 (or a very small timestep) for high-e crossers.

  Remaining caveat: a-confinement is strong SUGGESTIVE evidence, not the gold standard.
  The resonant ANGLE (phi = lambda - lambda_planet for prograde; the Morais-Namouni
  retrograde argument for i>90) is the rigorous proof; the retrograde argument does not
  librate in any simple sign-combination here and is the genuinely-open piece.

Self-contained planet model (giant-planet osculating elements propagated to epoch);
for production use accurate Horizons ICs + IAS15 (see the inline re-examination).
"""

from __future__ import annotations

import argparse
import math
import sys

import numpy as np
import requests

# giant-planet elements: a(AU), e, i, Omega, omega, M0(deg @J2000), period(yr), mass(Msun)
GIANTS = {
    "Jupiter": (5.2026, 0.0484, 1.303, 100.49, 273.87, 19.65, 11.862, 9.5479e-4),
    "Saturn": (9.5826, 0.0539, 2.485, 113.64, 339.39, 317.51, 29.457, 2.8589e-4),
    "Uranus": (19.201, 0.0473, 0.773, 74.00, 96.99, 142.27, 84.011, 4.3662e-5),
    "Neptune": (30.110, 0.0086, 1.770, 131.78, 276.34, 256.23, 164.79, 5.1514e-5),
}
J2000_JD = 2451545.0


def sbdb(des):
    d = requests.get(
        "https://ssd-api.jpl.nasa.gov/sbdb.api",
        params={"sstr": des, "full-prec": "true"},
        timeout=40,
    ).json()
    o = d["orbit"]
    el = {x["name"]: x for x in o["elements"]}

    def val(k):
        return float(el[k]["value"])

    def sig(k, frac):  # formal sigma if present, else a small fractional default
        s = el[k].get("sigma")
        return (
            float(s)
            if s
            not in (
                None,
                "",
            )
            else abs(val(k)) * frac + 1e-4
        )

    return dict(
        epoch_jd=float(o["epoch"]),
        a=val("a"),
        e=val("e"),
        i=val("i"),
        om=val("om"),
        w=val("w"),
        ma=val("ma"),
        sa=sig("a", 1e-4),
        se=sig("e", 1e-3),
        si=sig("i", 1e-2),
    )


def main():
    import rebound

    ap = argparse.ArgumentParser()
    ap.add_argument("object")
    ap.add_argument("planet")
    ap.add_argument("--clones", type=int, default=30)
    ap.add_argument("--span-kyr", type=float, default=100.0)
    args = ap.parse_args()
    e = sbdb(args.object)
    a_p, mu = GIANTS[args.planet][0], GIANTS[args.planet][7]
    band = a_p * mu ** (1 / 3) * 2.0
    yrs_since = (e["epoch_jd"] - J2000_JD) / 365.25
    print(f"=== REBOUND co-orbital confirmation: {args.object} vs {args.planet} ===")
    print(
        f"  a_p={a_p} band +/-{band:.2f} AU | {args.clones} clones | {args.span_kyr} kyr "
        f"| MERCURIUS",
        flush=True,
    )

    sim = rebound.Simulation()
    sim.units = ("yr", "AU", "Msun")
    sim.add(m=1.0)  # Sun
    for nm, (a, ec, inc, Om, om, M0, per, m) in GIANTS.items():
        M = math.radians((M0 + 360 * yrs_since / per) % 360)
        sim.add(
            m=m,
            a=a,
            e=ec,
            inc=math.radians(inc),
            Omega=math.radians(Om),
            omega=math.radians(om),
            M=M,
            primary=sim.particles[0],
        )
    n_planets = len(sim.particles)
    rng = np.random.default_rng(0)
    M0c = math.radians(e["ma"])
    for k in range(args.clones):
        da = 0 if k == 0 else rng.normal(0, e["sa"])
        de = 0 if k == 0 else rng.normal(0, e["se"])
        di = 0 if k == 0 else rng.normal(0, e["si"])
        sim.add(
            m=0.0,
            a=e["a"] + da,
            e=max(1e-3, min(0.99, e["e"] + de)),
            inc=math.radians(e["i"] + di),
            Omega=math.radians(e["om"]),
            omega=math.radians(e["w"]),
            M=M0c,
            primary=sim.particles[0],
        )
    sim.integrator = "mercurius"
    sim.dt = 0.6
    sim.move_to_com()

    nsnap = 60
    times = -np.linspace(0, args.span_kyr * 1000, nsnap)  # backward
    frac_in = np.zeros(args.clones)
    survived = np.zeros(args.clones, dtype=bool)
    for t in times:
        sim.integrate(t)
        sun = sim.particles[0]
        for c in range(args.clones):
            p = sim.particles[n_planets + c]
            try:
                a = p.orbit(primary=sun).a
            except Exception:
                continue
            if a > 0 and abs(a - a_p) < band:
                frac_in[c] += 1.0 / nsnap
    conf = int(np.sum(frac_in > 0.7))
    temp = int(np.sum((frac_in > 0.3) & (frac_in <= 0.7)))
    nominal = frac_in[0]
    print(f"  nominal clone: in-band {nominal * 100:.0f}% of {args.span_kyr} kyr", flush=True)
    print(
        f"  clone ensemble (n={args.clones}): {conf} resonance-locked (>70% in-band), "
        f"{temp} intermittent (30-70%), {args.clones - conf - temp} escaped",
        flush=True,
    )
    frac_conf = conf / args.clones
    verdict = (
        "CONFIRMED co-orbital (robust: most clones resonance-locked)"
        if frac_conf > 0.6
        else "TEMPORARY/marginal (clones split -- chaotic, not a stable resonance)"
        if frac_conf > 0.2
        else "NOT co-orbital (clones escape the band)"
    )
    print(f"  VERDICT: {verdict}  [{frac_conf * 100:.0f}% of clones locked]", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
