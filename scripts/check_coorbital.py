"""Confirm (or refute) a 1:1 co-orbital resonance candidate by integration.

A similar semimajor axis is necessary but NOT sufficient for co-orbital status. The
robust, convention-free signature (valid for retrograde resonances too, cf. Morais &
Namouni 2017): the object's semimajor axis LIBRATES around the planet's -- the
resonance protects it from drifting away. A non-resonant look-alike's a wanders out
of the resonance band (or the object suffers a close encounter and is ejected).

  python scripts/check_coorbital.py "2006 BZ8" Saturn

Integrates the candidate as a test particle under the giant planets and reports
whether a stays bounded near the planet (resonant) or escapes (just a look-alike).
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
from ariadne.dynamics.secular import (
    YEAR_S,
    add_test_particles,
    build_system,
    elements_to_state,
)
from ariadne.dynamics.secular_fast import integrate_fast_elements  # noqa: E402

PLANET_A = {"Jupiter": 5.2026, "Saturn": 9.5826, "Uranus": 19.201, "Neptune": 30.110}
# co-orbital half-width ~ a * (m_planet/M_sun)^(1/3) (order of magnitude)
PLANET_MU = {"Jupiter": 9.54e-4, "Saturn": 2.86e-4, "Uranus": 4.37e-5, "Neptune": 5.15e-5}


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
    ap.add_argument("object")
    ap.add_argument("planet")
    ap.add_argument("--span-kyr", type=float, default=100.0)
    ap.add_argument("--dt-days", type=float, default=20.0)
    args = ap.parse_args()
    e = sbdb_elements(args.object)
    a_p = PLANET_A[args.planet]
    width = a_p * PLANET_MU[args.planet] ** (1 / 3) * 2.0  # ~ resonance half-width
    from astropy.time import Time

    epoch = Time(e["epoch_jd"], format="jd", scale="tdb").utc.isot
    print(
        f"=== Co-orbital test: {args.object} vs {args.planet} (a_p={a_p}, "
        f"band +/-{width:.2f} AU) ==="
    )
    print(
        f"  object now: a={e['a']:.3f} (|a-a_p|={abs(e['a'] - a_p):.3f}), "
        f"e={e['e']:.2f}, i={e['i']:.1f} deg | {args.span_kyr} kyr backward",
        flush=True,
    )

    sys_ = build_system(epoch)  # GIANTS
    nu = M_to_nu(e["ma"], e["e"])
    sys_ = add_test_particles(
        sys_, [elements_to_state(e["a"], e["e"], e["i"], e["om"], e["w"], nu)]
    )
    dt = -args.dt_days * 86400.0
    n_steps = int(args.span_kyr * 1000 * YEAR_S / abs(dt))
    out = integrate_fast_elements(sys_, dt, n_steps, n_snap=200)

    a_hist = np.array([snap[0]["a_au"] for snap in out["elements"]])
    i_hist = np.array([snap[0]["i_deg"] for snap in out["elements"]])
    t = out["times_yr"]
    in_band = np.abs(a_hist - a_p) < width
    frac_in = float(in_band.mean())
    print(
        f"  a over {args.span_kyr} kyr: min={a_hist.min():.2f} max={a_hist.max():.2f} "
        f"(planet a={a_p})"
    )
    print(f"  fraction of time within the co-orbital band: {frac_in * 100:.0f}%")
    print(
        f"  inclination stays retrograde (i>90): {float((i_hist > 90).mean()) * 100:.0f}% of the time"
        if e["i"] > 90
        else f"  i range: {i_hist.min():.0f}-{i_hist.max():.0f} deg"
    )
    for k in range(0, len(t), max(1, len(t) // 8)):
        print(
            f"    t={t[k]:+8.0f} yr: a={a_hist[k]:.2f}  i={i_hist[k]:.0f}  "
            f"{'IN band' if in_band[k] else 'OUT'}"
        )
    if frac_in > 0.8:
        v = "CO-ORBITAL CONFIRMED: a librates within the resonance band (resonance-protected)"
    elif frac_in > 0.4:
        v = "INTERMITTENT/temporary co-orbital: a enters and leaves the band"
    else:
        v = "NOT co-orbital: a escapes the band -> just a transient semimajor-axis look-alike"
    print(f"  VERDICT: {v}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
