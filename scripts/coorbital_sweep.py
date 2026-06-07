"""Census of resonance-confirmed exotic (high-inclination / retrograde) co-orbitals
of the giant planets.

A similar semimajor axis is not enough; this integrates each exotic candidate and
keeps only those whose a LIBRATES within the planet's co-orbital band (resonance-
protected) -- the genuine captured-object class (e.g. Ka'epaoka'awela / 2015 BZ509).
Validates against that known object, then reports every confirmed co-orbital.

  python scripts/coorbital_sweep.py --max 15 --span-kyr 50
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
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

PL = {
    "Jupiter": (5.2026, 9.54e-4),
    "Saturn": (9.5826, 2.86e-4),
    "Uranus": (19.201, 4.37e-5),
    "Neptune": (30.110, 5.15e-5),
}
OUT = ROOT / "data" / "coorbital_census.json"


def M_to_nu(M, e):
    M = math.radians(M % 360)
    E = M
    for _ in range(80):
        E -= (E - e * math.sin(E) - M) / (1 - e * math.cos(E))
    return math.degrees(
        2 * math.atan2(math.sqrt(1 + e) * math.sin(E / 2), math.sqrt(1 - e) * math.cos(E / 2))
    )


def libration_status(nm, planet, el, span_kyr, dt_days):
    """el = dict(a,e,i,om,w,ma,epoch_jd) straight from the catalogue query."""
    from astropy.time import Time

    a_p, mu = PL[planet]
    band = a_p * mu ** (1 / 3) * 2.0
    epoch = Time(el["epoch_jd"], format="jd", scale="tdb").utc.isot
    sys_ = build_system(epoch)
    nu = M_to_nu(el["ma"], el["e"])
    sys_ = add_test_particles(
        sys_, [elements_to_state(el["a"], el["e"], el["i"], el["om"], el["w"], nu)]
    )
    dt = -dt_days * 86400.0
    out = integrate_fast_elements(sys_, dt, int(span_kyr * 1000 * YEAR_S / abs(dt)), n_snap=160)
    a_hist = np.array([s[0]["a_au"] for s in out["elements"]])
    i_hist = np.array([s[0]["i_deg"] for s in out["elements"]])
    frac = float((np.abs(a_hist - a_p) < band).mean())
    verdict = "CONFIRMED" if frac > 0.8 else "TEMPORARY" if frac > 0.4 else "look-alike"
    return {
        "object": nm,
        "planet": planet,
        "a": el["a"],
        "e": el["e"],
        "i": el["i"],
        "frac_in_band": round(frac, 2),
        "retrograde": el["i"] > 90,
        "a_min": round(float(a_hist.min()), 2),
        "a_max": round(float(a_hist.max()), 2),
        "verdict": verdict,
    }


def candidates():
    """Exotic (i>60) objects inside a giant-planet co-orbital band, with full
    elements pulled directly from the catalogue query (no per-object re-fetch)."""
    d = requests.get(
        "https://ssd-api.jpl.nasa.gov/sbdb_query.api"
        "?fields=full_name,a,e,i,om,w,ma,epoch"
        '&sb-cdata={"AND":["i|GT|60","a|GT|4"]}',
        timeout=120,
    ).json()
    out = []
    for r in d.get("data", []):
        try:
            el = dict(
                a=float(r[1]),
                e=float(r[2]),
                i=float(r[3]),
                om=float(r[4]),
                w=float(r[5]),
                ma=float(r[6]),
                epoch_jd=float(r[7]),
            )
        except (TypeError, ValueError):
            continue
        for pn, (pa, mu) in PL.items():
            if abs(el["a"] - pa) < pa * mu ** (1 / 3) * 2.0:
                out.append((abs(el["a"] - pa), r[0].strip(), pn, el))
    out.sort(key=lambda x: x[0])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max", type=int, default=15)
    ap.add_argument("--span-kyr", type=float, default=50.0)
    ap.add_argument("--dt-days", type=float, default=15.0)
    args = ap.parse_args()
    cands = candidates()
    print(
        f"=== Co-orbital confirmation sweep: {len(cands)} exotic candidates, "
        f"running closest {args.max} ===",
        flush=True,
    )
    results = []
    for da, nm, pn, el in cands[: args.max]:
        t0 = time.time()
        try:
            r = libration_status(nm, pn, el, args.span_kyr, args.dt_days)
        except Exception as ex:
            print(f"  {nm[:24]:24s} {pn:8s}: FAIL {type(ex).__name__}: {ex}", flush=True)
            continue
        results.append(r)
        tag = "RETRO" if r["retrograde"] else "hi-i"
        print(
            f"  {nm[:24]:24s} {pn:8s} i={el['i']:3.0f} {tag:5s}: a {r['a_min']}-{r['a_max']} "
            f"in-band {r['frac_in_band'] * 100:.0f}% -> {r['verdict']}  ({time.time() - t0:.0f}s)",
            flush=True,
        )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(results, indent=2))
    conf = [r for r in results if r["verdict"] == "CONFIRMED"]
    print(
        f"\n  CONFIRMED co-orbitals: {len(conf)} ({sum(r['retrograde'] for r in conf)} retrograde)",
        flush=True,
    )
    for r in conf:
        print(
            f"    {r['object']}  {r['planet']}  i={r['i']:.0f}deg  "
            f"{'RETROGRADE' if r['retrograde'] else 'high-i'}",
            flush=True,
        )
    print(f"  census -> {OUT}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
