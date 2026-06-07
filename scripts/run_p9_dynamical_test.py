"""Secular dynamical test: do the REAL extreme TNOs keep their clustering under the
giant planets alone, or does it take a Planet Nine to preserve it? Integrates the
real orbital elements ~1.5 Gyr for the null (giants only) and a few P9 candidates,
recording the clustering resultant R(varpi), R(Omega) vs time. Writes JSON.

  python scripts/run_p9_dynamical_test.py --a-min 250 --span-gyr 1.5
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from ariadne.discovery.frontier.planet_nine import (
    P9Orbit,
    derive_orientation,
    secular_dispersal_test,
)
from ariadne.discovery.frontier.tno_clustering import extreme_population, fetch_tnos  # noqa: E402

OUT = ROOT / "data" / "p9_dynamical_test.json"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a-min", type=float, default=250.0)
    ap.add_argument("--span-gyr", type=float, default=1.5)
    ap.add_argument("--dt-yr", type=float, default=1.5e5)
    args = ap.parse_args()
    ext = extreme_population(fetch_tnos(), a_min=args.a_min, q_min=30)
    o = derive_orientation(ext)
    om = (o["varpi_p9_deg"] - o["plane_node_deg"]) % 360
    cases = {
        "null_giants_only": None,
        "P9_from_data": P9Orbit(460, 0.25, 16, o["plane_node_deg"], om, 6),
        "P9_literature": P9Orbit(460, 0.25, 16, 100.0, 150.0, 6),
        "P9_massive_far": P9Orbit(700, 0.4, 20, o["plane_node_deg"], om, 10),
    }
    print(
        f"=== P9 dynamical test: {len(ext)} real TNOs (a>={args.a_min}), {args.span_gyr} Gyr ===",
        flush=True,
    )
    results = {"n_tno": len(ext), "a_min": args.a_min, "span_gyr": args.span_gyr, "cases": {}}
    for name, p9 in cases.items():
        t0 = time.time()
        r = secular_dispersal_test(ext, p9, span_gyr=args.span_gyr, dt_yr=args.dt_yr)
        results["cases"][name] = r
        dr_v = r["R_varpi"][-1] - r["R_varpi"][0]
        dr_o = r["R_Omega"][-1] - r["R_Omega"][0]
        print(
            f"  {name:20s} ({time.time() - t0:.0f}s): "
            f"R_varpi {r['R_varpi'][0]:.2f}->{r['R_varpi'][-1]:.2f} (d={dr_v:+.2f})  "
            f"R_Omega {r['R_Omega'][0]:.2f}->{r['R_Omega'][-1]:.2f} (d={dr_o:+.2f})  "
            f"da/dt={r['da_au_per_yr_max']:.1e}",
            flush=True,
        )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(results, indent=2))
    print(f"  wrote {OUT}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
