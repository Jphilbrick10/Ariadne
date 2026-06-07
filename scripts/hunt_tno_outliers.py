"""Extreme-TNO orbital-orientation outlier hunt (the Planet Nine clustering
signal), scored with the coherence engine.

  python scripts/hunt_tno_outliers.py
  python scripts/hunt_tno_outliers.py --a-min 250 --json out.json

Reports the clustering significance honestly across multiple a-cutoffs (the signal
is cutoff-dependent -- that IS the live debate), the coherence-selected cleanest
subsample, and the most orientation-anomalous orbits. Real JPL data; no downloads
beyond a small element table.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from ariadne.discovery.frontier.tno_clustering import outlier_hunt  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a-min", type=float, default=150.0)
    ap.add_argument("--q-min", type=float, default=30.0)
    ap.add_argument("--json", type=str, default="")
    ap.add_argument("--no-cache", action="store_true")
    args = ap.parse_args()

    rep = outlier_hunt(a_min=args.a_min, q_min=args.q_min, use_cache=not args.no_cache)
    if "error" in rep:
        print("  error:", rep["error"])
        return 1
    print("=== Extreme-TNO orientation outlier hunt ===")
    print(
        f"  {rep['n_total']} TNOs -> {rep['n_extreme']} extreme "
        f"(a>={rep['a_min']}, q>={rep['q_min']})"
    )
    for axis in ("Omega", "varpi"):
        s = rep[axis]
        sig = "SIGNIFICANT" if s["rayleigh_p"] < 0.05 else "not significant"
        print(
            f"  {axis:6s}: mean={s['mean_deg']:.0f} deg  R={s['R']:.3f}  "
            f"Rayleigh p={s['rayleigh_p']:.4f}  ({sig}, n={s['n']})"
        )
    print("  clustering vs a-cutoff (Omega) -- honest cutoff dependence:")
    for row in rep["Omega_select"]:
        print(
            f"    a>={row['a_min']:>3}: n={row['n']:>2}  R={row['R']:.3f}  "
            f"p={row['rayleigh_p']:.4f}  mean_E={row['mean_energy']:.3f}"
        )
    print("  most orientation-anomalous (least coherent) orbits:")
    for r in rep["most_anomalous"][:8]:
        print(
            f"    {r['name'][:30]:30s}  a={r['a']:>6.0f}  q={r['q']:>4.0f}  "
            f"Omega={r['Omega']:>5.0f}  E={r['energy']:.2f}"
        )
    if args.json:
        Path(args.json).write_text(json.dumps(rep, indent=2))
        print(f"  wrote {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
