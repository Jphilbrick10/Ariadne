"""Calibrate the ZTF/LSST anomaly basins against a real labeled ALeRCE sample.

Fetches the highest-probability objects of each class from the public ALeRCE
broker, computes their light-curve features, and fits robust (median + MAD) basins
saved to data/ztf_class_basins.json. The physics priors set which axes matter; this
sets the SCALES to the real survey so known objects score as known while genuine
novelties still cohere with nothing.

  python scripts/calibrate_ztf.py --per-class 25
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from ariadne.discovery.frontier import ztf_anomaly as Z  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-class", type=int, default=25)
    ap.add_argument("--classes", nargs="*", default=["RRL", "CEP", "E", "LPV", "QSO", "SNIa"])
    args = ap.parse_args()
    print(f"=== Calibrating ZTF basins from ALeRCE ({args.classes}) ===", flush=True)
    basins = Z.calibrate_from_alerce(tuple(args.classes), per_class=args.per_class)
    if not basins:
        print("  no basins fit (network? ALeRCE down?) -- physics priors remain", flush=True)
        return 1
    for k, b in basins.items():
        mu = " ".join(f"{a}={b['mu'][a]:.2f}" for a in b["mu"])
        print(f"  {k:6s} (n={b.get('n')}): {mu}", flush=True)
    print(f"  saved -> {Z.CALIBRATED_BASINS_PATH}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
