"""Transit search + coherence vetting of TESS light curves.

  python scripts/hunt_tess_candidates.py "Pi Mensae" "WASP-18" "TIC 307210830"
  python scripts/hunt_tess_candidates.py --targets-file my_tics.txt

For each target: pull the free TESS light curve (MAST), run Box Least Squares, and
score the candidate's coherence with a real planet. Surfaces planet-candidates and
flags likely eclipsing-binary / systematic false positives -- the vetting step that
is the actual bottleneck. Confirmation still needs follow-up we do not have; this
ranks what is worth that follow-up.
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from ariadne.discovery.frontier.tess_vetting import (  # noqa: E402
    fetch_tess_lightcurve,
    search_and_vet,
)

LEDGER = ROOT / "data" / "frontier_tess_ledger.jsonl"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("targets", nargs="*", help="TIC ids or resolvable names")
    ap.add_argument("--targets-file", type=str, default="")
    ap.add_argument("--min-period", type=float, default=0.5)
    ap.add_argument("--max-period", type=float, default=14.0)
    args = ap.parse_args()
    targets = list(args.targets)
    if args.targets_file:
        targets += [
            ln.strip()
            for ln in Path(args.targets_file).read_text().splitlines()
            if ln.strip() and not ln.startswith("#")
        ]
    if not targets:
        targets = ["Pi Mensae", "WASP-18"]  # demo defaults (real confirmed planets)

    print(f"=== TESS transit search + coherence vetting ({len(targets)} targets) ===")
    rows = []
    for tgt in targets:
        t, f = fetch_tess_lightcurve(tgt)
        if t is None:
            print(f"  {tgt:22s}  no light curve")
            continue
        c = search_and_vet(t, f, min_period=args.min_period, max_period=args.max_period)
        if c is None:
            print(f"  {tgt:22s}  no candidate")
            continue
        fe = c.features
        print(
            f"  {tgt:22s}  P={c.period:7.3f}d depth={fe.get('depth', 0) * 100:6.3f}% "
            f"SNR={fe.get('snr', 0):6.1f}  coh={c.coherence:.3f}  {c.verdict}"
        )
        rows.append(
            {
                "target": tgt,
                "period": c.period,
                "depth": fe.get("depth"),
                "snr": fe.get("snr"),
                "coherence": c.coherence,
                "verdict": c.verdict,
                "utc": datetime.now(timezone.utc).isoformat(),
            }
        )
    cand = [r for r in rows if r["verdict"] == "planet-candidate"]
    print(f"\n  {len(cand)}/{len(rows)} passed coherence vetting as planet-candidates")
    if rows:
        LEDGER.parent.mkdir(parents=True, exist_ok=True)
        with open(LEDGER, "a") as fh:
            for r in rows:
                fh.write(json.dumps(r) + "\n")
        print(f"  appended to {LEDGER}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
