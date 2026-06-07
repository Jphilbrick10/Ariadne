"""Physics-coherence novelty triage of ZTF light curves (the LSST firehose dry run).

  python scripts/triage_ztf_anomalies.py ZTF18abbuive ZTF20acobvxk
  python scripts/triage_ztf_anomalies.py --alerce-class SNIa --n 30

Scores each object's incoherence with the best-matching known class. High score =
coheres with nothing known = worth a human's attention. Uses the calibrated basins
(run scripts/calibrate_ztf.py first) if present, else physics priors. This is a
ranking/triage score, not a classifier or a discovery claim.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from ariadne.discovery.frontier import ztf_anomaly as Z  # noqa: E402

LEDGER = ROOT / "data" / "frontier_ztf_anomalies.jsonl"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("oids", nargs="*", help="ZTF object ids")
    ap.add_argument(
        "--alerce-class", type=str, default="", help="instead pull N objects of this ALeRCE class"
    )
    ap.add_argument("--n", type=int, default=20)
    args = ap.parse_args()
    oids = list(args.oids)
    if args.alerce_class:
        oids += Z.alerce_oids(args.alerce_class, n=args.n)
    if not oids:
        print("  give ZTF oids or --alerce-class")
        return 1

    using = "calibrated" if Z.CALIBRATED_BASINS_PATH.exists() else "physics-prior"
    print(f"=== ZTF novelty triage ({len(oids)} objects, {using} basins) ===")
    rows = []
    for oid in oids:
        t, m, e = Z.fetch_ztf_lightcurve(oid)
        if t is None:
            print(f"  {oid:16s}  no light curve")
            continue
        r = Z.score_lightcurve(t, m, e)
        fe = r["features"]
        flag = "  <== ANOMALOUS" if r["verdict"] == "anomalous" else ""
        print(
            f"  {oid:16s}  P={fe.get('period', float('nan')):8.3f} "
            f"amp={fe.get('amplitude', 0):.2f}  score={r['score']:5.2f}  {r['verdict']}{flag}"
        )
        rows.append(
            {
                "oid": oid,
                "score": r["score"],
                "verdict": r["verdict"],
                "best_class": r["best_class"],
                "period": fe.get("period"),
                "amplitude": fe.get("amplitude"),
            }
        )
    rows.sort(key=lambda r: -r["score"])
    anom = [r for r in rows if r["verdict"] == "anomalous"]
    print(f"\n  {len(anom)}/{len(rows)} flagged anomalous (top of the triage queue)")
    if rows:
        LEDGER.parent.mkdir(parents=True, exist_ok=True)
        with open(LEDGER, "a") as fh:
            for r in rows:
                fh.write(json.dumps(r) + "\n")
        print(f"  appended to {LEDGER}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
