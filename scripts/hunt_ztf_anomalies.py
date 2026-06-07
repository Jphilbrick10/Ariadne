"""Run the physics-coherence novelty detector on a broad live batch of ZTF objects
and surface the most anomalous -- the genuinely-novel / mis-fitting light curves.

Pulls a diverse sample across ALeRCE classes (so the batch isn't pre-sorted), scores
each light curve's incoherence with the best-matching known class, and ranks by
novelty. High score = coheres with nothing known = worth a human's attention. This
is the dress rehearsal for triaging the LSST alert stream.

  python scripts/hunt_ztf_anomalies.py --per-class 25
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from ariadne.discovery.frontier import ztf_anomaly as Z  # noqa: E402

LEDGER = ROOT / "data" / "ztf_anomaly_hunt.jsonl"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-class", type=int, default=25)
    ap.add_argument(
        "--classes",
        nargs="*",
        default=[
            "SNIa",
            "SNII",
            "QSO",
            "RRL",
            "CEP",
            "E",
            "LPV",
            "YSO",
            "Blazar",
            "CV/Nova",
            "AGN",
            "SLSN",
        ],
    )
    args = ap.parse_args()
    using = "calibrated" if Z.CALIBRATED_BASINS_PATH.exists() else "physics-prior"
    print(
        f"=== ZTF anomaly hunt: {len(args.classes)} classes x {args.per_class}, {using} basins ===",
        flush=True,
    )
    rows = []
    for cls in args.classes:
        oids = Z.alerce_oids(cls, n=args.per_class)
        for oid in oids:
            t, m, e = Z.fetch_ztf_lightcurve(oid)
            if t is None:
                continue
            r = Z.score_lightcurve(t, m, e)
            if r["best_class"] is None:
                continue
            rows.append(
                {
                    "oid": oid,
                    "alerce_class": cls,
                    "score": round(r["score"], 2),
                    "best_class": r["best_class"],
                    "verdict": r["verdict"],
                    "period": r["features"].get("period"),
                    "amplitude": round(r["features"].get("amplitude", 0), 2),
                }
            )
        print(
            f"  scanned {cls}: {sum(1 for x in rows if x['alerce_class'] == cls)} scored",
            flush=True,
        )
    rows.sort(key=lambda r: -r["score"])
    anom = [r for r in rows if r["verdict"] == "anomalous"]
    print(
        f"\n  {len(rows)} objects scored; {len(anom)} flagged anomalous "
        f"(score > {Z.ANOMALY_TAU}). Top of the triage queue:"
    )
    for r in rows[:15]:
        flag = "  <== ANOMALOUS" if r["verdict"] == "anomalous" else ""
        print(
            f"    {r['oid']:16s} ALeRCE={r['alerce_class']:8s} score={r['score']:5.2f} "
            f"(best-fit {r['best_class']}){flag}"
        )
    LEDGER.parent.mkdir(parents=True, exist_ok=True)
    with open(LEDGER, "w") as fh:
        for r in rows:
            fh.write(json.dumps(r) + "\n")
    print(f"  ledger -> {LEDGER}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
