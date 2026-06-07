"""Scan every tractable solar-system population for anomalies and write a ranked
report. Pulls each dynamical class from JPL (cached), runs the full anomaly suite
(clustering+null, plane warp, D-criterion family+null, exotic-origin scoring), and
ranks populations by how anomalous they are -- honestly, with chance-excess nulls.

  python scripts/scan_everything.py
  python scripts/scan_everything.py --no-cache      # force fresh JPL pulls

Writes data/anomaly_scan_report.md and data/anomaly_candidates.jsonl.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from ariadne.discovery.frontier.frontier_scanner import POPULATIONS, scan_all  # noqa: E402

REPORT = ROOT / "data" / "anomaly_scan_report.md"
LEDGER = ROOT / "data" / "anomaly_candidates.jsonl"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-cache", action="store_true")
    args = ap.parse_args()
    t0 = time.time()
    print(f"=== Scanning {len(POPULATIONS)} solar-system populations ===", flush=True)
    res = scan_all(use_cache=not args.no_cache, log=lambda m: print("  " + m, flush=True))

    # rank populations by number of (significant) flags, then exotic richness
    ranked = sorted(
        res.values(), key=lambda r: (-len(r.get("flags", [])), -r.get("exotic_frac", 0))
    )
    lines = [
        "# Solar-system anomaly scan",
        "",
        f"{len(res)} populations scanned in {time.time() - t0:.0f}s. "
        "Every flag carries a chance-excess null; unflagged = consistent with chance.",
        "",
    ]
    cand = []
    for r in ranked:
        if "note" in r and r["n"] < 8:
            continue
        flags = r.get("flags", [])
        head = f"## {r['name']}  (n={r['n']}, a_med={r.get('a_median', '?')} AU)"
        lines.append(head)
        lines.append(
            "**FLAGS: " + ("; ".join(flags) if flags else "none (consistent with chance)") + "**"
        )
        nd, ap = r.get("node", {}), r.get("apsidal", {})
        lines.append(
            f"- node clustering R={nd.get('R')} (p={nd.get('p_value')}); "
            f"apsidal R={ap.get('R')} (p={ap.get('p_value')})"
        )
        if "plane_warp_deg" in r:
            lines.append(f"- plane warp from invariable plane: {r['plane_warp_deg']} deg")
        fam = r.get("family", {})
        if "excess_sigma" in fam:
            lines.append(
                f"- D-criterion close pairs: {fam['observed_pairs']} obs, "
                f"excess {fam['excess_sigma']} sigma (p={fam['p_value']})"
            )
            for p in r.get("top_pairs", [])[:3]:
                lines.append(f"    - candidate pair D={p['D']}: {p['a']} <-> {p['b']}")
        if r.get("n_unbound"):
            lines.append(f"- UNBOUND/interstellar candidates (e>=1): {r['n_unbound']}")
        if r.get("top_exotic"):
            top = r["top_exotic"][0]
            lines.append(
                f"- exotic-origin objects: {r['n_exotic']} "
                f"({r['exotic_frac'] * 100:.0f}%); top: {top['name']} {top['flags']}"
            )
        lines.append("")
        # ledger: flagged populations + their top candidates
        if flags:
            for p in r.get("top_pairs", []):
                cand.append({"type": "orbit_pair", "population": r["name"], **p})
            for ex in r.get("top_exotic", [])[:5]:
                cand.append({"type": "exotic_object", "population": r["name"], **ex})

    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(lines), encoding="utf-8")
    with open(LEDGER, "w") as fh:
        for c in cand:
            fh.write(json.dumps(c) + "\n")

    flagged = [r["name"] for r in ranked if r.get("flags")]
    print(f"\n  {len(flagged)} populations flagged: {flagged}", flush=True)
    print(f"  report -> {REPORT}", flush=True)
    print(f"  {len(cand)} candidates -> {LEDGER}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
