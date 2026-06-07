"""Faint-object scan -- shift-and-stack synthetic tracking (the KBMOD-style deep /
TNO / distant regime).

Co-adds repeated exposures of the SAME sky patch (same CCD across a field's
exposures) along a grid of rate/position-angle hypotheses. A mover too faint to
see in any single frame accumulates signal-to-noise in the stack at its true
rate. This reaches objects below single-frame detection -- the regime the
point-source extractor cannot touch.

  python scripts/scan_faint_movers.py --data-dir data/decam_discovery_field

Honest scope: shift-stack has real chance false positives (noise peaks aligning),
so we keep only candidates with multi-frame CONSENSUS (signal present in >=N
frames). It needs repeated exposures of the same pointing. Heavy (a rate/PA grid
per CCD stack) -> background job.
"""

from __future__ import annotations

import argparse
import glob
import json
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
from run_auto_discovery import DATA, LEDGER, is_dup, load_ledger, write_watch  # noqa: E402

from ariadne.discovery.imaging.decam_instcal import load_decam_instcal  # noqa: E402
from ariadne.discovery.imaging.gpu_shift_stack import gpu_synthetic_tracking  # noqa: E402


def scan_field(
    data_dir, n_ccd=60, rate_min=0.5, rate_max=30.0, snr=5.0, min_epochs=3, min_consensus=3
):
    files = sorted(glob.glob(str(Path(data_dir) / "*.fits.fz")))
    by_ccd = defaultdict(list)  # ccd index -> [(science, wcs, mjd), ...]
    for f in files:
        try:
            inst = load_decam_instcal(f, read_dqm=False)
        except Exception as ex:
            print(f"  skip {Path(f).name}: {type(ex).__name__}", flush=True)
            continue
        for i, c in enumerate(inst.ccds[:n_ccd]):
            if c.wcs is not None:
                by_ccd[i].append((np.asarray(c.science, float), c.wcs, inst.mjd))
    cands = []
    for i, stack in sorted(by_ccd.items()):
        if len(stack) < min_epochs:
            continue
        imgs = [s[0] for s in stack]
        wcss = [s[1] for s in stack]
        mjds = [s[2] for s in stack]
        try:
            cc = gpu_synthetic_tracking(
                imgs,
                wcss,
                mjds,
                rate_min_arcsec_hr=rate_min,
                rate_max_arcsec_hr=rate_max,
                snr_threshold=snr,
                pixscale_arcsec=0.263,
            )
        except Exception:
            continue
        kept = [x for x in cc if x.consensus_count >= min_consensus]
        for x in kept:
            cands.append(
                {
                    "ra_deg": float(x.ra_deg),
                    "dec_deg": float(x.dec_deg),
                    "rate_arcsec_hr": float(x.rate_arcsec_hr),
                    "pa_deg": float(x.pa_deg),
                    "snr": round(float(x.stacked_snr), 1),
                    "n_images": int(x.n_images),
                    "consensus": int(x.consensus_count),
                    "ccd": i,
                }
            )
        if kept:
            print(f"  CCD {i}: {len(stack)} epochs -> {len(kept)} faint candidate(s)", flush=True)
    return cands


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="data/decam_discovery_field")
    ap.add_argument("--n-ccd", type=int, default=60)
    ap.add_argument("--rate-min", type=float, default=0.5)
    ap.add_argument(
        "--rate-max", type=float, default=30.0, help='"/hr; distant objects move slowly'
    )
    ap.add_argument("--snr", type=float, default=5.0)
    ap.add_argument(
        "--min-consensus",
        type=int,
        default=3,
        help="require signal in >=N frames (rejects chance noise-peak stacks)",
    )
    args = ap.parse_args()
    field_id = Path(args.data_dir).name
    print(f"=== FAINT-OBJECT (shift-stack) scan: {field_id} ===", flush=True)
    cands = scan_field(
        args.data_dir,
        args.n_ccd,
        args.rate_min,
        args.rate_max,
        args.snr,
        min_consensus=args.min_consensus,
    )
    cands.sort(key=lambda c: -c["snr"])
    print(f"\n  {len(cands)} faint shift-stack candidates (consensus>={args.min_consensus})")
    for c in cands[:15]:
        print(
            f'    ({c["ra_deg"]:.4f},{c["dec_deg"]:+.4f}) rate={c["rate_arcsec_hr"]:.1f}"/hr '
            f"SNR={c['snr']} consensus={c['consensus']}/{c['n_images']} CCD{c['ccd']}"
        )

    ledger = load_ledger()
    now = datetime.now(timezone.utc).isoformat()
    new = 0
    for c in cands:
        e = dict(c)
        e["nights"] = []
        e["source"] = "shift-stack-faint"
        e["field_id"] = field_id
        e["first_seen_utc"] = now
        e["above_floor"] = False
        if not is_dup(e, ledger):
            with open(LEDGER, "a") as fh:
                fh.write(json.dumps(e) + "\n")
            ledger.append(e)
            new += 1
    print(f"\n  {new} new faint candidates added to the ledger")
    state = (
        json.loads((DATA / "auto_discovery_state.json").read_text())
        if (DATA / "auto_discovery_state.json").exists()
        else {"runs": []}
    )
    state["runs"].append(
        {
            "utc": now,
            "field": f"{field_id} (faint shift-stack)",
            "recovered": 0,
            "recoverable": 0,
            "candidates": len(cands),
            "floor": "-",
            "above_floor": False,
            "new_in_ledger": new,
        }
    )
    (DATA / "auto_discovery_state.json").write_text(json.dumps(state, indent=2))
    write_watch(state, load_ledger())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
