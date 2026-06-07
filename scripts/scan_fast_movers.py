"""Fast-mover scan -- the 'Oumuamua / interstellar / fast-NEO regime.

A fast object TRAILS within a single exposure; the point-source rate-linker (capped
~120"/hr) cannot reach it. detect_streaks (Hough transform) finds the trail and
streak_to_endpoints converts it to a timed sky position + rate -- with NO upper
rate limit. This scans a field's exposures for streaks, reports fast-mover
candidates, and persists them to the shared discovery ledger.

  python scripts/scan_fast_movers.py --data-dir data/decam_discovery_field

Honest scope: this is the DETECTION capability for the fast regime. Confirming an
interstellar origin still needs a multi-night/observation arc -> hyperbolic orbit
(e>1), which classify_mover then flags. A single trail flags a *candidate fast
mover*; follow-up establishes the orbit.
"""

from __future__ import annotations

import argparse
import glob
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
from run_auto_discovery import DATA, LEDGER, is_dup, load_ledger, write_watch  # noqa: E402

from ariadne.discovery.imaging.decam_instcal import load_decam_instcal  # noqa: E402
from ariadne.discovery.imaging.streak_tracklets import streak_to_endpoints  # noqa: E402
from ariadne.discovery.imaging.streaks import classify_streak, detect_streaks  # noqa: E402


def scan_field(data_dir, n_ccd=60, min_length_px=20.0, sigma=4.0, max_per_ccd=10):
    cands = []
    files = sorted(glob.glob(str(Path(data_dir) / "*.fits.fz")))
    for f in files:
        try:
            inst = load_decam_instcal(f, read_dqm=False)
        except Exception as ex:
            print(f"  skip {Path(f).name}: {type(ex).__name__}", flush=True)
            continue
        exps = inst.exptime_s or 30.0
        nstreak = 0
        for c in inst.ccds[:n_ccd]:
            if c.wcs is None:
                continue
            try:
                streaks = detect_streaks(
                    np.asarray(c.science, float),
                    sigma_threshold=sigma,
                    min_length_px=min_length_px,
                    subtract_stars=True,
                    max_streaks=max_per_ccd,
                )
            except Exception:
                continue
            diag = float(np.hypot(*np.asarray(c.science).shape))
            for s in streaks:
                # reject satellites / cosmic rays / artifacts; keep only trails that
                # are PSF-thin at an asteroid-plausible rate (the real fast movers)
                cl = classify_streak(
                    s, exposure_seconds=exps, pixel_scale_arcsec=0.263, frame_diagonal_px=diag
                )
                if not cl.get("is_asteroid_candidate"):
                    continue
                try:
                    ep = streak_to_endpoints(s, c.wcs, inst.mjd, exps, zeropoint_mag=c.magzero)
                except Exception:
                    continue
                cands.append(
                    {
                        "ra_deg": float(ep.ra_start),
                        "dec_deg": float(ep.dec_start),
                        "rate_arcsec_hr": float(ep.rate_arcsec_hr),
                        "pa_deg": float(ep.pa_deg),
                        "mag": round(float(ep.mag), 2) if ep.mag == ep.mag else None,
                        "length_px": round(float(s.length_px), 1),
                        "streak_label": cl.get("label"),
                        "exposure": Path(f).stem,
                    }
                )
                nstreak += 1
        print(f"  {Path(f).name}: {nstreak} streak candidate(s)", flush=True)
    return cands


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="data/decam_discovery_field")
    ap.add_argument("--n-ccd", type=int, default=60)
    ap.add_argument(
        "--min-length-px",
        type=float,
        default=20.0,
        help="reject trails shorter than this (px); lower = more sensitive + noisier",
    )
    ap.add_argument("--sigma", type=float, default=4.0)
    args = ap.parse_args()
    field_id = Path(args.data_dir).name
    print(f"=== FAST-MOVER (streak) scan: {field_id} ===", flush=True)
    cands = scan_field(args.data_dir, args.n_ccd, args.min_length_px, args.sigma)
    cands.sort(key=lambda c: -c["length_px"])
    print(f"\n  {len(cands)} fast-mover (streak) candidates total")
    for c in cands[:15]:
        print(
            f'    ({c["ra_deg"]:.4f},{c["dec_deg"]:+.4f}) rate={c["rate_arcsec_hr"]:.0f}"/hr '
            f"PA={c['pa_deg']:.0f} len={c['length_px']}px mag~{c['mag']}"
        )

    # persist new fast-mover candidates to the shared ledger (dedup by position+rate)
    ledger = load_ledger()
    now = datetime.now(timezone.utc).isoformat()
    new = 0
    for c in cands:
        e = dict(c)
        e["nights"] = []
        e["source"] = "streak-fast-mover"
        e["field_id"] = field_id
        e["first_seen_utc"] = now
        e["above_floor"] = False
        if not is_dup(e, ledger):
            with open(LEDGER, "a") as fh:
                fh.write(json.dumps(e) + "\n")
            ledger.append(e)
            new += 1
    print(f"\n  {new} new fast-mover candidates added to the ledger")
    state = (
        json.loads((DATA / "auto_discovery_state.json").read_text())
        if (DATA / "auto_discovery_state.json").exists()
        else {"runs": []}
    )
    state["runs"].append(
        {
            "utc": now,
            "field": f"{field_id} (fast-mover scan)",
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
