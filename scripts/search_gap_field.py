"""Search an under-documented sky region for moving objects with archival DECam.

Targets the genuinely under-searched zones (galactic plane near the ecliptic, far
south) where surveys avoid the Milky Way's stellar confusion. Finds a same-pointing
multi-epoch set in the public NOIRLab DECam archive, downloads it, and runs the
validated faint-mover arsenal (GPU shift-and-stack + multi-frame consensus +
coherence vetting) -- the regime that reaches below single-frame detection.

  python scripts/search_gap_field.py --ra 269.6 --dec -23.4 --tag sgr --max-epochs 5

Honest scope: the galactic plane is the HARDEST sky (wall-to-wall stars), so expect
a high raw false-positive rate that consensus + vetting must beat down; few/no clean
candidates is the likely + honest outcome. This is exactly why the zone is
under-documented. Heavy (multi-GB downloads) -> background job.
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
from ariadne.discovery.imaging.noirlab_sia2 import (
    download_decam_exposure,
    query_decam_exposures,
)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ra", type=float, required=True)
    ap.add_argument("--dec", type=float, required=True)
    ap.add_argument("--tag", type=str, default="gap")
    ap.add_argument("--band", type=str, default="r")
    ap.add_argument("--max-epochs", type=int, default=5)
    ap.add_argument("--n-ccd", type=int, default=8)
    ap.add_argument("--pointing-tol", type=float, default=0.25, help="deg, same-pointing")
    args = ap.parse_args()

    print(f"=== Gap-field search: RA={args.ra} Dec={args.dec} ({args.tag}) ===", flush=True)
    recs = query_decam_exposures(
        args.ra,
        args.dec,
        radius_deg=0.4,
        band=args.band,
        proc_type="instcal",
        max_results=60,
        timeout_s=90,
    )
    # keep only the dominant same-pointing group (so CCD-i is the same sky each epoch)
    grp = defaultdict(list)
    for r in recs:
        grp[(round(r.ra_center * 4) / 4, round(r.dec_center * 4) / 4)].append(r)
    if not grp:
        print("  no exposures found")
        return 1
    pointing, rs = max(grp.items(), key=lambda kv: len(set(round(x.obs_mjd) for x in kv[1])))
    rs = sorted(rs, key=lambda r: r.obs_mjd)[: args.max_epochs]
    nights = sorted(set(round(r.obs_mjd) for r in rs))
    print(
        f"  pointing {pointing}: {len(rs)} exposures over {len(nights)} nights "
        f"(MJD {nights}, baseline {nights[-1] - nights[0]} d)",
        flush=True,
    )
    # CADENCE QUALITY GATE: moving-object detection needs epochs SPREAD in time.
    # A clump of near-simultaneous frames + one distant frame is degenerate -- a
    # stationary star is present in every near-simultaneous frame and fools the
    # consensus filter, and the lone far frame yields only spurious rate-alignments.
    mjds = sorted(r.obs_mjd for r in rs)
    intra = [mjds[i + 1] - mjds[i] for i in range(len(mjds) - 1)]
    well_spaced = sum(1 for g in intra if g > 0.02)  # >~30 min between epochs
    if len(nights) < 3 or well_spaced < 2:
        print(
            f"  WARNING: degenerate cadence ({well_spaced} well-spaced gaps, "
            f"{len(nights)} nights). This field cannot cleanly separate movers from "
            "stationary stars; any 'candidates' are confusion. A real search needs "
            ">=3 nights spaced hours-to-days apart -- rare in these gaps (the data "
            "there was taken for other science). Reporting, but treat as untrustworthy.",
            flush=True,
        )

    data_dir = ROOT / "data" / f"gap_{args.tag}"
    paths = []
    for i, r in enumerate(rs):
        print(f"  downloading {i + 1}/{len(rs)}: {r.obs_id} (MJD {r.obs_mjd:.2f})...", flush=True)
        p = download_decam_exposure(r, data_dir, timeout_s=900)
        if p:
            paths.append(p)
            print(f"    -> {p.name} ({p.stat().st_size / 1e6:.0f} MB)", flush=True)
        else:
            print("    download failed", flush=True)
    if len(paths) < 3:
        print(f"  only {len(paths)} exposures downloaded; need >=3 for shift-stack")
        return 1

    from scan_faint_movers import scan_field

    print(
        f"  running GPU shift-and-stack on {len(paths)} epochs, {args.n_ccd} CCDs "
        "(galactic plane: confusion-limited)...",
        flush=True,
    )
    try:
        cands = scan_field(str(data_dir), n_ccd=args.n_ccd, snr=6.0, min_consensus=len(paths) - 1)
    except Exception as ex:
        import traceback

        traceback.print_exc()
        print(f"  scan errored ({type(ex).__name__}); reporting nothing", flush=True)
        return 1
    cands.sort(key=lambda c: -c["snr"])
    print(
        f"\n  {len(cands)} faint shift-stack candidates (consensus>={len(paths) - 1}):", flush=True
    )
    for c in cands[:15]:
        print(
            f'    ({c["ra_deg"]:.4f},{c["dec_deg"]:+.4f}) rate={c["rate_arcsec_hr"]:.1f}"/hr '
            f"SNR={c['snr']} consensus={c['consensus']}/{c['n_images']} CCD{c['ccd']}",
            flush=True,
        )
    print(
        "\n  (galactic-center confusion makes clean detections unlikely; consensus +"
        " vetting reject most. Honest first-pass of the least-searched sky.)",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
