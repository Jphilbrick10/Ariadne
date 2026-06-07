"""GPU-accelerated shift-and-stack synthetic tracking (torch/CUDA).

The CPU `synthetic_tracking.fast_synthetic_tracking` is correct but bilinear-shifts
N full frames per rate/PA hypothesis on the CPU -- minutes per CCD. This does the
identical shift+coadd on the GPU (validated to match the CPU coadd to bilinear
rounding, ~1e-3), keeping the shifted stack on the device and transferring only the
small coadd per hypothesis (and the stack only when a peak exists, which is rare).
Falls back to the CPU implementation when no CUDA device is present, so callers can
always use it.

Speedup measured ~10x on a 2046x2046 x12 stack (RTX 5080); the win grows with frame
size and stack depth.
"""

from __future__ import annotations

import numpy as np

from .synthetic_tracking import (
    SyntheticCandidate,
    _per_image_signal_consensus,
    fast_synthetic_tracking,
    find_peaks_in_coadd,
    predicted_shift,
)

try:
    import torch

    HAVE_GPU = torch.cuda.is_available()
except Exception:  # pragma: no cover
    torch = None
    HAVE_GPU = False


def _gpu_coadd(img_t, base_gx, base_gy, shifts, nx, ny):
    """Shift each image by (dx,dy) on GPU and coadd (nanmean). Returns
    (coadd_np, coverage_np, stack_t) with the shifted stack left on the GPU."""
    grids, masks = [], []
    for dx, dy in shifts:
        sx = base_gx + dx
        sy = base_gy + dy
        masks.append((sx >= 0) & (sx <= nx - 1) & (sy >= 0) & (sy <= ny - 1))
        grids.append(torch.stack([2 * sx / (nx - 1) - 1, 2 * sy / (ny - 1) - 1], dim=-1))
    grid = torch.stack(grids)
    m = torch.stack(masks)
    out = torch.nn.functional.grid_sample(
        img_t.unsqueeze(1), grid, mode="bilinear", align_corners=True, padding_mode="zeros"
    ).squeeze(1)
    nan = torch.tensor(float("nan"), device=img_t.device)
    cov = m.sum(0)
    coadd = torch.where(
        cov > 0, torch.where(m, out, torch.zeros_like(out)).sum(0) / cov.clamp(min=1), nan
    )
    stack = torch.where(m, out, nan)
    return coadd.detach().cpu().numpy(), cov.detach().cpu().numpy().astype(np.int32), stack


def gpu_synthetic_tracking(
    images,
    wcs_list,
    image_mjds,
    *,
    t_ref_mjd=None,
    rate_min_arcsec_hr=0.5,
    rate_max_arcsec_hr=30.0,
    n_rates=24,
    n_pa=12,
    snr_threshold=5.0,
    pixscale_arcsec=1.0,
    n_top_per_hypothesis=5,
):
    """GPU shift-and-stack tracking. Same signature + output as
    fast_synthetic_tracking; uses CUDA when available, else falls back to CPU."""
    if not images:
        return []
    if not HAVE_GPU:
        return fast_synthetic_tracking(
            images,
            wcs_list,
            image_mjds,
            t_ref_mjd=t_ref_mjd,
            rate_min_arcsec_hr=rate_min_arcsec_hr,
            rate_max_arcsec_hr=rate_max_arcsec_hr,
            n_rates=n_rates,
            n_pa=n_pa,
            snr_threshold=snr_threshold,
            pixscale_arcsec=pixscale_arcsec,
            n_top_per_hypothesis=n_top_per_hypothesis,
        )

    if t_ref_mjd is None:
        t_ref_mjd = float(np.median(image_mjds))
    ny, nx = images[0].shape
    dev = "cuda"
    img_t = torch.as_tensor(
        np.stack([np.asarray(im, float) for im in images]), dtype=torch.float32, device=dev
    )
    ys = torch.arange(ny, device=dev, dtype=torch.float32)
    xs = torch.arange(nx, device=dev, dtype=torch.float32)
    base_gy, base_gx = torch.meshgrid(ys, xs, indexing="ij")

    rates = np.geomspace(rate_min_arcsec_hr, rate_max_arcsec_hr, n_rates)  # match CPU grid
    pas = np.linspace(0.0, 360.0, n_pa, endpoint=False)
    min_consensus = max(3, int(0.6 * len(images)))
    cands = []
    seen = []
    for rate in rates:
        for pa in pas:
            shifts = [
                predicted_shift(m, t_ref_mjd, float(rate), float(pa), pixscale_arcsec)
                for m in image_mjds
            ]
            coadd, cov, stack_t = _gpu_coadd(img_t, base_gx, base_gy, shifts, nx, ny)
            peaks = find_peaks_in_coadd(
                coadd,
                snr_threshold=snr_threshold,
                aperture_radius=3,
                min_separation_pix=8,
                coverage=cov,
                min_coverage=max(2, len(images) // 2),
            )
            if not peaks:
                continue
            peaks = peaks[:n_top_per_hypothesis]
            stack = stack_t.detach().cpu().numpy()  # transfer only when peaks (rare)
            for p in peaks:
                consensus = _per_image_signal_consensus(
                    stack, int(p["x"]), int(p["y"]), aperture_radius=3, signal_z_threshold=1.0
                )
                if consensus < min_consensus:
                    continue
                dup = any(
                    abs(p["x"] - sx) < 8 and abs(p["y"] - sy) < 8 and abs(rate - sr) < 1e-6
                    for sx, sy, sr in seen
                )
                if dup:
                    continue
                seen.append((p["x"], p["y"], rate))
                try:
                    ra, dec = wcs_list[0].pixel_to_world_values(p["x"], p["y"])
                    ra = float(ra) % 360.0
                    dec = float(dec)
                except Exception:
                    ra = dec = 0.0
                cands.append(
                    SyntheticCandidate(
                        ra_deg=ra,
                        dec_deg=dec,
                        rate_arcsec_hr=float(rate),
                        pa_deg=float(pa),
                        stacked_snr=float(p.get("snr", 0.0)),
                        n_images=len(images),
                        consensus_count=int(consensus),
                    )
                )
    return cands
