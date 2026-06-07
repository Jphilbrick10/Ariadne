"""GPU shift-and-stack must match the CPU implementation: identical coadd to
bilinear rounding, and it recovers the same planted moving object. Skips cleanly
when no CUDA device is present (the function falls back to CPU there anyway)."""

from __future__ import annotations

import math

import numpy as np
import pytest

torch = pytest.importorskip("torch")
if not torch.cuda.is_available():
    pytest.skip("no CUDA device", allow_module_level=True)


def _gaussian_image(npix, x, y, flux, sigma=1.5, bg=100, bg_sigma=8, seed=0):
    rng = np.random.default_rng(seed)
    im = rng.normal(bg, bg_sigma, (npix, npix))
    yy, xx = np.mgrid[0:npix, 0:npix]
    return im + flux * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2 * sigma**2))


def _scene():
    from astropy.wcs import WCS

    npix, ps = 200, 1.0
    mjds = [60450.0, 60450.083, 60453.0, 60453.083, 60456.0, 60456.083]
    tref = mjds[len(mjds) // 2]

    def mkw():
        w = WCS(naxis=2)
        w.wcs.crpix = [npix / 2, npix / 2]
        w.wcs.crval = [180.0, 20.0]
        w.wcs.cd = [[-ps / 3600, 0], [0, ps / 3600]]
        w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        return w

    imgs, wl = [], []
    for m in mjds:
        dt = (m - tref) * 24.0
        mo = 1.0 * dt / ps
        x = npix / 2 - mo * math.sin(math.radians(45))
        y = npix / 2 + mo * math.cos(math.radians(45))
        imgs.append(_gaussian_image(npix, x, y, 1500, seed=int(m * 1000)))
        wl.append(mkw())
    return imgs, wl, mjds, tref, ps


def test_gpu_coadd_matches_cpu():
    from ariadne.discovery.imaging.gpu_shift_stack import _gpu_coadd
    from ariadne.discovery.imaging.synthetic_tracking import predicted_shift, shift_and_stack

    imgs, _, mjds, tref, ps = _scene()
    shifts = [predicted_shift(m, tref, 1.0, 45.0, ps) for m in mjds]
    cc, ccov, _ = shift_and_stack(imgs, shifts, return_coverage=True, return_stack=True)
    img_t = torch.as_tensor(np.stack(imgs), dtype=torch.float32, device="cuda")
    ys = torch.arange(200, device="cuda", dtype=torch.float32)
    xs = torch.arange(200, device="cuda", dtype=torch.float32)
    bgy, bgx = torch.meshgrid(ys, xs, indexing="ij")
    gc, gcov, _ = _gpu_coadd(img_t, bgx, bgy, shifts, 200, 200)
    m = np.isfinite(cc) & np.isfinite(gc)
    assert np.max(np.abs(cc[m] - gc[m])) < 0.05  # bilinear rounding
    assert np.array_equal(ccov, gcov)


def test_gpu_recovers_moving_object_like_cpu():
    from ariadne.discovery.imaging.gpu_shift_stack import gpu_synthetic_tracking
    from ariadne.discovery.imaging.synthetic_tracking import fast_synthetic_tracking

    imgs, wl, mjds, tref, ps = _scene()
    kw = dict(
        t_ref_mjd=tref,
        rate_min_arcsec_hr=0.5,
        rate_max_arcsec_hr=3.0,
        n_rates=8,
        n_pa=8,
        snr_threshold=5.0,
        pixscale_arcsec=ps,
        n_top_per_hypothesis=3,
    )
    cpu = fast_synthetic_tracking(imgs, wl, mjds, **kw)
    gpu = gpu_synthetic_tracking(imgs, wl, mjds, **kw)
    assert cpu and gpu
    assert 0.7 <= gpu[0].rate_arcsec_hr <= 1.5
    pa_diff = abs(gpu[0].pa_deg - 45.0)
    pa_diff = min(pa_diff, 360 - pa_diff)
    assert pa_diff < 50
    # same best rate as CPU (same geomspace grid point)
    assert abs(gpu[0].rate_arcsec_hr - cpu[0].rate_arcsec_hr) < 1e-6
