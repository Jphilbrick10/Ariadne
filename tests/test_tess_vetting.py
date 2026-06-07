"""TESS transit search + coherence vetting. The search must recover an injected
period; the vetting must pass real-planet shapes and reject eclipsing-binary
impostors (deep, secondary eclipse, odd-even, grazing V). Network-free -- all on
synthetic light curves, so the MAST fetch is never hit."""

from __future__ import annotations

import numpy as np
import pytest

from ariadne.discovery.frontier.tess_vetting import (
    search_and_vet,
    search_transits,
    transit_incoherence,
    vet_features,
)


def _lc(
    period,
    depth,
    t0=1.0,
    dur=0.1,
    noise=0.0008,
    secondary=0.0,
    vshape=False,
    oddeven=0.0,
    days=27.0,
    cad=0.02,
    seed=0,
):
    rng = np.random.default_rng(seed)
    t = np.arange(0, days, cad)
    f = np.ones_like(t) + rng.normal(0, noise, t.size)
    ph = ((t - t0) / period + 0.5) % 1.0 - 0.5
    half = 0.5 * dur / period
    intr = np.abs(ph) < half
    if vshape:
        f[intr] -= depth * (1 - np.abs(ph[intr]) / half)
    else:
        prof = np.clip((half - np.abs(ph[intr])) / (0.25 * half), 0, 1)
        f[intr] -= depth * prof
    if oddeven > 0:
        ep = np.round((t - t0) / period).astype(int)
        f[intr & (ep % 2 == 1)] -= depth * oddeven
    if secondary > 0:
        f[np.abs(np.abs(ph) - 0.5) < half] -= depth * secondary
    return t, f


def test_bls_recovers_injected_period():
    t, f = _lc(3.5, 0.01)
    c = search_transits(t, f, min_period=1.0, max_period=10.0)
    assert c is not None
    assert abs(c.period - 3.5) < 0.05 or abs(c.period - 1.75) < 0.05  # P or P/2


@pytest.mark.parametrize(
    "name,kw",
    [
        ("clean_1pct", dict(period=3.5, depth=0.01)),
        ("shallow_0p3pct", dict(period=6.2, depth=0.003, dur=0.12)),
        ("jupiter_1p2pct", dict(period=8.1, depth=0.012, dur=0.15)),
    ],
)
def test_real_planet_shapes_pass(name, kw):
    t, f = _lc(**kw)
    c = search_and_vet(t, f, min_period=1.0, max_period=10.0)
    assert c.verdict == "planet-candidate", (name, c.coherence, c.features)
    assert c.coherence > 0.6


@pytest.mark.parametrize(
    "name,kw",
    [
        ("deep_v_secondary", dict(period=3.5, depth=0.18, vshape=True, secondary=0.5)),
        ("secondary_eclipse", dict(period=4.1, depth=0.05, secondary=0.6)),
        ("odd_even_blend", dict(period=3.5, depth=0.04, oddeven=0.7)),
        ("grazing_v", dict(period=2.8, depth=0.03, vshape=True)),
        ("deep_8pct", dict(period=5.0, depth=0.08)),
    ],
)
def test_eclipsing_binary_impostors_rejected(name, kw):
    t, f = _lc(**kw)
    c = search_and_vet(t, f, min_period=1.0, max_period=10.0)
    assert c.verdict == "likely-false-positive", (name, c.coherence, c.features)
    assert c.coherence < 0.42


def test_depth_penalty_is_one_sided():
    # a very shallow transit must incur NO depth penalty (planets can be tiny)
    E_shallow, sec = transit_incoherence({"depth": 0.0005, "snr": 30})
    assert sec["depth_excess"] == 0.0
    # a deep one must incur a large depth penalty
    E_deep, sec2 = transit_incoherence({"depth": 0.09, "snr": 30})
    assert sec2["depth_excess"] > 5.0


def test_shallow_low_snr_planet_not_overrejected():
    # a shallow, noisy real-planet shape must survive: shape features are unreliable
    # at low SNR, so they are down-weighted (the completeness fix).
    t, f = _lc(5.0, 0.0018, dur=0.12, noise=0.0030, seed=7)  # ~0.18% depth, low SNR
    c = search_and_vet(t, f, min_period=1.0, max_period=10.0)
    assert c.verdict == "planet-candidate", (c.coherence, c.features.get("snr"))


def test_flatten_removes_stellar_trend_keeps_transit():
    # a planet riding on slow stellar variability: flattening must remove the trend
    # (preserving the transit) so the candidate still vets as a planet
    from ariadne.discovery.frontier.tess_vetting import flatten_lightcurve

    t, f = _lc(4.0, 0.01, dur=0.1, noise=0.0008)
    f = f * (1.0 + 0.03 * np.sin(2 * np.pi * t / 9.0))  # 3% slow stellar wave
    tf, ff = flatten_lightcurve(t, f)
    assert np.std(ff) < np.std(f / np.median(f))  # trend reduced
    c = search_and_vet(t, f, min_period=1.0, max_period=10.0)
    assert c.verdict == "planet-candidate", (c.coherence,)


def test_period_alias_recovered():
    # plant a planet at P=2.0; even if BLS reports the 2x alias, vetting the half/
    # double period must recover a planet-candidate (real planets fold like EBs at 2x)
    t, f = _lc(2.0, 0.01, dur=0.08)
    c = search_and_vet(t, f, min_period=0.5, max_period=8.0)
    assert c.verdict == "planet-candidate"
    assert abs(c.period - 2.0) < 0.05 or abs(c.period - 1.0) < 0.05 or abs(c.period - 4.0) < 0.05


def test_secondary_significance_gate():
    # a clean shallow planet should not register a noise-driven secondary/odd-even
    t, f = _lc(6.2, 0.003, dur=0.12)
    c = search_and_vet(t, f, min_period=1.0, max_period=10.0)
    assert c.features["secondary_ratio"] == 0.0
    assert c.features["odd_even"] == 0.0
