"""ZTF/LSST physics-coherence novelty scoring. Known variable/transient classes
must cohere with a basin (low score); genuinely novel light curves -- physically
impossible feature combinations -- must cohere with nothing (high score). Network-
free: all synthetic, the ALeRCE fetch is never hit here."""

from __future__ import annotations

import numpy as np
import pytest

from ariadne.discovery.frontier.ztf_anomaly import (
    ANOMALY_TAU,
    CALIBRATED_BASINS_PATH,
    KNOWN_CLASS_BASINS,
    anomaly_score,
    lightcurve_features,
    load_basins,
    score_lightcurve,
)


def _t(n, span, seed=0):
    return np.sort(np.random.default_rng(seed).uniform(0, span, n))


def _rr_lyrae(seed=0):
    r = np.random.default_rng(seed)
    t = _t(150, 180, seed)
    return t, 15 + 0.8 * ((t / 0.55) % 1 - 0.5) + 0.05 * r.normal(size=t.size)


def _cepheid(seed=1):
    r = np.random.default_rng(seed)
    t = _t(150, 300, seed)
    return t, 14 + 0.8 * np.sin(2 * np.pi * t / 7.0) + 0.05 * r.normal(size=t.size)


def _agn(seed=2):
    r = np.random.default_rng(seed)
    t = _t(200, 400, seed)
    return t, 15 + np.cumsum(r.normal(0, 0.03, t.size))


def _eclipsing(seed=3):
    r = np.random.default_rng(seed)
    t = _t(150, 180, seed)
    ph = (t / 1.8) % 1
    m = 15 + np.where(np.abs(ph - 0.5) < 0.06, 0.5, 0.0) + np.where(np.abs(ph) < 0.06, 0.4, 0.0)
    return t, m + 0.04 * r.normal(size=t.size)


def test_lomb_scargle_recovers_short_period():
    # the periodogram MUST reach sub-day periods despite slow median cadence
    t, m = _rr_lyrae()
    f = lightcurve_features(t, m)
    assert abs(f["period"] - 0.55) < 0.02 or abs(f["period"] - 0.55 / 2) < 0.02
    assert f["ls_power"] > 0.3


@pytest.mark.parametrize("gen", [_rr_lyrae, _cepheid, _agn, _eclipsing])
def test_known_classes_score_low(gen):
    # synthetic curves use textbook amplitudes, so test the MECHANISM against the
    # physics priors they were built for; the calibrated (real-ZTF) basins are
    # validated on real data separately (see docs frontier scorecard).
    t, m = gen()
    r = score_lightcurve(t, m, basins=KNOWN_CLASS_BASINS)
    assert r["score"] < ANOMALY_TAU, (r["best_class"], r["score"])
    assert r["verdict"].startswith("known:")


def test_calibrated_basins_present_and_load():
    # the calibrated model ships committed; it must load and cover the core classes
    assert CALIBRATED_BASINS_PATH.exists()
    b = load_basins(force=True)
    assert {"RRL", "CEP", "QSO"} <= set(b)
    for cls, basin in b.items():
        assert basin["mu"] and basin["sig"]  # every class constrains >=1 axis


@pytest.mark.parametrize(
    "name,gen",
    [
        # 4-mag amplitude at an RR-Lyrae period: no pulsator gets that deep
        (
            "huge_amp_periodic",
            lambda: (
                lambda t: (
                    t,
                    15
                    + 4.0 * ((t / 0.5) % 1)
                    + 0.05 * np.random.default_rng(4).normal(size=t.size),
                )
            )(_t(200, 200, 4)),
        ),
        # strong periodicity at 0.08 d with 3-mag amplitude: fits nothing
        (
            "ultrafast_strong",
            lambda: (
                lambda t: (
                    t,
                    15
                    + 1.5 * np.sin(2 * np.pi * t / 0.08)
                    + 0.04 * np.random.default_rng(5).normal(size=t.size),
                )
            )(_t(250, 150, 5)),
        ),
    ],
)
def test_genuine_novelties_score_high(name, gen):
    t, m = gen()
    r = score_lightcurve(t, m)
    assert r["score"] > ANOMALY_TAU, (name, r["best_class"], r["score"])
    assert r["verdict"] == "anomalous"


def test_anomaly_score_structure():
    r = score_lightcurve(*_cepheid())
    assert set(r) >= {"score", "best_class", "verdict", "per_class", "features"}
    assert r["best_class"] in r["per_class"]
    assert r["per_class"][r["best_class"]] == round(r["score"], 2)


def test_insufficient_data_is_flagged():
    r = anomaly_score({})
    assert r["verdict"] == "insufficient-data" and r["score"] == float("inf")
