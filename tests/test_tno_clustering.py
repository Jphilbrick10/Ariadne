"""Extreme-TNO orientation outlier hunt: the circular statistics, the coherence
energy, and the outlier ranking must behave correctly. Network-free -- everything
runs on synthetic populations passed in, so the JPL fetch is never hit here."""

from __future__ import annotations

import math

import numpy as np

from ariadne.discovery.frontier.tno_clustering import (
    TNO,
    _ecliptic_lon_to_radec,
    circular_stats,
    extreme_population,
    orientation_energy,
    outlier_hunt,
    predict_planet_nine,
    rank_by_coherence,
    select_coherent_subsample,
)


def _mk(name, a, Omega, omega=0.0, q=40.0, e=0.8, inc=20.0):
    return TNO(name=name, a=a, e=e, inc=inc, Omega=Omega, omega=omega, q=q)


def test_circular_stats_uniform_vs_clustered():
    # uniform ring -> R ~ 0, not significant
    uni = np.linspace(0, 360, 60, endpoint=False)
    su = circular_stats(uni)
    assert su["R"] < 0.05 and su["rayleigh_p"] > 0.5
    # tight cluster near 130 deg -> high R, tiny p
    clus = 130 + np.random.default_rng(0).normal(0, 12, 40)
    sc = circular_stats(clus)
    assert sc["R"] > 0.7 and sc["rayleigh_p"] < 1e-3
    assert abs((sc["mean_deg"] - 130 + 180) % 360 - 180) < 8


def test_orientation_energy_bounds():
    assert orientation_energy(130, 130) == 0.0  # aligned -> 0
    assert abs(orientation_energy(130, 310) - 2.0) < 1e-9  # anti-aligned -> 2
    assert abs(orientation_energy(220, 130) - 1.0) < 1e-9  # orthogonal -> 1


def test_outlier_ranked_first():
    # 20 objects clustered at Omega=130, one planted anti-aligned at 310
    pop = [_mk(f"c{i}", a=300, Omega=130 + (i - 10)) for i in range(20)]
    pop.append(_mk("ODDBALL", a=300, Omega=310))
    ranked = rank_by_coherence(pop, axis="Omega")
    assert ranked[0]["name"] == "ODDBALL"  # least coherent first
    assert ranked[0]["energy"] > ranked[-1]["energy"]
    assert ranked[-1]["coherence"] > ranked[0]["coherence"]


def test_selector_prefers_cleaner_subsample():
    # low-a objects scattered, high-a objects tightly clustered -> selector should
    # pick a high a_min as the most coherent subsample
    rng = np.random.default_rng(1)
    scattered = [_mk(f"lo{i}", a=160 + i, Omega=float(rng.uniform(0, 360))) for i in range(30)]
    tight = [_mk(f"hi{i}", a=320 + i, Omega=130 + float(rng.normal(0, 8))) for i in range(20)]
    out = select_coherent_subsample(scattered + tight, axis="Omega", a_grid=(150, 250, 320))
    assert out["selected"]["a_min"] >= 250
    assert out["selected"]["mean_energy"] < out["table"][0]["mean_energy"]


def test_extreme_population_filter():
    pool = [
        _mk("keep", a=200, Omega=130, q=40),
        _mk("near", a=200, Omega=130, q=20),
        _mk("small", a=100, Omega=130, q=40),
        TNO("nan", float("nan"), 0, 0, 0, 0, 0),
    ]
    ext = extreme_population(pool, a_min=150, q_min=30)
    assert [t.name for t in ext] == ["keep"]


def test_outlier_hunt_offline():
    pop = [_mk(f"c{i}", a=300, Omega=130 + (i - 10)) for i in range(20)]
    rep = outlier_hunt(tnos=pop)
    assert rep["n_extreme"] == 20
    assert "Omega" in rep and rep["Omega"]["R"] > 0.5
    assert len(rep["most_anomalous"]) > 0


def test_planet_nine_is_apsidally_anti_aligned():
    # TNO perihelia clustered at varpi ~ 60 (omega=60, Omega=0) -> P9 perihelion
    # must point to the opposite side (~240), and P9 aphelion back near 60.
    pop = [_mk(f"c{i}", a=400, Omega=0.0, omega=60 + (i - 10)) for i in range(20)]
    p = predict_planet_nine(tnos=pop, a_min=300)
    assert (
        abs((p["predicted_orientation"]["longitude_of_perihelion_deg"] - 240 + 180) % 360 - 180) < 8
    )
    # mass/distance are explicitly NOT derived from this geometry
    assert "from_dynamics_NOT_this_geometry" in p


def test_ecliptic_to_radec_anchors():
    ra, dec = _ecliptic_lon_to_radec(0.0, 0.0)  # vernal equinox
    assert abs(ra) < 1e-6 and abs(dec) < 1e-6
    ra2, dec2 = _ecliptic_lon_to_radec(90.0, 0.0)  # summer solstice
    assert abs(ra2 - 90.0) < 1.0 and abs(dec2 - 23.4) < 0.5
