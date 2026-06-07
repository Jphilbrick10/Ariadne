"""Planet Nine inference from real data: the geometry (orbit -> sky), the magnitude,
the orientation derivation, and the open-minded multi-population anomaly scan. All
network-free and kernel-free here (the secular dynamical test needs DE440 + minutes,
so it is exercised by scripts/run_p9_dynamical_test.py, not in the unit suite)."""

from __future__ import annotations

import math

import numpy as np

from ariadne.discovery.frontier.planet_nine import (
    P9Orbit,
    apparent_magnitude,
    derive_orientation,
    ecliptic_to_radec,
    orbit_xyz,
    outer_system_anomaly_scan,
    plane_warp,
    planetary_perturbation_constraint,
    selection_bias_null_test,
    sky_position,
)
from ariadne.discovery.frontier.tno_clustering import TNO


def _tno(name, a, e=0.5, inc=15.0, Om=130.0, om=60.0, q=None):
    q = a * (1 - e) if q is None else q
    return TNO(name=name, a=a, e=e, inc=inc, Omega=Om, omega=om, q=q)


def test_orbit_geometry_aphelion_distance():
    o = P9Orbit(a_au=460, e=0.25, i_deg=16, Omega_deg=130, omega_deg=68)
    apo = orbit_xyz(o, 180.0)  # aphelion
    _, _, d = ecliptic_to_radec(apo)
    assert abs(d - o.Q_au) < 1.0  # aphelion distance = a(1+e)
    peri = orbit_xyz(o, 0.0)
    _, _, dp = ecliptic_to_radec(peri)
    assert abs(dp - o.q_au) < 1.0


def test_ecliptic_to_radec_vernal_equinox():
    ra, dec, _ = ecliptic_to_radec(np.array([1.0, 0.0, 0.0]))
    assert abs(ra) < 1e-6 and abs(dec) < 1e-6


def test_sky_position_peaks_at_aphelion():
    o = P9Orbit(a_au=500, e=0.4, i_deg=20, Omega_deg=100, omega_deg=150)
    sp = sky_position(o)
    # time-weighted most-likely distance is near aphelion (Kepler 2nd law)
    assert sp["aphelion_distance_AU"] > o.a_au
    assert "constellation_hint" in sp


def test_apparent_magnitude_fainter_when_farther():
    near = apparent_magnitude(400)
    far = apparent_magnitude(800)
    assert far > near  # bigger distance -> fainter (larger mag)
    assert 18 < near < 30  # physically sane for a super-Earth


def test_derive_orientation_anti_aligned():
    # TNO perihelia clustered at varpi ~ 60 -> P9 perihelion at ~240
    ext = [_tno(f"c{i}", a=300, Om=0.0, om=60 + (i - 10)) for i in range(20)]
    o = derive_orientation(ext)
    diff = abs((o["varpi_p9_deg"] - 240 + 180) % 360 - 180)
    assert diff < 8


def test_anomaly_scan_flags_detached_objects():
    # 30 regular scattered objects + 5 detached (q>40, a>150 -> Neptune can't make)
    pop = [_tno(f"reg{i}", a=60 + i, e=0.5) for i in range(30)]
    pop += [_tno(f"det{i}", a=300 + 50 * i, e=0.85, q=45 + i) for i in range(5)]
    r = outer_system_anomaly_scan(tnos=pop, include_comets=False)
    assert r["detached_anomaly"]["n_decoupled_q_gt_40_a_gt_150"] == 5
    assert r["detached_anomaly"]["max_perihelion_AU"] >= 45
    assert "a_bins" in r and r["n_tno"] == 35


def test_anomaly_scan_counts_retrograde():
    pop = [_tno(f"pro{i}", a=100 + i, inc=10) for i in range(10)]
    pop += [_tno(f"retro{i}", a=100 + i, inc=120) for i in range(7)]
    r = outer_system_anomaly_scan(tnos=pop, include_comets=False)
    assert r["retrograde_polar"]["n_i_gt_90"] == 7


def test_planetary_perturbation_constraint():
    # a 6 M_E perturber at 500 AU drifts Saturn by ~km over the Cassini era,
    # above the ~100 m ranging precision -> it constrains; farther out it weakens
    near = planetary_perturbation_constraint(6, 500)
    assert near["per_planet"]["Saturn"]["drift_m"] > 300
    assert near["saturn_constrains"] is True
    far = planetary_perturbation_constraint(6, 1500)  # ~3x farther -> ~27x weaker
    assert far["per_planet"]["Saturn"]["drift_m"] < near["per_planet"]["Saturn"]["drift_m"]


def test_plane_warp_detects_tilt():
    # a plane tilted 30 deg from the ecliptic is far from the ~1.6 deg invariable plane
    ext = [_tno(f"p{k}", a=300, inc=30, Om=130, om=float(k)) for k in range(20)]
    w = plane_warp(ext)
    assert w["mean_plane_tilt_deg"] > 20
    assert w["angle_from_invariable_deg"] > 20


def test_selection_bias_null_flags_real_cluster():
    # a tightly node-clustered population must SURVIVE the galactic-plane bias null
    rng = np.random.default_rng(1)
    ext = [
        _tno(f"c{k}", a=300, inc=20, Om=130 + float(rng.normal(0, 5)), om=float(k * 37 % 360))
        for k in range(30)
    ]
    r = selection_bias_null_test(ext, n_trials=20, pool=200, seed=1)
    assert r["Omega"]["observed_R"] > 0.8
    assert r["Omega"]["survives_bias"] is True and r["Omega"]["p_value"] < 0.05
