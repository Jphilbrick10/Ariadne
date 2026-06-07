"""Automated frontier scanner: per-population anomaly suite + significance.
Network-free (synthetic populations; the JPL fetch is never hit)."""

from __future__ import annotations

import numpy as np

from ariadne.discovery.frontier.exotic_orbit_hunt import Orb
from ariadne.discovery.frontier.frontier_scanner import (
    POPULATIONS,
    _family_sig,
    _orientation_sig,
    scan_population,
)


def _rand_pop(n, seed=0, a_lo=3, a_hi=30):
    rng = np.random.default_rng(seed)
    return [
        Orb(
            f"o{k}",
            float(rng.uniform(a_lo, a_hi)),
            float(rng.uniform(0, 0.6)),
            float(rng.uniform(2, 40)),
            float(rng.uniform(0, 360)),
            float(rng.uniform(0, 360)),
            0.0,
        )
        for k in range(n)
    ]


def test_orientation_sig_uniform_vs_clustered():
    rng = np.random.default_rng(0)
    uni = _orientation_sig(rng.uniform(0, 360, 80))
    assert uni["p_value"] > 0.05
    clus = 130 + rng.normal(0, 8, 60)
    cl = _orientation_sig(clus)
    assert cl["R"] > 0.7 and cl["p_value"] < 0.01 and cl["significant"]


def test_family_sig_detects_planted_family():
    pop = _rand_pop(40, seed=1)
    for k in range(4):
        pop.append(Orb(f"fam{k}", 8.0 + 0.01 * k, 0.30, 12.0, 200.0, 90.0, 5.6))
    s = _family_sig(pop, n_null=150, seed=2)
    assert s["observed_pairs"] >= 6 and s["excess_sigma"] > 3


def test_scan_population_structure_and_flags():
    # a population with a clear node cluster must be flagged
    rng = np.random.default_rng(3)
    pop = [
        Orb(
            f"c{k}",
            300 + k,
            0.8,
            18.0,
            130 + float(rng.normal(0, 6)),
            float(rng.uniform(0, 360)),
            45.0,
        )
        for k in range(40)
    ]
    r = scan_population("clustered", pop)
    assert r["n"] == 40 and "node" in r and "apsidal" in r
    assert any("node clustering" in f for f in r["flags"])
    assert "plane_warp_deg" in r  # distant -> warp computed


def test_scan_population_exotic_count():
    pop = _rand_pop(20, seed=4) + [Orb("RETRO", 9.6, 0.8, 165.0, 100.0, 60.0, 1.9)]
    r = scan_population("exo", pop)
    assert r["n_exotic"] >= 1
    assert r["top_exotic"][0]["name"] == "RETRO"


def test_population_catalogue_nonempty():
    assert len(POPULATIONS) >= 10
    assert all(isinstance(q, str) for q in POPULATIONS.values())
