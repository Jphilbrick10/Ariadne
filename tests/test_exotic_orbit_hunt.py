"""Exotic-orbit hunt: the D-criterion, the chance-excess significance, and the
exotic-origin scoring. Network-free (synthetic orbits)."""

from __future__ import annotations

import numpy as np

from ariadne.discovery.frontier.exotic_orbit_hunt import (
    Orb,
    dsh_matrix,
    exotic_score,
    find_pairs,
    hunt,
    pair_significance,
)


def _o(name, a, e=0.5, i=10.0, Om=100.0, om=60.0, q=None):
    return Orb(name, a, e, i, Om, om, a * (1 - e) if q is None else q)


def test_dsh_zero_for_identical():
    objs = [_o("A", 5.0), _o("B", 5.0)]
    D = dsh_matrix(objs)
    assert D[0, 1] < 1e-6


def test_find_pairs_flags_a_twin():
    rng = np.random.default_rng(0)
    pop = [
        _o(
            f"r{k}",
            a=float(rng.uniform(3, 30)),
            i=float(rng.uniform(5, 40)),
            Om=float(rng.uniform(0, 360)),
            om=float(rng.uniform(0, 360)),
        )
        for k in range(40)
    ]
    pop += [
        _o("TWIN_A", 8.0, 0.30, 12.0, 200.0, 90.0),
        _o("TWIN_B", 8.01, 0.30, 12.05, 200.5, 90.3),
    ]  # near-identical q,e,i,orient
    pairs = find_pairs(pop, 0.05)
    names = {(p["a"], p["b"]) for p in pairs}
    assert ("TWIN_A", "TWIN_B") in names or ("TWIN_B", "TWIN_A") in names


def test_significance_excess_for_real_cluster():
    rng = np.random.default_rng(1)
    pop = [
        _o(
            f"r{k}",
            a=float(rng.uniform(3, 30)),
            i=float(rng.uniform(5, 40)),
            Om=float(rng.uniform(0, 360)),
            om=float(rng.uniform(0, 360)),
        )
        for k in range(40)
    ]
    # plant a tight 4-member family (shared orbit incl. orientation)
    for k in range(4):
        pop.append(_o(f"fam{k}", 8.0 + 0.01 * k, 0.30, 12.0, 200.0, 90.0))
    s = pair_significance(pop, threshold=0.05, n_null=200, seed=2)
    assert s["observed_pairs"] >= 6  # 4 members -> 6 pairs
    assert s["excess_sigma"] > 3.0


def test_exotic_score_flags():
    retro = exotic_score(_o("R", 9.6, 0.8, i=165))
    assert "RETROGRADE" in retro["flags"] and "CO-ORBITAL/Saturn" in retro["flags"]
    assert retro["score"] > 3
    iso = exotic_score(_o("I", 20, 1.05, i=120))
    assert "UNBOUND/ISO" in iso["flags"] and "RETROGRADE" in iso["flags"]
    normal = exotic_score(_o("N", 3.0, 0.1, i=5))
    assert normal["flags"] == []


def test_hunt_offline():
    rng = np.random.default_rng(3)
    pop = [
        _o(
            f"r{k}",
            a=float(rng.uniform(3, 30)),
            i=float(rng.uniform(5, 40)),
            Om=float(rng.uniform(0, 360)),
            om=float(rng.uniform(0, 360)),
        )
        for k in range(30)
    ]
    pop.append(_o("RETRO", 9.6, 0.8, i=165))
    out = hunt(pop, label="test", n_null=100)
    assert out["n"] == 31 and "pair_significance" in out
    assert out["top_exotic"][0]["name"] == "RETRO"
