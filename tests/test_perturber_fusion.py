"""Multi-channel coherence fusion: combine every perturber signature into one
Equation-of-One estimate and classify each channel. Network-free (synthetic TNOs)."""

from __future__ import annotations

import numpy as np

from ariadne.discovery.frontier.perturber_fusion import (
    extract_channels,
    fuse,
    fusion_significance,
    run,
)
from ariadne.discovery.frontier.tno_clustering import TNO


def _tno(name, a, Om, om=60.0, e=0.8, inc=18.0, q=45.0):
    return TNO(name=name, a=a, e=e, inc=inc, Omega=Om, omega=om, q=q)


def _clustered(n=30, node=130.0, spread=8.0, seed=0):
    rng = np.random.default_rng(seed)
    return [
        _tno(
            f"c{k}",
            a=300 + k,
            Om=node + float(rng.normal(0, spread)),
            om=float(rng.uniform(0, 360)),
        )
        for k in range(n)
    ]


def test_fuse_recovers_node_and_anti_aligns():
    ext = _clustered(node=130.0)
    out = fuse(extract_channels(ext))
    o = out["fused_orbit"]
    assert abs((o["Omega_deg"] - 130 + 180) % 360 - 180) < 12  # node recovered
    # perihelion fused to be anti-aligned with the TNO apsidal mean
    ch = extract_channels(ext)
    assert abs((o["varpi_deg"] - (ch.apsidal_deg + 180) + 180) % 360 - 180) < 1


def test_node_is_evidential_apsidal_permissive():
    out = fuse(extract_channels(_clustered()))
    assert "node" in out["evidential_channels"]
    assert out["channel_verdict"]["node"].startswith("POSITIVE")
    assert out["channel_verdict"]["apsidal"].startswith("PERMISSIVE")  # default p=0.085


def test_exclusion_channels_present():
    out = fuse(extract_channels(_clustered()))
    for k in ("cassini", "comet", "dynamical"):
        assert k in out["channel_verdict"]
    # fused mass/distance must not be Cassini-excluded
    assert "CONFLICT" not in out["channel_verdict"]["cassini"]


def test_significance_runs_and_is_high_for_strong_cluster():
    ext = _clustered(spread=5.0)  # tight cluster
    s = fusion_significance(extract_channels(ext), n_null=300)
    assert s["combined_statistic"] > 0.8
    assert s["isotropic_null_p"] < 0.05


def test_run_offline():
    out = run(tnos=_clustered(node=100.0))
    assert out["n_extreme"] == 30
    assert abs((out["fused_orbit"]["Omega_deg"] - 100 + 180) % 360 - 180) < 12
    assert out["n_evidential"] >= 1
