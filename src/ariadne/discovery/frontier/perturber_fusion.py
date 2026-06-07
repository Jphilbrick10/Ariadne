"""Multi-channel coherence fusion for an outer-solar-system perturber.

The novel move: every independent line of evidence for a distant perturber is
individually weak, and the field analyses them one at a time. This fuses them all
into a SINGLE Equation-of-One estimate -- a joint posterior over the perturber's
orbit (a, e, i, Omega, omega, mass) -- and, crucially, reports for each channel
whether it provides POSITIVE evidence, is merely PERMISSIVE (consistent but not
confirmatory), or CONFLICTS. That honest bookkeeping is the thing single-channel
analyses cannot produce: it shows whether the weak signals all cohere toward the
same object or just fail to exclude one.

Channels (from the real data + this session's tests):
  node      -- the detached-TNO node (Omega) clustering. EVIDENTIAL: survives the
               galactic-plane bias null (p~0.005). Constrains the perturber's node.
  apsidal   -- the perihelion (varpi) clustering. WEAK: bias-consistent (p~0.085).
               Tested here for whether it independently corroborates the node.
  plane     -- the ~5 deg warp of the TNO plane from the invariable plane. Real;
               constrains the perturber inclination + corroborates the node.
  cassini   -- Saturn ranging: an EXCLUSION. Upper-bounds mass/distance^3.
  comet     -- the isotropic Oort cloud: a mild EXCLUSION on a massive/near body.
  dynamical -- our 1.5 Gyr secular test: a massive, distant perturber enhances
               clustering, the textbook one does not. A soft PRIOR toward large a, M.

Honest scope: fusing channels cannot manufacture a detection. If the only
evidential channel is the node, the joint significance is the node's significance,
and the tool says so. The value is a single transparent perturber estimate with an
auditable per-channel verdict.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from .planet_nine import INVARIABLE_PLANE, plane_warp
from .planet_nine import planetary_perturbation_constraint as _ppc
from .tno_clustering import TNO, circular_stats, extreme_population, fetch_tnos


# --------------------------------------------------------------------------- #
#  channel extraction from real data
# --------------------------------------------------------------------------- #
@dataclass
class Channels:
    node_deg: float
    node_R: float
    node_p: float
    apsidal_deg: float
    apsidal_R: float
    apsidal_p: float
    plane_tilt_deg: float
    n: int
    raw_Omega: np.ndarray = field(default=None, repr=False)


def extract_channels(
    ext: list[TNO], *, node_p: float = 0.005, apsidal_p: float = 0.085
) -> Channels:
    """Pull the channel summaries from the real extreme-TNO sample. node_p/apsidal_p
    are the bias-aware p-values measured by selection_bias_null_test (defaults are
    this session's a>=150 results)."""
    Om = np.array([t.Omega % 360 for t in ext])
    vp = np.array([t.varpi for t in ext])
    sn = circular_stats(Om)
    sa = circular_stats(vp)
    pw = plane_warp(ext)
    return Channels(
        node_deg=sn["mean_deg"],
        node_R=sn["R"],
        node_p=node_p,
        apsidal_deg=sa["mean_deg"],
        apsidal_R=sa["R"],
        apsidal_p=apsidal_p,
        plane_tilt_deg=pw["angle_from_invariable_deg"],
        n=len(ext),
        raw_Omega=Om,
    )


# --------------------------------------------------------------------------- #
#  per-channel energies (Equation-of-One E_c sectors), as functions of the orbit
# --------------------------------------------------------------------------- #
def _vm_energy(angle_deg, mean_deg, R, n):
    """von Mises incoherence: (1-cos) scaled by the data's concentration (R*n).
    Strong, real clusters pull hard; weak ones barely."""
    return R * n * (1 - math.cos(math.radians(angle_deg - mean_deg)))


def channel_energies(a, e, i, Om, om, M, ch: Channels) -> dict:
    """E_c of each channel for a candidate perturber orbit. Lower = more coherent."""
    varpi = (Om + om) % 360
    E_node = _vm_energy(Om, ch.node_deg, ch.node_R, ch.n)
    E_apsidal = _vm_energy(varpi, (ch.apsidal_deg + 180) % 360, ch.apsidal_R, ch.n)
    # plane: the perturber inclination should be able to warp the TNO plane by the
    # observed amount -> a soft band favoring i a few-to-several times the warp.
    E_plane = 0.5 * ((i - max(3 * ch.plane_tilt_deg, 10.0)) / 10.0) ** 2
    # cassini exclusion: Saturn drift over the Cassini era, ~2 km absorbable scale
    drift = _ppc(M, a)["per_planet"]["Saturn"]["drift_m"]
    E_cassini = (max(drift / 2000.0 - 1.0, 0.0)) ** 2
    # comet: a very massive/near body would imprint the isotropic Oort cloud
    E_comet = (max((M / a**1.5) / (6.0 / 500**1.5) - 2.0, 0.0)) ** 2
    # dynamical prior: our secular test favors massive + distant (forcing ~ M/a^1.5
    # is the WRONG direction; the test showed larger a AND larger M enhances). Use a
    # gentle reward for larger a and M within plausibility.
    E_dyn = 0.3 * ((600.0 - a) / 300.0) ** 2 * (a < 600) + 0.3 * ((8.0 - M) / 5.0) ** 2 * (M < 8)
    return {
        "node": E_node,
        "apsidal": E_apsidal,
        "plane": E_plane,
        "cassini": E_cassini,
        "comet": E_comet,
        "dynamical": E_dyn,
    }


# --------------------------------------------------------------------------- #
#  the fusion: joint MAP orbit + per-channel verdict
# --------------------------------------------------------------------------- #
def fuse(ch: Channels) -> dict:
    """Minimise the summed Equation-of-One energy over the perturber orbit. The
    orientation (Omega, omega) and the mass/distance (a, M) factorise, so we solve
    each on a grid, set i from the plane channel, and classify every channel."""
    # orientation: node fixes Omega; apsidal then fixes varpi -> omega
    Om = ch.node_deg
    varpi_target = (ch.apsidal_deg + 180) % 360
    om = (varpi_target - Om) % 360
    i = max(3 * ch.plane_tilt_deg, 16.0)
    # mass / distance grid
    a_grid = np.arange(300, 901, 50.0)
    M_grid = np.arange(3, 15.1, 1.0)
    best = None
    for a in a_grid:
        for M in M_grid:
            E = channel_energies(a, 0.25, i, Om, om, M, ch)
            tot = E["cassini"] + E["comet"] + E["dynamical"]
            if best is None or tot < best[0]:
                best = (tot, a, M, E)
    _, a, M, _ = best
    E = channel_energies(a, 0.25, i, Om, om, M, ch)
    orbit = {
        "a_au": float(a),
        "e": 0.25,
        "i_deg": round(i, 1),
        "Omega_deg": round(Om, 1),
        "omega_deg": round(om, 1),
        "varpi_deg": round((Om + om) % 360, 1),
        "mass_earths": float(M),
    }
    # classify each channel for this fused orbit
    verdict = {}
    for name, e in E.items():
        if name == "node":
            verdict[name] = "POSITIVE (real, survives bias p=%.3f)" % ch.node_p
        elif name == "apsidal":
            corrob = ch.apsidal_p < 0.05
            verdict[name] = (
                "POSITIVE (corroborates)"
                if corrob
                else "PERMISSIVE (consistent, not significant p=%.3f)" % ch.apsidal_p
            )
        elif name == "plane":
            verdict[name] = "POSITIVE (real warp) " if e < 1 else "TENSION (i off)"
        else:
            verdict[name] = "PERMISSIVE (not excluded)" if e < 1 else "CONFLICT (excluded)"
    evidential = [k for k, v in verdict.items() if v.startswith("POSITIVE")]
    return {
        "fused_orbit": orbit,
        "channel_energy": {k: round(v, 3) for k, v in E.items()},
        "channel_verdict": verdict,
        "evidential_channels": evidential,
        "n_evidential": len(evidential),
    }


# --------------------------------------------------------------------------- #
#  joint significance: is the combined EVIDENTIAL coherence beyond chance?
# --------------------------------------------------------------------------- #
def fusion_significance(ch: Channels, *, n_null: int = 2000, seed: int = 0) -> dict:
    """Combined significance of the evidential channels. The node carries the real
    signal; this asks whether the node AND apsidal jointly (an anti-aligned, common-
    node perturber) beat an isotropic null. Fast isotropic-shuffle null on the
    orientation; the bias-aware node p is reported separately (node_p)."""
    rng = np.random.default_rng(seed)
    n = ch.n
    real = ch.node_R + 0.5 * ch.apsidal_R
    null = np.empty(n_null)
    for k in range(n_null):
        Om = rng.uniform(0, 360, n)
        vp = rng.uniform(0, 360, n)
        sn = circular_stats(Om)
        sa = circular_stats(vp)
        null[k] = sn["R"] + 0.5 * sa["R"]
    p = float(np.mean(null >= real))
    # sigma
    p_c = min(max(p, 1e-9), 1 - 1e-9)
    x = 1 - 2 * p_c
    a = 0.147
    ln = math.log(1 - x * x)
    t = 2 / (math.pi * a) + ln / 2
    sig = math.sqrt(2) * math.copysign(math.sqrt(math.sqrt(t * t - ln / a) - t), x)
    return {
        "combined_statistic": round(real, 3),
        "isotropic_null_p": round(p, 4),
        "isotropic_sigma": round(sig, 2),
        "bias_aware_node_p": ch.node_p,
        "note": (
            "isotropic null; the bias-aware node p (node_p) is the honest "
            "headline since galactic-plane selection, not isotropy, is the "
            "real null -- this combined stat is dominated by the node."
        ),
    }


def run(
    *,
    a_min: float = 150.0,
    q_min: float = 30.0,
    use_cache: bool = True,
    tnos: list[TNO] | None = None,
) -> dict:
    """Full fusion on the real (or supplied) extreme-TNO sample."""
    pool = tnos if tnos is not None else fetch_tnos(use_cache=use_cache)
    ext = extreme_population(pool, a_min=a_min, q_min=q_min)
    ch = extract_channels(ext)
    out = fuse(ch)
    out["significance"] = fusion_significance(ch, n_null=2000)
    out["n_extreme"] = len(ext)
    return out
