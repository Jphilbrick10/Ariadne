"""Black-hole discriminator built from the Coherence/Equation-of-One framework.
Checks the framework geometry (coherence horizon r_s, tau_c -> 0, force ratio ~1 at
orbital distance) and the Equation-of-One model selection (planet vs black hole)."""

from __future__ import annotations

import math

from ariadne.discovery.frontier.black_hole_discriminator import (
    AU,
    M_EARTH,
    M_SUN,
    Observation,
    coherence_force_ratio,
    coherence_time_ratio,
    discriminate,
    min_coherence_time,
    schwarzschild_radius,
)


def test_coherence_horizon_radius():
    # a 5 Earth-mass black hole has a ~few-cm horizon; the Sun ~3 km
    rs = schwarzschild_radius(5 * M_EARTH)
    assert 0.03 < rs < 0.08  # ~4-5 cm
    assert abs(schwarzschild_radius(M_SUN) - 2950) < 100  # ~2.95 km


def test_coherence_time_zero_at_horizon_one_far():
    m = 6 * M_EARTH
    rs = schwarzschild_radius(m)
    assert coherence_time_ratio(rs, m) == 0.0  # tau_c -> 0 at horizon
    assert coherence_time_ratio(rs * 0.5, m) == 0.0  # inside: undefined -> 0
    assert coherence_time_ratio(1e6 * rs, m) > 0.999  # ~1 far away


def test_orbit_cannot_distinguish():
    # the framework's own prediction: coherence force = Newton at orbital distance
    ratio = coherence_force_ratio(500 * AU, 6 * M_EARTH)
    assert abs(ratio - 1.0) < 1e-10


def test_min_tau_c_separates_bh_from_planet():
    m = 6 * M_EARTH
    assert min_coherence_time(m, as_black_hole=True) == 0.0
    planet = min_coherence_time(m, as_black_hole=False, body_radius_m=1.0e7)
    assert planet > 0.999  # finite, coherent matter


def test_unconfirmed_mass_is_inconclusive():
    r = discriminate(
        6,
        500,
        Observation(
            optical_V_limit=24, optical_covered=0.9, detected_optical=False, mass_confirmed=False
        ),
    )
    assert "inconclusive" in r["verdict"]


def test_optical_detection_favors_planet():
    r = discriminate(
        6,
        500,
        Observation(
            optical_V_limit=24, optical_covered=0.9, detected_optical=True, mass_confirmed=True
        ),
    )
    assert r["P_black_hole"] < 0.05
    assert "planet favored" in r["verdict"]


def test_faint_planet_is_indistinguishable():
    # at 1500 AU a 6 M_E planet is V~27, below a V=24 survey -> darkness uninformative
    r = discriminate(
        6,
        1500,
        Observation(
            optical_V_limit=24, optical_covered=0.99, detected_optical=False, mass_confirmed=True
        ),
    )
    assert "indistinguishable" in r["verdict"]


def test_bh_probability_rises_with_coverage():
    def p(cov):
        return discriminate(
            6,
            500,
            Observation(
                optical_V_limit=24, optical_covered=cov, detected_optical=False, mass_confirmed=True
            ),
        )["P_black_hole"]

    assert p(0.99) > p(0.9) > p(0.5)  # more coverage -> stronger dark-mass case


def test_significance_is_prior_dominated():
    # the honest finding: even at near-complete coverage, the call depends on the prior
    lo = discriminate(
        6,
        500,
        Observation(
            optical_V_limit=24, optical_covered=0.999, detected_optical=False, mass_confirmed=True
        ),
        pbh_prior=0.001,
    )
    hi = discriminate(
        6,
        500,
        Observation(
            optical_V_limit=24, optical_covered=0.999, detected_optical=False, mass_confirmed=True
        ),
        pbh_prior=0.5,
    )
    assert hi["P_black_hole"] > 0.9 and lo["P_black_hole"] < 0.1
