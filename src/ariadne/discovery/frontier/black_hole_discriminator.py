"""Is a distant gravitating object a black hole? -- a discriminator built from the
Coherence / Equation-of-One framework's OWN black-hole physics.

The framework (One; Unity-Black Holes Outline; THE COMPLETE MATHEMATICAL CHAIN)
models a black hole as a coherence-core object whose defining feature is a
COHERENCE HORIZON: the coherence-time field reaches zero on a closed surface,

    tau_c(r) / tau_inf = sqrt( -g_tt ) = sqrt( 1 - r_s/r ),   r_s = 2GM/c^2

with the recursive collapse stress  R(r) = (1/2) f'(r)/f(r) = d ln(tau_c)/dr -> inf
at r_s. Ordinary matter (a planet) has a HARD SURFACE at R_body >> r_s, so its
tau_c stays finite (~1) everywhere and R is finite. The framework is explicit that
a black hole is otherwise GR-IDENTICAL: all coherence-core corrections scale as
(L/M)^2 ~ 10^-88 for planetary masses, and the coherence force correction at orbital
distances (Form 30: F_coh/F_Newton = 1/(1 - r_s/r)) is ~10^-13. So GRAVITY AND THE
ORBIT CANNOT TELL THEM APART -- which this module verifies numerically.

The ONE thing tau_c -> 0 changes that is observable: light cannot leave a coherence
horizon, so a black hole is ELECTROMAGNETICALLY DARK despite carrying gravitational
mass, while a finite-tau_c planet REFLECTS sunlight and EMITS heat. The discriminator
therefore scores the hypotheses PLANET vs BLACK-HOLE against the real electromagnetic
(non-)detections using the Equation-of-One coherence selector, and returns a
significance.

Honest scope: this tells you BH-vs-planet ONLY IF a perturber's gravitational mass is
independently established AND the predicted sky region has deep multi-band coverage.
It does not, by itself, prove a perturber exists.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

# physical constants (SI)
G = 6.674e-11
C = 2.998e8
M_SUN = 1.989e30
M_EARTH = 5.972e24
R_EARTH = 6.371e6
AU = 1.495978707e11
SIGMA_SB = 5.670e-8
L_SUN = 3.828e26


# --------------------------------------------------------------------------- #
#  framework black-hole geometry (verbatim from the Coherence/EoO documents)
# --------------------------------------------------------------------------- #
def schwarzschild_radius(mass_kg: float) -> float:
    """Coherence horizon radius r_s = 2GM/c^2 (m). The surface where tau_c -> 0."""
    return 2 * G * mass_kg / C**2


def coherence_time_ratio(r_m: float, mass_kg: float) -> float:
    """tau_c(r)/tau_inf = sqrt(1 - r_s/r) (the framework's load-bearing identity).
    Zero at the coherence horizon, ~1 far away. Undefined (returns 0) inside r_s."""
    rs = schwarzschild_radius(mass_kg)
    x = 1.0 - rs / r_m
    return math.sqrt(x) if x > 0 else 0.0


def collapse_stress(r_m: float, mass_kg: float) -> float:
    """Recursive collapse stress R(r) = (1/2) f'/f for Schwarzschild f=1-r_s/r:
    R = r_s / (2 r (r - r_s)). Diverges at the coherence horizon -> the BH signature."""
    rs = schwarzschild_radius(mass_kg)
    denom = 2 * r_m * (r_m - rs)
    return rs / denom if denom != 0 else math.inf


def coherence_force_ratio(r_m: float, mass_kg: float) -> float:
    """Form 30: F_coh/F_Newton = 1/(1 - r_s/r). At orbital distances this is
    1 + ~1e-13, i.e. a black hole and a planet pull on the TNOs identically --
    the proof that dynamics cannot discriminate them."""
    rs = schwarzschild_radius(mass_kg)
    return 1.0 / (1.0 - rs / r_m) if r_m > rs else math.inf


def min_coherence_time(
    mass_kg: float, *, as_black_hole: bool, body_radius_m: float | None = None
) -> float:
    """The MINIMUM of tau_c/tau_inf the object reaches -- the framework discriminator.
    Black hole: tau_c -> 0 (coherence horizon). Planet: tau_c at its hard surface
    R_body >> r_s, so ~1 (finite, coherent matter)."""
    if as_black_hole:
        return 0.0
    R = body_radius_m if body_radius_m else R_EARTH
    return coherence_time_ratio(R, mass_kg)


# --------------------------------------------------------------------------- #
#  observable predictions of each hypothesis
# --------------------------------------------------------------------------- #
def planet_radius(mass_kg: float) -> float:
    """Mass-radius relation for a rocky/ice super-Earth (R ~ M^0.27, R_E units)."""
    return R_EARTH * (mass_kg / M_EARTH) ** 0.27


def reflected_V_mag(mass_kg: float, dist_au: float, *, albedo: float = 0.4) -> float:
    """Apparent V magnitude of sunlight reflected by a planet at heliocentric ~=
    geocentric distance dist_au. H + 5 log10(r*d), H from radius+albedo."""
    diam_km = 2 * planet_radius(mass_kg) / 1000.0
    H = 5 * math.log10(1329.0 / (diam_km * math.sqrt(albedo)))
    return H + 5 * math.log10(max(dist_au, 1e-3) * max(dist_au - 1.0, 1e-3))


def equilibrium_temp(dist_au: float, *, albedo: float = 0.4, internal_K: float = 30.0) -> float:
    """Planet equilibrium temperature (K) at dist_au, with a floor from residual
    internal heat (a super-Earth stays ~30 K even in the dark)."""
    T_eq = 278.0 * (1 - albedo) ** 0.25 / math.sqrt(dist_au)
    return (T_eq**4 + internal_K**4) ** 0.25


def thermal_peak_um(temp_K: float) -> float:
    """Wien peak wavelength (micron). WISE bands: W1 3.4, W2 4.6, W3 12, W4 22 um."""
    return 2898.0 / temp_K


def accretion_luminosity(
    mass_kg: float, dist_au: float, *, n_ism_cm3: float = 0.1, v_kms: float = 1.0
) -> float:
    """Bondi accretion luminosity of a PBH from the tenuous outer-solar-system
    medium (W). Vanishingly small -> a PBH gives essentially no positive signal
    either; the discriminator is darkness, not a glow."""
    rho = n_ism_cm3 * 1.67e-24 * 1e6  # kg/m^3 (protons)
    v = v_kms * 1000.0
    mdot = 4 * math.pi * (G * mass_kg) ** 2 * rho / v**3
    return 0.1 * mdot * C**2  # 10% radiative efficiency


@dataclass
class Hypothesis:
    name: str
    is_black_hole: bool
    reflected_V: float  # apparent V mag (inf if dark)
    thermal_peak_um: float
    thermal_temp_K: float
    accretion_W: float
    min_tau_c: float  # framework coherence-horizon discriminator
    prior_weight: float = 1.0  # E_q sector: astrophysical prior


def build_hypotheses(mass_earths: float, dist_au: float, *, pbh_prior: float = 0.02) -> dict:
    """The PLANET and BLACK-HOLE hypotheses for a perturber of given mass/distance,
    with their framework + observable predictions. pbh_prior = a-priori odds a
    planetary-mass dark object is a PBH (planetary-mass PBHs are a small allowed
    fraction of dark matter; default conservative 2%)."""
    m = mass_earths * M_EARTH
    T = equilibrium_temp(dist_au)
    planet = Hypothesis(
        "planet",
        False,
        reflected_V_mag(m, dist_au),
        thermal_peak_um(T),
        T,
        0.0,
        min_coherence_time(m, as_black_hole=False, body_radius_m=planet_radius(m)),
        prior_weight=1.0,
    )
    bh = Hypothesis(
        "black_hole",
        True,
        math.inf,  # dark: no reflection
        0.0,
        0.0,
        accretion_luminosity(m, dist_au),
        min_coherence_time(m, as_black_hole=True),
        prior_weight=pbh_prior,
    )
    return {
        "planet": planet,
        "black_hole": bh,
        "r_s_m": schwarzschild_radius(m),
        "planet_radius_m": planet_radius(m),
    }


# --------------------------------------------------------------------------- #
#  Equation-of-One model selection: which hypothesis coheres with the data?
# --------------------------------------------------------------------------- #
@dataclass
class Observation:
    """The real electromagnetic (non-)detection in the predicted region."""

    optical_V_limit: float  # surveys reach this depth with NO detection
    optical_covered: float  # fraction of the predicted region searched (0..1)
    detected_optical: bool = False
    ir_W4_limit_Jy: float = 0.0  # WISE 22um limit (0 = not used)
    detected_ir: bool = False
    mass_confirmed: bool = False  # is the gravitational mass independently real?


def _nondetect_likelihood(
    predicted_V: float, V_limit: float, covered: float, *, p_miss_searched: float = 0.02
) -> float:
    """P(no optical detection | hypothesis). A reflected source brighter than the
    survey limit is missed only if it sits in the UNSEARCHED gap (1-covered) or via
    a small false-negative rate where searched; a source fainter than the limit is
    undetectable everywhere, so non-detection is uninformative (likelihood ~1)."""
    if predicted_V >= V_limit:
        return 1.0  # too faint to see anyway
    return (1.0 - covered) + covered * p_miss_searched


def discriminate(
    mass_earths: float, dist_au: float, obs: Observation, *, pbh_prior: float = 0.02
) -> dict:
    """Score PLANET vs BLACK-HOLE with the Equation-of-One coherence selector and
    return the posterior + a significance. Proper-Bayesian form: coherence(h) =
    prior(h) * P(data|h) = exp(-(E_q + E_c)/2), with E_q = -2 ln prior and
    E_c = -2 ln P(data|h). Significance = the Gaussian sigma of the posterior."""
    h = build_hypotheses(mass_earths, dist_au, pbh_prior=pbh_prior)
    planet, bh = h["planet"], h["black_hole"]
    energies = {}
    for hyp in (planet, bh):
        if obs.detected_optical:
            # a reflective planet explains a detection; a coherence horizon (dark) does not
            L = 0.85 if hyp is planet else 1e-4
        else:
            L = (
                _nondetect_likelihood(hyp.reflected_V, obs.optical_V_limit, obs.optical_covered)
                if hyp is planet
                else 1.0
            )  # a BH is always optically dark
        Ec = -2.0 * math.log(max(L, 1e-9))
        Eq = -2.0 * math.log(max(hyp.prior_weight, 1e-9))  # prior as the E_q sector
        energies[hyp.name] = {
            "E_c": round(Ec, 3),
            "E_q": round(Eq, 3),
            "E_total": Ec + Eq,
            "coherence": math.exp(-0.5 * (Ec + Eq)),
            "likelihood_nondetect": round(L, 5),
        }
    s = sum(v["coherence"] for v in energies.values()) or 1.0
    post = {k: v["coherence"] / s for k, v in energies.items()}
    p_bh = post["black_hole"]
    # significance: Gaussian sigma equivalent of the posterior probability
    sig = _p_to_sigma(p_bh)
    optical_informative = planet.reflected_V < obs.optical_V_limit
    return {
        "mass_earths": mass_earths,
        "dist_au": dist_au,
        "r_s_cm": round(h["r_s_m"] * 100, 2),
        "planet_radius_km": round(h["planet_radius_m"] / 1000, 0),
        "planet_predicted_V": round(planet.reflected_V, 1),
        "planet_min_tau_c": round(planet.min_tau_c, 6),
        "bh_min_tau_c": bh.min_tau_c,
        "coherence_force_ratio_at_dist": coherence_force_ratio(dist_au * AU, mass_earths * M_EARTH),
        "energies": energies,
        "posterior": {k: round(v, 4) for k, v in post.items()},
        "P_black_hole": round(p_bh, 4),
        "significance_sigma": round(sig, 2),
        "verdict": _verdict(p_bh, obs, optical_informative),
    }


def _p_to_sigma(p: float) -> float:
    """Two-sided Gaussian sigma for a probability p (p=0.997 -> ~3 sigma)."""
    p = min(max(p, 1e-9), 1 - 1e-9)
    # inverse error function via Winitzki approximation
    x = 2 * p - 1
    a = 0.147
    ln = math.log(1 - x * x)
    t = 2 / (math.pi * a) + ln / 2
    return math.sqrt(2) * math.copysign(math.sqrt(math.sqrt(t * t - ln / a) - t), x)


def _verdict(p_bh: float, obs: Observation, optical_informative: bool = True) -> str:
    if not obs.mass_confirmed:
        return (
            "inconclusive: no independently-confirmed gravitational mass -- "
            "cannot claim a dark object until a perturber is established"
        )
    if obs.detected_optical:
        return "planet favored: a reflected optical source IS present (finite tau_c)"
    if not optical_informative:
        return (
            "indistinguishable: a planet this cold/distant is itself below the "
            "survey limit, so optical darkness does not favor a coherence horizon"
        )
    if obs.optical_covered < 0.5:
        return "inconclusive: predicted region not deeply searched (coverage < 50%)"
    if p_bh > 0.997:
        return "black-hole / dark-mass favored at >3 sigma (gravitationally real, optically absent)"
    if p_bh > 0.95:
        return "black-hole / dark-mass favored (~2 sigma)"
    if p_bh < 0.05:
        return "planet favored: a reflected source should be present in the searched region but is absent only in the gap"
    return "ambiguous: optical non-detection is real but prior-dominated (need higher coverage or a PBH prior)"
