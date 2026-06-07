"""Hunt for genetically-related orbit pairs/clusters and exotic-origin objects in
under-mapped solar-system populations.

Two objects on nearly identical orbits can share an origin: a recent breakup, a
rotational-fission pair, or a group captured together. The Southworth-Hawkins
D-criterion measures orbital similarity; pairs below a small D are candidate
associations. The honest hard part is significance -- low-inclination objects near a
planet cluster naturally -- so we compare the observed pair count to a null that
shuffles the orientation angles (destroying real apsidal/nodal coincidences while
keeping the a/e/i distribution). Only an EXCESS over that null is interesting.

Also scores each object for an exotic origin (retrograde, polar, co-orbital, unbound)
-- the captured-interstellar and Oort-return candidates.

Honest scope: these are CANDIDATE associations on osculating elements (noisier than
proper elements), to be checked against known families and refined with proper
elements + backward integration. The tool flags; it does not confirm.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass

import numpy as np


def _base_designation(name: str) -> str:
    """Strip a fragment suffix so split-comet pieces share a parent id:
    'C/1882 R1-C (Great September comet)' -> 'C/1882 R1'. Used to exclude
    same-parent fragments (a known, not novel, source of identical orbits)."""
    n = name.split("(")[0].strip()
    return re.sub(r"-[A-Z]{1,2}$", "", n).strip()


def same_parent(a: str, b: str) -> bool:
    base = _base_designation(a)
    return bool(base) and base == _base_designation(b)


PLANET_A = {
    "Mercury": 0.39,
    "Venus": 0.72,
    "Earth": 1.0,
    "Mars": 1.52,
    "Jupiter": 5.20,
    "Saturn": 9.54,
    "Uranus": 19.2,
    "Neptune": 30.1,
}


@dataclass
class Orb:
    name: str
    a: float
    e: float
    inc: float
    Omega: float
    omega: float
    q: float


# --------------------------------------------------------------------------- #
#  Southworth-Hawkins D-criterion (vectorized)
# --------------------------------------------------------------------------- #
def dsh_matrix(objs: list[Orb]) -> np.ndarray:
    """Pairwise Southworth-Hawkins D for a list of orbits. q in AU."""
    e = np.array([o.e for o in objs])
    q = np.array([o.q for o in objs])
    i = np.radians([o.inc for o in objs])
    Om = np.radians([o.Omega for o in objs])
    om = np.radians([o.omega for o in objs])
    iA, iB = i[:, None], i[None, :]
    OA, OB = Om[:, None], Om[None, :]
    cosI = np.clip(np.cos(iA) * np.cos(iB) + np.sin(iA) * np.sin(iB) * np.cos(OB - OA), -1, 1)
    I = np.arccos(cosI)
    dO = OB - OA
    t = np.clip(np.cos((iA + iB) / 2) * np.sin(dO / 2) / np.maximum(np.cos(I / 2), 1e-9), -1, 1)
    P = (om[None, :] - om[:, None]) + 2 * np.arcsin(t)
    P = np.where(np.abs(dO) > np.pi, (om[None, :] - om[:, None]) - 2 * np.arcsin(t), P)
    D = np.sqrt(
        (e[None, :] - e[:, None]) ** 2
        + (q[None, :] - q[:, None]) ** 2
        + (2 * np.sin(I / 2)) ** 2
        + ((e[:, None] + e[None, :]) / 2) ** 2 * (2 * np.sin(P / 2)) ** 2
    )
    return D


def _count_below(D: np.ndarray, thr: float) -> int:
    iu = np.triu_indices(D.shape[0], k=1)
    return int((D[iu] < thr).sum())


def find_pairs(
    objs: list[Orb], threshold: float = 0.05, *, exclude_same_parent: bool = True
) -> list[dict]:
    """Candidate genetic pairs (D below threshold), closest first. By default drops
    same-parent fragments (e.g. split-comet pieces), which are a known -- not novel
    -- source of identical orbits, and tags sequential-designation pairs (likely
    same-survey discovery, a selection artifact rather than a confirmed pair)."""
    D = dsh_matrix(objs)
    out = []
    iu, ju = np.triu_indices(len(objs), k=1)
    for a, b in zip(iu, ju):
        if D[a, b] >= threshold:
            continue
        na, nb = objs[a].name, objs[b].name
        if exclude_same_parent and same_parent(na, nb):
            continue
        seq = _base_designation(na)[:8] == _base_designation(nb)[:8] and _base_designation(
            na
        ) != _base_designation(nb)
        out.append(
            {
                "D": round(float(D[a, b]), 4),
                "a": na,
                "b": nb,
                "orb_a": objs[a],
                "orb_b": objs[b],
                "same_survey_suspect": bool(seq),
            }
        )
    out.sort(key=lambda r: r["D"])
    return out


def pair_significance(
    objs: list[Orb], *, threshold: float = 0.05, n_null: int = 200, seed: int = 0
) -> dict:
    """Is the observed close-pair count an EXCESS over chance? Null shuffles the
    orientation angles (Omega, omega) independently, keeping the a/e/i/q marginals,
    so genuine apsidal/nodal coincidences are destroyed but the population's natural
    clumpiness is preserved."""
    obs = _count_below(dsh_matrix(objs), threshold)
    rng = np.random.default_rng(seed)
    Om = np.array([o.Omega for o in objs])
    om = np.array([o.omega for o in objs])
    null = np.empty(n_null)
    base = [Orb(o.name, o.a, o.e, o.inc, 0, 0, o.q) for o in objs]
    for k in range(n_null):
        po, pw = rng.permutation(Om), rng.permutation(om)
        for j, o in enumerate(base):
            o.Omega = po[j]
            o.omega = pw[j]
        null[k] = _count_below(dsh_matrix(base), threshold)
    mean, std = float(null.mean()), float(null.std() or 1e-9)
    return {
        "observed_pairs": obs,
        "null_mean": round(mean, 1),
        "null_std": round(std, 2),
        "excess_sigma": round((obs - mean) / std, 2),
        "p_value": round(float(np.mean(null >= obs)), 4),
    }


# --------------------------------------------------------------------------- #
#  exotic-origin scoring (captured interstellar / Oort-return candidates)
# --------------------------------------------------------------------------- #
def exotic_score(o: Orb) -> dict:
    """Flag and score an object for a non-disk origin. Retrograde and polar orbits
    cannot form in the prograde protoplanetary disk; co-orbital + retrograde is the
    captured-interstellar signature (cf. 2015 BZ509)."""
    flags = []
    score = 0.0
    if o.inc > 90:
        flags.append("RETROGRADE")
        score += (o.inc - 90) / 90 * 2
    if 80 < o.inc < 100:
        flags.append("POLAR")
        score += 1.0
    if o.e >= 1.0:
        flags.append("UNBOUND/ISO")
        score += 3.0
    if o.inc > 60:
        for pn, pa in PLANET_A.items():
            if abs(o.a - pa) / pa < 0.03:
                flags.append(f"CO-ORBITAL/{pn}")
                score += 2.0
    return {
        "name": o.name,
        "score": round(score, 2),
        "flags": flags,
        "a": o.a,
        "e": o.e,
        "inc": o.inc,
    }


def hunt(
    objs: list[Orb], *, label: str = "", pair_threshold: float = 0.05, n_null: int = 200
) -> dict:
    """Full exotic-orbit hunt on a population: candidate pairs + their significance,
    and the ranked exotic-origin candidates."""
    pairs = find_pairs(objs, pair_threshold)
    sig = pair_significance(objs, threshold=pair_threshold, n_null=n_null)
    exotic = sorted((exotic_score(o) for o in objs), key=lambda r: -r["score"])
    exotic = [e for e in exotic if e["flags"]]
    return {
        "label": label,
        "n": len(objs),
        "pairs": pairs[:12],
        "pair_significance": sig,
        "n_exotic": len(exotic),
        "top_exotic": exotic[:12],
    }
