# Frontier capabilities: where the coherence selector has a real, honest shot

Three capabilities aimed beyond plain solar-system moving-object detection, each at
a place where the bottleneck is **selection / vetting / anomaly** (the regime the
Equation-of-One engine actually wins at) rather than raw detection (mature tools own
that), and each on **free public data**. The honest edge throughout: a coherence
energy needs no training set, so it ranks candidates by agreement with *known
physics* -- the low-data / novel / anomaly regime where supervised ML starves.

These are **triage and ranking tools**, not discovery claims. None of them confirms
anything; confirmation needs follow-up we do not have. They decide what is *worth*
that follow-up. What they are validated to do is below, with the honest limits.

---

## 1. Extreme-TNO orbital-orientation outlier hunt  (`frontier.tno_clustering`)

The distant solar system's real, debated anomaly: the most extreme trans-Neptunian
objects appear to have orbits clustered in orientation, which -- if real and not an
observational-selection artifact -- points to a distant perturber ("Planet Nine").

**Data:** all 6,043 TNOs with orbital elements from JPL SBDB (a small free table).

**What it does:** reproduces the clustering signal from scratch, scores each orbit's
circular coherence energy `E_c = 1 - cos(angle - mean)` (the von Mises form of the
S_One coherence term), and ranks the least-coherent (most anomalous) orbits.

**Validated result (real JPL data):**

| a-cutoff | n  | Omega R | Rayleigh p | significant? |
|----------|----|---------|------------|--------------|
| a >= 150 | 74 | 0.310   | 0.0007     | yes          |
| a >= 250 | 37 | 0.218   | 0.17       | no           |
| a >= 400 | 16 | 0.182   | 0.60       | no           |

This faithfully reproduces the **live scientific debate**: the original 2016 claim
used a >= 250 AU on a tiny sample; with today's larger sample that "clean" detached
population is consistent with uniform. The tool reports the cutoff-dependence
transparently instead of cherry-picking. Famous oddballs (Leleakuhonua / "The
Goblin", 2017 OF201) fall out correctly as the most orientation-anomalous orbits.

**Honest limit:** this reproduces known science and adds a coherence-based scoring /
outlier tool. It does not discover Planet Nine (that needs deep imaging we lack), and
the signal is sensitive to observational selection effects we do not model.

---

## 2. TESS transit search + coherence vetting  (`frontier.tess_vetting`)

Finding a periodic dip is commodity (Box Least Squares). The bottleneck is deciding
whether it is a planet or an impostor (eclipsing binary, blend, systematic). A real
planet's transit is a tightly constrained physical object -- a coherence basin; we
score the S_One `E_total` (a SUM of physical penalty sectors, not a mean, so one
unambiguous disqualifier like a 9% depth is not diluted) and label it.

**Data:** free TESS light curves from MAST (via `lightkurve`).

**Validated -- synthetic separation (8/8):** real-planet shapes pass (coherence
0.98-1.00); all eclipsing-binary impostors rejected (coherence 0.00-0.16): deep
eclipses, secondary eclipses, odd-even depth differences, and grazing V-shapes,
including when BLS locks onto the half-period alias.

**Validated -- real TESS data:**

| target | true | recovered period | recovered depth | coherence | verdict |
|--------|------|------------------|-----------------|-----------|---------|
| Pi Mensae c | P=6.27d, ~0.031% | **6.271d** | **0.030%** | 1.000 | planet-candidate |
| WASP-18 b   | P=0.941d, ~0.9%  | **0.942d** | 1.01%      | 0.912 | planet-candidate |

A real, very shallow confirmed super-Earth and a hot Jupiter, both recovered at the
correct period and depth and vetted correctly.

**Honest limit:** surfaces and ranks candidates; does not confirm planets (needs
radial-velocity / imaging follow-up). The basin is a first-pass physical prior.

---

## 3. ZTF/LSST physics-coherence novelty triage  (`frontier.ztf_anomaly`)

LSST will issue ~10 million alerts a night. The open problem is triage: which few
objects are novel enough to deserve a telescope. Trained classifiers push a
never-seen phenomenon into whichever known bin fits least badly; a coherence energy
inverts that -- novelty is incoherence with the **best-matching** known class, so
something that coheres with nothing floats to the top.

**Data:** free public ZTF light curves via the ALeRCE broker (no auth).

**Calibration (honest):** the physics sets which axes matter (period, amplitude,
periodicity strength, asymmetry) and the class topology; the *scales* are calibrated
to a real labeled ALeRCE sample (`scripts/calibrate_ztf.py`, 6 classes x ~20 objects,
robust median + MAD), because textbook amplitudes (~0.8 mag) do not match real ZTF
amplitudes (~1.9 mag for RRab). The calibrated basins ship committed.

**Validated both directions:**
- **Real known objects recognized as known: 14/15** with the calibrated basins
  (RR Lyrae at real ~1.9-mag amplitudes, LPVs at P~195d, QSOs low-amplitude).
- **Genuine novelties still flagged anomalous** (score 7.9-12.9 vs threshold 3.0):
  physically impossible combinations (4-mag amplitude at an RR-Lyrae period;
  3-mag-amplitude strong periodicity at 0.08 d) -- calibration did not become
  permissive enough to absorb them. The synthetic-mechanism test runs against the
  physics priors deterministically.

**Honest limit:** a triage / ranking score, not a classifier or a discovery. The
calibration validation above is in-sample (top-probability objects); a fully
held-out evaluation is the next step. High score means "worth a human's attention,"
nothing more.

---

## Running them

```
python scripts/hunt_tno_outliers.py
python scripts/hunt_tess_candidates.py "Pi Mensae" "WASP-18"
python scripts/calibrate_ztf.py --per-class 20      # (re)fit the ZTF basins
python scripts/triage_ztf_anomalies.py --alerce-class SNIa --n 20
```

Tests (network-free, all synthetic/offline): `tests/test_tno_clustering.py`,
`tests/test_tess_vetting.py`, `tests/test_ztf_anomaly.py` -- 27 tests.
