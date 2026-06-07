# Characterization Engine — "spot → analyze → identify"

A unified engine that takes any change-detection — a **moving** object or a
**varying/stationary** source — and returns a structured verdict: a probability
distribution over the object taxonomy, the derived physical properties with
**calibrated uncertainty**, an honest confidence, and the single most useful
next observation. Built 2026-06-02 on top of the validated detection +
sub-arcsec ephemeris foundation.

## Modules (`src/ariadne/discovery/imaging/`)

| Module | Role |
|---|---|
| `characterize.py` | The spine. `characterize_mover()`, `characterize_variable()`, `characterize(candidate)` dispatcher → `ObjectDossier`. |
| `color.py` | Multi-band color → taxonomy (Ivezić a* for S/C, i−z for V-type, red-vs-neutral for TNOs). |
| `light_curve.py` | Lomb-Scargle period + false-alarm prob + folded shape → variable-star typing; rotation period for movers. |
| `transient_detection.py` | Engine-2 spotter: difference-image residuals clustered across epochs → variable / transient / slow-mover, with a fractional-change real/bogus filter. |
| `classifier.py` | `identify(candidate)` — top-level fusion → composite label + ranked posterior + confidence + evidence + next step. |

## What it produces

**Mover** → dynamical class (from rate→distance or a fitted orbit), surface type
(from color), activity (coma/tail from morphology), size with a 16–84% CI
(Monte-Carlo-propagated from the distance uncertainty), known/new, and — if the
orbit is hyperbolic — an interstellar-object flag.

**Variable** → period (or "aperiodic"), amplitude, folded shape, and a typed
posterior (eclipsing binary / RR Lyrae / Cepheid / Mira / microlensing / nova /
transient), each with the next observation that would disambiguate it.

The design rule: **never assert a hard label from thin data.** Confidence
reflects how much the data actually constrains the answer; `next_step` says how
to sharpen it (e.g. "a second night yields a real orbit," "a g-band frame gives
the taxonomy").

## Honest status (read this)

These are **correct, tested, literature-grounded first-pass implementations** —
the right architecture, real algorithms (not stubs), 20 passing tests including
period recovery, color classification, and variable/transient separation from a
star field. They are **not yet** validated against real catalogs and machine-learning benchmarks.
The remaining road, stated plainly:

1. **Validate against real labeled catalogs**, not just synthetic + literature
   thresholds: color taxonomy vs a real asteroid spectral-type set; variable
   typing vs VSX/AAVSO; the real/bogus filter vs real difference-image stamps.
2. **Replace heuristics with trained classifiers** where production pipelines do
   (real/bogus CNN; light-curve type RNN/RF).
3. **Make the fusion a rigorous Bayesian/NPE posterior** with full uncertainty
   propagation (currently: structured rules + MC size-CI + erf-integrated
   distance class-probs — calibrated on the axes that are validated, not yet a
   single learned posterior).

Each upgrade is a staged, validate-against-truth effort — the same discipline
that took detection recall from 31% → 91%. The foundation is built to extend
into all of them without changing callers.
