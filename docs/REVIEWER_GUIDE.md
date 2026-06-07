# Reviewer Guide

This guide is for scientists, reviewers, and maintainers who want to decide
whether Ariadne's public claims are supported by runnable evidence. It is not a
marketing page; it is a map from claims to commands, artifacts, and known
limits.

## Quick Verdict Path

From a fresh checkout:

```bash
python -m pip install --upgrade pip
pip install -e ".[dev,docs]"
PYTHONPATH=src python scripts/reviewer_quickcheck.py
```

The script writes `results/reviewer_quickcheck_report.json` with per-gate
commands, status, durations, and output tails. The equivalent manual commands
are:

```bash
ruff check src/ariadne --select E9,F63,F7,F82,B023,B904
pytest -m "not slow" -q
sphinx-build -b html -W --keep-going docs/sphinx docs/sphinx/_build/html
python scripts/build_closure_report.py --fail-on-critical
```

Expected release evidence for `1.0.0rc2`:

- fast suite: 1048 passed, 5 skipped, 30 deselected
- production lint sanity: passed
- Sphinx docs with warnings as errors: passed
- closure ledger: status `complete`, readiness `1.000000`, critical failures `0`

The exact launch checklist is tracked in
[`PUBLIC_LAUNCH_AUDIT.md`](PUBLIC_LAUNCH_AUDIT.md).

## Claim-To-Evidence Map

| Public claim | Primary evidence | Reproducibility path | What it does not prove |
|---|---|---|---|
| Ariadne is an installable Python research package | `pyproject.toml`, `src/ariadne/`, CI workflow | `python -m build`; `python -m twine check dist/*`; wheel install smoke | Does not imply operational certification |
| CR3BP/Lagrange/orbit-family machinery is numerically stable | `src/ariadne/validate/stage*.py`, tests, docs tables | `pytest -m "not slow" -q`; targeted stage scripts under `src/ariadne/validate/` | Does not replace full mission-ops tools |
| Dynamics are cross-checked against external references | `docs/WHITE_PAPER.md`, `docs/VALIDATION_RESULTS.md`, `data/benchmarks/closure/` | Review validation docs and rerun stage scripts with required optional tools/data | GMAT availability and external tool versions can change |
| Discovery pipeline can recover known TNO-style cases | `examples/04_tno_orbit_fit.py`, discovery tests, `data/benchmarks/real_corpus_mpc_500/` | `ariadne discover 90377`; discovery-focused pytest slices | Does not establish new-object discovery from a blind live survey |
| Inference benchmarks include adversarial and calibration checks | `data/benchmarks/external_inference/`, `data/benchmarks/inference_adversarial_builtin/` | Inspect metrics, reliability outputs, and benchmark scripts | Built-in CI proxies are not a substitute for a large labelled alert stream |
| Real-data discovery results are reported honestly | `docs/VALIDATION_RESULTS.md`, `docs/REAL_BENCHMARKS.md`, `docs/FRONTIER_FINDINGS.md` | Compare stored benchmark artifacts to the documented null findings | Null results do not prove no objects exist in deeper or unqueried surveys |
| Solar/interplanetary navigation tools are exploratory research engines | `docs/SOLAR_SYSTEM_NAVIGATOR.md`, `data/benchmarks/solar_navigator_benchmark/`, `src/ariadne/interplanetary/` | Run examples and benchmark scripts for specific scenarios | Does not certify globally optimal mission design across all constraints |
| Public release hygiene has been checked | `docs/PUBLIC_LAUNCH_AUDIT.md`, `CITATION.cff`, `SECURITY.md`, `LICENSE` | Re-run the audit commands and inspect git-tracked files | Full-repo cosmetic lint remains a known non-gate |

## Deep Review Path

A deeper review should separate three levels of evidence:

1. **Unit and integration behavior**: run the full fast test suite and inspect
   `tests/` for edge cases and adversarial coverage.
2. **Numerical and scientific validation**: run selected validation stages from
   `src/ariadne/validate/`, then compare results with `docs/WHITE_PAPER.md`,
   `docs/VALIDATION_RESULTS.md`, and stored benchmark artifacts.
3. **External replication**: rerun cross-checks that depend on external tools or
   live catalogs, documenting exact tool versions, catalog epochs, and network
   acquisition dates.

Recommended focused commands:

```bash
pytest tests/test_integrators.py tests/test_packaging.py -q
pytest tests/test_external_corpora.py tests/test_discovery_complete.py -q
PYTHONPATH=src python scripts/reviewer_quickcheck.py --profile smoke
PYTHONPATH=src python -m ariadne.validate.stage24
PYTHONPATH=src python -m ariadne.validate.stage32
```

On Windows PowerShell, use:

```powershell
$env:PYTHONPATH = "src"
python -m ariadne.validate.stage24
python -m ariadne.validate.stage32
```

## Known Reviewer Attack Points

These are legitimate questions to ask during review:

- The package is **source-available under PolyForm Noncommercial**, not OSI
  open-source. Public text should use that wording.
- The repository has a full-lint backlog outside the production correctness
  gate. The gated lint subset catches fatal Python errors, loop-closure bugs,
  and exception-chain issues in production source.
- Some benchmark artifacts are dated. Network-backed results should be rerun
  before making new scientific claims.
- Built-in inference benchmarks are serious but not equivalent to a large blind
  ZTF/LSST/NOIRLab labelled-corpus competition.
- Operational mission design still belongs in GMAT, Monte, Copernicus, or a
  comparable certified workflow.

## What Would Raise Confidence Further

The strongest next evidence would be:

- a frozen, citable labelled discovery corpus with exact acquisition scripts
- a paper-style reproduction notebook for each headline table
- independent reruns by a third party
- archived GMAT/JPL/cross-tool inputs and outputs with version pins
- a CI job for the reviewer guide's quick verdict path

Those are evidence upgrades, not substitutes for the existing tests.
