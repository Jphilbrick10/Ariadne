# Public Launch Audit

Audit date: 2026-06-05
Target version: 1.0.0rc2

This is the reviewer-facing launch checklist for Ariadne. It records what was
actually run before public release and keeps the scope scientifically honest.
For a claim-by-claim evidence map, see [`REVIEWER_GUIDE.md`](REVIEWER_GUIDE.md).

## Verified Gates

| Gate | Command | Result |
|---|---|---|
| One-command reviewer quickcheck | `python scripts\reviewer_quickcheck.py --report results\reviewer_quickcheck_report.json` | passed; lint, fast tests, Sphinx docs, closure ledger |
| Frozen reviewer evidence manifest | `python scripts\build_reviewer_evidence_manifest.py` | complete; 21 hashed evidence files |
| Fast offline test suite | `pytest -m "not slow" -q` | 1050 passed, 5 skipped, 30 deselected |
| Production lint sanity | `ruff check src\ariadne --select E9,F63,F7,F82,B023,B904` | passed |
| Closure proof ledger | `python scripts/build_closure_report.py --fail-on-critical` | status complete, readiness 1.000000, critical failures 0 |
| Stage 32 backend selector validation | `$env:PYTHONPATH='src'; python -m ariadne.validate.stage32` | all gates pass; selector matched measured parallel/GPU winners |
| Package build | `python -m build` | wheel and sdist built |
| Distribution metadata | `python -m twine check dist\*` | wheel passed, sdist passed |
| Wheel install smoke | `pip install --force-reinstall --no-deps dist\ariadne_astro-1.0.0rc2-py3-none-any.whl`; import/system smoke | passed |
| Sphinx docs | `sphinx-build -b html -W --keep-going docs\sphinx docs\sphinx\_build\html` | passed |
| Focused IOD/discovery docs regression | `pytest tests\test_image_pipeline_tier1.py tests\test_discovery_coverage_boost.py tests\test_discovery_complete.py -q` | 126 passed |

## Hygiene Checks

- No tracked file larger than 1 MiB was found with `git ls-files`.
- Raw survey imagery, SPICE kernels, local databases, build products, root logs,
  and frontier/anomaly scan caches are ignored.
- The secret scan found only test placeholders and environment-variable based
  configuration examples; no real token was detected.
- README/docs are UTF-8 clean: no replacement characters were detected.
- `docs/sphinx/_build` is ignored and not tracked.

## Scientific Scope

Ariadne is a research toolkit, not an operational mission-planning authority or
an MPC discovery service. The validation artifacts should be read as evidence
for reproducible algorithms and benchmark behavior, not as a claim that every
mission design or survey-discovery edge case has been certified for operations.

The strongest public claims are backed by tracked artifacts in:

- `data/benchmarks/closure/`
- `data/benchmarks/engine_improvement/`
- `data/benchmarks/external_inference/`
- `data/benchmarks/real_decam_discovery/`
- `data/benchmarks/real_corpus_alerce_ztf_asteroid_2026/`
- `docs/VALIDATION_RESULTS.md`
- `docs/COHERENCE_REAL_DATA_SCORECARD.md`

## Non-Gates To Avoid Overclaiming

- Full-repository `ruff check .` is not currently a release gate. The release
  does gate on fatal Python diagnostics, loop-closure correctness, and exception
  chaining in production source via the `Production lint sanity` check above.
- Network-backed benchmarks depend on live third-party services and catalog
  versions. Treat stored benchmark artifacts as dated evidence, and rerun the
  acquisition scripts when making new scientific claims.
- The frozen ALeRCE/ZTF alert corpus is a public broker-classified replay
  corpus. It improves external-data provenance, but it is not a private
  NOIRLab/LSST labelled benchmark.
- Commercial use is not open-source licensed; the project is source-available
  under PolyForm Noncommercial plus a separate commercial license path.

## Recommended Release Procedure

1. Start from a clean worktree.
2. Run the verified gates above.
3. Regenerate figures or benchmark artifacts only when their source scripts and
   provenance are also updated.
4. Tag the release after the closure report and `CITATION.cff` version agree.
5. Attach the built wheel/sdist or publish from the exact tagged commit.
