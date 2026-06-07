# Ariadne documentation

Start with the top-level [README](../README.md) for install and a 30-second tour, and the
rendered [Sphinx docs](sphinx/) (also on ReadTheDocs). The markdown documents here are
grouped by audience.

## Start here (users)
- [`WHITE_PAPER.md`](WHITE_PAPER.md) — capstone write-up: methods, full validation table, results, honest limitations
- [`FRONTIER_CAPABILITIES.md`](FRONTIER_CAPABILITIES.md) — the outer-solar-system / exoplanet / anomaly toolkit
- [`SOLAR_SYSTEM_NAVIGATOR.md`](SOLAR_SYSTEM_NAVIGATOR.md) — route search across Lambert / flyby / moon-tour engines
- [`CISLUNAR_MISSION_ARCHITECT.md`](CISLUNAR_MISSION_ARCHITECT.md) — cislunar mission-architecture ranking
- [`LSST_READINESS.md`](LSST_READINESS.md) — how the pipeline maps onto the Rubin/LSST stream
- [`OPERATIONS.md`](OPERATIONS.md) — running the discovery tooling

## Validation & results
- [`REVIEWER_GUIDE.md`](REVIEWER_GUIDE.md) - claim-to-evidence map and quick verdict path for skeptical review
- [`VALIDATION_RESULTS.md`](VALIDATION_RESULTS.md) — cross-validation against GMAT, REBOUND, independent libraries
- [`REAL_BENCHMARKS.md`](REAL_BENCHMARKS.md) — benchmarks on real data
- [`SENSITIVITY.md`](SENSITIVITY.md) — sensitivity / parameter studies
- [`COHERENCE_REAL_DATA_SCORECARD.md`](COHERENCE_REAL_DATA_SCORECARD.md) — coherence-selector A/B results on real data
- [`CAPABILITY_SUMMARY.md`](CAPABILITY_SUMMARY.md) — what the engine is and is not
- [`FRONTIER_FINDINGS.md`](FRONTIER_FINDINGS.md) — the honest Planet Nine / anomaly-hunt findings

## Internal / development notes
These are working notes kept for transparency and reproducibility, not polished user guides.
- [`PUBLIC_LAUNCH_AUDIT.md`](PUBLIC_LAUNCH_AUDIT.md) — public-release readiness audit
- [`PRODUCTION_DISCOVERY_PIPELINE.md`](PRODUCTION_DISCOVERY_PIPELINE.md) — the nightly pipeline internals
- [`OPERATIONAL_DISCOVERY.md`](OPERATIONAL_DISCOVERY.md) — operator path for building a labeled corpus
- [`CHARACTERIZATION_ENGINE.md`](CHARACTERIZATION_ENGINE.md) — object-characterization engine notes

## Design reference
- [`../MASTER_PLAN.md`](../MASTER_PLAN.md) — the full design bible: vision, prior art, CR3BP foundations, the 45-stage validation roadmap. Source-file docstrings reference its sections (for example `MASTER_PLAN.md §3.10`).
