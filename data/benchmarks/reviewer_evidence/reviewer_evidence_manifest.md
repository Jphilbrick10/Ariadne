# Reviewer Evidence Manifest

- status: `complete`
- file count: `21`
- certificate hash: `061a00885093742bd47734c6d8c58ad3d8093c5551b9e8c30a451b5031edd501`

## Reproduction Commands

| Name | Command | Purpose |
|---|---|---|
| `reviewer_quickcheck` | `PYTHONPATH=src python scripts/reviewer_quickcheck.py` | Run production lint, fast tests, strict docs, and closure ledger. |
| `closure_ledger` | `python scripts/build_closure_report.py --fail-on-critical` | Regenerate closure evidence and fail on critical residuals. |
| `stage24_packaging` | `PYTHONPATH=src python -m ariadne.validate.stage24` | Verify installability, license metadata, CI, and core imports. |
| `stage32_selector` | `PYTHONPATH=src python -m ariadne.validate.stage32` | Verify ensemble backend faithfulness and selector calibration. |

## Frozen Evidence Files

| Role | Path | Bytes | SHA-256 |
|---|---|---:|---|
| `public_entrypoint` | `README.md` | 13754 | `c3ebe3c53312736c8a4d8bda99d3a5b3a59deb4a9ea994f0ad3f065a19939e19` |
| `review_map` | `docs/REVIEWER_GUIDE.md` | 6862 | `a218459113eef11ffdda58d7f524e813882cc7cd9990c8915dbdba56c1247363` |
| `release_audit` | `docs/PUBLIC_LAUNCH_AUDIT.md` | 4097 | `35d223d7b08176426b974d3c504046cbe44d359cdd91b59e9c21cb50f3592261` |
| `scientific_narrative` | `docs/WHITE_PAPER.md` | 13621 | `a0bba8ffb26bce5aa76d029ab5173ebcbab635c4ae9640ddd6f93a0e73bc0360` |
| `validation_report` | `docs/VALIDATION_RESULTS.md` | 46459 | `654f8bb1fc7bb5ead33474fe8b12db18160727950d776d41503dc3937300d652` |
| `real_data_report` | `docs/REAL_BENCHMARKS.md` | 2886 | `310a34eb7b4c4e3db2cd09feaea7f676260a180c6955da755ff968eedee3a009` |
| `closure_ledger` | `data/benchmarks/closure/closure_report.json` | 17038 | `52b637fbce22fc5161489aec49c36923c22516607c366b698e1ca664593434d5` |
| `closure_summary` | `data/benchmarks/closure/closure_report.md` | 3113 | `f6351ad0f909c3a82b7a76ef2d00581f997eac3cd0c33b7b908b15d9901057c2` |
| `frozen_labelled_corpus_manifest` | `data/benchmarks/real_corpus_mpc_500/corpus_manifest.json` | 2510 | `c7076711bb43d708f012d03bee2d72593279bd9f89c36385c96791ae67c29803` |
| `frozen_labelled_cases` | `data/benchmarks/real_corpus_mpc_500/labelled_cases.jsonl` | 362093 | `7830bd7b19d0d78e1aae21ec5805b18fb0b1f69bc302a589b0ac94ddc5bc6a97` |
| `real_corpus_metrics` | `data/benchmarks/real_corpus_mpc_500/benchmark/metrics.json` | 163445 | `067d59dc51d06dac0058ed28195955d5a7a225ab93b4e60acf8e2994a194c15a` |
| `real_corpus_reliability` | `data/benchmarks/real_corpus_mpc_500/benchmark/reliability_curve.csv` | 69 | `925e913a8598741d555eb8864e60969872590ca90c3d3f6f3c3e8414f2a1c33b` |
| `external_ztf_alert_corpus_manifest` | `data/benchmarks/real_corpus_alerce_ztf_asteroid_2026/corpus_manifest.json` | 1773 | `2034e950ad8190324009fc9e7f6da8679caeb6a6211499c9e10b0f87f8de61e0` |
| `external_ztf_alerts` | `data/benchmarks/real_corpus_alerce_ztf_asteroid_2026/alerts.jsonl` | 247577 | `453af4622f357167315b03f40e55a1da1aa3d74734aebd181c284a86ca21bbd7` |
| `external_ztf_replay_manifest` | `data/benchmarks/real_corpus_alerce_ztf_asteroid_2026/replay_manifest.json` | 425 | `84518aa922f80814a20e002fdbc61b01f5f63154f0b113e9742619d5f11cdd2f` |
| `external_ztf_replay_provenance` | `data/benchmarks/real_corpus_alerce_ztf_asteroid_2026/provenance.jsonl` | 170401 | `98a275edde8f01855ad8aa677fa827aad55031936a948adefee7958632056ef3` |
| `external_inference_metrics` | `data/benchmarks/external_inference/metrics.json` | 753918 | `4636c503a52365406c5a7697be9e5d7c5df08f31b313c24a69db2ecf5f2e36dd` |
| `adversarial_metrics` | `data/benchmarks/inference_adversarial_builtin/metrics.json` | 84725 | `b5eb626b97c86b5feb466f1b570d26bea6e5806c3590c061f5e04c6fcaba4b2e` |
| `artifact_integrity_manifest` | `data/benchmarks/artifact_integrity/artifact_manifest.json` | 14300 | `d333646b74e10414cc6951ce18fb6b49a46d1092144cbc1c531ed06caaeded03` |
| `citation_metadata` | `CITATION.cff` | 548 | `e079661d8e9b27754efd6b4df0d57e761741527ade68a9be6eeb38b0bc0751bd` |
| `license` | `LICENSE` | 4672 | `772505897bc3cc6f375d69d309d488577b90741d1dbb91e566b026bfaf0fd912` |

## Limits

- Source-available under PolyForm Noncommercial, not OSI open-source.
- Frozen MPC corpus is useful evidence, not a blind NOIRLab/ZTF/LSST competition.
- Frozen ALeRCE/ZTF corpus is a public broker-classified alert replay corpus, not a ground-truth labelled LSST/NOIRLab benchmark.
- Network-backed external benchmarks should be rerun before new scientific claims.
- Operational mission design still requires certified workflows such as GMAT/Monte/Copernicus.
