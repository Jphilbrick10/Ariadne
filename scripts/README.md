# `scripts/` - runnable tools and reproducible pipelines

These are operator-facing scripts: discovery runs, calibration, validation, and
benchmark/proof generation. They are not part of the importable `ariadne` package
(that lives in [`src/ariadne/`](../src/ariadne/)); they orchestrate it.

**Run pattern:** `PYTHONPATH=src python scripts/<name>.py [args]`. Every script has a
module docstring at the top explaining what it does and what data or credentials it
needs. Scripts that touch external archives (MPC, ALeRCE, MAST, NOIRLab) require
network access and, for NOIRLab Data Lab, `DATALAB_USER` / `DATALAB_PASS`.

## Mission design & navigation
- `architect_cislunar_mission.py` - rank Earth-Moon mission architectures
- `navigate_solar_system.py` - search/rank solar-system routes (Lambert, flyby, moon tour)
- `build_solar_transfer_atlas.py` - whole-solar-system Lambert corridor atlas
- `map_solar_system_3d.py` - interactive 3D solar-system map (plotly HTML)
- `render_asteroid_frame.py` - render a mission/asteroid plate

## TNO & outer-solar-system frontier
- `hunt_tno_outliers.py` - extreme-TNO orbital-clustering outlier hunt
- `hunt_tess_candidates.py` - TESS transit search + coherence vetting
- `hunt_ztf_anomalies.py` / `triage_ztf_anomalies.py` - ZTF/ALeRCE novelty triage
- `run_p9_dynamical_test.py` - Planet Nine secular/dynamical test
- `gaia_encounters.py` - Gaia close-stellar-encounter search
- `scan_everything.py`, `scan_faint_movers.py`, `scan_fast_movers.py` - broad anomaly sweeps

## Co-orbital & pair analysis
- `check_coorbital.py`, `coorbital_sweep.py` - 1:1 resonance screening / census
- `coorbital_rebound.py` - REBOUND clone-ensemble a-confinement (use IAS15 for high-e crossers)
- `coorbital_resonant_angle.py` - resonant-angle libration test (validated on Achilles to L4)
- `confirm_pair.py` - backward N-body asteroid-pair confirm/refute

## Survey discovery pipeline (nightly / DECam / NSC / DES)
- `run_full_discovery_night.py` - end-to-end nightly discovery run
- `run_auto_discovery.py`, `run_catalog_survey.py` - automated survey drivers
- `run_nsc_discovery.py`, `run_des_discovery.py` - NOIRLab Source Catalog / DES runs
- `run_decam_e2e.py`, `run_real_decam_e2e.py`, `run_real_decam_phase_a.py` - DECam end-to-end
- `run_coherence_field_discovery.py`, `run_full_image_pipeline.py` - pipeline variants
- `process_decam_night.py`, `search_gap_field.py` - per-night processing / under-documented sky
- `fetch_discovery_field.py`, `fetch_nsc_field.py`, `find_discovery_field.py`, `extract_field_to_npz.py` - field acquisition

## Calibration
- `calibrate_confidence.py` - posterior confidence (temperature/ECE) on real orbits
- `calibrate_mover.py`, `calibrate_vet.py`, `calibrate_ztf.py` - selector/basin calibration
- `fit_calibration_large_corpus.py`, `train_neural_orbit_prior.py` - large-corpus fits

## Validation (real data)
- `validate_mover_real.py`, `validate_coherence_vet.py`, `validate_coherence_vet_real.py` - selector A/B vs baselines
- `validate_lightcurve.py` - light-curve analyzer vs labeled VSX/ZTF
- `validate_repeated.py` - repeated draws with confidence intervals
- `validate_robustness.py` - graceful-degradation stress test
- `validate_chain_quality_coherence.py` - chain-quality checks
- `run_real_recovery.py`, `run_real_backtest.py`, `run_real_benchmarks.py` - real-data recovery/backtests
- `fetch_and_validate_jpl.py`, `acquire_real_labelled_corpus.py` - ground-truth/corpus acquisition

## Benchmarks & proof artifacts
- `reviewer_quickcheck.py` - one-command local reviewer audit; writes a JSON evidence report
- `build_reviewer_evidence_manifest.py` - frozen hash manifest for public review evidence
- `run_discovery_benchmark.py`, `benchmark_engine_improvement.py`, `benchmark_image_pipeline.py` - benchmark drivers
- `benchmark_solar_navigator.py`, `benchmark_transport_admissibility.py` - navigation/transport benchmarks
- `build_closure_report.py`, `build_dream_run.py`, `build_artifact_manifest.py`, `build_external_inference_benchmark.py` - proof/closure artifacts
- `compare_replay_manifests.py`, `compare_solar_navigator_benchmarks.py` - manifest/benchmark diffs
- `promote_navigator_routes.py`, `bootstrap_scheduler_from_benchmark.py` - promotion/scheduling

## Imaging & detection diagnostics
- `diagnose_2015_field_knowns.py`, `diagnose_detection_fix.py`, `diagnose_match_distance.py`, `diagnose_pixel_truth.py`, `diagnose_recall_vs_mag.py` - detection diagnostics
- `test_difference_completeness.py`, `test_psf_matched_completeness.py`, `test_rate_floor.py`, `test_single_snapshot_rate.py` - completeness/rate studies
- `measure_velocity_coherent_linking.py` - coherent-linking velocity study

## Console & tooling
- `dashboard.py` - local discovery console (launch via `start_dashboard.bat` on Windows)
- `cockpit_worker.py` - background worker for the console
- `fix_notebook_ids.py` - notebook metadata maintenance
