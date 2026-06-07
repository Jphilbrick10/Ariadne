# Rubin / LSST readiness for the frontier toolkit

This session repeatedly reached the same wall: re-mining existing catalogs and the
poorly-cadenced archival gaps is exhausted. The genuine discovery frontier is the
Vera Rubin Observatory's LSST -- now coming online, surveying the entire southern sky
(including the galactic plane the old surveys avoided) with, for the first time, a
cadence built for moving-object and transient detection. Our validated, skeptical
pipeline is exactly the triage + vetting layer that flood will need. This documents
how each tool maps onto LSST data, and the concrete steps to switch it on.

Honest status: LSST is not yet delivering a public science stream we can pull, so
this is an integration *plan*, not a running adapter. What makes it credible is that
**most of the tools already run on the LSST precursor (ZTF) through the same brokers
LSST will use** -- so the code path is exercised today and extends with a source swap.

## The LSST data products and which tool consumes each

| LSST product | What it is | Our tool |
|---|---|---|
| **Alerts (DIASources)** | ~10 million difference-image detections per night, served via community brokers (ALeRCE, Lasair, Fink, ANTARES) | `ztf_anomaly` -- physics-coherence novelty triage. **Already runs on ZTF via ALeRCE**; ALeRCE will ingest LSST, so the same `fetch -> features -> score` path applies with the broker source swapped. |
| **SSObjects / SSSource** | the linked moving-object catalogue (tracklets joined into orbits) | `exotic_orbit_hunt` (D-criterion families, capture candidates), `planet_nine`/`tno_clustering` (clustering, plane warp, bias nulls) -- run directly on the orbital-element table, as we do now on JPL SBDB. |
| **Nightly visits (calexp/diffim)** | the cadenced images of each field | the moving-object arsenal (`scan_faint_movers` GPU shift-stack, `scan_fast_movers` streak, the rate-linker + coherence vetting) + `search_gap_field` -- the cadence gate we added is satisfied by LSST's design (repeat visits over nights). |
| **Forced photometry / light curves** | long light curves per object | `tess_vetting`-style transit search is TESS-specific, but `ztf_anomaly`'s feature + novelty scoring applies to any light curve. |

## Why the precursor work transfers

- **Brokers are the abstraction layer.** ALeRCE already serves ZTF and is a designated
  LSST broker. `ztf_anomaly.fetch_ztf_lightcurve` + `alerce_oids` hit ALeRCE today; the
  same calls return LSST objects once LSST flows. The calibrated class basins
  (`ztf_class_basins.json`) recalibrate on LSST samples with the same
  `calibrate_from_alerce` routine.
- **The orbital-element tools are catalogue-agnostic.** They consume (a, e, i, Omega,
  omega, ...) rows; LSST's SSObjects table is the same shape as the JPL SBDB queries we
  already run, just deeper and southern.
- **The moving-object arsenal is validated on real DECam** (the same camera lineage and
  CP product format as Rubin's pipeline outputs), and the cadence-quality gate ensures
  it only trusts properly-spaced multi-epoch data -- which LSST provides by design.

## Switch-on checklist (when the LSST stream is public)

1. Point `ztf_anomaly` at the LSST broker endpoint (ALeRCE/Lasair), recalibrate the
   class basins on a labeled LSST sample, and **add the resonant-angle and richer
   features identified this session** before trusting novelty flags at scale.
2. Run `exotic_orbit_hunt` + the clustering/bias-null suite on the LSST SSObjects table
   as it grows (the southern + galactic-plane objects that are missing today).
3. Run the GPU shift-stack arsenal on LSST deep-drilling fields and the galactic-plane
   visits -- the under-documented sky we could not search for lack of cadenced data.
4. **Carry the discipline forward:** every flag gets a chance-excess or dynamical null,
   every confirmation gets a validated diagnostic (the co-orbital retraction this
   session is the cautionary example -- the resonant-angle diagnostic must be validated
   against a known object before any co-orbital claim).

## Known gaps to close before the stream (from this session's honest findings)

- **Co-orbital confirmation:** replace the a-band proxy with validated resonant-angle
  libration (REBOUND + Horizons machinery is built; the diagnostic is the missing
  piece). Do not claim co-orbitals until it confirms Ka'epaoka'awela.
- **ZTF/LSST novelty:** complete the class basis (all major classes) + richer features
  (added eta, frac_bright; the full ALeRCE feature set is the target) so anomalies are
  genuine novelties, not unmodeled known classes.
- **TESS-style transit completeness:** the SNR-confidence + period-alias fixes recovered
  shallow planets (TOI-700 d, LHS 3844 b); stellar-activity contamination (GJ 1132)
  remains a hard case needing better detrending.

The durable product of all this work is the validated, self-skeptical pipeline. It is
ready to be pointed at Rubin, with the precursor (ZTF/DECam/JPL) paths exercised today
and the specific refinements needed before the firehose clearly itemized above.
