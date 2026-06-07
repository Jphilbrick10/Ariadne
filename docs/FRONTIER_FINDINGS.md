# Frontier anomaly-hunt: the toolkit and the honest findings

A complete, reusable anomaly-hunting pipeline for the outer solar system and beyond,
built and validated on real public data, plus the honest record of what it found.
The governing discipline throughout: every candidate carries a chance-excess or
dynamical null, so the pipeline resists manufacturing discoveries and reports clean
negatives plainly.

## The toolkit (all in `src/ariadne/discovery/frontier/` + `scripts/`)

| Module / script | What it does |
|---|---|
| `tno_clustering` | Extreme-TNO orbital-orientation clustering + outlier scoring; Planet Nine signal |
| `planet_nine` | Rigorous clustering (Rayleigh+Kuiper), plane warp, sky map, magnitude, **secular dynamical test**, **selection-bias null**, **planetary-perturbation (Cassini) constraint**, multi-population scan |
| `perturber_fusion` | Multi-channel Equation-of-One fusion: one perturber estimate from all signals, with per-channel POSITIVE/PERMISSIVE/CONFLICT verdicts |
| `black_hole_discriminator` | Built from the Coherence framework's coherence-horizon physics; planet-vs-black-hole model selection |
| `exotic_orbit_hunt` | D-criterion genetic-pair search (fragment-aware) + capture-candidate scoring, with a chance-excess null |
| `frontier_scanner` + `scan_everything.py` | Automated anomaly suite over 14 populations (~98k objects), honest nulls |
| `confirm_pair.py` | Backward N-body integration to confirm/refute asteroid pairs |
| `check_coorbital.py` / `coorbital_sweep.py` | a-libration test for 1:1 co-orbital resonance; census of exotic co-orbitals |
| `ztf_anomaly` + `hunt_ztf_anomalies.py` | Physics-coherence novelty triage for the ZTF/LSST stream (calibrated on ALeRCE) |
| `tess_vetting` + `hunt_tess_candidates.py` | Transit search + coherence vetting of TESS light curves |
| `search_gap_field.py` | Mine under-documented sky (galactic plane / far south) with archival DECam |
| `map_solar_system_3d.py` | Interactive 3D map with the perturber placed |

## The Planet Nine verdict (the central investigation)

Run on the real extreme-TNO data with every test we could devise:

- **The textbook evidence does not survive.** The perihelion (apsidal) clustering --
  *the* signature of a single shepherding planet -- is statistically consistent with
  galactic-plane selection bias (our bias Monte Carlo: p=0.085; the OSSOS survey
  reached the same conclusion). A planet does not maintain the clustering better than
  the giant planets do alone (secular test), and the clustering oscillates on its own
  over ~100 Myr (backward integration).
- **Two anomalies do survive every test:** the distant objects share a common orbital
  plane warped ~3-5 deg from the invariable plane (bias-tested, p=0.005), and objects
  exist on orbits the known planets cannot have made (perihelia to 81 AU).
- **The fusion of all channels** lands on a single best estimate (a~600 AU, ~8 Earth
  masses, node ~131 deg, anti-aligned) but is honest that only the node + plane are
  POSITIVE evidence (and they are one underlying signal); the apsidal, Cassini, comet,
  and dynamical channels are PERMISSIVE, not confirmatory. Joint significance is the
  node's ~2.8 sigma -- fusion does not manufacture more.
- **Black hole vs planet:** the Coherence framework's own physics says a black hole is
  gravitationally identical to a planet (we verified the coherence-force ratio = 1 at
  orbital distance); only electromagnetic darkness could distinguish them, and that is
  prior-dominated + needs near-complete deep coverage. Neither favored nor excluded.

**Bottom line:** most likely there is no textbook Planet Nine; the real anomaly (the
warped plane + the detached objects) is as well explained by a more distant/massive
perturber, a long-departed stellar flyby, or collective dynamics. Rubin/LSST settles
it within ~2 years.

## The discovery sweep (what scanning everything actually found)

- **Scanned ~98,000 objects** across 14 populations. After honest nulls, the inner
  populations' "clustering" is known selection + secular-resonance + Lagrange-cloud
  structure (not anomalies); the comet "1000-sigma families" are split-comet fragments
  + the real Kreutz sungrazer family (known).
- **Every candidate pair refuted:** the Atira and Centaur pairs diverge chaotically
  within ~10 kyr (no traceable common origin); the scattered-disk "pair" was four
  objects from the *same survey field on the same night* (selection artifact).
- **Exotic candidates flagged** (genuinely retrograde, near a giant planet's
  semimajor axis): 2012 TL139 (Neptune), 2006 BZ8 (Saturn), and the known retrograde
  Jupiter co-orbital 514107 Ka'epaoka'awela. These are the captured-object class and
  worth real follow-up.
- **RETRACTION (honesty correction):** an earlier draft "confirmed" 2012 TL139 as a
  retrograde Neptune co-orbital via a semimajor-axis-band libration test. On rigorous
  re-examination with REBOUND (MERCURIUS) + orbit-uncertainty clones + accurate
  JPL-Horizons initial conditions, that diagnostic FAILS VALIDATION: it does not
  confirm even the known co-orbital Ka'epaoka'awela (~44% in-band). The a-band proxy
  is inadequate -- a co-orbital's semimajor axis can swing widely while the RESONANT
  ANGLE librates. So the co-orbital confirmations are **withdrawn**; a trustworthy
  result needs resonant-angle libration analysis (subtle for retrograde resonances),
  which is identified as the missing piece. The integrator + clone machinery is sound;
  the diagnostic is not yet correct. This is the discipline working as intended --
  the rigorous follow-up caught an over-claim before it stood.
- **Resolution part 1 -- prograde (validated):** resonant-angle libration
  (phi = lambda - lambda_planet) on REBOUND + Horizons VALIDATES for the prograde case:
  588 Achilles (a textbook L4 Jupiter Trojan) librates in 6% of the circle, centered at
  61 deg = the L4 point exactly. Prograde co-orbital confirmation is trustworthy.
- **Resolution part 2 -- the retraction was itself a numerical artifact (RESTORED):**
  the a-band failures (2006 BZ8 "escaping to 31 AU", Ka'epaoka'awela "44% in-band") were
  traced to a too-COARSE integration timestep mishandling the deep planet-crossing
  encounters -- not a flaw in the concept. With the adaptive IAS15 integrator (30 kyr),
  the validation object Ka'epaoka'awela is 100% a-confined (5.08-5.36 around Jupiter's
  5.20), AND both candidates are a-confined: 2012 TL139 98% (a 29-37 around Neptune),
  2006 BZ8 100% (a 9.3-9.8 around Saturn). So 2012 TL139 and 2006 BZ8 are RESTORED as
  likely (temporary) retrograde co-orbitals -- the captured-object class.
- **Remaining honest caveats:** a-confinement under accurate integration is strong
  SUGGESTIVE evidence, not the gold standard. The retrograde resonant ANGLE still does
  not librate in any simple sign-combination (the Morais-Namouni retrograde argument is
  the genuinely-open piece), and 30 kyr is short vs the Myr-scale stability a publication
  claim needs. These objects are also likely already cataloged. Net: validated prograde
  tool; retrograde candidates restored as likely co-orbitals on accurate-integration
  a-confinement (validated on Ka'epaoka'awela), with the resonant-angle proof still open.
  The lesson cuts both ways -- retract when a method fails validation, RESTORE when the
  failure turns out to be a numerical artifact.

## Searching the under-documented sky

- Quantified the gaps: known TNOs are ~30% deficient in the galactic plane (\|b\|<15)
  and near-absent at Dec<-40.
- Built the full archival-search pipeline (find field -> download -> shift-stack ->
  vet) and ran it on the Sgr galactic center. Result: confusion artifacts only, and a
  key lesson -- the gaps lack the *well-cadenced* multi-epoch data moving-object search
  needs (the archival data there was taken for other science). Added a cadence-quality
  gate so the tool refuses to trust degenerate cadences.

## The honest meta-conclusion

Re-mining the solar-system catalog is exhausted: it is picked over, and the genuine
gaps lack the cadenced data to search. No new solar-system discovery emerged -- but
the toolkit refused to manufacture false ones, confirmed a real exotic object, and
reproduced the live Planet Nine debate from scratch. The durable result is the
validated, honest pipeline itself, ready to point at the Rubin/LSST stream -- the
firehose of fresh, properly-cadenced, all-southern-sky data now coming online -- where
its triage and vetting are exactly what the flood will need.
