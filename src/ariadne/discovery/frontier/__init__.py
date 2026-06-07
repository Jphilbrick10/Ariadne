"""Frontier targets -- where Ariadne's coherence selector has a real, honest shot
beyond plain solar-system moving-object detection.

Three capabilities, each aimed at a place where the bottleneck is SELECTION /
VETTING / ANOMALY (the regime the Equation-of-One engine actually wins at) rather
than raw detection (which mature tools already own), and each on free public data:

  tno_clustering  -- orbital-orientation outlier hunt in the extreme-TNO
                     population (the Planet Nine clustering signal), JPL SBDB.
  tess_vetting    -- transit search + coherence vetting of the candidates a
                     standard pipeline discards, free TESS light curves (MAST).
  ztf_anomaly     -- physics-coherence novelty score over transient light curves
                     so genuinely-new objects float to the top, public ZTF.

The honest edge across all three: a coherence energy needs no training set, so it
ranks candidates by agreement with KNOWN PHYSICS -- exactly the low-data / novel /
anomaly regime where supervised ML starves.
"""
