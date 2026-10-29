"""
nmrmatch.reflib
================

Utilities for building and curating the reference NMR peak library that
`nmrmatch.matching` matches experimental peaks against: fetching/scraping
candidate peak tables, pulling BMRB's metabolomics collection in bulk,
filtering a broad library down to a target organism's known metabolites,
and merging multiple reference libraries into one CSV. These are the
scripts used once, offline, to build `reference_peaks*.csv` -- not part
of the per-sample matching pipeline itself.
"""
