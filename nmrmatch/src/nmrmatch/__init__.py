"""
nmrmatch
========

Reference-library matching, confidence scoring, and fold-change
computation for 1D NMR metabolomics. Companion package to `nmrpick`,
which handles baseline correction, peak detection, and group-comparison
statistics.

The most commonly used functions are re-exported here so callers can do:

    from nmrmatch import score_compound, match_sample, match_all_samples

rather than reaching into the submodules directly.
"""

from .scoring import (
    gaussian_match_weight,
    candidate_matches,
    score_compound,
    gaussian_match_weight_2d,
    candidate_matches_2d,
    score_compound_2d,
)
from .matching import match_sample, match_all_samples

__all__ = [
    "gaussian_match_weight",
    "candidate_matches",
    "score_compound",
    "gaussian_match_weight_2d",
    "candidate_matches_2d",
    "score_compound_2d",
    "match_sample",
    "match_all_samples",
]

__version__ = "0.1.0"
