"""
peak_matcher.py (thin shim)

The real implementation now lives in the installed `nmrmatch` package
(nmrmatch.matching), built from this exact file when the project was
packaged for pip/GitHub installation -- see ../nmrmatch/. This file is
kept only so `python peak_matcher.py ...` keeps working unchanged for
anyone with existing muscle memory or automation, without this
directory also carrying its own, possibly-drifting copy of the matching
logic (utils.py's score_compound() through a local `from utils import
score_compound`, which this shim no longer needs).

Prefer going forward:
    nmrmatch-match ...                        (installed console script)
    from nmrmatch import match_sample, match_all_samples
"""
from nmrmatch.matching import (  # noqa: F401
    DEFAULT_TOLERANCE_PPM,
    match_sample,
    match_all_samples,
    main,
)

if __name__ == "__main__":
    main()
