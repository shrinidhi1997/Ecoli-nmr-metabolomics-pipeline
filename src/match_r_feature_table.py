"""
match_r_feature_table.py (thin shim)

The real implementation now lives in the installed `nmrmatch` package
(nmrmatch.cli.match_feature_table), built from this exact file when the
project was packaged for pip/GitHub installation -- see ../nmrmatch/.
This file is kept only so `python match_r_feature_table.py ...` keeps
working unchanged for anyone with existing muscle memory or automation,
without this directory also carrying its own, possibly-drifting copy of
the reshaping/matching logic (the local `from peak_matcher import
match_all_samples`, which this shim no longer needs).

Prefer going forward:
    nmrmatch-match-feature-table ...          (installed console script)
    from nmrmatch.cli.match_feature_table import wide_to_long, parse_ppm_label
"""
from nmrmatch.cli.match_feature_table import (  # noqa: F401
    PPM_LABEL_RE,
    parse_ppm_label,
    wide_to_long,
    main,
)

if __name__ == "__main__":
    main()
