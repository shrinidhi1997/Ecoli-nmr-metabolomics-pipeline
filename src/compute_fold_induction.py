"""
compute_fold_induction.py (thin shim)

The real implementation now lives in the installed `nmrmatch` package
(nmrmatch.fold_change), built from this exact file when the project was
packaged for pip/GitHub installation -- see ../nmrmatch/. This file is
kept only so `python compute_fold_induction.py ...` keeps working
unchanged for anyone with existing muscle memory or automation.

Prefer going forward:
    nmrmatch-fold-change ...                  (installed console script)
    from nmrmatch.fold_change import compute_fold_induction, parse_group_arg
"""
from nmrmatch.fold_change import (  # noqa: F401
    parse_group_arg,
    compute_fold_induction,
    main,
)

if __name__ == "__main__":
    main()
