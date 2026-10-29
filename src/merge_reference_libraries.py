"""
merge_reference_libraries.py (thin shim)

The real implementation now lives in the installed `nmrmatch` package
(nmrmatch.reflib.merge), built from this exact file when the project
was packaged for pip/GitHub installation -- see ../nmrmatch/. This file
is kept only so `python merge_reference_libraries.py ...` keeps working
unchanged for anyone with existing muscle memory or automation.

Prefer going forward:
    nmrmatch-merge-libraries ...              (installed console script)
"""
from nmrmatch.reflib.merge import main  # noqa: F401

if __name__ == "__main__":
    main()
