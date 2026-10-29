"""
filter_reference_to_ecoli.py (thin shim)

The real implementation now lives in the installed `nmrmatch` package
(nmrmatch.reflib.filter_ecoli), built from this exact file when the
project was packaged for pip/GitHub installation -- see ../nmrmatch/.
This file is kept only so `python filter_reference_to_ecoli.py ...`
keeps working unchanged for anyone with existing muscle memory or
automation.

Prefer going forward:
    nmrmatch-filter-ecoli ...                 (installed console script)
"""
from nmrmatch.reflib.filter_ecoli import main  # noqa: F401

if __name__ == "__main__":
    main()
