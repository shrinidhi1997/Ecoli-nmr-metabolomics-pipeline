"""
build_bmrb_library.py (thin shim)

The real implementation now lives in the installed `nmrmatch` package
(nmrmatch.reflib.build_bmrb), built from this exact file when the
project was packaged for pip/GitHub installation -- see ../nmrmatch/.
This file is kept only so `python build_bmrb_library.py ...` keeps
working unchanged for anyone with existing muscle memory or automation.

Prefer going forward:
    nmrmatch-build-bmrb ...                   (installed console script)
"""
from nmrmatch.reflib.build_bmrb import main  # noqa: F401

if __name__ == "__main__":
    main()
