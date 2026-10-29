"""
fetch_nmr_reference.py (thin shim)

The real implementation now lives in the installed `nmrmatch` package
(nmrmatch.reflib.fetch), built from this exact file when the project
was packaged for pip/GitHub installation -- see ../nmrmatch/. This file
is kept only so `python fetch_nmr_reference.py ...` keeps working
unchanged for anyone with existing muscle memory or automation, and so
`scrape_ecmdb_compound_page`/`build_reference_csv` are still importable
from this filename if anything (a notebook, a one-off script) does so.

See nmrmatch.reflib.fetch's own module docstring for the full context
on why this file exists (ECMDB's bulk export has no shift data; this is
option B, a scrape template; option C, build_bmrb_library.py, is what
this pipeline's shipped reference library actually comes from).

Prefer going forward:
    nmrmatch-fetch-reference ...              (installed console script)
    from nmrmatch.reflib.fetch import build_reference_csv
"""
from nmrmatch.reflib.fetch import (  # noqa: F401
    HEADERS,
    scrape_ecmdb_compound_page,
    build_reference_csv,
    main,
)

if __name__ == "__main__":
    main()
