"""
fetch_nmr_reference.py

WHY THIS FILE EXISTS (read this before running anything)
----------------------------------------------------------------
ECMDB's bulk "All metabolite information" download (see build_library.py)
gives you compound identity/structure data for every E. coli metabolite,
but it does NOT include the actual NMR chemical-shift peak tables. The
field documentation for the bulk export lists only structural/taxonomic
fields -- name, formula, SMILES, InChI, classification, etc. No shift or
spectrum fields. The actual peak data lives behind ECMDB's interactive
spectral tools (1D NMR search, 2D NMR search) and per-compound spectrum
pages, not a single downloadable file.

So there are three honest ways to get from "my peak list" to "candidate E.
coli compounds":

OPTION A -- fastest, no code, right now
  Go to https://ecmdb.ca/spectra/nmr/one_d/search/new , choose 1H or 13C,
  paste your peak list (one shift per line), set a tolerance, and search.
  This queries ECMDB's own curated spectral library directly. For a
  handful of spectra, just use this. (Note: ECMDB has been reported
  unreachable from some networks -- if so, this option isn't available to
  you and the BMRB-based library (option C) is the path forward.)

OPTION B -- manual/semi-manual export from ECMDB
  Open each compound's page (https://ecmdb.ca/compounds/<ECMDB_ID>), copy
  its NMR peak table(s) into a CSV. Realistic for a target list of tens of
  compounds, tedious for thousands.

OPTION C -- scriptable / batch, the one actually used by this pipeline
  build_bmrb_library.py pulls real experimental NMR peak data in bulk from
  BMRB's metabolomics collection (a different, independent database from
  ECMDB) via their documented REST API. This is what reference_peaks_bmrb_
  library.csv in this project actually comes from. See that file for
  details -- it covers several hundred well-characterized compounds with
  real experimental spectra, not all ~3,755 ECMDB compounds, but it is a
  genuinely scriptable and verified path, unlike scraping ECMDB pages.

reference_peaks.csv schema (what peak_matcher.py expects):
    compound_id, compound_name, nucleus, shift_ppm, intensity, multiplicity
    ECMDB00042, L-Alanine, 1H, 1.47, 100, d
    ECMDB00042, L-Alanine, 1H, 3.78, 60, q
    ...
`intensity` and `multiplicity` may be left blank if unknown; the matcher
degrades gracefully (it just skips the intensity-correlation bonus).

The scrape_ecmdb_compound_page() function below is kept as a template for
option B, in case you want to automate scraping a small target list of
ECMDB compound pages yourself. It has NOT been verified against a live
page in this environment -- treat it as a starting point: run it on ONE
compound first, print what it finds, and fix the BeautifulSoup selectors
to match what you actually see before trusting a batch run.
"""

from __future__ import annotations
import time
from pathlib import Path

import pandas as pd
import requests
from bs4 import BeautifulSoup

HEADERS = {"User-Agent": "research-script/1.0 (contact: set-your-email-here)"}


def scrape_ecmdb_compound_page(ecmdb_id: str, nucleus: str = "1H",
                                delay_s: float = 1.0) -> pd.DataFrame:
    """
    Best-effort fetch of one ECMDB compound page's NMR peak table.
    See the module docstring: verify selectors before batch use.
    """
    url = f"https://ecmdb.ca/compounds/{ecmdb_id}"
    resp = requests.get(url, headers=HEADERS, timeout=30)
    resp.raise_for_status()
    time.sleep(delay_s)  # be polite; don't hammer a public research DB

    soup = BeautifulSoup(resp.text, "html.parser")
    rows = []

    spectrum_links = [
        a["href"] for a in soup.find_all("a", href=True)
        if "/spectra/" in a["href"] and nucleus.lower() in a["href"].lower()
    ]

    for link in spectrum_links:
        spec_url = link if link.startswith("http") else f"https://ecmdb.ca{link}"
        s_resp = requests.get(spec_url, headers=HEADERS, timeout=30)
        time.sleep(delay_s)
        if s_resp.status_code != 200:
            continue
        s_soup = BeautifulSoup(s_resp.text, "html.parser")

        for table in s_soup.find_all("table"):
            headers_text = [th.get_text(strip=True).lower() for th in table.find_all("th")]
            if not any("shift" in h or "ppm" in h for h in headers_text):
                continue
            for tr in table.find_all("tr")[1:]:
                cells = [td.get_text(strip=True) for td in tr.find_all("td")]
                if len(cells) >= 1:
                    try:
                        shift = float(cells[0])
                    except ValueError:
                        continue
                    intensity = cells[1] if len(cells) > 1 else None
                    multiplicity = cells[2] if len(cells) > 2 else None
                    rows.append({
                        "compound_id": ecmdb_id,
                        "nucleus": nucleus,
                        "shift_ppm": shift,
                        "intensity": intensity,
                        "multiplicity": multiplicity,
                    })

    return pd.DataFrame(rows)


def build_reference_csv(ecmdb_ids: list[str], out_path: Path,
                         nucleus: str = "1H") -> pd.DataFrame:
    """Loop scrape_ecmdb_compound_page over a target compound list."""
    frames = []
    for cid in ecmdb_ids:
        try:
            df = scrape_ecmdb_compound_page(cid, nucleus=nucleus)
            if df.empty:
                print(f"  {cid}: no peaks parsed (check selectors, see module docstring)")
            else:
                print(f"  {cid}: {len(df)} peaks")
            frames.append(df)
        except requests.RequestException as e:
            print(f"  {cid}: request failed ({e})")
    out = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    out.to_csv(out_path, index=False)
    print(f"Wrote {len(out)} rows to {out_path}")
    return out


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ecmdb-ids", nargs="+", required=True,
                     help="e.g. --ecmdb-ids ECMDB00042 ECMDB00061")
    ap.add_argument("--nucleus", default="1H", choices=["1H", "13C"])
    ap.add_argument("--out", default="../output/reference_peaks.csv")
    args = ap.parse_args()
    build_reference_csv(args.ecmdb_ids, Path(args.out), nucleus=args.nucleus)


if __name__ == "__main__":
    main()
