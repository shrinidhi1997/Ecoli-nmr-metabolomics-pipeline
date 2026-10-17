"""
merge_reference_libraries.py

Combines the two independently-built, real-data reference peak libraries
(build_bmrb_library.py's BMRB pull and build_hmdb_library.py's HMDB
peaklist-archive pull, both already confirmed against ECMDB) into ONE
reference table for peak_matcher.py to search against, so a real sample
spectrum gets checked against compounds from both sources in a single
run instead of two separate, smaller searches.

Why this is a *merge*, not a plain concat:
  * The two builders use different column sets. BMRB's output has
    chebi_id/pubchem_id but no intensity (BMRB's own _Peak_char loop
    doesn't carry an intensity value in the real data seen so far); HMDB's
    has intensity but no chebi_id/pubchem_id. peak_matcher.py only
    strictly requires compound_id/compound_name/nucleus/shift_ppm and
    treats everything else (including intensity) as optional, so this
    harmonizes both onto one shared column set rather than leaving ragged
    columns or silently dropping data either side has.
  * The two builders also used different names for the same real flag:
    BMRB's is "is_ecoli_metabolite", HMDB's is "is_ecoli". This renames
    HMDB's to match so the combined file has ONE flag column, not two
    that mean the same thing.
  * A "source" column ("BMRB" or "HMDB") is added. The same real compound
    can legitimately show up under both databases' own accession scheme
    (e.g. a bmseXXXXXX row and an HMDBXXXXXXX row for the same molecule) --
    peak_matcher.py already scores two different compound_ids
    independently, so nothing breaks, but "source" lets you tell, after
    matching, which database a given hit's peaks actually came from.

Does NOT deduplicate rows across sources by compound identity (e.g. by
name or InChIKey): the same real compound can carry a different,
independently-measured peak set in each database, and collapsing those
would discard real corroborating evidence rather than combine it.
peak_matcher.py already scores each (compound_id, nucleus) group
separately, so a compound appearing under both a BMRB and an HMDB
accession just shows up as two candidate rows in match_results.csv --
recognizable by compound_name/source, not a bug.

Run:
    python merge_reference_libraries.py \
        --bmrb ../output/reference_peaks_bmrb_library.csv \
        --hmdb ../output/reference_peaks_hmdb_library.csv \
        --out ../output/reference_peaks_combined.csv

Then point peak_matcher.py at the combined file:
    python peak_matcher.py --reference ../output/reference_peaks_combined.csv \
        --experimental ... --out ../output/match_results.csv
"""
from __future__ import annotations
import argparse
from pathlib import Path

import pandas as pd

# compound_id/compound_name/nucleus/shift_ppm are peak_matcher.py's hard
# requirement; the rest are optional extras either builder can supply.
OUTPUT_COLUMNS = [
    "compound_id", "compound_name", "nucleus", "shift_ppm", "intensity",
    "kegg_id", "chebi_id", "pubchem_id", "inchikey",
    "is_ecoli_metabolite", "source",
]

# Explicit dtype for a column this builder's output never has (e.g. HMDB
# rows have no chebi_id/pubchem_id; BMRB rows have no intensity). Without
# this, filling a missing column with bare pd.NA leaves its dtype for
# pandas to guess from an all-NA column, which pandas itself now warns is
# changing in a future version (seen live: a real FutureWarning on a real
# 4,109-compound merge) -- giving every column an explicit dtype up front
# avoids relying on that guess at all, now or after a pandas upgrade.
OUTPUT_DTYPES = {
    "compound_id": "object", "compound_name": "object", "nucleus": "object",
    "shift_ppm": "float64", "intensity": "float64",
    "kegg_id": "object", "chebi_id": "object", "pubchem_id": "object",
    "inchikey": "object", "is_ecoli_metabolite": "object", "source": "object",
}


def load_bmrb(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    df["source"] = "BMRB"
    return df


def load_hmdb(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    df["source"] = "HMDB"
    if "is_ecoli" in df.columns:
        df = df.rename(columns={"is_ecoli": "is_ecoli_metabolite"})
    return df


def merge_libraries(bmrb_df: pd.DataFrame, hmdb_df: pd.DataFrame) -> pd.DataFrame:
    """Harmonizes both frames onto OUTPUT_COLUMNS (any column one side
    lacks is filled with NA rather than dropped) and stacks them, keeping
    every row from both -- see module docstring for why no dedup happens
    here."""
    frames = []
    for df in (bmrb_df, hmdb_df):
        df = df.copy()
        for col in OUTPUT_COLUMNS:
            if col not in df.columns:
                dtype = OUTPUT_DTYPES[col]
                # pd.NA itself can't be cast into a float64 array (real
                # error hit while testing this fix); float columns need
                # float("nan") instead, object columns take None either way.
                fill_value = float("nan") if dtype == "float64" else None
                df[col] = pd.Series([fill_value] * len(df), dtype=dtype, index=df.index)
        frames.append(df[OUTPUT_COLUMNS])
    return pd.concat(frames, ignore_index=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bmrb", required=True,
                     help="reference_peaks_bmrb_library.csv from build_bmrb_library.py")
    ap.add_argument("--hmdb", required=True,
                     help="reference_peaks_hmdb_library.csv from build_hmdb_library.py -- "
                          "build it with --include-unconfirmed for the biggest combined "
                          "database; is_ecoli_metabolite in the combined output still "
                          "tells you which rows are confirmed E. coli compounds")
    ap.add_argument("--out", default="../output/reference_peaks_combined.csv")
    args = ap.parse_args()

    bmrb_df = load_bmrb(args.bmrb)
    hmdb_df = load_hmdb(args.hmdb)
    combined = merge_libraries(bmrb_df, hmdb_df)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(args.out, index=False)

    n_bmrb_compounds = bmrb_df["compound_id"].nunique()
    n_hmdb_compounds = hmdb_df["compound_id"].nunique()
    n_confirmed = combined.loc[
        combined["is_ecoli_metabolite"] == True, "compound_id"  # noqa: E712
    ].nunique()
    print(f"Wrote {len(combined)} peaks across "
          f"{combined['compound_id'].nunique()} compounds "
          f"({n_bmrb_compounds} from BMRB + {n_hmdb_compounds} from HMDB) to {args.out}")
    print(f"{n_confirmed} of those compounds are confirmed E. coli metabolites "
          f"(is_ecoli_metabolite == True).")


if __name__ == "__main__":
    main()
