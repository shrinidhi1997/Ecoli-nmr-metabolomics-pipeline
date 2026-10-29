"""
nmrmatch.reflib.filter_ecoli
(originally filter_reference_to_ecoli.py)

Filters reference_peaks_combined.csv (or the raw BMRB/HMDB library CSVs
individually) down to rows confirmed as real E. coli metabolites, using
the is_ecoli_metabolite column build_bmrb_library.py / build_hmdb_library.py
already compute from an exact-ID cross-reference against ECMDB
(compound_library.csv, built by build_library.py from a real ecmdb.json).

Why this, instead of matching against ASICS's bundled compound list:
ASICS's 191-compound library isn't organism-specific, and cross-referencing
it against BMRB/HMDB would need fuzzy name matching (the same real compound
is often named differently across databases -- "D-Glucose" vs
"alpha-D-Glucose", "Glutamate" vs "L-Glutamic acid" -- already confirmed
unreliable elsewhere in this project, see build_kegg_ecoli_list.py's
normalize_compound_name()). ECMDB cross-referencing is exact-ID (KEGG
Compound ID, ChEBI ID, PubChem Compound ID, or InChIKey -- any one
matching is enough, see build_bmrb_library.py's is_ecoli_by_ecmdb()), and
ECMDB is itself curated specifically for E. coli K-12 MG1655, so this is
the more authoritative filter already sitting in your output/ folder, not
a new thing to build.

WHY THIS WAS PROBABLY MISSING FROM YOUR LAST RUN
is_ecoli_metabolite is only populated as True/False if build_bmrb_library.py
was run with --ecmdb-library ../output/compound_library.csv and/or
build_hmdb_library.py with --compound-library ../output/compound_library.csv.
Without that flag, the column is either absent or every value is blank
(is_ecoli_by_ecmdb() returns None, meaning "not checked", not "confirmed
absent" -- see that function's own docstring). The lignin/plant-flavonoid
flood you saw is consistent with that cross-reference never having been
run, or having been run but not yet filtered on.

This script treats a blank/NaN is_ecoli_metabolite as "not checked" and
reports it separately from "checked and False" -- rather than silently
dropping both into the same bucket -- since those are different situations
needing different fixes (re-run the builders with the ECMDB flag, vs. this
compound genuinely isn't E. coli-linked).

Run:
    python filter_reference_to_ecoli.py \\
        --reference ../output/reference_peaks_combined.csv \\
        --out ../output/reference_peaks_combined_ecoli_only.csv

Then point match_r_feature_table.py / peak_matcher.py at the filtered file
instead of the full combined one.
"""
from __future__ import annotations
import argparse
from pathlib import Path

import pandas as pd

TRUE_STRINGS = {"true", "1", "1.0", "yes"}
FALSE_STRINGS = {"false", "0", "0.0", "no"}


def to_tri_state(value) -> bool | None:
    """CSV round-tripping turns a real Python bool column into a mix of
    actual bools, string "True"/"False" (what pandas writes a bool Series
    as), and blank/NaN for untested rows -- normalizing all of that here so
    filtering doesn't silently miscount any of the three states.
    Returns True, False, or None (not checked)."""
    if pd.isna(value):
        return None
    if isinstance(value, bool):
        return value
    s = str(value).strip().lower()
    if s in TRUE_STRINGS:
        return True
    if s in FALSE_STRINGS:
        return False
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reference", required=True,
                     help="reference_peaks_combined.csv (or a single-source "
                          "reference_peaks_bmrb_library.csv / "
                          "reference_peaks_hmdb_library.csv)")
    ap.add_argument("--out", default="../output/reference_peaks_combined_ecoli_only.csv")
    args = ap.parse_args()

    df = pd.read_csv(args.reference)
    if "is_ecoli_metabolite" not in df.columns:
        raise SystemExit(
            "No is_ecoli_metabolite column in this file at all. Re-run "
            "build_bmrb_library.py with --ecmdb-library ../output/compound_library.csv "
            "and/or build_hmdb_library.py with --compound-library "
            "../output/compound_library.csv, then re-run merge_reference_libraries.py, "
            "before filtering."
        )

    state = df["is_ecoli_metabolite"].map(to_tri_state)
    n_total_compounds = df["compound_id"].nunique()
    n_true_compounds = df.loc[state == True, "compound_id"].nunique()  # noqa: E712
    n_false_compounds = df.loc[state == False, "compound_id"].nunique()  # noqa: E712
    n_unchecked_compounds = df.loc[state.isna(), "compound_id"].nunique()

    print(f"{n_total_compounds} total compounds in {args.reference}:")
    print(f"  {n_true_compounds} confirmed E. coli metabolites (kept)")
    print(f"  {n_false_compounds} checked against ECMDB and NOT confirmed (dropped)")
    print(f"  {n_unchecked_compounds} never checked against ECMDB at all (dropped, "
          f"but these are a DIFFERENT situation -- see note below)")

    if n_true_compounds == 0:
        raise SystemExit(
            "\nZero compounds are confirmed E. coli metabolites -- filtering would "
            "produce an empty, useless reference library. This almost always means "
            "is_ecoli_metabolite was never populated (every row is 'not checked', "
            "not 'checked and false'). Re-run build_bmrb_library.py / "
            "build_hmdb_library.py with the ECMDB cross-reference flag (see this "
            "script's module docstring), then merge_reference_libraries.py again, "
            "before filtering."
        )

    if n_unchecked_compounds > 0:
        print(f"\nNote: {n_unchecked_compounds} compound(s) were never cross-referenced "
              f"against ECMDB and are being dropped along with the confirmed-not-E.coli "
              f"ones. If that number is large relative to {n_total_compounds}, the ECMDB "
              f"cross-reference may only be partially run -- re-running the builders with "
              f"the ECMDB flag could recover some of them as confirmed matches instead of "
              f"losing them here.")

    filtered = df.loc[state == True].copy()  # noqa: E712
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    filtered.to_csv(args.out, index=False)
    print(f"\nWrote {len(filtered)} peak rows across {n_true_compounds} confirmed "
          f"E. coli compounds to {args.out}")


if __name__ == "__main__":
    main()
