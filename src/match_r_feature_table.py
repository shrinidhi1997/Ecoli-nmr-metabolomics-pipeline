"""
match_r_feature_table.py

Runs peak_matcher.py's real multi-peak-evidence compound matching
against the R pipeline's merged peak list, instead of
annotate_r_feature_table.py's single-bin lookup.

Why this replaces annotate_r_feature_table.py for identification
-------------------------------------------------------------------
annotate_r_feature_table.py labels each ppm bin with every reference
compound that has ANY one peak within tolerance of that bin alone.
Against a broad, multi-thousand-compound BMRB+HMDB library, that
floods every bin with dozens of compounds that happen to share a peak
near that one shift (plant flavonoids, lignin breakdown products,
drugs -- nothing plausible in an E. coli lysate), because a single
coincidental peak carries no real evidence either way.

peak_matcher.py's match_sample()/score_compound() instead scores a
compound by how many of ITS OWN characteristic peaks show up together
in the sample (coverage = n_matched / n_ref), with min_ref_peaks=2 by
default specifically to reject single-peak coincidences -- already
justified in that file's own docstring against a ~3,300-compound
library test. A compound only ranks highly when several of its real
peaks line up with real peaks in your spectrum, not one lucky
coincidence.

The only gap is shape: the R pipeline writes a WIDE table (one row
per ppm feature/"bin", one column per sample --
write_metaboanalyst_csv()'s own output), and match_all_samples()
expects a LONG table (one row per sample per detected peak). This
script does that reshape, then runs the real matcher.

Run (from the project's src/ folder):
    python match_r_feature_table.py \
        --feature-table ../output/control_vs_sample/peak_list_control_vs_sample.csv \
        --reference ../output/reference_peaks_combined.csv \
        --min-ref-peaks 2 \
        --top-n 10 \
        --out ../output/control_vs_sample/compound_matches_control_vs_sample.csv

Works against either the un-annotated feature table
(peak_list_control_vs_sample.csv) or the _annotated.csv file
annotate_r_feature_table.py wrote -- PPM_LABEL_RE below only reads the
leading ppm_X.XXXX off the Sample label and ignores anything glued on
after it, so an already-annotated file's compound-name text doesn't
need stripping first.
"""
from __future__ import annotations
import argparse
import re
from pathlib import Path

import pandas as pd

from peak_matcher import match_all_samples

# Same label shape write_metaboanalyst_csv() always writes:
# paste0("ppm_", sprintf("%.4f", merged$ppm)). No trailing $ here on
# purpose: annotate_r_feature_table.py's output glues " (compound1,
# compound2, ...)" onto this same prefix, and the ppm value itself is
# unchanged by that, so matching just the leading number lets this
# script work directly off an already-annotated file too, not only
# the original un-annotated one.
PPM_LABEL_RE = re.compile(r"^ppm_(-?\d+\.\d+)")


def parse_ppm_label(label: str) -> float | None:
    m = PPM_LABEL_RE.match(str(label).strip())
    return float(m.group(1)) if m else None


def wide_to_long(feature_table_path: str, nucleus: str = "1H") -> pd.DataFrame:
    """Reshapes write_metaboanalyst_csv()'s output (row 0 Sample names,
    row 1 Label/group, rows 2+ one 'ppm_X.XXXX' feature per row, blank
    where a sample has no value for that feature) into one row per
    (sample_id, detected peak) -- the shape match_all_samples()
    expects. A blank/NA cell means no peak was detected there for
    that sample, so it's dropped rather than treated as a real
    zero-intensity peak."""
    raw = pd.read_csv(feature_table_path, header=None, dtype=str)
    header_row = raw.iloc[0]
    body = raw.iloc[2:].reset_index(drop=True)
    body.columns = header_row.values
    sample_cols = [c for c in body.columns if c != "Sample"]

    rows = []
    skipped_labels = 0
    for _, row in body.iterrows():
        ppm = parse_ppm_label(row["Sample"])
        if ppm is None:
            skipped_labels += 1
            continue
        for s in sample_cols:
            val = row[s]
            if pd.isna(val) or str(val).strip() == "":
                continue
            rows.append({"sample_id": s, "nucleus": nucleus, "shift_ppm": ppm,
                         "intensity": float(val)})
    if skipped_labels:
        print(f"Skipped {skipped_labels} row(s) whose Sample label didn't start with "
              f"'ppm_X.XXXX' at all -- check --feature-table is really a "
              f"write_metaboanalyst_csv()-shaped file (annotated or not).")
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--feature-table", required=True,
                     help="write_metaboanalyst_csv()'s output -- annotated or not, "
                          "both work (see module docstring)")
    ap.add_argument("--reference", required=True,
                     help="reference_peaks*.csv -- compound_id, compound_name, "
                          "nucleus, shift_ppm columns")
    ap.add_argument("--nucleus", default="1H")
    ap.add_argument("--tolerance", type=float, default=None,
                     help="ppm tolerance; default is nucleus-specific (1H=0.008) "
                          "if not set -- see peak_matcher.py's DEFAULT_TOLERANCE_PPM")
    ap.add_argument("--top-n", type=int, default=10,
                     help="max ranked candidate compounds to keep per sample -- if you're "
                          "also passing --min-confidence, set this generously high (e.g. "
                          "500) so top-n itself never becomes the thing truncating the "
                          "list before the confidence floor gets a say")
    ap.add_argument("--min-ref-peaks", type=int, default=2,
                     help="minimum independent reference peaks a compound needs "
                          "to be considered (default 2, see peak_matcher.py's "
                          "own docstring for why)")
    ap.add_argument("--min-confidence", type=float, default=None,
                     help="keep every compound with confidence >= this value (0-1) "
                          "instead of only whatever lands in the top-n by rank -- e.g. "
                          "--min-confidence 0.88 to see every reasonably-supported match, "
                          "not just the highest-confidence handful. confidence blends "
                          "coverage, evidence_score, shift_score, and peak_prominence; "
                          "see utils.score_compound()'s docstring for the exact weights "
                          "and the rationale. Not set by default.")
    ap.add_argument("--out", default="../output/compound_matches.csv")
    args = ap.parse_args()

    experimental = wide_to_long(args.feature_table, args.nucleus)
    if experimental.empty:
        raise SystemExit("No usable peak rows found -- check --feature-table is the "
                          "UN-annotated write_metaboanalyst_csv() output.")

    reference = pd.read_csv(args.reference)
    missing_ref = {"compound_id", "compound_name", "nucleus", "shift_ppm"} - set(reference.columns)
    if missing_ref:
        raise SystemExit(f"--reference file missing columns: {missing_ref}")

    print(f"{len(experimental)} peak rows across {experimental['sample_id'].nunique()} "
          f"sample(s), reshaped from {args.feature_table}")

    results = match_all_samples(experimental, reference, args.tolerance, args.top_n,
                                 args.min_ref_peaks, args.min_confidence)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(args.out, index=False)
    print(f"Wrote {len(results)} ranked candidate-compound rows to {args.out}")

    if not results.empty:
        print("\nTop match per sample:")
        top1 = (results.sort_values(["sample_id", "confidence"], ascending=[True, False])
                        .groupby("sample_id").head(1))
        print(top1.to_string(index=False))

    # How many UNIQUE compounds you'd keep at a few confidence floors -- a
    # SEPARATE, uncapped, unfiltered pass (top_n huge, no min_confidence),
    # so this always reflects the TRUE full distribution, not just whatever
    # is left after your own --top-n/--min-confidence already trimmed it.
    # Printed every run so you can see where 0.80-0.95 actually fall without
    # re-running the matcher repeatedly to feel it out.
    full_ranked = match_all_samples(experimental, reference, args.tolerance,
                                     top_n=1_000_000, min_ref_peaks=args.min_ref_peaks,
                                     min_confidence=None)
    if not full_ranked.empty:
        print("\nUnique compounds by confidence floor (full distribution, "
              "ignoring --top-n/--min-confidence):")
        for floor in (0.95, 0.90, 0.88, 0.85, 0.80, 0.70):
            n_unique = full_ranked.loc[full_ranked["confidence"] >= floor, "compound_id"].nunique()
            print(f"  confidence >= {floor:.2f}: {n_unique} unique compound(s)")

    if results.empty:
        print("\nNo compound met min_ref_peaks and tolerance -- try --min-ref-peaks 1 "
              "or a wider --tolerance if that's too strict for your data.")


if __name__ == "__main__":
    main()
