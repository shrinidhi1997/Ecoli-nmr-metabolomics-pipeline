"""
compute_fold_induction.py

Turns match_r_feature_table.py's per-sample compound-match output into a
per-compound fold-change table between two growth-condition groups
(e.g. control vs sample), using the matched_intensity_mean column utils.py's
score_compound() now reports (the mean of a sample's own experimental
intensities at that compound's matched peaks -- see that function's
docstring for exactly what it is and isn't).

Why this is a SEPARATE script rather than folded into match_r_feature_table.py
-------------------------------------------------------------------------------
Matching (is this compound present, and how confidently) and quantifying
(how much of it is there, relative to another condition) are different
questions with different failure modes -- a compound can be identified with
high confidence in every sample and still have no reliable fold-change if
it wasn't both detected AND kept inside --top-n in enough samples of both
groups to compare. Keeping this separate means you can re-run just this
cheap step after changing group assignments or which match file to use,
without re-running the matcher.

What this does NOT do (read before trusting a number):
  * It is NOT true qNMR absolute quantification. matched_intensity_mean is
    whatever your experimental "intensity" column held upstream -- if that
    was find_shape_based_peaks()'s raw peak height, a fold-change here can
    be distorted by linewidth differences between samples (shimming,
    viscosity, temperature) even when the real molar amount is unchanged.
    If you ran refine_peaks_lorentzian() and fed its "area" column through
    as "intensity" into merge_peaks_reference_guided() before building the
    feature table, this fold-change is the textbook-correct qNMR quantity
    instead. Which one you actually have depends on your R run -- this
    script can't tell the difference, so that choice is still on you.
  * A compound only gets a fold-change here if match_r_feature_table.py's
    output has an intensity value for it in at least one sample of EACH
    group. The output's n_<group_name>_detected columns (named from
    whatever --group names you pass in) tell you how many of each
    group's samples actually contributed -- a compound built from 1 sample
    per group is much weaker evidence than one built from every replicate,
    even though both get a number. Filter on those two columns yourself
    before trusting any single row.
  * Averaging matched_intensity_mean is a mean over whichever reference
    peaks happened to match in that one sample -- it is NOT guaranteed to
    be the same physical peak across every sample if coverage differs
    sample to sample (e.g. compound X matches 2 of its 3 reference peaks in
    one sample but all 3 in another). coverage_min/coverage_max per group
    are reported so you can see when that's happened, rather than silently
    averaging over mismatched evidence.

Run (replace the sample names below with your own; they must match
sample_group_map in config_control_vs_sample.R exactly):
    python compute_fold_induction.py \\
        --matches ../output/control_vs_sample/compound_matches_control_vs_sample.csv \\
        --group control:sample_01,sample_02,sample_03 \\
        --group sample:sample_04,sample_05,sample_06 \\
        --out ../output/control_vs_sample/fold_induction_control_vs_sample.csv
"""
from __future__ import annotations
import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def parse_group_arg(spec: str) -> tuple[str, list[str]]:
    """'control:sample_01,sample_02' -> ('control', ['sample_01', 'sample_02'])"""
    if ":" not in spec:
        raise SystemExit(f"--group '{spec}' must be 'group_name:sample1,sample2,...'")
    name, rest = spec.split(":", 1)
    samples = [s.strip() for s in rest.split(",") if s.strip()]
    if not samples:
        raise SystemExit(f"--group '{spec}' listed no sample IDs")
    return name.strip(), samples


def compute_fold_induction(matches: pd.DataFrame, group_a: tuple[str, list[str]],
                            group_b: tuple[str, list[str]],
                            min_samples_per_group: int = 1) -> pd.DataFrame:
    a_name, a_samples = group_a
    b_name, b_samples = group_b

    missing = (set(a_samples) | set(b_samples)) - set(matches["sample_id"].unique())
    if missing:
        print(f"Warning: {len(missing)} sample(s) named in --group aren't in the match file "
              f"at all (typo, or that sample had zero matches): {sorted(missing)}")

    if "matched_intensity_mean" not in matches.columns:
        raise SystemExit(
            "No matched_intensity_mean column in --matches. Re-run match_r_feature_table.py "
            "with the updated peak_matcher.py/utils.py (score_compound() now returns this "
            "field) -- an older match file won't have it."
        )

    rows = []
    for (compound_id, compound_name), grp in matches.groupby(["compound_id", "compound_name"]):
        a_rows = grp[grp["sample_id"].isin(a_samples)]
        b_rows = grp[grp["sample_id"].isin(b_samples)]
        a_vals = a_rows["matched_intensity_mean"].dropna()
        b_vals = b_rows["matched_intensity_mean"].dropna()

        if len(a_vals) < min_samples_per_group or len(b_vals) < min_samples_per_group:
            continue

        a_mean = float(a_vals.mean())
        b_mean = float(b_vals.mean())

        if a_mean <= 0 or b_mean <= 0:
            fold_change = None
            log2fc = None
        else:
            fold_change = b_mean / a_mean
            log2fc = float(np.log2(fold_change))

        rows.append({
            "compound_id": compound_id,
            "compound_name": compound_name,
            f"{a_name}_mean_intensity": round(a_mean, 6),
            f"{b_name}_mean_intensity": round(b_mean, 6),
            f"n_{a_name}_detected": len(a_vals),
            f"n_{a_name}_total": len(a_samples),
            f"n_{b_name}_detected": len(b_vals),
            f"n_{b_name}_total": len(b_samples),
            "coverage_min": round(float(grp["coverage"].min()), 4) if "coverage" in grp else None,
            "coverage_max": round(float(grp["coverage"].max()), 4) if "coverage" in grp else None,
            "fold_change": round(fold_change, 4) if fold_change is not None else None,
            "log2_fold_change": round(log2fc, 4) if log2fc is not None else None,
        })

    df = pd.DataFrame(rows)
    if not df.empty:
        # most differential first -- biggest |log2FC| is what a reader
        # actually wants to see first; undefined fold-changes (None) sort last
        df["_abs_log2fc"] = df["log2_fold_change"].abs()
        df = df.sort_values("_abs_log2fc", ascending=False, na_position="last").drop(columns="_abs_log2fc")
        df = df.reset_index(drop=True)
    return df


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--matches", required=True,
                     help="match_r_feature_table.py's output -- must have been run with "
                          "the updated peak_matcher.py so matched_intensity_mean is present")
    ap.add_argument("--group", action="append", required=True, dest="groups",
                     help="'group_name:sample1,sample2,...' -- pass exactly twice, once per "
                          "group (e.g. --group control:... --group sample:...)")
    ap.add_argument("--min-samples-per-group", type=int, default=1,
                     help="a compound needs at least this many samples WITH a matched "
                          "intensity in EACH group to get a fold-change at all (default 1 -- "
                          "raise to 2 or 3 if you only trust a fold-change backed by most/all "
                          "of a group's replicates)")
    ap.add_argument("--out", default="../output/fold_induction.csv")
    args = ap.parse_args()

    if len(args.groups) != 2:
        raise SystemExit(f"Expected exactly 2 --group arguments, got {len(args.groups)}")
    group_a = parse_group_arg(args.groups[0])
    group_b = parse_group_arg(args.groups[1])

    matches = pd.read_csv(args.matches)
    result = compute_fold_induction(matches, group_a, group_b, args.min_samples_per_group)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.out, index=False)
    print(f"Wrote {len(result)} compound(s) with a fold-change between "
          f"'{group_a[0]}' and '{group_b[0]}' to {args.out}")
    if not result.empty:
        print(f"\nTop 10 by |log2 fold change| ({group_b[0]} / {group_a[0]}):")
        print(result.head(10).to_string(index=False))
    else:
        print("\nNo compound had matched_intensity_mean in enough samples of both groups. "
              "Check --matches actually has the matched_intensity_mean column populated "
              "(needs real intensity values in the experimental peak table, not blanks), "
              "and that your --group sample IDs match the match file exactly.")


if __name__ == "__main__":
    main()
