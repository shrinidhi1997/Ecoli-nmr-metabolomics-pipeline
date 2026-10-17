"""
hsqc_matcher.py
OPTIONAL. Matches one or many experimental 2D [1H,13C]-HSQC cross-peak
lists against an HSQC reference library (built by
build_bmrb_library.py --hsqc) and ranks candidate compounds by confidence.

This is a completely separate module from peak_matcher.py -- nothing here
is imported by the 1D pipeline, and 1D-only use of this project never
touches this file.

Inputs
------
reference_hsqc.csv    compound_id, compound_name, h_ppm, c_ppm
                       -> build_bmrb_library.py --hsqc (not shipped in
                       reference_library/ by default -- run that command
                       yourself to generate it, same as the 1D BMRB/HMDB
                       builders when you want to regenerate rather than
                       use the shipped CSVs)

experimental_hsqc_peaks.csv   sample_id, h_ppm, c_ppm, intensity (optional)
                              one row per picked 2D cross-peak
                              -> hsqc_analysis.py

Output
------
match_results_hsqc.csv   sample_id, compound_id, compound_name, n_ref,
                          n_exp, n_matched, coverage, purity, shift_score,
                          confidence, sorted by sample_id then confidence

Run:
    python hsqc_matcher.py \\
        --reference ../output/reference_peaks_hsqc_bmrb_library.csv \\
        --experimental ../output/experimental_hsqc_peaks.csv \\
        --out ../output/match_results_hsqc.csv
"""

from __future__ import annotations
import argparse
from pathlib import Path

import pandas as pd

from utils import score_compound_2d

DEFAULT_TOLERANCE_H_PPM = 0.03
DEFAULT_TOLERANCE_C_PPM = 0.4


def match_sample_2d(sample_peaks: pd.DataFrame, reference: pd.DataFrame,
                     tolerance_h_ppm: float = DEFAULT_TOLERANCE_H_PPM,
                     tolerance_c_ppm: float = DEFAULT_TOLERANCE_C_PPM,
                     min_ref_peaks: int = 2) -> pd.DataFrame:
    """
    sample_peaks: rows for ONE sample, columns h_ppm, c_ppm.
    reference: full HSQC reference table (all compounds).
    Returns one scored row per candidate compound, best first.

    min_ref_peaks defaults to 2 for the same reason as peak_matcher.py's
    1D matcher: a compound with only one reference cross-peak gets
    coverage=1.0 from any single coincidental hit, which is not real
    corroborating evidence. Pass min_ref_peaks=1 to allow single-cross-peak
    compounds through deliberately.
    """
    if sample_peaks.empty:
        return pd.DataFrame()

    results = []
    for (compound_id, compound_name), ref_group in reference.groupby(
        ["compound_id", "compound_name"]
    ):
        if len(ref_group) < min_ref_peaks:
            continue

        scores = score_compound_2d(
            experimental_h=sample_peaks["h_ppm"].to_numpy(),
            experimental_c=sample_peaks["c_ppm"].to_numpy(),
            reference_h=ref_group["h_ppm"].to_numpy(),
            reference_c=ref_group["c_ppm"].to_numpy(),
            tolerance_h_ppm=tolerance_h_ppm,
            tolerance_c_ppm=tolerance_c_ppm,
        )
        if scores["n_matched"] == 0:
            continue
        results.append({
            "compound_id": compound_id,
            "compound_name": compound_name,
            **scores,
        })

    df = pd.DataFrame(results)
    if not df.empty:
        df = df.sort_values("confidence", ascending=False).reset_index(drop=True)
    return df


def match_all_samples_2d(experimental: pd.DataFrame, reference: pd.DataFrame,
                          tolerance_h_ppm: float = DEFAULT_TOLERANCE_H_PPM,
                          tolerance_c_ppm: float = DEFAULT_TOLERANCE_C_PPM,
                          top_n: int = 10, min_ref_peaks: int = 2) -> pd.DataFrame:
    all_rows = []
    for sample_id, sample_peaks in experimental.groupby("sample_id"):
        ranked = match_sample_2d(sample_peaks, reference, tolerance_h_ppm,
                                  tolerance_c_ppm, min_ref_peaks)
        if ranked.empty:
            continue
        ranked.insert(0, "sample_id", sample_id)
        all_rows.append(ranked.head(top_n))
    return pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--experimental", required=True)
    ap.add_argument("--tolerance-h", type=float, default=DEFAULT_TOLERANCE_H_PPM)
    ap.add_argument("--tolerance-c", type=float, default=DEFAULT_TOLERANCE_C_PPM)
    ap.add_argument("--top-n", type=int, default=10)
    ap.add_argument("--min-ref-peaks", type=int, default=2,
                     help="minimum independent reference cross-peaks a compound needs "
                          "to be considered (default 2, same reasoning as peak_matcher.py)")
    ap.add_argument("--out", default="../output/match_results_hsqc.csv")
    args = ap.parse_args()

    reference = pd.read_csv(args.reference)
    experimental = pd.read_csv(args.experimental)

    required_ref_cols = {"compound_id", "compound_name", "h_ppm", "c_ppm"}
    required_exp_cols = {"sample_id", "h_ppm", "c_ppm"}
    missing_ref = required_ref_cols - set(reference.columns)
    missing_exp = required_exp_cols - set(experimental.columns)
    if missing_ref:
        raise SystemExit(f"reference file missing columns: {missing_ref}")
    if missing_exp:
        raise SystemExit(f"experimental file missing columns: {missing_exp}")

    results = match_all_samples_2d(experimental, reference, args.tolerance_h,
                                    args.tolerance_c, args.top_n, args.min_ref_peaks)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(args.out, index=False)
    print(f"Wrote {len(results)} ranked HSQC hits to {args.out}")
    if not results.empty:
        print(results.groupby("sample_id").head(3).to_string(index=False))


if __name__ == "__main__":
    main()
