"""
nmrmatch.matching
(originally peak_matcher.py)
Matches one or many experimental NMR peak lists against a reference
compound peak library and ranks candidate compound IDs by confidence.

Inputs
------
reference_peaks.csv   compound_id, compound_name, nucleus, shift_ppm,
                       intensity (optional), multiplicity (optional)
                       -> see fetch_nmr_reference.py / build_bmrb_library.py
                       for how to build this

experimental_peaks.csv  sample_id, nucleus, shift_ppm, intensity (optional)
                         one row per detected peak; a sample_id can repeat
                         across many rows (one row per peak)

Output
------
match_results.csv   sample_id, compound_id, compound_name, n_ref, n_exp,
                     n_matched, coverage, purity, shift_score,
                     intensity_corr, confidence
                     sorted by sample_id then confidence (desc)

Run:
    python peak_matcher.py \
        --reference ../output/reference_peaks.csv \
        --experimental ../example_data/experimental_peaks_example.csv \
        --tolerance 0.03 \
        --out ../output/match_results.csv
"""

from __future__ import annotations
import argparse
from pathlib import Path

import pandas as pd

from .scoring import score_compound

DEFAULT_TOLERANCE_PPM = {"1H": 0.03, "13C": 0.5}


def match_sample(sample_peaks: pd.DataFrame, reference: pd.DataFrame,
                  tolerance_ppm: float | None, min_ref_peaks: int = 2) -> pd.DataFrame:
    """
    sample_peaks: rows for ONE sample, columns nucleus, shift_ppm, [intensity]
    reference: full reference table (all compounds, all nuclei)
    Returns one scored row per candidate compound, best first.

    min_ref_peaks defaults to 2, not 1: a compound with only one reference
    peak gets coverage=1.0 from any single coincidental match, regardless
    of how many other compounds in a large library could equally well
    explain that one peak -- confirmed in practice against a ~3,300-
    compound library, where every top hit was a different 1-peak compound
    explaining under 3% of the total spectrum each. Requiring at least two
    independent reference peaks means a match only scores well when
    multiple lines land together, real corroborating evidence rather than
    one lucky coincidence. Pass min_ref_peaks=1 to disable this if you
    specifically want single-peak compounds considered (legitimate for
    compounds that only have one real proton environment, e.g. formate).
    """
    results = []
    for (compound_id, compound_name, nucleus), ref_group in reference.groupby(
        ["compound_id", "compound_name", "nucleus"]
    ):
        if len(ref_group) < min_ref_peaks:
            continue
        exp_group = sample_peaks[sample_peaks["nucleus"] == nucleus]
        if exp_group.empty:
            continue

        tol = tolerance_ppm or DEFAULT_TOLERANCE_PPM.get(nucleus, 0.03)
        exp_intensity = exp_group["intensity"].to_numpy() if "intensity" in exp_group and exp_group["intensity"].notna().all() else None
        ref_intensity = ref_group["intensity"].to_numpy() if "intensity" in ref_group and ref_group["intensity"].notna().all() else None

        scores = score_compound(
            experimental_ppm=exp_group["shift_ppm"].to_numpy(),
            experimental_intensity=exp_intensity,
            reference_ppm=ref_group["shift_ppm"].to_numpy(),
            reference_intensity=ref_intensity,
            tolerance_ppm=tol,
        )
        if scores["n_matched"] == 0:
            continue
        results.append({
            "compound_id": compound_id,
            "compound_name": compound_name,
            "nucleus": nucleus,
            **scores,
        })

    df = pd.DataFrame(results)
    if not df.empty:
        df = df.sort_values("confidence", ascending=False).reset_index(drop=True)
    return df


def match_all_samples(experimental: pd.DataFrame, reference: pd.DataFrame,
                       tolerance_ppm: float | None, top_n: int = 10,
                       min_ref_peaks: int = 2,
                       min_confidence: float | None = None) -> pd.DataFrame:
    """
    min_confidence, when set, is applied AFTER the top_n rank cutoff: keep
    only the rows (among that sample's top_n ranked candidates) whose
    confidence is >= min_confidence. This is why match_r_feature_table.py's
    --min-confidence help text says to set --top-n generously high (e.g.
    500) alongside it -- with the default top_n=10, a real compound that
    clears min_confidence could still get cut here first, before
    min_confidence ever gets a say, if 10+ other compounds outrank it.
    """
    all_rows = []
    for sample_id, sample_peaks in experimental.groupby("sample_id"):
        ranked = match_sample(sample_peaks, reference, tolerance_ppm, min_ref_peaks)
        ranked = ranked.head(top_n)
        if min_confidence is not None and not ranked.empty:
            ranked = ranked[ranked["confidence"] >= min_confidence]
        ranked.insert(0, "sample_id", sample_id)
        all_rows.append(ranked)
    return pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--experimental", required=True)
    ap.add_argument("--tolerance", type=float, default=None,
                     help="ppm tolerance; default is nucleus-specific "
                          "(1H=0.03, 13C=0.5) if not set")
    ap.add_argument("--top-n", type=int, default=10,
                     help="max candidate compounds to keep per sample")
    ap.add_argument("--min-ref-peaks", type=int, default=2,
                     help="minimum independent reference peaks a compound needs to be "
                          "considered (default 2 -- filters out trivial single-peak "
                          "matches, which get artificially high confidence from any one "
                          "coincidental hit; pass 1 to allow single-peak compounds)")
    ap.add_argument("--out", default="../output/match_results.csv")
    args = ap.parse_args()

    reference = pd.read_csv(args.reference)
    experimental = pd.read_csv(args.experimental)

    required_ref_cols = {"compound_id", "compound_name", "nucleus", "shift_ppm"}
    required_exp_cols = {"sample_id", "nucleus", "shift_ppm"}
    missing_ref = required_ref_cols - set(reference.columns)
    missing_exp = required_exp_cols - set(experimental.columns)
    if missing_ref:
        raise SystemExit(f"reference file missing columns: {missing_ref}")
    if missing_exp:
        raise SystemExit(f"experimental file missing columns: {missing_exp}")

    results = match_all_samples(experimental, reference, args.tolerance, args.top_n,
                                 args.min_ref_peaks)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(args.out, index=False)
    print(f"Wrote {len(results)} ranked hits to {args.out}")
    if not results.empty:
        print(results.groupby("sample_id").head(3).to_string(index=False))


if __name__ == "__main__":
    main()
