"""
analyze_standard_dilution_linearity.py

Checks whether one compound's matched peak intensity scales linearly with
known concentration across a standard dilution series -- the basic sanity
check you'd want before trusting this pipeline's peak picking and matching
on real standards at different concentrations.

This is a GENERIC script, not alanine-specific: it only needs
match_r_feature_table.py's own output (for --matches) plus a small
sample_id -> concentration lookup CSV (for --concentrations). It is meant
to be reused as-is once you have real standards data, not just for the
simulated demo it ships with.

Input
-----
--matches         match_r_feature_table.py's output CSV. Needs
                   sample_id, compound_id, compound_name,
                   matched_intensity_mean, n_distinct_matched, coverage,
                   confidence columns (all present in that script's
                   normal output).
--concentrations  a CSV with sample_id, concentration_mM columns -- ground
                   truth, not something this pipeline measures itself.
--compound-id     which compound_id to check (must appear in --matches).

What this does NOT do: it is not a replacement for compute_fold_induction.py
(which compares two discrete groups) or for a real calibration-curve LOD/LOQ
analysis. It is a quick, honest linearity check: ordinary least-squares
fit of matched_intensity_mean against concentration, with R-squared,
reported alongside each point's own n_distinct_matched/coverage so a drop
in coverage at low concentration (common, and not a bug -- see
compute_fold_induction.py's own docstring on why coverage can differ
sample to sample) is visible next to the number, not hidden behind it.

Run:
    python analyze_standard_dilution_linearity.py \
        --matches ../example_results/1d/compound_matches_alanine_dilution.csv \
        --concentrations ../example_data/1d/concentrations.csv \
        --compound-id HMDB0000161 \
        --out ../example_results/1d/alanine_dilution_linearity.csv
"""
from __future__ import annotations
import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def compute_linearity(matches: pd.DataFrame, concentrations: pd.DataFrame,
                       compound_id: str) -> tuple[pd.DataFrame, dict]:
    sub = matches[matches["compound_id"] == compound_id].copy()
    if sub.empty:
        raise SystemExit(f"compound_id '{compound_id}' not found in --matches at all.")

    merged = sub.merge(concentrations, on="sample_id", how="inner")
    missing = set(sub["sample_id"]) - set(merged["sample_id"])
    if missing:
        print(f"Warning: {len(missing)} sample(s) in --matches have no row in "
              f"--concentrations, dropped: {sorted(missing)}")

    merged = merged.dropna(subset=["matched_intensity_mean", "concentration_mM"])
    merged = merged.sort_values("concentration_mM").reset_index(drop=True)

    cols = ["sample_id", "concentration_mM", "matched_intensity_mean",
            "n_distinct_matched", "coverage", "confidence"]
    result = merged[cols].copy()

    stats = {"n_points": len(merged), "slope": None, "intercept": None,
              "r_squared": None, "pearson_r": None}
    if len(merged) >= 2:
        x = merged["concentration_mM"].to_numpy(dtype=float)
        y = merged["matched_intensity_mean"].to_numpy(dtype=float)
        slope, intercept = np.polyfit(x, y, 1)
        y_pred = slope * x + intercept
        ss_res = float(np.sum((y - y_pred) ** 2))
        ss_tot = float(np.sum((y - y.mean()) ** 2))
        r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
        pearson_r = float(np.corrcoef(x, y)[0, 1]) if np.std(x) > 0 and np.std(y) > 0 else float("nan")
        stats.update(slope=round(float(slope), 6), intercept=round(float(intercept), 6),
                      r_squared=round(r_squared, 6), pearson_r=round(pearson_r, 6))
    return result, stats


def plot_linearity(result: pd.DataFrame, stats: dict, compound_name: str, out_png: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6, 5))
    x = result["concentration_mM"].to_numpy(dtype=float)
    y = result["matched_intensity_mean"].to_numpy(dtype=float)
    ax.scatter(x, y, color="#2f6690", zorder=5)
    for _, row in result.iterrows():
        ax.annotate(f"{row['n_distinct_matched']}/6 peaks", (row["concentration_mM"], row["matched_intensity_mean"]),
                     textcoords="offset points", xytext=(6, 6), fontsize=8, color="#555")

    if stats["slope"] is not None:
        xs = np.linspace(0, x.max() * 1.05, 50)
        ax.plot(xs, stats["slope"] * xs + stats["intercept"], color="#d9534f",
                linestyle="--", label=f"fit: R²={stats['r_squared']:.3f}")
        ax.legend()

    ax.set_xlabel("concentration (mM)")
    ax.set_ylabel("matched_intensity_mean (a.u.)")
    ax.set_title(f"{compound_name} -- standard dilution linearity")
    fig.tight_layout()
    fig.savefig(out_png, dpi=130)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--matches", required=True)
    ap.add_argument("--concentrations", required=True)
    ap.add_argument("--compound-id", required=True)
    ap.add_argument("--out", default="../example_results/dilution_linearity.csv")
    args = ap.parse_args()

    matches = pd.read_csv(args.matches)
    concentrations = pd.read_csv(args.concentrations)
    result, stats = compute_linearity(matches, concentrations, args.compound_id)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(out_path, index=False)
    print(f"Wrote {len(result)} concentration point(s) to {out_path}")
    print(result.to_string(index=False))
    print(f"\nLinear fit: slope={stats['slope']}, intercept={stats['intercept']}, "
          f"R-squared={stats['r_squared']}, Pearson r={stats['pearson_r']}")
    if result["n_distinct_matched"].nunique() > 1:
        print("\nNote: n_distinct_matched differs across concentrations -- the lower "
              "point(s) are averaged over fewer matched peaks than the higher ones. "
              "This is a real, honest difference in coverage, not a computation error; "
              "see compute_fold_induction.py's own docstring for why that happens and "
              "why it can bias a naive average. Treat the fit above as a first pass, "
              "and compare it against a version restricted to peaks matched at every "
              "concentration if you need a strictly apples-to-apples slope.")

    compound_name = matches.loc[matches["compound_id"] == args.compound_id, "compound_name"].iloc[0]
    plot_linearity(result, stats, compound_name, out_path.with_suffix(".png"))
    print(f"Wrote {out_path.with_suffix('.png')}")


if __name__ == "__main__":
    main()
