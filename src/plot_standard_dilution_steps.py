"""
plot_standard_dilution_steps.py

The "many graphs for each step, how it's performing" figures for a real
N-point standard dilution series, reading only the CSVs
pipeline_standard_dilution_core.R and match_r_feature_table.py already
wrote to ../example_results/<dataset>_real_dilution/ -- this script does
no peak picking or matching of its own, it only visualizes what the real
pipeline already produced. Generalized (dataset name + sample list passed
in) so the same code runs on haird_real_dilution, swepf_real_dilution, or
any future standard series.

Run from this src/ folder:
    python plot_standard_dilution_steps.py haird HAIRD-30-1H HAIRD-60-1H HAIRD-120-1H
    python plot_standard_dilution_steps.py swepf SWEPF-30-1H SWEPF-60-1H SWEPF-120-1H
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PALETTE = ["#2f6690", "#d9534f", "#5a8f3c", "#8e5aa3"]


def load_sample_csvs(res_dir: Path, sample_ids: list[str], suffix: str) -> dict[str, pd.DataFrame]:
    return {s: pd.read_csv(res_dir / f"{s}_{suffix}.csv") for s in sample_ids}


def step1_raw_overlay(res_dir: Path, sample_ids: list[str], colors: dict):
    raw = load_sample_csvs(res_dir, sample_ids, "raw_spectrum")
    fig, axes = plt.subplots(2, 1, figsize=(10, 8))
    for s in sample_ids:
        d = raw[s]
        axes[0].plot(d["ppm"], d["intensity"], color=colors[s], label=s, linewidth=0.6)
    axes[0].set_xlim(11, 0)
    axes[0].set_xlabel("ppm")
    axes[0].set_ylabel("raw intensity (a.u.)")
    axes[0].set_title("Step 0/1: raw imported spectra, 0-11 ppm\n(the reporting range -- real 1H spectra don't have anything past this)")
    axes[0].legend(fontsize=8)

    for s in sample_ids:
        d = raw[s]
        mask = (d["ppm"] >= 0.5) & (d["ppm"] <= 4.5)
        axes[1].plot(d.loc[mask, "ppm"], d.loc[mask, "intensity"], color=colors[s], label=s, linewidth=0.8)
    axes[1].set_xlim(4.5, 0.5)
    axes[1].set_xlabel("ppm")
    axes[1].set_ylabel("raw intensity (a.u.)")
    axes[1].set_title("Same spectra, zoomed to the aliphatic/amino-acid region (0.5-4.5 ppm)")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    out = res_dir / "step1_raw_spectra_overlay.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"Wrote {out}")


def step2_baseline_example(res_dir: Path, sample_id: str):
    raw = pd.read_csv(res_dir / f"{sample_id}_raw_spectrum.csv")
    base = pd.read_csv(res_dir / f"{sample_id}_baseline.csv")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for ax, (lo, hi) in zip(axes, [(0.5, 4.5), (6.5, 9.5)]):
        mask = (raw["ppm"] >= lo) & (raw["ppm"] <= hi)
        ax.plot(raw.loc[mask, "ppm"], raw.loc[mask, "intensity"], color="#333333", linewidth=0.6, label="raw")
        ax.plot(base.loc[mask, "ppm"], base.loc[mask, "baseline"], color="#d9534f", linewidth=1.2, label="airPLS baseline")
        ax.set_xlim(hi, lo)
        ax.set_xlabel("ppm")
        ax.legend(fontsize=8)
    axes[0].set_ylabel("intensity (a.u.)")
    fig.suptitle(f"Step 1: airPLS baseline fit on real data ({sample_id})")
    fig.tight_layout()
    out = res_dir / "step2_baseline_fit_example.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"Wrote {out}")


def step3_detected_peaks(res_dir: Path, sample_ids: list[str], colors: dict):
    corrected = load_sample_csvs(res_dir, sample_ids, "corrected_spectrum")
    peaks = load_sample_csvs(res_dir, sample_ids, "detected_peaks")
    fig, axes = plt.subplots(len(sample_ids), 1, figsize=(10, 3.3 * len(sample_ids)), sharex=True)
    axes = np.atleast_1d(axes)
    for ax, s in zip(axes, sample_ids):
        d, p = corrected[s], peaks[s]
        mask = (d["ppm"] >= 0.5) & (d["ppm"] <= 4.5)
        ax.plot(d.loc[mask, "ppm"], d.loc[mask, "corrected"], color=colors[s], linewidth=0.6)
        pm = (p["ppm"] >= 0.5) & (p["ppm"] <= 4.5)
        ax.scatter(p.loc[pm, "ppm"], p.loc[pm, "intensity"], color="black", s=14, zorder=5, label=f"{pm.sum()} peaks (shape detector)")
        ax.set_xlim(4.5, 0.5)
        ax.set_ylabel("corrected (a.u.)")
        ax.set_title(s)
        ax.legend(fontsize=8, loc="upper right")
    axes[-1].set_xlabel("ppm")
    fig.suptitle("Step 2: baseline-corrected spectra with detected peaks (aliphatic region, recommended shape detector)")
    fig.tight_layout()
    out = res_dir / "step3_detected_peaks_aliphatic.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"Wrote {out}")


def step4_detector_comparison(res_dir: Path):
    df = pd.read_csv(res_dir / "detector_comparison.csv")
    x = np.arange(len(df))
    width = 0.35
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.bar(x - width / 2, df["n_shape_only"], width, label="shape detector only (used for the feature table)", color="#2f6690")
    ax.bar(x + width / 2, df["n_combined_with_deriv"], width, label="+ second-derivative pass (NOT used)", color="#d9534f")
    ax.set_xticks(x)
    ax.set_xticklabels(df["sample_id"])
    ax.set_ylabel("peaks detected")
    ax.set_title("Step 2: real data reproduces the second-derivative detector's\nknown noise sensitivity -- not used downstream")
    ax.legend(fontsize=8)
    for i, (a, b) in enumerate(zip(df["n_shape_only"], df["n_combined_with_deriv"])):
        ax.annotate(str(a), (i - width / 2, a), textcoords="offset points", xytext=(0, 4), ha="center", fontsize=8)
        ax.annotate(str(b), (i + width / 2, b), textcoords="offset points", xytext=(0, 4), ha="center", fontsize=8)
    fig.tight_layout()
    out = res_dir / "step4_detector_comparison.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"Wrote {out}")


def step5_merged_feature_heatmap(res_dir: Path, sample_ids: list[str], feature_table_csv: str):
    # write_metaboanalyst_csv()'s own shape: first column "Sample" holds
    # "Label" (group names) in row 0, then one "ppm_<value>" row per
    # merged feature; the sample columns hold that feature's intensity,
    # blank where that sample had no peak within the merge tolerance.
    raw = pd.read_csv(res_dir / feature_table_csv)
    feat = raw[raw["Sample"] != "Label"].copy()
    feat["ppm"] = feat["Sample"].str.replace("ppm_", "", regex=False).astype(float)
    feat = feat.sort_values("ppm", ascending=False).reset_index(drop=True)
    mat = feat[sample_ids].apply(pd.to_numeric, errors="coerce").to_numpy()

    fig, ax = plt.subplots(figsize=(6, 10))
    plot_mat = np.log10(np.where(mat > 0, mat, np.nan))
    im = ax.imshow(plot_mat, aspect="auto", cmap="viridis",
                    extent=[0, len(sample_ids), len(feat), 0])
    ax.set_xticks(np.arange(len(sample_ids)) + 0.5)
    ax.set_xticklabels(sample_ids, rotation=20, ha="right")
    step = max(1, len(feat) // 40)
    ytick_idx = np.arange(0, len(feat), step)
    ax.set_yticks(ytick_idx + 0.5)
    ax.set_yticklabels([f"{feat.loc[i, 'ppm']:.3f}" for i in ytick_idx], fontsize=7)
    ax.set_ylabel("ppm (merged feature)")
    ax.set_title(f"Step 3: merged feature table\n({len(feat)} features x {len(sample_ids)} samples, log10 intensity, white = not detected)")
    cbar = fig.colorbar(im, ax=ax, fraction=0.04)
    cbar.set_label("log10(intensity)")
    fig.tight_layout()
    out = res_dir / "step5_merged_feature_heatmap.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    n_full = int((mat > 0).all(axis=1).sum())
    print(f"Wrote {out} ({n_full}/{len(feat)} features detected in all {len(sample_ids)} samples)")


def step6_confidence_distribution(res_dir: Path, sample_ids: list[str], colors: dict, matches_csv: str):
    df = pd.read_csv(res_dir / matches_csv)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    axes[0].hist(df["confidence"], bins=30, color="#2f6690")
    axes[0].set_xlabel("confidence")
    axes[0].set_ylabel("candidate-compound rows")
    axes[0].set_title(f"Step 4: confidence distribution,\nall candidates, all {len(sample_ids)} samples")

    for s in sample_ids:
        sub = df[df["sample_id"] == s]
        axes[1].hist(sub["confidence"], bins=20, alpha=0.5, label=s, color=colors[s])
    axes[1].set_xlabel("confidence")
    axes[1].set_title("Same, split by sample")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    out = res_dir / "step6_confidence_distribution.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"Wrote {out}")


def step7_top_matches_bar(res_dir: Path, sample_ids: list[str], matches_csv: str, top_n: int = 15):
    df = pd.read_csv(res_dir / matches_csv)
    counts = df.groupby("compound_id")["sample_id"].nunique()
    in_all = counts[counts == len(sample_ids)].index
    sub = df[df["compound_id"].isin(in_all)]
    agg = (sub.groupby(["compound_id", "compound_name"])["confidence"].mean()
           .sort_values(ascending=False).head(top_n))
    fig, ax = plt.subplots(figsize=(8, 6))
    labels = [f"{name}\n({cid})" for cid, name in agg.index]
    ax.barh(range(len(agg)), agg.values, color="#2f6690")
    ax.set_yticks(range(len(agg)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel(f"mean confidence across all {len(sample_ids)} samples")
    ax.set_title(f"Step 4: top {min(top_n, len(agg))} compounds matched in ALL {len(sample_ids)} real samples\n(fixed evidence_score scorer)")
    fig.tight_layout()
    out = res_dir / "step7_top_matches_all_samples.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"Wrote {out}")
    return agg


if __name__ == "__main__":
    dataset = sys.argv[1]
    sample_ids = sys.argv[2:]
    res_dir = Path(f"../example_results/{dataset}_real_dilution")
    colors = {s: PALETTE[i % len(PALETTE)] for i, s in enumerate(sample_ids)}
    feature_table_csv = f"peak_list_{dataset}_real_dilution.csv"
    matches_csv = f"compound_matches_{dataset}_real_dilution.csv"

    step1_raw_overlay(res_dir, sample_ids, colors)
    step2_baseline_example(res_dir, sample_ids[-1])
    step3_detected_peaks(res_dir, sample_ids, colors)
    step4_detector_comparison(res_dir)
    step5_merged_feature_heatmap(res_dir, sample_ids, feature_table_csv)
    step6_confidence_distribution(res_dir, sample_ids, colors, matches_csv)
    step7_top_matches_bar(res_dir, sample_ids, matches_csv)
