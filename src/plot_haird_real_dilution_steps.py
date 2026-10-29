"""
plot_haird_real_dilution_steps.py

Builds the "many graphs for each step, how it's performing" figures for
the three real HAIRD standard samples (HAIRD-30/60/120-1H), reading only
the CSVs run_pipeline_haird_real_dilution.R and match_r_feature_table.py
already wrote to ../example_results/haird_real_dilution/ -- this script
does no peak picking or matching of its own, it only visualizes what the
real pipeline already produced.

Run from this src/ folder:
    python plot_haird_real_dilution_steps.py
"""
from __future__ import annotations
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RES = Path("../example_results/haird_real_dilution")
SAMPLE_IDS = ["HAIRD-30-1H", "HAIRD-60-1H", "HAIRD-120-1H"]
COLORS = {"HAIRD-30-1H": "#2f6690", "HAIRD-60-1H": "#d9534f", "HAIRD-120-1H": "#5a8f3c"}


def load_sample_csvs(suffix: str) -> dict[str, pd.DataFrame]:
    return {s: pd.read_csv(RES / f"{s}_{suffix}.csv") for s in SAMPLE_IDS}


def step1_raw_overlay():
    raw = load_sample_csvs("raw_spectrum")
    fig, axes = plt.subplots(2, 1, figsize=(10, 8))
    for s in SAMPLE_IDS:
        d = raw[s]
        axes[0].plot(d["ppm"], d["intensity"], color=COLORS[s], label=s, linewidth=0.6)
    axes[0].set_xlim(16, -4)
    axes[0].set_xlabel("ppm")
    axes[0].set_ylabel("raw intensity (a.u.)")
    axes[0].set_title("Step 0/1: raw imported spectra, full range")
    axes[0].legend(fontsize=8)

    for s in SAMPLE_IDS:
        d = raw[s]
        mask = (d["ppm"] >= 0.5) & (d["ppm"] <= 4.5)
        axes[1].plot(d.loc[mask, "ppm"], d.loc[mask, "intensity"], color=COLORS[s], label=s, linewidth=0.8)
    axes[1].set_xlim(4.5, 0.5)
    axes[1].set_xlabel("ppm")
    axes[1].set_ylabel("raw intensity (a.u.)")
    axes[1].set_title("Same spectra, zoomed to the aliphatic/amino-acid region (0.5-4.5 ppm)")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(RES / "step1_raw_spectra_overlay.png", dpi=130)
    plt.close(fig)
    print("Wrote step1_raw_spectra_overlay.png")


def step2_baseline_example(sample_id: str = "HAIRD-120-1H"):
    raw = pd.read_csv(RES / f"{sample_id}_raw_spectrum.csv")
    base = pd.read_csv(RES / f"{sample_id}_baseline.csv")
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
    fig.savefig(RES / "step2_baseline_fit_example.png", dpi=130)
    plt.close(fig)
    print("Wrote step2_baseline_fit_example.png")


def step3_detected_peaks():
    corrected = load_sample_csvs("corrected_spectrum")
    peaks = load_sample_csvs("detected_peaks")
    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
    for ax, s in zip(axes, SAMPLE_IDS):
        d, p = corrected[s], peaks[s]
        mask = (d["ppm"] >= 0.5) & (d["ppm"] <= 4.5)
        ax.plot(d.loc[mask, "ppm"], d.loc[mask, "corrected"], color=COLORS[s], linewidth=0.6)
        pm = (p["ppm"] >= 0.5) & (p["ppm"] <= 4.5)
        ax.scatter(p.loc[pm, "ppm"], p.loc[pm, "intensity"], color="black", s=14, zorder=5, label=f"{pm.sum()} peaks (shape detector)")
        ax.set_xlim(4.5, 0.5)
        ax.set_ylabel("corrected (a.u.)")
        ax.set_title(s)
        ax.legend(fontsize=8, loc="upper right")
    axes[-1].set_xlabel("ppm")
    fig.suptitle("Step 2: baseline-corrected spectra with detected peaks (aliphatic region, recommended shape detector)")
    fig.tight_layout()
    fig.savefig(RES / "step3_detected_peaks_aliphatic.png", dpi=130)
    plt.close(fig)
    print("Wrote step3_detected_peaks_aliphatic.png")


def step4_detector_comparison():
    df = pd.read_csv(RES / "detector_comparison.csv")
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
    fig.savefig(RES / "step4_detector_comparison.png", dpi=130)
    plt.close(fig)
    print("Wrote step4_detector_comparison.png")


def step5_merged_feature_heatmap():
    # write_metaboanalyst_csv()'s own shape: first column "Sample" holds
    # "Label" (group names) in row 0, then one "ppm_<value>" row per
    # merged feature; the sample columns hold that feature's intensity,
    # blank where that sample had no peak within the merge tolerance.
    raw = pd.read_csv(RES / "peak_list_haird_real_dilution.csv")
    feat = raw[raw["Sample"] != "Label"].copy()
    feat["ppm"] = feat["Sample"].str.replace("ppm_", "", regex=False).astype(float)
    feat = feat.sort_values("ppm", ascending=False).reset_index(drop=True)
    mat = feat[SAMPLE_IDS].apply(pd.to_numeric, errors="coerce").to_numpy()

    fig, ax = plt.subplots(figsize=(6, 10))
    plot_mat = np.log10(np.where(mat > 0, mat, np.nan))
    im = ax.imshow(plot_mat, aspect="auto", cmap="viridis",
                    extent=[0, len(SAMPLE_IDS), len(feat), 0])
    ax.set_xticks(np.arange(len(SAMPLE_IDS)) + 0.5)
    ax.set_xticklabels(SAMPLE_IDS, rotation=20, ha="right")
    step = max(1, len(feat) // 40)
    ytick_idx = np.arange(0, len(feat), step)
    ax.set_yticks(ytick_idx + 0.5)
    ax.set_yticklabels([f"{feat.loc[i, 'ppm']:.3f}" for i in ytick_idx], fontsize=7)
    ax.set_ylabel("ppm (merged feature)")
    ax.set_title(f"Step 3: merged feature table\n({len(feat)} features x {len(SAMPLE_IDS)} samples, log10 intensity, white = not detected)")
    cbar = fig.colorbar(im, ax=ax, fraction=0.04)
    cbar.set_label("log10(intensity)")
    fig.tight_layout()
    fig.savefig(RES / "step5_merged_feature_heatmap.png", dpi=130)
    plt.close(fig)
    n_full = int((mat > 0).all(axis=1).sum())
    print(f"Wrote step5_merged_feature_heatmap.png ({n_full}/{len(feat)} features detected in all 3 samples)")


def step6_confidence_distribution():
    df = pd.read_csv(RES / "compound_matches_haird_real_dilution.csv")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    axes[0].hist(df["confidence"], bins=30, color="#2f6690")
    axes[0].set_xlabel("confidence")
    axes[0].set_ylabel("candidate-compound rows")
    axes[0].set_title("Step 4: confidence distribution,\nall candidates, all 3 samples")

    for s in SAMPLE_IDS:
        sub = df[df["sample_id"] == s]
        axes[1].hist(sub["confidence"], bins=20, alpha=0.5, label=s, color=COLORS[s])
    axes[1].set_xlabel("confidence")
    axes[1].set_title("Same, split by sample")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(RES / "step6_confidence_distribution.png", dpi=130)
    plt.close(fig)
    print("Wrote step6_confidence_distribution.png")


def step7_top_matches_bar(top_n: int = 15):
    df = pd.read_csv(RES / "compound_matches_haird_real_dilution.csv")
    counts = df.groupby("compound_id")["sample_id"].nunique()
    in_all = counts[counts == 3].index
    sub = df[df["compound_id"].isin(in_all)]
    agg = (sub.groupby(["compound_id", "compound_name"])["confidence"].mean()
           .sort_values(ascending=False).head(top_n))
    fig, ax = plt.subplots(figsize=(8, 6))
    labels = [f"{name}\n({cid})" for cid, name in agg.index]
    ax.barh(range(len(agg)), agg.values, color="#2f6690")
    ax.set_yticks(range(len(agg)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("mean confidence across all 3 samples")
    ax.set_title(f"Step 4: top {top_n} compounds matched in ALL 3 real samples\n(fixed evidence_score scorer)")
    fig.tight_layout()
    fig.savefig(RES / "step7_top_matches_all_samples.png", dpi=130)
    plt.close(fig)
    print("Wrote step7_top_matches_all_samples.png")
    return agg


if __name__ == "__main__":
    step1_raw_overlay()
    step2_baseline_example()
    step3_detected_peaks()
    step4_detector_comparison()
    step5_merged_feature_heatmap()
    step6_confidence_distribution()
    step7_top_matches_bar()
