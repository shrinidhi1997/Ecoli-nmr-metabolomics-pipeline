"""
plot_standard_dilution_detail.py

Two extra figure types beyond plot_haird_real_dilution_steps.py's step1-7:

1. <dataset>_drift_correction_full.png
   Raw spectrum + fitted airPLS baseline overlaid, FULL ppm range (not a
   couple of small windows) -- so the drift (and the baseline tracking
   it) is visible everywhere at once, next to the corrected result.

2. <dataset>_peak_zoom_panels.png
   A tight zoom (+/-0.03 ppm) around each of a sample's biggest detected
   peaks, corrected spectrum + the exact picked point marked -- so the
   actual detection decision is visible at the resolution of the real
   line shape, not just as a dot on a regional overview plot.

Generalized over dataset name / sample list / results folder so the same
script runs on haird_real_dilution and swepf_real_dilution (or any future
one) without copying code.

Run from this src/ folder:
    python plot_standard_dilution_detail.py haird HAIRD-30-1H HAIRD-60-1H HAIRD-120-1H
    python plot_standard_dilution_detail.py swepf SWEPF-30-1H SWEPF-60-1H SWEPF-120-1H
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


def drift_correction_full(dataset: str, sample_ids: list[str], res_dir: Path):
    # The water peak (~4.7-5.3 ppm) is 100-1000x taller than everything
    # else in a real 1H spectrum -- on a linear axis sized to fit it, the
    # actual drift in the rest of the spectrum (the thing this figure is
    # supposed to show) flattens to an invisible line at the bottom. So
    # the y-axis here is capped to real signal outside the water/DSS
    # exclusion zones (same ranges find_shape_based_peaks() itself
    # excludes by default), and the water peak is deliberately clipped --
    # noted in the title, not hidden.
    EXCLUDE = [(3.9, 5.3), (-0.05, 0.05)]

    def _non_excluded_max(ppm, y):
        mask = np.ones(len(ppm), dtype=bool)
        for lo, hi in EXCLUDE:
            mask &= ~((ppm >= lo) & (ppm <= hi))
        return np.nanmax(y[mask]) if mask.any() else np.nanmax(y)

    fig, axes = plt.subplots(len(sample_ids), 2, figsize=(12, 3.2 * len(sample_ids)))
    if len(sample_ids) == 1:
        axes = axes.reshape(1, 2)
    for row, s in enumerate(sample_ids):
        raw = pd.read_csv(res_dir / f"{s}_raw_spectrum.csv")
        base = pd.read_csv(res_dir / f"{s}_baseline.csv")
        corr = pd.read_csv(res_dir / f"{s}_corrected_spectrum.csv")

        ymax_raw = _non_excluded_max(raw["ppm"].to_numpy(), raw["intensity"].to_numpy()) * 1.4
        ymax_corr = _non_excluded_max(corr["ppm"].to_numpy(), corr["corrected"].to_numpy()) * 1.4

        ax = axes[row, 0]
        ax.plot(raw["ppm"], raw["intensity"], color="#333333", linewidth=0.5, label="raw")
        ax.plot(base["ppm"], base["baseline"], color="#d9534f", linewidth=1.2, label="big-bin baseline")
        ax.set_xlim(11, 0)
        ax.set_ylim(-ymax_raw * 0.15, ymax_raw)
        ax.set_ylabel("intensity (a.u.)")
        ax.set_title(f"{s}: raw + fitted baseline, 0-11 ppm reporting range\n(water peak clipped so real drift is visible)")
        ax.legend(fontsize=7, loc="upper right")

        ax2 = axes[row, 1]
        ax2.plot(corr["ppm"], corr["corrected"], color="#2f6690", linewidth=0.5)
        ax2.axhline(0, color="#999999", linewidth=0.6, linestyle="--")
        ax2.set_xlim(11, 0)
        ax2.set_ylim(-ymax_corr * 0.15, ymax_corr)
        ax2.set_title(f"{s}: after baseline subtraction\n(same y-scale cap, water clipped)")
    axes[-1, 0].set_xlabel("ppm")
    axes[-1, 1].set_xlabel("ppm")
    fig.suptitle(f"{dataset.upper()}: drift correction across the whole spectrum, not just a zoomed window")
    fig.tight_layout()
    out = res_dir / f"{dataset}_drift_correction_full.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    print(f"Wrote {out}")


def peak_zoom_panels(dataset: str, sample_ids: list[str], res_dir: Path, n_peaks: int = 6, half_win_ppm: float = 0.03):
    for s in sample_ids:
        corr = pd.read_csv(res_dir / f"{s}_corrected_spectrum.csv")
        peaks = pd.read_csv(res_dir / f"{s}_detected_peaks.csv")
        top = peaks.nlargest(n_peaks, "intensity").sort_values("ppm", ascending=False).reset_index(drop=True)

        ncols = 3
        nrows = int(np.ceil(len(top) / ncols))
        fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3 * nrows))
        axes = np.atleast_2d(axes)
        for i, (_, pk) in enumerate(top.iterrows()):
            ax = axes[i // ncols, i % ncols]
            lo, hi = pk["ppm"] - half_win_ppm, pk["ppm"] + half_win_ppm
            mask = (corr["ppm"] >= lo) & (corr["ppm"] <= hi)
            ax.plot(corr.loc[mask, "ppm"], corr.loc[mask, "corrected"], color="#2f6690", linewidth=1.0)
            ax.scatter([pk["ppm"]], [pk["intensity"]], color="#d9534f", zorder=5, s=40, marker="x")
            ax.set_xlim(hi, lo)
            ax.set_title(f"{pk['ppm']:.4f} ppm", fontsize=9)
            ax.tick_params(labelsize=7)
        for j in range(len(top), nrows * ncols):
            axes[j // ncols, j % ncols].axis("off")
        fig.suptitle(f"{s}: zoomed view of its {len(top)} largest detected peaks\n(red x = exact point the shape detector picked)")
        fig.tight_layout()
        out = res_dir / f"{s}_peak_zoom_panels.png"
        fig.savefig(out, dpi=130)
        plt.close(fig)
        print(f"Wrote {out}")


if __name__ == "__main__":
    dataset = sys.argv[1]
    sample_ids = sys.argv[2:]
    res_dir = Path(f"../example_results/{dataset}_real_dilution")
    drift_correction_full(dataset, sample_ids, res_dir)
    peak_zoom_panels(dataset, sample_ids, res_dir)
