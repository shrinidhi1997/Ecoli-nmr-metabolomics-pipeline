"""
plot_multiplet_overlap_benchmark.py

Draws the MultiNMRFit-style accuracy/precision panels for this project's
own pipeline, from run_multiplet_overlap_benchmark.R's output CSVs:
  <res_dir>/accuracy_precision_summary.csv
  <res_dir>/resolution_rate.csv

Three figures, not one -- each answers a different question:

  1. multiplet_benchmark_accuracy_precision.png
     The direct analogue of the paper's own figure: for each parameter
     (x, J, I, lw), BIAS (mean fitted-minus-true, left column) and
     PRECISION (SD of that error across repeated noise draws, right
     column) vs. S/N, one line per Hz separation. Solid = triplet,
     dashed = quartet. Uses only the "production_bound" setting (this
     project's real, shipped refine_peaks_lorentzian() defaults) -- this
     is what the pipeline actually does today, not a best-case tuning.

  2. multiplet_benchmark_resolution_rate.png
     The main result of this benchmark: the fraction of
     draws where every expected line of a multiplet came back as its own
     distinct fitted peak (not merged with a neighbour, not missed
     entirely), vs. S/N, one line per Hz separation, triplet and quartet
     as separate panels. Bias/precision above are only meaningful where
     this is high -- a "great" bias number averaged over a cell that
     resolved on 2 of 10 draws is an artifact of survivorship, not
     accuracy, which is why this script ships as a SEPARATE, equally
     prominent figure rather than a footnote.

  3. multiplet_benchmark_linewidth_bound_effect.png
     Linewidth bias specifically, production_bound vs. widened_bound
     side by side. refine_peaks_lorentzian()'s linewidth floor
     (0.001 ppm) was tuned for this project's real data, where it's
     irrelevant (real linewidths are always well above it); a
     synthetic line narrower than that floor gets its width
     measurement artificially inflated towards the floor, not a real
     detector limitation -- this figure shows exactly how much of the
     "production" linewidth bias is that artifact vs. genuine fit error.

Color: hz_sep uses this project's validated 3-slot categorical order
(blue=15, orange=20, aqua=30 Hz) -- fixed across every panel in every
figure here, never re-cycled per subplot, so the same color always means
the same separation.

Run from this src/ folder (after run_multiplet_overlap_benchmark.R):
    python plot_multiplet_overlap_benchmark.py
"""
from __future__ import annotations
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

# Fixed categorical order -- validated 3-slot subset of this project's
# default palette (see the dataviz skill's references/palette.md): worst
# adjacent/all-pairs CVD Delta-E and normal-vision Delta-E both clear the
# floors for exactly 3 series. Never reassign these per-panel.
HZ_COLORS = {15: "#2a78d6", 20: "#eb6834", 30: "#1baf7a"}
MULTIPLET_STYLE = {"triplet": {"linestyle": "-", "marker": "o"},
                    "quartet": {"linestyle": "--", "marker": "s"}}

PARAM_LABELS = {
    "x_hz": ("chemical shift x", "Hz"),
    "J_hz": ("coupling constant J (derived)", "Hz"),
    "I_area": ("intensity (fitted area)", "a.u."),
    "lw_hz": ("linewidth lw (FWHM)", "Hz"),
}


def _style_axis(ax):
    ax.set_xscale("log")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.axhline(0, color="#b5b4ad", linewidth=0.8, zorder=0)
    ax.grid(True, which="major", axis="both", color="#e8e7e2", linewidth=0.6, zorder=0)


def plot_accuracy_precision(summary: pd.DataFrame, res_dir: Path, bound_setting: str = "production_bound"):
    sub = summary[summary["bound_setting"] == bound_setting]
    params = ["x_hz", "J_hz", "I_area", "lw_hz"]
    fig, axes = plt.subplots(len(params), 2, figsize=(10, 3.1 * len(params)))

    for row, param in enumerate(params):
        for col, metric in enumerate(["bias", "sd"]):
            ax = axes[row, col]
            pdat = sub[sub["parameter"] == param]
            for hz_sep, color in HZ_COLORS.items():
                for mtype, style in MULTIPLET_STYLE.items():
                    d = pdat[(pdat["hz_sep"] == hz_sep) & (pdat["multiplet_type"] == mtype)]
                    d = d.sort_values("snr_level")
                    if d.empty:
                        continue
                    ax.plot(d["snr_level"], d[metric], color=color, linewidth=2,
                             marker=style["marker"], markersize=5, linestyle=style["linestyle"])
            _style_axis(ax)
            label, unit = PARAM_LABELS[param]
            metric_name = "bias (fit - true)" if metric == "bias" else "precision (SD across draws)"
            ax.set_title(f"{label} -- {metric_name}", fontsize=10)
            ax.set_ylabel(unit, fontsize=9)
            if row == len(params) - 1:
                ax.set_xlabel("S/N (log scale)", fontsize=9)

    # legend: color = Hz separation, linestyle/marker = multiplet type --
    # identity never carried by color alone (dataviz check 6).
    sep_handles = [mlines.Line2D([], [], color=c, linewidth=2, label=f"{hz} Hz separation")
                   for hz, c in HZ_COLORS.items()]
    type_handles = [mlines.Line2D([], [], color="#52514e", linewidth=2,
                                   linestyle=s["linestyle"], marker=s["marker"], markersize=5,
                                   label=mtype)
                    for mtype, s in MULTIPLET_STYLE.items()]
    fig.legend(handles=sep_handles + type_handles, loc="lower center", ncol=5,
               bbox_to_anchor=(0.5, 0.0), fontsize=9, frameon=False)
    fig.suptitle(f"This pipeline's accuracy & precision vs. known synthetic ground truth\n"
                 f"(refine_peaks_lorentzian(), {bound_setting} -- see script header)", fontsize=11)
    fig.tight_layout(rect=[0, 0.035, 1, 0.965])
    out = res_dir / "multiplet_benchmark_accuracy_precision.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"Wrote {out}")


def plot_resolution_rate(resolution_rate: pd.DataFrame, res_dir: Path, bound_setting: str = "production_bound"):
    sub = resolution_rate[resolution_rate["bound_setting"] == bound_setting]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.6), sharey=True)
    for ax, mtype in zip(axes, ["triplet", "quartet"]):
        d_all = sub[sub["multiplet_type"] == mtype]
        for hz_sep, color in HZ_COLORS.items():
            d = d_all[d_all["hz_sep"] == hz_sep].sort_values("snr_level")
            if d.empty:
                continue
            ax.plot(d["snr_level"], d["resolution_rate"], color=color, linewidth=2.2,
                    marker="o", markersize=6)
        ax.set_xscale("log")
        ax.set_ylim(-0.05, 1.05)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(True, which="major", color="#e8e7e2", linewidth=0.6, zorder=0)
        ax.set_title(mtype, fontsize=11)
        ax.set_xlabel("S/N (log scale)", fontsize=9)
    axes[0].set_ylabel("fraction of draws fully resolved\n(every expected line its own fitted peak)", fontsize=9)

    sep_handles = [mlines.Line2D([], [], color=c, linewidth=2.2, marker="o", markersize=6,
                                  label=f"{hz} Hz separation") for hz, c in HZ_COLORS.items()]
    fig.legend(handles=sep_handles, loc="lower center", ncol=3, bbox_to_anchor=(0.5, 0.0),
               fontsize=9, frameon=False)
    fig.suptitle("Where this pipeline starts losing distinct peaks\n"
                 "(fraction of synthetic draws where the triplet/quartet fully resolved)", fontsize=11)
    fig.tight_layout(rect=[0, 0.1, 1, 0.9])
    out = res_dir / "multiplet_benchmark_resolution_rate.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"Wrote {out}")


def plot_linewidth_bound_effect(summary: pd.DataFrame, res_dir: Path):
    sub = summary[summary["parameter"] == "lw_hz"]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    for ax, bound_setting in zip(axes, ["production_bound", "widened_bound"]):
        d_all = sub[sub["bound_setting"] == bound_setting]
        for hz_sep, color in HZ_COLORS.items():
            for mtype, style in MULTIPLET_STYLE.items():
                d = d_all[(d_all["hz_sep"] == hz_sep) & (d_all["multiplet_type"] == mtype)].sort_values("snr_level")
                if d.empty:
                    continue
                ax.plot(d["snr_level"], d["bias"], color=color, linewidth=2,
                         marker=style["marker"], markersize=5, linestyle=style["linestyle"])
        _style_axis(ax)
        ax.set_title(bound_setting.replace("_", " "), fontsize=11)
        ax.set_xlabel("S/N (log scale)", fontsize=9)
    axes[0].set_ylabel("linewidth bias, Hz FWHM (fit - true)", fontsize=9)

    sep_handles = [mlines.Line2D([], [], color=c, linewidth=2, label=f"{hz} Hz separation")
                   for hz, c in HZ_COLORS.items()]
    type_handles = [mlines.Line2D([], [], color="#52514e", linewidth=2,
                                   linestyle=s["linestyle"], marker=s["marker"], markersize=5, label=mtype)
                    for mtype, s in MULTIPLET_STYLE.items()]
    fig.legend(handles=sep_handles + type_handles, loc="lower center", ncol=5,
               bbox_to_anchor=(0.5, 0.0), fontsize=9, frameon=False)
    fig.suptitle("How much of the linewidth bias is refine_peaks_lorentzian()'s\n"
                 "production floor (0.001 ppm), not the fit itself", fontsize=11)
    fig.tight_layout(rect=[0, 0.09, 1, 0.9])
    out = res_dir / "multiplet_benchmark_linewidth_bound_effect.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"Wrote {out}")


if __name__ == "__main__":
    res_dir = Path("../example_results/multiplet_benchmark")
    summary = pd.read_csv(res_dir / "accuracy_precision_summary.csv")
    resolution_rate = pd.read_csv(res_dir / "resolution_rate.csv")

    plot_accuracy_precision(summary, res_dir)
    plot_resolution_rate(resolution_rate, res_dir)
    plot_linewidth_bound_effect(summary, res_dir)
