"""
plot_multiplet_benchmark_heatmap.py

Modelled on the accuracy/precision figures used in published NMR
fitting benchmarks (heatmaps of chemical shift, coupling constant,
intensity and linewidth vs. S/N, plus an example-fit panel split by
quartet/triplet and Hz separation). This script renders that visual
grammar --
heatmaps, not line plots -- for THIS project's own pipeline, from the
already-existing multiplet-overlap benchmark
(simulate_multiplet_overlap_benchmark.R -> run_multiplet_overlap_
benchmark.R), reading its per-line output rather than its pre-aggregated
bias/sd-in-Hz summary, because the reference figure's cells are
dimensionless percentages, and that conversion has to happen here.

WHAT "ACCURACY (%)" AND "PRECISION (%)" MEAN HERE -- READ THIS FIRST:
The reference figure does not state its own normalization, and the
earlier line-plot script (plot_multiplet_overlap_benchmark.py) reports
bias/sd in each parameter's native unit (Hz, a.u. area) specifically to
avoid inventing one. This script has to pick a normalization to get a
single "%" scale across four different-unit parameters, so it uses, per
matched line, the fitted-minus-true error divided by that line's own
TRUE value of the most natural reference quantity:
    chemical shift x  -> error (Hz) / true linewidth (Hz FWHM)
    coupling J         -> error (Hz) / true J (Hz)
    intensity I        -> error (area) / true area
    linewidth lw       -> error (Hz) / true linewidth (Hz FWHM)
Chemical-shift error is expressed relative to linewidth (not to the
absolute ppm position, which is arbitrary and would make every cell
read near 0%) -- a standard way to ask "how many linewidths off was the
fit." Per (hz_sep, S/N, multiplet_type, parameter) cell:
    Accuracy (%)  = |mean(relative error)|   -- systematic bias
    Precision (%) = SD(relative error)       -- scatter across draws
matching this project's own existing bias/sd convention
(run_multiplet_overlap_benchmark.R's agg_err()), just rescaled to %.
These numbers are therefore NOT directly comparable to the original
paper's own percentages (different fitting method -- independent
per-line Lorentzians here vs. a shared constrained multiplet model
there -- and almost certainly a different normalization choice); they
are this pipeline's own accuracy/precision, in the same visual form.

J is only defined for a multiplet where EVERY line resolved to its own
fitted peak in that draw (see run_multiplet_overlap_benchmark.R's own
j_derived flag) -- cells with zero resolved draws are plotted as grey
"n/a" rather than silently dropped. Quartets at 15/20 Hz separation
never resolve at any S/N with this pipeline's production bound (see
resolution_rate.csv) -- the heatmap shows that honestly as a grey J
column, not a guessed number.

Reads (production_bound setting only -- this project's real, shipped
defaults, same choice the line-plot script makes):
    ../example_results/multiplet_benchmark/per_line_results.csv
    ../example_results/multiplet_benchmark/j_derived_results.csv
    ../example_data/multiplet_benchmark/spectra_long.csv   (figure 2 panel C only)

Writes two figures to ../example_results/multiplet_benchmark/:
  1. multiplet_benchmark_heatmap_overall.png
     Overall view: one
     Accuracy (%) heatmap and one Precision (%) heatmap, S/N (2.5-50)
     down the rows, the four parameters across the columns -- pooled
     across both multiplet types and all three Hz separations, i.e.
     "how does this pipeline do overall, as a function of S/N."

  2. multiplet_benchmark_heatmap_by_separation.png
     Per-separation view: one row
     block per Hz separation (30, 20, 15 -- i.e. decreasing, so
     increasing overlap top to bottom), each with an Accuracy heatmap
     and a Precision heatmap split into quartet | triplet column groups
     (x, J, I, lw sub-columns), at a 4-value S/N subset (2.5, 5, 10, 20)
     for legibility, plus a third panel showing one real example draw at
     S/N=10 for that separation: the raw synthetic points, this
     pipeline's total fitted curve, and each matched line's own fitted
     Lorentzian shaded underneath (quartet lines in one color, triplet
     lines in another) -- so the heatmap numbers have a picture to go
     with them, same as the reference figure's panel C.

Run from this src/ folder (after simulate_multiplet_overlap_benchmark.R
and run_multiplet_overlap_benchmark.R have (re)written their outputs):
    python plot_multiplet_benchmark_heatmap.py
"""
from __future__ import annotations
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

RES_DIR = Path("../example_results/multiplet_benchmark")
DATA_DIR = Path("../example_data/multiplet_benchmark")
BOUND_SETTING = "production_bound"  # this project's real, shipped default

PARAM_ORDER = ["x_hz", "J_hz", "I_area", "lw_hz"]
PARAM_LONG = {"x_hz": "chemical shift", "J_hz": "coupling constant",
              "I_area": "intensity", "lw_hz": "linewidth"}
PARAM_SHORT = {"x_hz": "x", "J_hz": "J", "I_area": "I", "lw_hz": "lw"}

MULTIPLET_COLOR = {"quartet": "#2a78d6", "triplet": "#eb6834"}  # fig-2 panel C fill colors

CMAP = matplotlib.colormaps["YlOrRd"]
CMAP_BAD = matplotlib.colormaps["YlOrRd"].copy()
CMAP_BAD.set_bad("#d9d9d4")


def lorentzian(x, x0, height, hwhm):
    return height / (1 + ((x - x0) / hwhm) ** 2)


def load_relative_errors() -> pd.DataFrame:
    """One row per matched line/derivation, with a `parameter` column and
    a `rel_pct` column already normalized per the header above."""
    per_line = pd.read_csv(RES_DIR / "per_line_results.csv")
    per_line = per_line[per_line["bound_setting"] == BOUND_SETTING]
    j_tbl = pd.read_csv(RES_DIR / "j_derived_results.csv")
    j_tbl = j_tbl[j_tbl["bound_setting"] == BOUND_SETTING]

    keep_cols = ["hz_sep", "snr_level", "multiplet_type"]
    rows = []

    d = per_line.copy()
    d["parameter"] = "x_hz"
    d["rel_pct"] = d["err_x_hz"] / d["true_lw_hz"] * 100
    rows.append(d[keep_cols + ["parameter", "rel_pct"]])

    d = per_line.copy()
    d["parameter"] = "I_area"
    d["rel_pct"] = d["err_I"] / d["true_area"] * 100
    rows.append(d[keep_cols + ["parameter", "rel_pct"]])

    d = per_line.copy()
    d["parameter"] = "lw_hz"
    d["rel_pct"] = d["err_lw_hz"] / d["true_lw_hz"] * 100
    rows.append(d[keep_cols + ["parameter", "rel_pct"]])

    d = j_tbl.copy()
    d["parameter"] = "J_hz"
    d["rel_pct"] = d["err_J_hz"] / d["true_J_hz"] * 100
    rows.append(d[keep_cols + ["parameter", "rel_pct"]])

    return pd.concat(rows, ignore_index=True)


def aggregate_cells(long_df: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)  # all-NaN groups (unresolved J)
        g = long_df.groupby(group_cols)["rel_pct"]
        out = g.agg(accuracy=lambda v: np.abs(np.nanmean(v)) if v.notna().any() else np.nan,
                     precision=lambda v: np.nanstd(v) if v.notna().any() else np.nan,
                     n=lambda v: int(v.notna().sum())).reset_index()
    return out


VMAX_PERCENTILE = 80  # see header note on why color scale is capped, not stretched to the max


def capped_vmax(values: np.ndarray) -> float:
    """80th-percentile color-scale cap, same spirit as this project's own
    CLIP_PCT convention elsewhere (fit_binned_baseline()'s callers, the
    unlabeled-stacked-spectra plots): a couple of catastrophic S/N=2.5
    cells (independent-Lorentzian fitting genuinely breaks down there --
    a real finding, not a bug) would otherwise stretch the whole
    colormap and flatten every other, more informative cell to pale
    yellow. Capped cells still show their true printed number and still
    render as the darkest, most-saturated color -- nothing is hidden,
    only the color gradient is kept usable for the cells that aren't
    off the charts."""
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return 1.0
    v = float(np.percentile(finite, VMAX_PERCENTILE))
    return v if v > 0 else float(np.nanmax(finite)) or 1.0


def draw_heatmap(ax, matrix, row_labels, col_labels, vmax, title, cbar_label):
    masked = np.ma.masked_invalid(matrix)
    im = ax.imshow(masked, cmap=CMAP_BAD, vmin=0, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(col_labels)))
    ax.set_xticklabels(col_labels, fontsize=8.5)
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=8.5)
    ax.set_title(title, fontsize=10.5, pad=8)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(length=0)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            v = matrix[i, j]
            txt = "n/a" if np.isnan(v) else f"{v:.1f}"
            lum = 0 if np.isnan(v) else (1 - v / vmax if vmax else 0)
            color = "#1a1a19" if (np.isnan(v) or lum > 0.45) else "#ffffff"
            ax.text(j, i, txt, ha="center", va="center", fontsize=7, color=color)
    return im


def plot_overall(long_df: pd.DataFrame, res_dir: Path):
    snr_order = sorted(long_df["snr_level"].unique())
    cells = aggregate_cells(long_df, ["snr_level", "parameter"])
    col_labels = [PARAM_LONG[p] for p in PARAM_ORDER]

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 5.6))
    for ax, metric, label in zip(axes, ["accuracy", "precision"],
                                  ["Accuracy (%)", "Precision (%)"]):
        mat = np.full((len(snr_order), len(PARAM_ORDER)), np.nan)
        for i, snr in enumerate(snr_order):
            for j, p in enumerate(PARAM_ORDER):
                r = cells[(cells.snr_level == snr) & (cells.parameter == p)]
                if len(r):
                    mat[i, j] = r[metric].iloc[0]
        vmax = capped_vmax(mat)
        row_labels = [f"{s:g}" for s in snr_order]
        im = draw_heatmap(ax, mat, row_labels, col_labels, vmax, label, label)
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label(label, fontsize=9)
    axes[0].set_ylabel("S/N", fontsize=10)
    fig.suptitle("This pipeline's accuracy & precision vs. known synthetic ground truth\n"
                 f"(pooled across quartet/triplet and 15/20/30 Hz separation, {BOUND_SETTING})",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    out = res_dir / "multiplet_benchmark_heatmap_overall.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Wrote {out}")


def build_example_panel(ax, hz_sep: int, snr_level: float):
    """One representative draw for this hz_sep at the given S/N: raw
    points, total fitted curve, and each matched line's own fitted
    Lorentzian shaded underneath. Uses the raw simulated intensity
    directly as "data" -- this generator injects no baseline drift (see
    simulate_multiplet_overlap_benchmark.R's own header), so raw and
    baseline-corrected are the same trace here."""
    per_line = pd.read_csv(RES_DIR / "per_line_results.csv")
    per_line = per_line[(per_line.bound_setting == BOUND_SETTING) &
                         (per_line.hz_sep == hz_sep) & (per_line.snr_level == snr_level)]
    gt = pd.read_csv(DATA_DIR / "ground_truth.csv")
    gt = gt[(gt.hz_sep == hz_sep) & (gt.snr_level == snr_level)]
    draw_id = int(sorted(gt.draw_id.unique())[0])

    sp = pd.read_csv(DATA_DIR / "spectra_long.csv")
    sp = sp[sp.draw_id == draw_id].sort_values("ppm")
    ax.plot(sp.ppm, sp.intensity, "o", color="#1a1a19", markersize=1.6,
            markerfacecolor="none", markeredgewidth=0.5, zorder=3)

    row = per_line[per_line.draw_id == draw_id]
    ppm_grid = sp.ppm.to_numpy()
    total = np.zeros_like(ppm_grid)
    for mtype in ("quartet", "triplet"):
        sub = row[(row.multiplet_type == mtype) & row.matched]
        comp_sum = np.zeros_like(ppm_grid)
        for _, ln in sub.iterrows():
            height = ln.fit_area / (np.pi * ln.fit_hwhm_ppm)
            comp = lorentzian(ppm_grid, ln.fit_ppm, height, ln.fit_hwhm_ppm)
            comp_sum += comp
        total += comp_sum
        if comp_sum.max() > 0:
            ax.fill_between(ppm_grid, comp_sum, color=MULTIPLET_COLOR[mtype], alpha=0.28,
                             linewidth=0, zorder=1, label=f"{mtype} (fitted)")
    ax.plot(ppm_grid, total, color="#2d2d2b", linewidth=1.1, zorder=2, label="total fit")
    ax.invert_xaxis()
    ax.set_xlabel("chemical shift (ppm)", fontsize=8)
    ax.set_yticks([])
    for sp_ in ("top", "right", "left"):
        ax.spines[sp_].set_visible(False)
    ax.set_title(f"example draw, S/N={snr_level:g}", fontsize=9)
    ax.tick_params(labelsize=7.5)
    ax.legend(fontsize=6.5, frameon=False, loc="upper left")


def plot_by_separation(long_df: pd.DataFrame, res_dir: Path):
    hz_sep_order = [30, 20, 15]
    snr_subset = [2.5, 5, 10, 20]
    col_labels = []
    for mtype in ("quartet", "triplet"):
        for p in PARAM_ORDER:
            col_labels.append(PARAM_SHORT[p])

    cells = aggregate_cells(long_df, ["hz_sep", "snr_level", "multiplet_type", "parameter"])
    # shared color scale per metric across all three row-blocks, so the
    # three rows are honestly comparable to each other (not each
    # auto-stretched to its own max).
    sub_all = cells[cells.snr_level.isin(snr_subset)]
    vmax_acc = capped_vmax(sub_all["accuracy"].to_numpy())
    vmax_prec = capped_vmax(sub_all["precision"].to_numpy())

    fig, axes = plt.subplots(3, 3, figsize=(17, 13.5),
                              gridspec_kw={"width_ratios": [1.3, 1.3, 1]})
    for r, hz_sep in enumerate(hz_sep_order):
        for col, (metric, vmax, label) in enumerate([
            ("accuracy", vmax_acc, "Accuracy (%)"), ("precision", vmax_prec, "Precision (%)")]):
            mat = np.full((len(snr_subset), 8), np.nan)
            for i, snr in enumerate(snr_subset):
                for j, mtype in enumerate(("quartet", "triplet")):
                    for k, p in enumerate(PARAM_ORDER):
                        rr = cells[(cells.hz_sep == hz_sep) & (cells.snr_level == snr) &
                                   (cells.multiplet_type == mtype) & (cells.parameter == p)]
                        if len(rr):
                            mat[i, j * 4 + k] = rr[metric].iloc[0]
            ax = axes[r, col]
            row_labels = [f"{s:g}" for s in snr_subset] if col == 0 else [""] * len(snr_subset)
            im = draw_heatmap(ax, mat, row_labels, col_labels, vmax, "", label)
            ax.axvline(3.5, color="white", linewidth=2)
            ax.text(1.5, -1.0, "quartet", ha="center", fontsize=9, fontweight="bold")
            ax.text(5.5, -1.0, "triplet", ha="center", fontsize=9, fontweight="bold")
            # one colorbar per (row, metric) rather than spanning all three
            # rows: a shared vmax already makes the three blocks
            # comparable, and a spanning colorbar in matplotlib landed on
            # top of the middle row's cell text -- this avoids that.
            cbar = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.03)
            cbar.ax.tick_params(labelsize=7)
            if r == 0:
                cbar.set_label(label, fontsize=8.5)
        axes[r, 0].set_ylabel(f"{hz_sep} Hz separation\nS/N", fontsize=9.5)
        build_example_panel(axes[r, 2], hz_sep, snr_level=10)

    fig.suptitle("Accuracy & precision by quartet/triplet separation (Hz) -- this pipeline, "
                 f"synthetic ground truth ({BOUND_SETTING})\n"
                 "grey \"n/a\" = that multiplet never resolved to one fitted peak per line "
                 "at that S/N, so J could not be derived",
                 fontsize=11)
    fig.tight_layout(rect=[0.03, 0, 0.97, 0.93])
    out = res_dir / "multiplet_benchmark_heatmap_by_separation.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Wrote {out}")


if __name__ == "__main__":
    long_df = load_relative_errors()
    plot_overall(long_df, RES_DIR)
    plot_by_separation(long_df, RES_DIR)
