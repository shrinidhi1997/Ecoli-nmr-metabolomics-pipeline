"""
hsqc_analysis.py
OPTIONAL 2D [1H,13C]-HSQC add-on. Nothing in the 1D pipeline
(baseline_peak_picking.R, run_peak_picking_control_vs_sample.R) imports
this file. Skip it entirely (don't run it) and every 1D workflow behaves
exactly as it did before HSQC support existed.

WHAT THIS DOES for one sample:
  1. reads the raw processed 2D Bruker HSQC spectrum (2rr + procs + proc2s)
  2. optionally applies a fixed ppm offset you supply (see OFFSET, below)
  3. picks 2D cross-peaks, excluding the water/solvent region on the 1H axis
  4. saves <sample_id>_hsqc_peaks.csv and <sample_id>_hsqc_plot.png

OFFSET: this repo's 1D side (baseline_peak_picking.R) does not compute or
export its own DSS/TMS calibration offset -- it assumes your Bruker
processing already referenced the spectrum correctly, and picks peaks off
that ppm axis as given. There is therefore no automatic offset to pull in
from the 1D step for HSQC either. --offset defaults to 0.0 (use the axis
as read) and exists only so you can supply a manual correction yourself if
you independently determine this sample's 2D data needs one, e.g. from
comparing a known cross-peak's expected position against where it actually
lands.

Run standalone (useful for testing against one HSQC folder without a
larger pipeline):
    python hsqc_analysis.py --sample-id sample_treatment_01 \\
        --path /path/to/your/staged/bruker/samples/sample_treatment_01-HSQC/pdata/1 \\
        --offset 0.0 --out-dir ../output

Requires the nmrglue package (pip install nmrglue), used only by this
file and peak_picking_2d.py -- the rest of this pipeline does not need it.
"""

from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from peak_picking_2d import read_bruker_2d, pick_peaks_2d  # noqa: E402


def process_hsqc_sample(sample_id: str, hsqc_path: str, offset: float | None,
                         out_dir: str | Path,
                         exclude_h_ppm_ranges: list[tuple[float, float]] | None = None,
                         threshold_multiplier: float = 6.0,
                         min_separation_h_ppm: float = 0.03,
                         min_separation_c_ppm: float = 0.4) -> pd.DataFrame:
    """
    Full per-sample 2D pipeline: read -> apply an optional manual ppm
    offset -> pick cross-peaks -> save CSV + plot. Returns the picked-peaks
    DataFrame (columns: sample_id, h_ppm, c_ppm, intensity).

    offset: a manual ppm shift on the 1H axis, or None/0.0 to use the axis
    as read (see this file's header for why there's no automatic value to
    pull in here).
    exclude_h_ppm_ranges: defaults to a single water-region cut, [(4.6, 5.0)],
    if not given; pass [] explicitly to disable exclusion entirely.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if exclude_h_ppm_ranges is None:
        exclude_h_ppm_ranges = [(4.6, 5.0)]

    h_scale, c_scale, intensity_2d = read_bruker_2d(hsqc_path)
    if offset:
        h_scale = h_scale + offset

    peaks = pick_peaks_2d(
        h_scale, c_scale, intensity_2d,
        threshold_multiplier=threshold_multiplier,
        min_separation_h_ppm=min_separation_h_ppm,
        min_separation_c_ppm=min_separation_c_ppm,
        exclude_h_ppm_ranges=exclude_h_ppm_ranges,
    )
    peaks.insert(0, "sample_id", sample_id)

    peaks.to_csv(out_dir / f"{sample_id}_hsqc_peaks.csv", index=False)

    _plot_hsqc(sample_id, h_scale, c_scale, intensity_2d, peaks,
                exclude_h_ppm_ranges, out_dir)

    return peaks


def _plot_hsqc(sample_id: str, h_scale: np.ndarray, c_scale: np.ndarray,
                intensity_2d: np.ndarray, peaks: pd.DataFrame,
                exclude_h_ppm_ranges: list[tuple[float, float]],
                out_dir: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 6))

    positive = intensity_2d[intensity_2d > 0]
    if positive.size:
        floor = np.percentile(positive, 80)
        ceiling = intensity_2d.max()
        levels = np.geomspace(max(floor, 1e-9), max(ceiling, floor * 1.01), 12)
        ax.contour(h_scale, c_scale, intensity_2d, levels=levels,
                   linewidths=0.4, colors="#3b6ea5")

    if not peaks.empty:
        ax.scatter(peaks["h_ppm"], peaks["c_ppm"], color="#d9534f", s=14,
                   zorder=5, label="picked cross-peaks")

    for lo, hi in exclude_h_ppm_ranges:
        ax.axvspan(lo, hi, color="grey", alpha=0.15)

    ax.set_xlim(h_scale.max(), h_scale.min())
    ax.set_ylim(c_scale.max(), c_scale.min())
    ax.set_xlabel("1H (ppm)")
    ax.set_ylabel("13C (ppm)")
    ax.set_title(f"{sample_id} -- HSQC")
    ax.legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(out_dir / f"{sample_id}_hsqc_plot.png", dpi=130)
    plt.close(fig)


def main():
    import argparse
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--path", required=True, help="HSQC pdata folder, e.g. .../pdata/1")
    ap.add_argument("--offset", type=float, default=0.0,
                     help="manual ppm shift on the 1H axis (0.0 = use the axis as read; "
                          "see this file's header for why there's no value to auto-derive)")
    ap.add_argument("--out-dir", default="../output")
    ap.add_argument("--threshold-multiplier", type=float, default=6.0)
    ap.add_argument("--min-separation-h-ppm", type=float, default=0.03)
    ap.add_argument("--min-separation-c-ppm", type=float, default=0.4)
    ap.add_argument("--exclude-h-ppm", nargs="*", default=["4.6,5.0"],
                     help="pairs 'lo,hi' to blank out on the 1H axis, e.g. 4.6,5.0 for water")
    args = ap.parse_args()

    exclude = [tuple(float(x) for x in pair.split(",")) for pair in args.exclude_h_ppm]
    peaks = process_hsqc_sample(
        args.sample_id, args.path, args.offset, args.out_dir,
        exclude_h_ppm_ranges=exclude,
        threshold_multiplier=args.threshold_multiplier,
        min_separation_h_ppm=args.min_separation_h_ppm,
        min_separation_c_ppm=args.min_separation_c_ppm,
    )
    print(f"Picked {len(peaks)} HSQC cross-peaks for {args.sample_id}")
    print(peaks.to_string(index=False))


if __name__ == "__main__":
    main()
