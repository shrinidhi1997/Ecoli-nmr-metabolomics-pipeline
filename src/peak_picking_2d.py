"""
peak_picking_2d.py
OPTIONAL. Reads a raw 2D [1H,13C]-HSQC Bruker spectrum and picks cross-
peaks. Nothing in the 1D pipeline (baseline_peak_picking.R,
run_peak_picking_control_vs_sample.R) imports or depends on this file --
it exists only for hsqc_analysis.py / hsqc_matcher.py, and skipping HSQC
entirely leaves the rest of the pipeline untouched.

This is a deliberately small, separate module rather than a shared file
with the 1D reader: the 1D side of this pipeline reads Bruker data through
R (ASICS's importSpectraBruker()/createSpectra()), not through this
Python/nmrglue path at all, so there is no 1D code here to share.

read_bruker_2d() and pick_peaks_2d() below are unchanged from an earlier,
all-Python version of this pipeline that read both 1D and 2D data through
nmrglue directly; that 1D path was later replaced by the current R/ASICS
pipeline, but the 2D reader never depended on it and still works standalone.
"""

from __future__ import annotations
from pathlib import Path

import numpy as np
import pandas as pd

import nmrglue as ng


def read_bruker_2d(pdata_dir: str | Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    pdata_dir: path to a 2D experiment's pdata folder (e.g. '.../pdata/1'
    containing 2rr + procs + proc2s). Returns (h_ppm_scale, c_ppm_scale,
    intensity_2d) where intensity_2d has shape (len(c_ppm_scale),
    len(h_ppm_scale)) -- i.e. rows are the indirect (13C) dimension and
    columns are the direct (1H) dimension, matching how nmrglue returns
    Bruker 2D data (dim 0 = indirect, dim 1 = direct) and how it is
    conventionally displayed (1H along x, 13C along y).
    """
    dic, data = ng.bruker.read_pdata(str(pdata_dir))
    udic = ng.bruker.guess_udic(dic, data)
    if udic["ndim"] != 2:
        raise ValueError(
            f"{pdata_dir}: read {udic['ndim']}D data, expected 2D. "
            "Check this pdata folder is really a 2D experiment (should contain "
            "2rr, not 1r)."
        )
    uc0 = ng.fileiobase.uc_from_udic(udic, dim=0)
    uc1 = ng.fileiobase.uc_from_udic(udic, dim=1)
    scale0 = np.asarray(uc0.ppm_scale(), dtype=float)  # indirect dim (rows)
    scale1 = np.asarray(uc1.ppm_scale(), dtype=float)  # direct dim (columns)

    label0 = str(udic[0].get("label", "")).upper()
    label1 = str(udic[1].get("label", "")).upper()
    data = np.asarray(data, dtype=float)

    # Normalize to (h_scale, c_scale, data[c, h]) regardless of which
    # physical dimension nmrglue put first, using the axis labels it
    # derived from the acqus/procs nucleus tags rather than assuming order.
    if "H" in label0 and "C" in label1:
        # dim0 = 1H (rows), dim1 = 13C (cols) -- transpose so rows=13C, cols=1H
        h_scale, c_scale = scale0, scale1
        data = data.T
    elif "C" in label0 and "H" in label1:
        c_scale, h_scale = scale0, scale1
        # already rows=13C(dim0), cols=1H(dim1) -- no transpose needed
    else:
        raise ValueError(
            f"{pdata_dir}: could not identify a 1H/13C dimension pair from "
            f"axis labels {label0!r}, {label1!r}. This reader is written "
            "specifically for 2D [1H,13C]-HSQC data."
        )
    return h_scale, c_scale, data


def pick_peaks_2d(h_scale: np.ndarray, c_scale: np.ndarray, intensity_2d: np.ndarray,
                   threshold_multiplier: float = 6.0,
                   min_separation_h_ppm: float = 0.03,
                   min_separation_c_ppm: float = 0.5,
                   exclude_h_ppm_ranges: list[tuple[float, float]] | None = None) -> pd.DataFrame:
    """
    Local-maximum peak picking for a 2D [1H,13C]-HSQC intensity array
    (rows=13C, cols=1H, as returned by read_bruker_2d). Optional and
    separate from the 1D picker in baseline_peak_picking.R -- HSQC
    processing never touches 1D results.

    Peak selection is two stages, mirroring the 1D picker's height +
    minimum-separation combination but adapted to 2D:
      1. A coarse candidate pass: any point that is the maximum inside a
         small fixed 3x3-pixel window AND clears threshold_multiplier *
         sigma (sigma a single global MAD-based noise estimate over the
         whole array). This is only a cheap prefilter to drop non-maxima
         and flat noise, not the real separation enforcement -- on a
         coarsely-digitized array, min_separation_h_ppm/min_separation_c_ppm
         can correspond to less than one pixel, which would make a maximum-
         filter window alone useless (every noisy point above threshold
         would pass as its own 1-pixel-window "local max").
      2. Greedy non-maximum suppression in actual ppm units: candidates are
         visited tallest-first, and a candidate is kept only if it does not
         fall within min_separation_h_ppm AND min_separation_c_ppm of an
         already-kept peak (a rectangular exclusion zone around each kept
         peak). This is what actually enforces "peaks must be this far
         apart" regardless of how finely the spectrum is digitized, the way
         the 1D picker's `distance` parameter does for scipy.find_peaks.

    (Unlike the 1D picker, this does not use a locally-adaptive noise
    baseline -- HSQC's dominant artifact is a water/solvent STRIPE along
    the 1H axis, which exclude_h_ppm_ranges handles directly, rather than a
    slowly-varying shoulder that would need a rolling baseline.)
    """
    from scipy.ndimage import maximum_filter

    h_scale = np.asarray(h_scale, dtype=float)
    c_scale = np.asarray(c_scale, dtype=float)
    intensity_2d = np.asarray(intensity_2d, dtype=float)
    if intensity_2d.shape != (len(c_scale), len(h_scale)):
        raise ValueError(
            f"intensity_2d shape {intensity_2d.shape} does not match "
            f"(len(c_scale), len(h_scale)) = ({len(c_scale)}, {len(h_scale)})"
        )

    flat = intensity_2d.ravel()
    mad = np.median(np.abs(flat - np.median(flat)))
    sigma = 1.4826 * mad if mad > 0 else np.std(flat)
    threshold = threshold_multiplier * sigma

    # Stage 1: cheap fixed-window prefilter (candidates only, not the real
    # separation logic -- see docstring).
    local_max = maximum_filter(intensity_2d, size=(3, 3), mode="nearest")
    is_candidate = (intensity_2d == local_max) & (intensity_2d > threshold)
    c_idx, h_idx = np.nonzero(is_candidate)

    candidates = pd.DataFrame({
        "h_ppm": h_scale[h_idx],
        "c_ppm": c_scale[c_idx],
        "intensity": intensity_2d[c_idx, h_idx],
    }).sort_values("intensity", ascending=False).reset_index(drop=True)

    if exclude_h_ppm_ranges:
        keep = np.ones(len(candidates), dtype=bool)
        for lo, hi in exclude_h_ppm_ranges:
            lo, hi = min(lo, hi), max(lo, hi)
            keep &= ~((candidates["h_ppm"] >= lo) & (candidates["h_ppm"] <= hi))
        candidates = candidates[keep].reset_index(drop=True)

    # Stage 2: greedy non-maximum suppression in real ppm units.
    kept_h: list[float] = []
    kept_c: list[float] = []
    kept_rows: list[int] = []
    for idx, row in candidates.iterrows():
        too_close = any(
            abs(row["h_ppm"] - kh) < min_separation_h_ppm
            and abs(row["c_ppm"] - kc) < min_separation_c_ppm
            for kh, kc in zip(kept_h, kept_c)
        )
        if not too_close:
            kept_h.append(row["h_ppm"])
            kept_c.append(row["c_ppm"])
            kept_rows.append(idx)

    peaks = candidates.loc[kept_rows].copy()
    peaks["h_ppm"] = np.round(peaks["h_ppm"], 4)
    peaks["c_ppm"] = np.round(peaks["c_ppm"], 3)
    return peaks.sort_values("intensity", ascending=False).reset_index(drop=True)
