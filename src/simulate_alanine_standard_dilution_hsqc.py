"""
simulate_alanine_standard_dilution_hsqc.py

The 2D [1H,13C]-HSQC counterpart to simulate_alanine_standard_dilution_1d.R
-- a SIMULATED (not real instrument data) standard dilution series for
L-alanine's two cross-peaks, run through the real, unmodified
pick_peaks_2d() (peak_picking_2d.py) and the real _plot_hsqc() plotting
code (hsqc_analysis.py), then matched against a small demo-only HSQC
reference library with the real hsqc_matcher.py matching code
(match_sample_2d / score_compound_2d).

Why this bypasses read_bruker_2d(): read_bruker_2d() only knows how to
parse a real Bruker pdata folder (2rr + procs + proc2s binary files).
Fabricating a convincing fake one is high-effort and would mimic a
proprietary binary format rather than test anything about this pipeline's
own peak-picking or matching logic. This script builds the (h_scale,
c_scale, intensity_2d) array that read_bruker_2d() would have returned,
directly in memory, and feeds it to pick_peaks_2d() exactly the way
hsqc_analysis.py's process_hsqc_sample() does after its own
read_bruker_2d() call. Skipping only that one Bruker-specific read step
is the same choice simulate_alanine_standard_dilution_1d.R makes on the
1D side (see that file's own header).

Cross-peak positions used below: h_ppm matches the 1H demo's own CH3/CH
multiplet centers (1.47, 3.78 ppm, from the real shifts in
reference_library/reference_peaks_combined_ecoli_only.csv). The two 13C
values (16.9 and 51.3 ppm) are typical published shifts for L-alanine's
methyl and alpha carbons in aqueous solution -- NOT pulled from a live
BMRB query (this environment has no network access to BMRB's API), and
not independently re-verified against a primary source here. Treat them
as representative, not authoritative, and check them yourself (e.g. via
build_bmrb_library.py --hsqc once you have API access) before relying on
this demo's exact numbers for anything beyond trying the pipeline out.

No HSQC reference library ships in reference_library/ (see README's
"Optional: 2D HSQC support" section) -- so this script also writes a
tiny, clearly-labeled demo-only one
(example_data/hsqc/reference_hsqc_demo_alanine.csv), just large enough
for hsqc_matcher.py to have something real to match these simulated
cross-peaks against. This is NOT a substitute for the real
build_bmrb_library.py --hsqc library.

Regenerate with (from this src/ folder):
    python simulate_alanine_standard_dilution_hsqc.py

Outputs:
  example_data/hsqc/reference_hsqc_demo_alanine.csv
      the 2-row demo-only HSQC reference library described above
  example_results/hsqc/<sample_id>_hsqc_peaks.csv, _hsqc_plot.png
      per concentration, in exactly the shape hsqc_analysis.py produces
      for real data
  example_results/hsqc/match_results_hsqc_alanine_dilution.csv
      hsqc_matcher.py's real matching output across all 5 concentrations
  example_results/hsqc/alanine_dilution_hsqc_sensitivity.csv/.png
      each cross-peak's own picked intensity vs. concentration -- see
      this file's bottom section for why this, not matched_intensity_mean,
      is what's tracked here (score_compound_2d has no intensity output
      at all; see utils.py's own docstring).
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from peak_picking_2d import pick_peaks_2d  # noqa: E402
from hsqc_analysis import _plot_hsqc  # noqa: E402
from hsqc_matcher import match_all_samples_2d  # noqa: E402

CONCENTRATIONS_MM = [0.5, 1, 2, 5, 10]
SAMPLE_IDS = [f"std_{str(c).rstrip('0').rstrip('.') if '.' in str(c) else str(c)}mM"
              for c in CONCENTRATIONS_MM]

# (h_ppm, c_ppm, protons) -- protons sets the relative scaling between the
# two cross-peaks, same reasoning as the 1D script (3 for CH3, 1 for CH).
CROSS_PEAKS = {
    "CH3": (1.47, 16.9, 3),
    "CH": (3.78, 51.3, 1),
}

H_AXIS = np.linspace(0.5, 4.5, 300)   # direct (1H) dimension
C_AXIS = np.linspace(10.0, 60.0, 130)  # indirect (13C) dimension, coarser


def build_hsqc_matrix(conc_mM: float, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    sigma_h, sigma_c = 0.015, 0.8
    k = 400.0

    intensity = np.zeros((len(C_AXIS), len(H_AXIS)))
    hh, cc = np.meshgrid(H_AXIS, C_AXIS)
    for h0, c0, protons in CROSS_PEAKS.values():
        height = protons * k * conc_mM
        intensity += height * np.exp(-((hh - h0) ** 2) / (2 * sigma_h ** 2)
                                      - ((cc - c0) ** 2) / (2 * sigma_c ** 2))
    intensity += rng.normal(0, 15, size=intensity.shape)
    return intensity


def main():
    example_data = Path("../example_data/hsqc")
    example_results = Path("../example_results/hsqc")
    example_data.mkdir(parents=True, exist_ok=True)
    example_results.mkdir(parents=True, exist_ok=True)

    # Demo-only HSQC reference library -- see this file's header.
    ref_rows = [
        {"compound_id": "DEMO_ALANINE", "compound_name": "L-Alanine (demo reference, not from build_bmrb_library.py)",
         "h_ppm": h_ppm, "c_ppm": c_ppm}
        for h_ppm, c_ppm, _ in CROSS_PEAKS.values()
    ]
    ref_df = pd.DataFrame(ref_rows)
    ref_path = example_data / "reference_hsqc_demo_alanine.csv"
    ref_df.to_csv(ref_path, index=False)
    print(f"Wrote {ref_path}")

    all_peaks = []
    sensitivity_rows = []
    for sample_id, conc in zip(SAMPLE_IDS, CONCENTRATIONS_MM):
        seed = 20261018 + CONCENTRATIONS_MM.index(conc)
        intensity_2d = build_hsqc_matrix(conc, seed)

        peaks = pick_peaks_2d(H_AXIS, C_AXIS, intensity_2d,
                               threshold_multiplier=6.0,
                               min_separation_h_ppm=0.05,
                               min_separation_c_ppm=2.0,
                               exclude_h_ppm_ranges=[])
        peaks.insert(0, "sample_id", sample_id)
        peaks_path = example_results / f"{sample_id}_hsqc_peaks.csv"
        peaks.to_csv(peaks_path, index=False)
        print(f"{sample_id} ({conc} mM): {len(peaks)} cross-peak(s) picked -> {peaks_path}")
        all_peaks.append(peaks)

        _plot_hsqc(sample_id, H_AXIS, C_AXIS, intensity_2d, peaks, [],
                   example_results)

        # Sensitivity tracking: for each KNOWN cross-peak, the nearest
        # picked peak's own intensity (not a matcher output -- see header).
        for label, (h0, c0, _protons) in CROSS_PEAKS.items():
            if peaks.empty:
                picked_intensity = None
            else:
                dist = np.hypot((peaks["h_ppm"] - h0) / 0.05, (peaks["c_ppm"] - c0) / 2.0)
                best = dist.idxmin()
                picked_intensity = peaks.loc[best, "intensity"] if dist[best] <= 1.0 else None
            sensitivity_rows.append({"sample_id": sample_id, "concentration_mM": conc,
                                      "cross_peak": label, "picked_intensity": picked_intensity})

    experimental = pd.concat(all_peaks, ignore_index=True)

    # Real matcher, real reference-library shape, just pointed at the
    # demo-only reference built above.
    match_results = match_all_samples_2d(experimental, ref_df, top_n=10, min_ref_peaks=2)
    match_path = example_results / "match_results_hsqc_alanine_dilution.csv"
    match_results.to_csv(match_path, index=False)
    print(f"\nWrote {len(match_results)} matched row(s) to {match_path}")
    if not match_results.empty:
        print(match_results.to_string(index=False))
    else:
        print("No concentration matched both demo reference cross-peaks at "
              "min_ref_peaks=2 -- see the per-concentration peak counts above.")

    sens_df = pd.DataFrame(sensitivity_rows)
    sens_path = example_results / "alanine_dilution_hsqc_sensitivity.csv"
    sens_df.to_csv(sens_path, index=False)
    print(f"\nWrote {sens_path}")
    print(sens_df.to_string(index=False))

    _plot_sensitivity(sens_df, example_results / "alanine_dilution_hsqc_sensitivity.png")
    print(f"Wrote {example_results / 'alanine_dilution_hsqc_sensitivity.png'}")


def _plot_sensitivity(sens_df: pd.DataFrame, out_png: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6, 5))
    for label, grp in sens_df.groupby("cross_peak"):
        grp = grp.sort_values("concentration_mM")
        ax.plot(grp["concentration_mM"], grp["picked_intensity"], marker="o", label=label)
    ax.set_xlabel("concentration (mM)")
    ax.set_ylabel("picked cross-peak intensity (a.u.)")
    ax.set_title("L-Alanine HSQC cross-peaks -- picked intensity vs. concentration")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_png, dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    main()
