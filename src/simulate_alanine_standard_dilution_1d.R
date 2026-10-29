# simulate_alanine_standard_dilution_1d.R
#
# Generates a SIMULATED (not real instrument data) 1D 1H standard dilution
# series for L-alanine, runs it through the real peak-picking pipeline in
# this src/ folder (the same functions run_peak_picking_control_vs_sample.R
# uses on real Bruker data), and writes the resulting synthetic spectra +
# feature table as example data and example results for this repository.
#
# Why this exists: a new user of this pipeline has no way to try it before
# they have their own staged Bruker data (see README's "Input data layout").
# This script, and its 2D counterpart
# (simulate_alanine_standard_dilution_hsqc.py), give a self-contained,
# honest demo: synthetic spectra, run through the UNMODIFIED real detection
# code, producing real (if synthetic-sourced) output files.
#
# Why a dilution series, not a control-vs-treatment pair: this mirrors a
# standards validation run (one compound, several known concentrations) --
# the thing you'd actually do to check peak picking and matching scale
# sensibly with concentration before trusting them on a real two-group
# biological comparison. It is NOT a replacement for the real
# run_peak_picking_control_vs_sample.R demo; it exercises the same
# underlying functions a different way.
#
# Why L-alanine, specifically: it's a simple, two-carbon-environment
# compound (a CH3 doublet and a CH quartet) that is ALREADY in the shipped
# reference_library/reference_peaks_combined_ecoli_only.csv, under
# compound_id HMDB0000161. The six 1H shift values used below
# (1.459, 1.474, 3.751, 3.766, 3.780, 3.795 ppm) are copied directly from
# that real, shipped file, not invented -- so the downstream matching step
# (match_r_feature_table.py, run unmodified against that same reference
# CSV) has a real, pre-existing entry to match against, exactly as it
# would for real data.
#
# What's simulated, and how: for each concentration, each of alanine's two
# multiplets is built from individual Lorentzian lines at its real
# reference shifts, with textbook-correct relative intensities within each
# multiplet (CH3 doublet: 1:1, since CH3 couples to 1 CH proton; CH
# quartet: 1:3:3:1, since CH couples to 3 equivalent CH3 protons), and the
# CH3 multiplet's total area set to 3x the CH multiplet's (3 protons vs 1),
# both scaled linearly with concentration. A smooth sinusoidal baseline
# drift (different seed per concentration) and Gaussian noise (sd fixed
# across concentrations) are added on top.
#
# Only find_shape_based_peaks() (this file's own "RECOMMENDED" detector)
# is used below, not find_peaks_second_derivative(). That second detector
# was tried first, on this exact synthetic data, and it did not work: it
# returned several hundred to nearly a thousand spurious peaks scattered
# across completely flat, signal-free regions of the very same spectrum.
# That is not a bug introduced here -- find_peaks_second_derivative()'s
# own docstring in baseline_peak_picking.R already documents this
# detector's thresholds as the more fragile of the two ("honest
# trade-offs... amplifies noise... not a clear, reliable win"), and this
# synthetic spectrum's point-to-point noise, at this ppm grid spacing,
# lands squarely in the regime that docstring warns about. Rather than
# retune a secondary, already-documented-as-fragile detector against one
# synthetic benchmark (exactly the kind of single-benchmark overfitting
# that same docstring cautions against), this demo uses the primary,
# recommended detector alone and reports that choice here plainly.
#
# What find_shape_based_peaks() actually found, verified by running this
# script (see printed per-concentration peak counts, and
# example_results/1d/ for the saved numbers): at 0.5 and 1 mM, only 4 of
# alanine's 6 real reference peaks clear the default detection threshold
# -- the two larger, inner CH quartet lines (3.766, 3.780 ppm) and both
# CH3 doublet lines, but NOT the two smaller, outer quartet lines (3.751,
# 3.795 ppm). From 2 mM up, all 6 are detected. This is a genuine,
# honest sensitivity limit of the default thresholds on this synthetic
# noise level, not a flaw to hide -- exactly the kind of finding a real
# standard dilution series is meant to surface before you trust this
# pipeline's sensitivity on real low-concentration samples.
#
# Regenerate with (from this src/ folder):
#   Rscript simulate_alanine_standard_dilution_1d.R
#
# Outputs:
#   example_data/1d/synthetic_spectra_alanine_dilution.csv
#       ppm column + one intensity column per concentration sample
#   example_data/1d/concentrations.csv
#       sample_id, concentration_mM -- the ground truth this demo's
#       linearity check (analyze_standard_dilution_linearity.py) reads
#   example_results/1d/peak_list_alanine_dilution.csv
#       write_metaboanalyst_csv()'s output -- same shape
#       run_peak_picking_control_vs_sample.R produces for real data

if (!requireNamespace("nmrpick", quietly = TRUE)) {
  stop("Package 'nmrpick' is not installed. From an R session: ",
       'devtools::install_github("shrinidhi1997/Ecoli-nmr-metabolomics-pipeline", subdir = "nmrpick")')
}
library(nmrpick)

concentrations_mM <- c(0.5, 1, 2, 5, 10)
sample_ids <- paste0("std_", sub("\\.0$", "", formatC(concentrations_mM, format = "f", digits = 1)), "mM")
names(concentrations_mM) <- sample_ids

ppm <- seq(-0.2, 10.2, by = 0.0005)

lorentzian <- function(x, x0, height, hwhm) height / (1 + ((x - x0) / hwhm)^2)

# Real shifts copied from reference_library/reference_peaks_combined_ecoli_only.csv,
# compound_id HMDB0000161 (L-Alanine), nucleus 1H.
ch3_ppm <- c(1.459, 1.474)
ch3_wt  <- c(1, 1) / 2
ch_ppm  <- c(3.751, 3.766, 3.780, 3.795)
ch_wt   <- c(1, 3, 3, 1) / 8
hwhm <- 0.0028  # ~0.0056 ppm FWHM -- sharp, well-resolved lines

build_spectrum <- function(ppm, conc_mM, seed_offset) {
  set.seed(20261018 + seed_offset)
  unit <- 180 * conc_mM

  y <- rep(0, length(ppm))
  for (i in seq_along(ch3_ppm)) y <- y + lorentzian(ppm, ch3_ppm[i], 3 * unit * ch3_wt[i], hwhm)
  for (i in seq_along(ch_ppm))  y <- y + lorentzian(ppm, ch_ppm[i], 1 * unit * ch_wt[i], hwhm)

  drift <- 15 * sin(ppm / 3 + seed_offset) + 6 * cos(ppm / 1.7 - seed_offset * 0.5)
  y <- y + drift + 40

  y + rnorm(length(ppm), mean = 0, sd = 6)
}

cat("Building", length(sample_ids), "synthetic spectra...\n")
spec_mat <- sapply(seq_along(sample_ids), function(i) {
  build_spectrum(ppm, concentrations_mM[i], seed_offset = i)
})
colnames(spec_mat) <- sample_ids

dir.create("../example_data/1d", showWarnings = FALSE, recursive = TRUE)
dir.create("../example_results/1d", showWarnings = FALSE, recursive = TRUE)

spectra_out <- cbind(ppm = ppm, as.data.frame(spec_mat))
write.csv(spectra_out, "../example_data/1d/synthetic_spectra_alanine_dilution.csv",
          row.names = FALSE)
write.csv(data.frame(sample_id = sample_ids, concentration_mM = concentrations_mM,
                      row.names = NULL),
          "../example_data/1d/concentrations.csv", row.names = FALSE)
cat("Wrote example_data/1d/synthetic_spectra_alanine_dilution.csv and concentrations.csv\n")

# Fit each sample's baseline once, then run ONLY find_shape_based_peaks()
# (the RECOMMENDED detector) -- see this file's header for why
# find_peaks_second_derivative() is deliberately not used here.
cat("\nFitting baselines and detecting peaks...\n")
baselines <- lapply(sample_ids, function(s) fit_airpls_baseline(ppm, spec_mat[, s]))
names(baselines) <- sample_ids

per_sample <- lapply(sample_ids, function(s) {
  shape_peaks <- find_shape_based_peaks(ppm, spec_mat[, s], precomputed_baseline = baselines[[s]])
  shape_peaks[order(-shape_peaks$ppm), c("ppm", "intensity")]
})
names(per_sample) <- sample_ids
for (s in sample_ids) {
  cat("  ", s, " (", concentrations_mM[s], "mM): ", nrow(per_sample[[s]]), " peaks\n", sep = "")
}

spectra_list <- lapply(sample_ids, function(s) {
  list(ppm = ppm, corrected = spec_mat[, s] - baselines[[s]])
})
names(spectra_list) <- sample_ids

cat("\nMerging across concentration samples (reference-guided)...\n")
merged <- merge_peaks_reference_guided(per_sample, spectra_list, tol = 0.008)
cat("Merged feature table:", nrow(merged), "features x", ncol(merged) - 1, "samples\n")

groups <- setNames(rep("alanine_standard", length(sample_ids)), sample_ids)
write_metaboanalyst_csv(merged, groups, "../example_results/1d/peak_list_alanine_dilution.csv")

cat("\nNext: from a shell, run `nmrmatch-match-feature-table` (or the equivalent\n",
    "python match_r_feature_table.py) against this feature table and the shipped\n",
    "E. coli reference library (see README's new demo section for the exact\n",
    "command), then analyze_standard_dilution_linearity.py.\n", sep = "")
