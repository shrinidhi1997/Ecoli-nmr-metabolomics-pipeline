# run_peak_picking_control_vs_sample.R
#
# Full custom-pipeline run for a two-group comparison, using
# find_shape_based_peaks()/refine_peaks_lorentzian()/merge_peaks_reference_guided()
# from this src/ folder instead of ASICS's own quantification. No ASICS()
# call anywhere in this script; ASICS is used only for importSpectraBruker()/
# createSpectra() to read raw Bruker spectra.
#
# All dataset-specific settings (where your data is, which sample goes in
# which group, detection thresholds) live in config_control_vs_sample.R, not
# here. Edit that file for a new dataset; this script should not need to
# change. See README.md for the input folder layout this expects.
#
# Run from this src/ folder:
#   source("run_peak_picking_control_vs_sample.R")
#
# Two detection knobs worth understanding before you change them in
# config_control_vs_sample.R:
#
#   INTENSITY = real peak AREA, not raw height, when use_lorentzian_
#      refinement is TRUE. With it FALSE (this pipeline's default),
#      "intensity" is the plain baseline-corrected HEIGHT at each
#      detected peak -- fast, no per-peak/per-cluster curve fitting at
#      all -- but a real accuracy cost: height is NOT proportional to
#      molar amount when linewidth differs between samples (shimming,
#      viscosity, temperature), area is. Set it to TRUE if a fold-change
#      result looks off and you want to rule out linewidth differences
#      as the cause.
#
#   The cross-SAMPLE alignment tolerance (merge_peaks_reference_guided's
#      `tol`, set from alignment_tol_ppm) is a DIFFERENT knob from the ppm
#      tolerance you'd raise for IDENTIFICATION (match_r_feature_table.py's
#      --tolerance, used in a later step). This one only controls how
#      close another sample's own peak has to be to call it "the same
#      feature" when building the merged table across your samples --
#      raise alignment_tol_ppm if your samples show more run-to-run drift
#      than the default; raise --tolerance on the Python side separately
#      for identification drift against the reference library.

source("config_control_vs_sample.R")
if (!requireNamespace("nmrpick", quietly = TRUE)) {
  stop("Package 'nmrpick' is not installed. From an R session: ",
       'devtools::install_github("shrinidhi1997/Ecoli-nmr-metabolomics-pipeline", subdir = "nmrpick")')
}
library(nmrpick)

library(ASICS)  # only for importSpectraBruker()/createSpectra() -- no ASICS() call below

dir.create(out_base, showWarnings = FALSE, recursive = TRUE)

if (!dir.exists(data_dir)) {
  stop("Can't find this folder: ", data_dir,
       "\nEdit data_dir in config_control_vs_sample.R, and see README.md's ",
       "\"Staging real Bruker data\" section for the junction/symlink commands ",
       "that build this staging folder from raw Bruker data.")
}

message("Importing ", length(sample_group_map), " samples from: ", data_dir)
my_spectra <- importSpectraBruker(data_dir)
message(ncol(my_spectra), " sample(s) imported: ", paste(colnames(my_spectra), collapse = ", "))

if (!all(colnames(my_spectra) %in% names(sample_group_map))) {
  stop("These imported sample names don't match sample_group_map's names: ",
       paste(setdiff(colnames(my_spectra), names(sample_group_map)), collapse = ", "))
}

spectra_obj <- createSpectra(my_spectra)
ppm <- getPpmGrid(spectra_obj)
spec_mat <- getSpectra(spectra_obj)
if (nrow(spec_mat) != length(ppm)) spec_mat <- t(spec_mat)
sample_names <- colnames(spec_mat)

# ============================================================
# STEP 1 -- fit each sample's baseline ONCE, up front, and reuse it
# everywhere downstream (both detectors in STEP 2, and Lorentzian
# refinement / the merge fallback check later).
#
# This used to run THREE separate times per sample: once silently inside
# find_shape_based_peaks(), once silently inside
# find_peaks_second_derivative(), and once again explicitly here for
# refinement/merging -- airPLS is the single most expensive step in this
# whole pipeline (an iterative reweighted fit), so with the second
# detector added, every sample was paying for it 3x for no reason, the
# exact same deterministic computation on the exact same data every
# time. This is the real fix for that (not a different, cheaper, less
# accurate baseline algorithm -- this IS airPLS, just run once instead
# of three times): ~3x less baseline-fitting work per sample, zero
# accuracy change, since it's the identical computation being reused
# instead of repeated.
#
# This also still picks up fit_airpls_baseline()'s current default
# (max_iter=4) -- changed from 15 after real testing found the OLD,
# fully-converged baseline runs systematically too low in crowded
# regions (it over-commits to a sparse, biased-low set of anchor points
# once there are very few genuinely peak-free stretches left).
# ============================================================
message("\nFitting baselines (once per sample, reused by both detectors below)...")
baselines <- lapply(sample_names, function(s) fit_airpls_baseline(ppm, spec_mat[, s]))
names(baselines) <- sample_names

# ============================================================
# STEP 2 -- detect peaks per sample at the lowered threshold, passing
# each sample's already-fit baseline in so neither detector refits it.
# ============================================================
message("\nDetecting peaks (prominence_mult=", prominence_mult,
        ", min_height_mult=", min_height_mult,
        ", plus second-derivative pass at prominence_mult=", deriv_prominence_mult,
        ", min_height_mult=", deriv_min_height_mult, ")...")
per_sample <- lapply(sample_names, function(s) {
  shape_peaks <- find_shape_based_peaks(ppm, spec_mat[, s],
                                         prominence_mult = prominence_mult,
                                         min_height_mult = min_height_mult,
                                         precomputed_baseline = baselines[[s]])
  deriv_peaks <- find_peaks_second_derivative(ppm, spec_mat[, s],
                                               prominence_mult = deriv_prominence_mult,
                                               min_height_mult = deriv_min_height_mult,
                                               precomputed_baseline = baselines[[s]])
  # Keep every shape-detector peak as-is (that detector's own height/width
  # handling is what refine_peaks_lorentzian() below expects). From the
  # derivative detector, keep ONLY the peaks that are genuinely NEW --
  # more than deriv_dedup_tol_ppm from anything the shape detector already
  # found -- so a peak both detectors agree on isn't counted twice.
  is_new <- sapply(deriv_peaks$ppm, function(p) {
    nrow(shape_peaks) == 0 || all(abs(shape_peaks$ppm - p) > deriv_dedup_tol_ppm)
  })
  new_from_deriv <- deriv_peaks[is_new, c("ppm", "intensity")]
  combined <- rbind(shape_peaks[, c("ppm", "intensity")], new_from_deriv)
  combined <- combined[order(-combined$ppm), ]
  rownames(combined) <- NULL
  attr(combined, "n_shape") <- nrow(shape_peaks)
  attr(combined, "n_new_from_deriv") <- nrow(new_from_deriv)
  combined
})
names(per_sample) <- sample_names
for (s in sample_names) {
  message("  ", s, ": ", nrow(per_sample[[s]]), " peaks (",
          attr(per_sample[[s]], "n_shape"), " shape-detector, +",
          attr(per_sample[[s]], "n_new_from_deriv"), " new from second-derivative pass)")
}

# ============================================================
# STEP 3 -- Lorentzian refinement, then carry AREA through as the
# table's "intensity" column (see note 2 above for why) -- SKIPPED when
# use_lorentzian_refinement is FALSE, see that flag's own comment above
# for the real speed-vs-accuracy trade-off this switches.
# ============================================================
if (use_lorentzian_refinement) {
  message("\nRefining peaks (Lorentzian fit for real area)...")
  per_sample_qnmr <- lapply(sample_names, function(s) {
    corrected <- spec_mat[, s] - baselines[[s]]
    refined <- refine_peaks_lorentzian(ppm, corrected, per_sample[[s]])
    n_failed <- sum(!refined$converged)
    if (n_failed > 0) {
      message("  ", s, ": ", n_failed, " peak(s) didn't converge to a clean Lorentzian fit -- ",
              "falling back to the originally detected height for just those, area for the rest.")
    }
    # area is NA exactly where converged is FALSE (refine_peaks_lorentzian's own
    # behaviour) -- fall back to the plain detected height there rather than
    # dropping the peak outright, so a hard-to-fit but real peak still gets a
    # usable (if less precise) intensity value instead of disappearing.
    area_or_height <- ifelse(is.na(refined$area), refined$intensity, refined$area)
    data.frame(ppm = refined$ppm, intensity = area_or_height)
  })
} else {
  message("\nSkipping Lorentzian refinement (use_lorentzian_refinement = FALSE) -- ",
          "intensity below is plain detected peak height, not fitted area.")
  per_sample_qnmr <- lapply(sample_names, function(s) {
    per_sample[[s]][, c("ppm", "intensity")]
  })
}
names(per_sample_qnmr) <- sample_names

# ============================================================
# STEP 4 -- merge across samples, reference-guided, with the
# baseline-corrected spectra available for the raw-signal fallback.
# ============================================================
message("\nMerging across samples (reference-guided, tol=", alignment_tol_ppm, " ppm)...")
spectra_list <- lapply(sample_names, function(s) {
  list(ppm = ppm, corrected = spec_mat[, s] - baselines[[s]])
})
names(spectra_list) <- sample_names

merged <- merge_peaks_reference_guided(per_sample_qnmr, spectra_list, tol = alignment_tol_ppm)
message("Merged feature table: ", nrow(merged), " features x ", ncol(merged) - 1, " samples")

# ============================================================
# STEP 5 -- write the feature table match_r_feature_table.py reads.
# ============================================================
out_file <- file.path(out_base, "peak_list_control_vs_sample.csv")
write_metaboanalyst_csv(merged, sample_group_map, out_file)

# Built from control_samples/test_samples (computed in config_control_vs_sample.R
# from your actual group_labels), not retyped and not a hardcoded "control"/
# "sample" string match, so this printed command stays correct whatever two
# group labels you chose there, and never goes stale when the sample set
# changes.
control_names_str <- paste(control_samples, collapse = ",")
test_names_str <- paste(test_samples, collapse = ",")

message("\nDone. Next, from a shell (nmrmatch installed via `pip install -e ./nmrmatch` ",
        "or from GitHub -- see nmrmatch/README.md):\n",
        "  nmrmatch-match-feature-table \\\n",
        "      --feature-table ", out_file, " \\\n",
        "      --reference ../output/reference_peaks_combined_ecoli_only.csv \\\n",
        "      --tolerance 0.008 \\\n",
        "      --top-n 500 \\\n",
        "      --min-confidence 0.88 \\\n",
        "      --out ../output/control_vs_sample/compound_matches_control_vs_sample.csv\n",
        "  (equivalently: python match_r_feature_table.py with the same flags, or ",
        "python -m nmrmatch.cli.match_feature_table)\n",
        "\nThen EITHER (a), the full R stats+plots script (recommended -- t-test/",
        "Wilcoxon/BH-FDR plus 8 figures, reads control_samples/test_samples from ",
        "config_control_vs_sample.R, same as this script):\n",
        "  source(\"run_stats_and_plots_control_vs_sample.R\")\n",
        "\nOR (b), just the plain fold-change numbers with no R session and no plots:\n",
        "  nmrmatch-fold-change \\\n",
        "      --matches ../output/control_vs_sample/compound_matches_control_vs_sample.csv \\\n",
        "      --group control:", control_names_str, " \\\n",
        "      --group sample:", test_names_str, " \\\n",
        "      --out ../output/control_vs_sample/fold_induction_control_vs_sample.csv")
