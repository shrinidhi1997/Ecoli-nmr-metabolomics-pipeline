# pipeline_standard_dilution_core.R
#
# The real, reusable body of run_pipeline_haird_real_dilution.R, pulled out
# into one function so the same code runs on any real N-point standard
# dilution series (HAIRD, SWEPF, or whatever comes next) without copying
# the whole script per dataset. See run_pipeline_haird_real_dilution.R's
# original header for the full rationale (direct Bruker reader standing in
# for ASICS, shape-detector-only feature table, etc.) -- that reasoning is
# unchanged here, just parameterized.
#
# Each per-dataset caller (run_pipeline_<name>_real_dilution.R) just sets
# dataset_name/data_root/sample_ids/label_value/out_dir and calls
# run_standard_dilution_pipeline(...).

if (!requireNamespace("nmrpick", quietly = TRUE)) {
  stop("Package 'nmrpick' is not installed. From an R session: ",
       'devtools::install_github("shrinidhi1997/Ecoli-nmr-metabolomics-pipeline", subdir = "nmrpick")')
}
library(nmrpick)

run_standard_dilution_pipeline <- function(dataset_name, data_root, sample_ids,
                                            label_value, out_dir,
                                            label_units_confirmed = FALSE,
                                            report_ppm_range = c(0, 11)) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # ---- STEP 0: import ----
  # Cropping the spectrum to a narrower ppm window BEFORE baseline
  # fitting was tried and reverted: it changes fit_airpls_baseline()'s
  # boundary behavior (its own endpoint-anchoring logic) enough to swing
  # detected peak counts around unpredictably (confirmed: HAIRD-60-1H
  # went from 182 to 27 peaks, not a small effect, and not one I
  # understood well enough to trust). Baseline fitting and peak
  # detection below still run on the FULL native digitized range, same
  # as the already-validated run. report_ppm_range instead filters
  # PEAKS/FEATURES AFTER detection, for reporting only -- see STEP 2b.
  message("Importing ", length(sample_ids), " real sample(s) from: ", data_root)
  raw <- lapply(sample_ids, function(s) read_bruker_1r(file.path(data_root, s)))
  names(raw) <- sample_ids
  for (s in sample_ids) {
    m <- raw[[s]]$meta
    message("  ", s, ": SI=", m$si, " points, ppm range [",
            round(min(raw[[s]]$ppm), 3), ", ", round(max(raw[[s]]$ppm), 3),
            "], SF=", round(m$sf, 6), " MHz")
  }
  for (s in sample_ids) {
    write.csv(data.frame(ppm = raw[[s]]$ppm, intensity = raw[[s]]$intensity),
              file.path(out_dir, paste0(s, "_raw_spectrum.csv")), row.names = FALSE)
  }
  message("Wrote raw spectra to ", out_dir)

  # ---- STEP 1: baseline (big-bin robust percentile, once per sample) ----
  # Uses fit_binned_baseline(), not fit_airpls_baseline(), as of this
  # version -- see that function's own header comment in
  # baseline_peak_picking.R for the full rationale and the numbers from
  # testing. Short version: airPLS's single global smoothness-penalized
  # solve let the water peak's hump bleed into distant, unrelated regions
  # (confirmed: baseline running up to ~20,000 a.u. above the real floor
  # several ppm away, and up to 89% above the real signal on the water
  # peak's own falling shoulder) -- exactly the "baseline going beyond
  # where it should / hallucinating" problem reported on real data. Big,
  # independent ppm bins with a low percentile and NO cross-bin smoothing
  # fix this by construction: nothing outside a bin can influence that
  # bin's own estimate, so there's no long-range (or even short-range)
  # bleed to begin with. Confirmed across all 6 real samples: far-field
  # overshoot dropped by roughly 5-10x, near-water overshoot dropped from
  # +60%..+89% to essentially 0 (-7%..+6%), with no new interpolation-edge
  # artifacts (checked: extra peaks found afterward are NOT clustered at
  # bin edges). Peak detection below still runs on the full native,
  # fine-grained spectrum exactly as before -- only the baseline that gets
  # subtracted from it changed, per the original "big bins for baseline
  # only, not for peak picking" ask.
  message("\nFitting baselines (big-bin robust percentile)...")
  baselines <- lapply(sample_ids, function(s) fit_binned_baseline(raw[[s]]$ppm, raw[[s]]$intensity))
  names(baselines) <- sample_ids
  for (s in sample_ids) {
    write.csv(data.frame(ppm = raw[[s]]$ppm, baseline = baselines[[s]]),
              file.path(out_dir, paste0(s, "_baseline.csv")), row.names = FALSE)
    write.csv(data.frame(ppm = raw[[s]]$ppm, corrected = raw[[s]]$intensity - baselines[[s]]),
              file.path(out_dir, paste0(s, "_corrected_spectrum.csv")), row.names = FALSE)
  }
  message("Wrote baselines and corrected spectra to ", out_dir)

  # ---- STEP 2: peak detection (shape detector; deriv pass for comparison only) ----
  prominence_mult <- 3
  min_height_mult <- 1.5
  deriv_prominence_mult <- 1.5
  deriv_min_height_mult <- 0.75
  deriv_dedup_tol_ppm <- 0.004
  alignment_tol_ppm <- 0.008

  message("\nDetecting peaks (prominence_mult=", prominence_mult,
          ", min_height_mult=", min_height_mult,
          ", plus second-derivative pass at prominence_mult=", deriv_prominence_mult,
          ", min_height_mult=", deriv_min_height_mult, ")...")
  shape_only <- lapply(sample_ids, function(s) {
    find_shape_based_peaks(raw[[s]]$ppm, raw[[s]]$intensity,
                            prominence_mult = prominence_mult,
                            min_height_mult = min_height_mult,
                            precomputed_baseline = baselines[[s]])[, c("ppm", "intensity")]
  })
  names(shape_only) <- sample_ids

  # Report/use only peaks within report_ppm_range (default 0-11 ppm,
  # the conventional small-molecule 1H window) -- baseline fitting and
  # detection above still ran on the full native range (see STEP 0's
  # note), this only trims what gets written out and merged. Any peak
  # dropped here is logged, not silently discarded.
  for (s in sample_ids) {
    before_n <- nrow(shape_only[[s]])
    keep <- shape_only[[s]]$ppm >= report_ppm_range[1] & shape_only[[s]]$ppm <= report_ppm_range[2]
    shape_only[[s]] <- shape_only[[s]][keep, , drop = FALSE]
    dropped <- before_n - nrow(shape_only[[s]])
    if (dropped > 0) {
      message("  ", s, ": dropped ", dropped, " peak(s) outside the ",
              report_ppm_range[1], "-", report_ppm_range[2], " ppm reporting range.")
    }
  }

  per_sample <- lapply(sample_ids, function(s) {
    shape_peaks <- shape_only[[s]]
    deriv_peaks <- find_peaks_second_derivative(raw[[s]]$ppm, raw[[s]]$intensity,
                                                 prominence_mult = deriv_prominence_mult,
                                                 min_height_mult = deriv_min_height_mult,
                                                 precomputed_baseline = baselines[[s]])
    is_new <- sapply(deriv_peaks$ppm, function(p) {
      nrow(shape_peaks) == 0 || all(abs(shape_peaks$ppm - p) > deriv_dedup_tol_ppm)
    })
    new_from_deriv <- deriv_peaks[is_new, c("ppm", "intensity")]
    combined <- rbind(shape_peaks[, c("ppm", "intensity")], new_from_deriv)
    # Same 0-11 ppm reporting filter as shape_only above, so the
    # detector-comparison chart compares both detectors' fragility
    # within the real spectrum, not conflated with out-of-range noise.
    combined <- combined[combined$ppm >= report_ppm_range[1] & combined$ppm <= report_ppm_range[2], , drop = FALSE]
    combined <- combined[order(-combined$ppm), ]
    rownames(combined) <- NULL
    attr(combined, "n_new_from_deriv") <- nrow(new_from_deriv)
    combined
  })
  names(per_sample) <- sample_ids

  detector_comparison <- data.frame(
    sample_id = sample_ids,
    n_shape_only = sapply(sample_ids, function(s) nrow(shape_only[[s]])),
    n_combined_with_deriv = sapply(sample_ids, function(s) nrow(per_sample[[s]]))
  )
  write.csv(detector_comparison, file.path(out_dir, "detector_comparison.csv"), row.names = FALSE)

  for (s in sample_ids) {
    message("  ", s, ": ", nrow(shape_only[[s]]), " peaks from the recommended shape detector alone; ",
            "+", attr(per_sample[[s]], "n_new_from_deriv"),
            " more if the second-derivative pass is also enabled (not used downstream -- see detector_comparison.csv).")
    write.csv(shape_only[[s]], file.path(out_dir, paste0(s, "_detected_peaks.csv")), row.names = FALSE)
  }
  message("Wrote per-sample detected peaks (shape-detector) and detector_comparison.csv to ", out_dir)

  # ---- STEP 3: merge across samples (reference-guided) ----
  message("\nMerging across samples (reference-guided, tol=", alignment_tol_ppm, " ppm)...")
  spectra_list <- lapply(sample_ids, function(s) list(ppm = raw[[s]]$ppm, corrected = raw[[s]]$intensity - baselines[[s]]))
  names(spectra_list) <- sample_ids
  merged <- merge_peaks_reference_guided(shape_only, spectra_list, tol = alignment_tol_ppm)
  message("Merged feature table: ", nrow(merged), " features x ", ncol(merged) - 1, " samples")

  # ---- STEP 4: write feature table + concentrations lookup ----
  groups <- setNames(rep(paste0(dataset_name, "_standard"), length(sample_ids)), sample_ids)
  out_file <- file.path(out_dir, paste0("peak_list_", dataset_name, "_real_dilution.csv"))
  write_metaboanalyst_csv(merged, groups, out_file)
  message("Wrote ", out_file)

  write.csv(data.frame(sample_id = sample_ids, concentration_mM = label_value[sample_ids]),
            file.path(out_dir, "concentrations.csv"), row.names = FALSE)
  units_note <- if (label_units_confirmed) {
    "units confirmed (final concentration)."
  } else {
    "UNITS NOT CONFIRMED -- the numeric label is taken directly from the sample's own folder name; treat the x-axis on any linearity plot as 'dilution label', not a confirmed concentration, until the unit is confirmed."
  }
  message("Wrote ", file.path(out_dir, "concentrations.csv"), " -- ", units_note)

  invisible(list(raw = raw, baselines = baselines, shape_only = shape_only,
                 detector_comparison = detector_comparison, merged = merged,
                 feature_table_path = out_file))
}
