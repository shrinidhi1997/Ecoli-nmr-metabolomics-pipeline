# run_multiplet_overlap_benchmark.R
#
# Runs THIS project's own, unmodified detection pipeline (fit_airpls_
# baseline() -> find_shape_based_peaks() -> refine_peaks_lorentzian(),
# the exact functions every other script in src/ uses on real data) on a
# synthetic triplet+quartet benchmark, and scores the result against
# KNOWN ground truth -- the same comparison MultiNMRFit's own paper makes
# (accuracy/precision of x, J, I, lw across controlled Hz separation and
# S/N), but for this pipeline instead of theirs.
#
# INPUT CONTRACT -- reads exactly two files, column names below required:
#   <in_dir>/ground_truth.csv
#     draw_id, hz_sep, snr_level, field_mhz, multiplet_type, line_index,
#     true_ppm, true_J_hz, true_lw_hz (FWHM, Hz), true_height, true_area
#   <in_dir>/spectra_long.csv
#     draw_id, ppm, intensity
# simulate_multiplet_overlap_benchmark.R writes exactly this; swap in a
# different generator (e.g. a more detailed one of your own) and this
# script does not need to change AT ALL, as long as it writes the same
# two files. intensity may already include noise and/or drift -- this
# script fits its own baseline per draw either way, same as it would on
# real data, so a nonzero simulated baseline is handled, not assumed away.
#
# WHY TWO HWHM BOUND SETTINGS, NOT ONE: refine_peaks_lorentzian()'s
# lower linewidth bound used to be hardcoded at 0.001 ppm -- fine for
# this project's real data (always >=0.007-0.01 ppm), but a real
# constraint if the benchmark's own true linewidth sits near or below
# it (1 Hz FWHM = 0.000833 ppm HWHM at 600.17 MHz, already below 0.001).
# Fitting with the old bound in that regime would measure the bound, not
# the pipeline's real resolving power. Now exposed as a parameter (see
# baseline_peak_picking.R's own comment on refine_peaks_lorentzian()), so
# this script reports BOTH: "production_bound" (the exact 0.001-0.05 ppm
# bound every other script in this project actually ships with) and
# "widened_bound" (a much looser floor, 0.00005 ppm, to see what the
# detector/fit can do without that constraint). Report both, don't just
# pick the flattering one -- the gap between them is itself a finding.
#
# MATCHING, one true line to one fitted peak: nearest-ppm, greedy,
# without replacement (a fitted peak claimed by one true line can't also
# be claimed by another), within a tolerance tied to the multiplet's own
# J -- tight enough that a match can't accidentally jump to the NEXT
# line over (adjacent true lines within one multiplet are exactly J Hz
# apart by construction). A true line with no fitted peak left within
# tolerance is UNMATCHED: either merged into a neighbour (the detector
# found one peak where two true lines sit) or not detected at all
# (below threshold). Either way it counts against that multiplet's
# "resolved" flag below -- the benchmark's headline metric: where the
# detector starts losing distinct peaks.
#
# J is NEVER fit directly -- this pipeline has no shared-multiplet model
# (each Lorentzian in refine_peaks_lorentzian()'s cluster fit is
# independent; see that function's own docstring). It is DERIVED, only
# when every expected line in a multiplet resolved to its own distinct
# fitted peak, as the mean consecutive spacing (in Hz) between that
# multiplet's SORTED fitted centers -- flagged j_derived = TRUE in the
# output so this is never mistaken for an independent measurement the
# way MultiNMRFit's own constrained-model J actually is.
#
# Run from this src/ folder (after simulate_multiplet_overlap_benchmark.R,
# or after the real generator has written the same two files):
#   Rscript run_multiplet_overlap_benchmark.R

if (!requireNamespace("nmrpick", quietly = TRUE)) {
  stop("Package 'nmrpick' is not installed. From an R session: ",
       'devtools::install_github("shrinidhi1997/Ecoli-nmr-metabolomics-pipeline", subdir = "nmrpick")')
}
library(nmrpick)

run_multiplet_overlap_benchmark <- function(
    in_dir  = "../example_data/multiplet_benchmark",
    out_dir = "../example_results/multiplet_benchmark",
    prominence_mult  = 4,
    min_height_mult  = 1.5,
    cluster_tol      = 0.03,
    match_tol_frac_of_J = 0.4   # match tolerance = this fraction of J, in ppm
) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  ground_truth <- read.csv(file.path(in_dir, "ground_truth.csv"))
  spectra_long <- read.csv(file.path(in_dir, "spectra_long.csv"))

  bound_settings <- list(
    production_bound = c(0.001, 0.05),   # this project's real, shipped default
    widened_bound    = c(0.00005, 0.05)  # floor relaxed, to see real resolving power
  )

  draw_ids <- sort(unique(ground_truth$draw_id))
  per_line_rows <- list()

  for (bound_name in names(bound_settings)) {
    hwhm_bounds <- bound_settings[[bound_name]]
    message("\n=== bound setting: ", bound_name, " (hwhm in [",
            hwhm_bounds[1], ", ", hwhm_bounds[2], "] ppm) ===")

    for (did in draw_ids) {
      gt <- ground_truth[ground_truth$draw_id == did, ]
      sp <- spectra_long[spectra_long$draw_id == did, ]
      sp <- sp[order(sp$ppm), ]
      ppm <- sp$ppm; intensity <- sp$intensity
      field_mhz <- gt$field_mhz[1]
      j_hz <- gt$true_J_hz[1]
      match_tol_ppm <- match_tol_frac_of_J * j_hz / field_mhz

      baseline <- fit_airpls_baseline(ppm, intensity)
      corrected <- intensity - baseline
      peaks <- find_shape_based_peaks(ppm, intensity, exclude = list(),
                                       prominence_mult = prominence_mult,
                                       min_height_mult = min_height_mult,
                                       precomputed_baseline = baseline)
      fitted <- refine_peaks_lorentzian(ppm, corrected, peaks,
                                         cluster_tol = cluster_tol,
                                         hwhm_bounds_ppm = hwhm_bounds)

      # greedy nearest-ppm matching, one true line <-> one fitted peak,
      # without replacement, across BOTH multiplets at once (so a
      # triplet line and a quartet line competing for the same fitted
      # peak -- the actual overlap failure mode -- is handled correctly:
      # whichever true line is closer claims it, the other stays
      # unmatched, exactly reflecting a real merge).
      gt_ord <- gt[order(gt$true_ppm), ]
      avail <- rep(TRUE, nrow(fitted))
      matched_fit_idx <- rep(NA_integer_, nrow(gt_ord))
      if (nrow(fitted) > 0) {
        dist_mat <- outer(gt_ord$true_ppm, fitted$ppm, function(a, b) abs(a - b))
        repeat {
          dist_mat_avail <- dist_mat
          dist_mat_avail[, !avail] <- Inf
          dist_mat_avail[!is.na(matched_fit_idx), ] <- Inf
          if (all(!is.finite(dist_mat_avail))) break
          best <- which(dist_mat_avail == min(dist_mat_avail), arr.ind = TRUE)[1, ]
          if (dist_mat_avail[best[1], best[2]] > match_tol_ppm) break
          matched_fit_idx[best[1]] <- best[2]
          avail[best[2]] <- FALSE
        }
      }

      for (i in seq_len(nrow(gt_ord))) {
        row <- gt_ord[i, ]
        fidx <- matched_fit_idx[i]
        matched <- !is.na(fidx)
        per_line_rows[[length(per_line_rows) + 1]] <- data.frame(
          bound_setting = bound_name, draw_id = did, hz_sep = row$hz_sep,
          snr_level = row$snr_level, multiplet_type = row$multiplet_type,
          line_index = row$line_index, true_ppm = row$true_ppm,
          true_lw_hz = row$true_lw_hz, true_area = row$true_area,
          matched = matched,
          fit_ppm = if (matched) fitted$ppm[fidx] else NA_real_,
          fit_area = if (matched) fitted$area[fidx] else NA_real_,
          fit_hwhm_ppm = if (matched) fitted$hwhm[fidx] else NA_real_,
          fit_converged = if (matched) fitted$converged[fidx] else NA
        )
      }
    }
  }

  per_line <- do.call(rbind, per_line_rows)
  per_line$field_mhz <- ground_truth$field_mhz[1]
  per_line$err_x_hz  <- (per_line$fit_ppm - per_line$true_ppm) * per_line$field_mhz
  per_line$err_I     <- per_line$fit_area - per_line$true_area
  per_line$err_lw_hz <- (2 * per_line$fit_hwhm_ppm * per_line$field_mhz) - per_line$true_lw_hz
  write.csv(per_line, file.path(out_dir, "per_line_results.csv"), row.names = FALSE)
  message("\nWrote ", nrow(per_line), " per-line result rows to ",
          file.path(out_dir, "per_line_results.csv"))

  # resolved: EVERY expected line of this multiplet, in this draw,
  # matched a DISTINCT fitted peak -- the "did it lose a distinct peak"
  # flag, at the multiplet level.
  resolved_tbl <- aggregate(matched ~ bound_setting + draw_id + hz_sep + snr_level + multiplet_type,
                             data = per_line, FUN = all)
  names(resolved_tbl)[names(resolved_tbl) == "matched"] <- "resolved"

  # J, derived ONLY for resolved multiplets, from sorted fitted centers.
  j_rows <- list()
  for (bound_name in names(bound_settings)) {
    sub_all <- per_line[per_line$bound_setting == bound_name, ]
    for (did in draw_ids) {
      for (mtype in c("triplet", "quartet")) {
        sub <- sub_all[sub_all$draw_id == did & sub_all$multiplet_type == mtype, ]
        if (nrow(sub) == 0) next
        resolved <- all(sub$matched)
        j_val <- NA_real_
        if (resolved) {
          centers <- sort(sub$fit_ppm)
          j_val <- mean(diff(centers)) * sub$field_mhz[1]
        }
        j_rows[[length(j_rows) + 1]] <- data.frame(
          bound_setting = bound_name, draw_id = did, hz_sep = sub$hz_sep[1],
          snr_level = sub$snr_level[1], multiplet_type = mtype,
          true_J_hz = 7,  # matches the generator's fixed J; override if your generator varies it
          j_derived_hz = j_val, j_derived = TRUE, resolved = resolved
        )
      }
    }
  }
  j_tbl <- do.call(rbind, j_rows)
  j_tbl$err_J_hz <- j_tbl$j_derived_hz - j_tbl$true_J_hz
  write.csv(j_tbl, file.path(out_dir, "j_derived_results.csv"), row.names = FALSE)
  message("Wrote ", nrow(j_tbl), " J-derivation row(s) (NA where not resolved) to ",
          file.path(out_dir, "j_derived_results.csv"))

  # ---- aggregate: bias (mean error) and precision (SD) per cell ----
  # bias_median alongside the mean-based bias/sd: found directly while
  # building this benchmark -- a handful of catastrophic high-noise fits
  # (e.g. linewidth fit to 30 Hz FWHM against a true 1 Hz line) can
  # dominate a MEAN enough to hide a real, separate effect that's
  # actually visible in the median (refine_peaks_lorentzian()'s
  # production linewidth floor: 74% of matched/resolved lines landed
  # EXACTLY on the 0.001 ppm wall in one real run of this benchmark, yet
  # production_bound vs. widened_bound's mean bias barely differed,
  # because a few extreme outliers swamped the average either way).
  # Report both rather than only the paper's mean-based convention.
  agg_err <- function(df, value_col) {
    stats::aggregate(df[[value_col]], by = list(bound_setting = df$bound_setting,
                                                 hz_sep = df$hz_sep, snr_level = df$snr_level,
                                                 multiplet_type = df$multiplet_type),
                      FUN = function(v) c(bias = mean(v, na.rm = TRUE), sd = stats::sd(v, na.rm = TRUE),
                                           bias_median = stats::median(v, na.rm = TRUE),
                                           n = sum(!is.na(v))))
  }
  make_long <- function(agg, param_name) {
    m <- as.matrix(agg$x)
    data.frame(bound_setting = agg$bound_setting, hz_sep = agg$hz_sep,
               snr_level = agg$snr_level, multiplet_type = agg$multiplet_type,
               parameter = param_name, bias = m[, "bias"], sd = m[, "sd"],
               bias_median = m[, "bias_median"], n = m[, "n"])
  }

  summary_rows <- rbind(
    make_long(agg_err(per_line, "err_x_hz"), "x_hz"),
    make_long(agg_err(per_line, "err_I"), "I_area"),
    make_long(agg_err(per_line, "err_lw_hz"), "lw_hz"),
    make_long(agg_err(j_tbl, "err_J_hz"), "J_hz")
  )

  resolution_rate <- stats::aggregate(resolved ~ bound_setting + hz_sep + snr_level + multiplet_type,
                                       data = resolved_tbl, FUN = mean)
  names(resolution_rate)[names(resolution_rate) == "resolved"] <- "resolution_rate"

  write.csv(summary_rows, file.path(out_dir, "accuracy_precision_summary.csv"), row.names = FALSE)
  write.csv(resolution_rate, file.path(out_dir, "resolution_rate.csv"), row.names = FALSE)
  message("Wrote ", file.path(out_dir, "accuracy_precision_summary.csv"), " and ",
          file.path(out_dir, "resolution_rate.csv"))

  invisible(list(per_line = per_line, j_tbl = j_tbl,
                  summary = summary_rows, resolution_rate = resolution_rate))
}

if (sys.nframe() == 0) {
  run_multiplet_overlap_benchmark()
}
