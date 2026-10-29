# peak_picking_functions.R
#
# Reduces a raw TopSpin peak-picking export (hundreds of points, many of
# them noise/baseline wiggle or the same multiplet line picked several
# times) down to a smaller set of real, distinct peaks, then aligns that
# reduced list across several samples into one feature table ready for
# the MetaboAnalyst stats script.
#
# Two functions:
#   pick_strong_peaks()         -- run once per sample's TopSpin export
#   merge_peaks_across_samples() -- combine the per-sample results

# ------------------------------------------------------------------
# pick_strong_peaks()
#
# How it filters:
#  1. Local baseline/noise floor: a TopSpin peak list picked with a
#     permissive threshold reports negative-going points too -- points
#     where the spectrum dipped below zero from imperfect phasing or
#     baseline correction. Those negative values are themselves a
#     direct, local measurement of how big the noise is in that part
#     of the spectrum. For each candidate peak, this looks at the
#     negative points within `noise_window` rows on either side and
#     takes `noise_mult` times their median absolute size as that
#     point's local noise floor. A positive point has to clear that
#     floor to count as a real peak, not noise.
#  2. Cluster-merge: many of the remaining points are the same peak's
#     top picked two or three times a fraction of a ppm apart. Points
#     within `cluster_tol` ppm of each other are treated as one peak,
#     keeping only the tallest. `cluster_tol` defaults to 0.004 ppm,
#     tighter than a real resolvable J-coupling (usually >=0.007 ppm
#     at 600 MHz), so genuinely separate multiplet lines are kept
#     apart rather than merged into one.
#
# Defaults (noise_mult = 1.5, cluster_tol = 0.004) were tuned against
# a real 839-point TopSpin export and landed at 275 peaks, inside the
# 200-300 target range. Raise noise_mult for fewer, more
# conservative peaks; lower it for more. Raise cluster_tol to merge
# more aggressively.
pick_strong_peaks <- function(ppm, intensity, noise_window = 15,
                               noise_mult = 1.5, cluster_tol = 0.004) {
  n <- length(ppm)
  stopifnot(length(intensity) == n, all(diff(ppm) <= 0))  # expects ppm descending, as TopSpin exports it

  noise_floor <- numeric(n)
  for (i in seq_len(n)) {
    lo <- max(1, i - noise_window)
    hi <- min(n, i + noise_window)
    neg_here <- intensity[lo:hi]
    neg_here <- neg_here[neg_here < 0]
    noise_floor[i] <- if (length(neg_here) >= 2) {
      noise_mult * stats::median(abs(neg_here))
    } else {
      # no negative points nearby (happens in the cleaner aliphatic
      # region) -- fall back to a low quantile of the local |intensity|
      noise_mult * stats::quantile(abs(intensity[lo:hi]), 0.25)
    }
  }

  keep <- intensity > 0 & intensity > noise_floor
  cand <- data.frame(ppm = ppm[keep], intensity = intensity[keep])
  cand <- cand[order(-cand$intensity), ]  # tallest first for greedy clustering

  final_rows <- data.frame(ppm = numeric(0), intensity = numeric(0))
  claimed <- rep(FALSE, nrow(cand))
  for (i in seq_len(nrow(cand))) {
    if (claimed[i]) next
    final_rows <- rbind(final_rows, cand[i, ])
    nearby <- which(abs(cand$ppm - cand$ppm[i]) <= cluster_tol)
    claimed[nearby] <- TRUE
  }
  final_rows <- final_rows[order(-final_rows$ppm), ]
  rownames(final_rows) <- NULL
  final_rows
}

# ------------------------------------------------------------------
# merge_peaks_across_samples()
#
# Takes a NAMED list of per-sample tables, each already run through
# pick_strong_peaks() (each a data.frame with ppm, intensity columns),
# and aligns them into one feature table: one row per distinct peak
# position, one column per sample. Peaks within `tol` ppm of each
# other across different samples are treated as the same feature,
# which allows for the small sample-to-sample shift that's normal
# even after good referencing. A sample with no matching peak in a
# given cluster gets NA for that feature -- a real peak legitimately
# absent in that sample, not a bug.
merge_peaks_across_samples <- function(table_list, tol = 0.01) {
  sample_names <- names(table_list)
  if (is.null(sample_names) || any(sample_names == "")) {
    stop("table_list must be a named list, e.g. list(Sample1 = tbl1, Sample2 = tbl2)")
  }

  pooled <- do.call(rbind, lapply(sample_names, function(s) {
    tbl <- table_list[[s]]
    data.frame(ppm = tbl$ppm, intensity = tbl$intensity, sample = s,
               stringsAsFactors = FALSE)
  }))
  pooled <- pooled[order(-pooled$intensity), ]

  claimed <- rep(FALSE, nrow(pooled))
  clusters <- list()
  for (i in seq_len(nrow(pooled))) {
    if (claimed[i]) next
    center_ppm <- pooled$ppm[i]
    nearby <- which(!claimed & abs(pooled$ppm - center_ppm) <= tol)
    clusters[[length(clusters) + 1]] <- pooled[nearby, ]
    claimed[nearby] <- TRUE
  }

  out <- data.frame(ppm = vapply(clusters, function(cl) cl$ppm[which.max(cl$intensity)], numeric(1)))
  for (s in sample_names) {
    out[[s]] <- vapply(clusters, function(cl) {
      hit <- cl[cl$sample == s, ]
      if (nrow(hit) == 0) NA_real_ else hit$intensity[which.max(hit$intensity)]
    }, numeric(1))
  }
  out <- out[order(-out$ppm), ]
  rownames(out) <- NULL
  out
}

# ------------------------------------------------------------------
# merge_peaks_reference_guided()  (RECOMMENDED over merge_peaks_across_samples())
#
# merge_peaks_across_samples() above pools every sample's peaks
# together and lets whichever point happens to be tallest define each
# cluster's center -- so the anchor ppm for a feature can drift
# depending on which sample happened to have the biggest peak there,
# and a sample with no matching peak is always recorded as NA even if
# it actually has real, if weaker, signal at that position.
#
# This version fixes both: it fixes the feature positions to ONE
# reference sample's own peak list (by default, whichever sample has
# the most picked peaks -- usually the cleanest/most complete one), so
# every feature's ppm is anchored consistently; and for any other
# sample with no picked peak at a given position, it checks that
# SAMPLE'S OWN raw (baseline-corrected) spectrum near that exact ppm
# before giving up -- a real signal that just fell under that sample's
# own peak-picking threshold is recorded, and only a position with
# genuinely no signal (flat, noise-only) is left as NA. That matters
# here specifically: a missing value should mean "this metabolite
# really isn't present/detectable in this sample," not "the detector
# happened to miss it by a hair."
#
# table_list:    same as merge_peaks_across_samples() -- a named list of
#                per-sample peak tables (ppm, intensity columns).
# spectra_list:  OPTIONAL, but needed for the fallback signal check.
#                A named list (same names as table_list) of
#                list(ppm = <full ppm grid>, corrected = <baseline-
#                corrected intensity>, noise = <local noise, optional>)
#                for each sample -- i.e. what you already have after
#                calling fit_airpls_baseline() and subtracting it.
#                Without this, a sample with no matching peak just
#                gets NA, same as the older function.
# reference:     which sample's peak list to use as the fixed set of
#                feature positions. Defaults to whichever sample has
#                the most peaks.
# tol:           how close (in ppm) another sample's own picked peak
#                has to be to a reference position to count as the
#                same feature.
# fallback_window_ppm, fallback_noise_mult: for the raw-signal check --
#                look within fallback_window_ppm of the reference
#                position, and only accept it as real signal (not NA)
#                if it clears fallback_noise_mult times that sample's
#                own local noise level there.
#
# NOTE: the fallback-signal lookup below calls .window_bounds(), which is
# DEFINED in baseline_peak_picking.R, not this file. Every script in this
# project sources both files, in the order peak_picking_functions.R then
# baseline_peak_picking.R (see baseline_peak_picking.R's own header) --
# that's fine even though this function is DEFINED first, because R only
# looks up .window_bounds() when merge_peaks_reference_guided() is
# actually CALLED, which always happens later, after both files are
# loaded. Don't call this function without also sourcing
# baseline_peak_picking.R first.
merge_peaks_reference_guided <- function(table_list, spectra_list = NULL,
                                          reference = NULL, tol = 0.01,
                                          fallback_window_ppm = 0.01,
                                          fallback_noise_mult = 2) {
  sample_names <- names(table_list)
  if (is.null(sample_names) || any(sample_names == "")) {
    stop("table_list must be a named list, e.g. list(Sample1 = tbl1, Sample2 = tbl2)")
  }
  if (is.null(reference)) {
    reference <- sample_names[which.max(sapply(table_list, nrow))]
  }
  ref_tbl <- table_list[[reference]][order(-table_list[[reference]]$intensity), ]

  out <- data.frame(ppm = ref_tbl$ppm)
  for (s in sample_names) {
    if (s == reference) { out[[s]] <- ref_tbl$intensity; next }
    tbl <- table_list[[s]]
    # Fallback spectrum, prepared ONCE per sample rather than re-scanned on
    # every feature. Two independent reasons to sort it here:
    #  1. Speed: the old code did abs(sp$ppm - center) <= fallback_window_ppm,
    #     a full scan of the WHOLE raw spectrum, for EVERY merged feature
    #     (nrow(out), 200+ in a real run) x EVERY non-reference sample --
    #     O(n_features x n_samples x spectrum_length). .window_bounds()
    #     (defined in baseline_peak_picking.R, sourced alongside this file
    #     in every script that calls this function) turns each lookup into
    #     an O(log n) binary search via findInterval() instead, but that
    #     only works on an ASCENDING-sorted ppm axis.
    #  2. Correctness: getPpmGrid()/the raw import never guaranteed an
    #     ascending ppm axis -- the rest of this codebase defensively
    #     re-sorts for exactly that reason (see baseline_peak_picking.R) --
    #     and stats::approx() below silently requires non-decreasing x
    #     anyway. So this sort was already implicitly assumed; doing it
    #     explicitly closes a latent correctness gap, not just a speed one.
    sp <- NULL
    if (!is.null(spectra_list) && !is.null(spectra_list[[s]])) {
      sp_raw <- spectra_list[[s]]
      sp_ord <- order(sp_raw$ppm)
      sp <- list(ppm = sp_raw$ppm[sp_ord], corrected = sp_raw$corrected[sp_ord])
      if (!is.null(sp_raw$noise)) sp$noise <- sp_raw$noise[sp_ord]
    }
    vals <- numeric(nrow(out))
    for (i in seq_len(nrow(out))) {
      center <- out$ppm[i]
      hit <- tbl[abs(tbl$ppm - center) <= tol, ]
      if (nrow(hit) > 0) {
        vals[i] <- hit$intensity[which.max(hit$intensity)]
      } else if (!is.null(sp)) {
        b <- .window_bounds(sp$ppm, center, fallback_window_ppm)
        if (b$hi >= b$lo) {
          local_val <- max(sp$corrected[b$lo:b$hi])
          noise_here <- if (!is.null(sp$noise)) {
            stats::approx(sp$ppm, sp$noise, xout = center, rule = 2)$y
          } else {
            NA_real_
          }
          if (!is.na(noise_here) && local_val > fallback_noise_mult * noise_here) {
            vals[i] <- local_val  # real evidence of signal, below the full peak-picking bar, recorded rather than discarded
          } else {
            vals[i] <- NA_real_   # genuinely no signal there -- a true missing feature, not a detector miss
          }
        } else {
          vals[i] <- NA_real_
        }
      } else {
        vals[i] <- NA_real_
      }
    }
    out[[s]] <- vals
  }
  out <- out[order(-out$ppm), ]
  rownames(out) <- NULL
  attr(out, "reference") <- reference
  out
}

# ------------------------------------------------------------------
# write_metaboanalyst_csv()
#
# Writes a merged feature table (ppm + one column per sample) out in
# MetaboAnalystR's "colu" + "disc" format: samples across the top,
# group label in the second row, one feature per row below. `groups`
# must be a named vector mapping each sample column name to its group
# label, e.g. c(SampleA = "control", SampleB = "sample").
write_metaboanalyst_csv <- function(merged, groups, out_file) {
  sample_cols <- setdiff(colnames(merged), "ppm")
  missing <- setdiff(sample_cols, names(groups))
  if (length(missing) > 0) {
    stop("No group given for: ", paste(missing, collapse = ", "))
  }

  feature_label <- paste0("ppm_", sprintf("%.4f", merged$ppm))
  body <- cbind(Sample = feature_label, merged[, sample_cols, drop = FALSE])

  header_row <- c("Sample", sample_cols)
  label_row  <- c("Label", unname(groups[sample_cols]))

  con <- file(out_file, "w")
  writeLines(paste(header_row, collapse = ","), con)
  writeLines(paste(label_row, collapse = ","), con)
  write.table(body, con, sep = ",", row.names = FALSE, col.names = FALSE,
              na = "", append = TRUE, quote = FALSE)
  close(con)
  message("Wrote ", nrow(merged), " features x ", length(sample_cols),
          " samples to ", normalizePath(out_file))
}
