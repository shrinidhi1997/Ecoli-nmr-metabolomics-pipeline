# baseline_peak_picking.R
#
# Picks peaks directly off the full, continuous spectra you already have
# in session from the import step (ppm, spec_mat), no ASICS quantification
# involved.
#
# Functions, in the order you'd normally use them:
#   fit_airpls_baseline()         -- adaptive baseline (replaces the old
#                                     fixed-percentile window baseline)
#   savgol_filter()                -- Savitzky-Golay smoothing, used only
#                                     to find peak/valley shape
#   find_shape_based_peaks()       -- RECOMMENDED peak picker. Valley-to-
#                                     valley prominence, now built on the
#                                     two functions above plus a
#                                     bias-corrected noise estimate and a
#                                     minimum-width check.
#   refine_peaks_lorentzian()      -- OPTIONAL. Fits a real NMR lineshape
#                                     (Lorentzian) to each detected peak
#                                     for a refined center/height/width/
#                                     area -- use when you need actual
#                                     quantification, not just a peak list.
#   find_baseline_corrected_peaks() -- the original point-by-point
#                                     detector, kept only for comparison.
#
# Source peak_picking_functions.R yourself before this file -- this
# script reuses merge_peaks_across_samples()/merge_peaks_reference_guided()
# and write_metaboanalyst_csv() from it:
#
#   source("peak_picking_functions.R")
#   source("baseline_peak_picking.R")

# ==================================================================
# .window_bounds() -- fast index-range lookup, used everywhere below that
# needs "every point within +/-half_win of each of many control points."
#
# The straightforward way to write that (and how every function below
# used to) is: for each control point, scan the WHOLE ppm array with a
# boolean comparison. That is genuinely O(control points x spectrum
# length) -- fine on a short test window, but a real full spectrum run
# at a real sampling density has both many more control points AND many
# more spectrum points than any synthetic benchmark in this project, and
# this exact pattern got called repeatedly (baseline fit, both
# detectors' own noise estimation, more than once in some of them) on
# every sample. Confirmed as the real, not hypothetical, slow part when
# a run on real 6-sample data stalled out specifically during baseline
# fitting.
#
# The fix: ppm is already sorted everywhere this is used (every caller
# does `ord <- order(ppm)` first), so "which index range falls within
# this window" is a binary search (findInterval()), not a full scan --
# O(log n) per control point instead of O(n), and R's findInterval()
# does all of them in one vectorized pass. Returns 1-based (lo, hi)
# INDEX pairs, one per center, instead of a boolean mask -- slicing
# y[lo[i]:hi[i]] then only touches the points actually inside that one
# window, not the whole array. An empty window (no points in range)
# comes back as hi[i] < lo[i]; every caller already checks for that.
# Verified to return EXACTLY the same point sets as the old brute-force
# scan on real test data before being wired into anything below.
.window_bounds <- function(x_sorted, centers, half_win) {
  lo <- findInterval(centers - half_win, x_sorted) + 1L
  hi <- findInterval(centers + half_win, x_sorted)
  list(lo = lo, hi = hi)
}

# ==================================================================
# PART 1 -- baseline: fit_airpls_baseline()
#
# Replaces the earlier "5th percentile within a sliding window" baseline
# with airPLS (adaptive iteratively reweighted penalized least squares --
# Zhang et al. 2010), the standard baseline-correction algorithm used in
# chromatography/spectroscopy. Rather than picking a fixed low quantile
# everywhere (which has to be tuned and can be wrong in regions with
# different peak density), airPLS FITS a smooth curve and then
# iteratively re-weights: points where the data sits below the current
# fit get a large weight (pull the curve down to them -- that's where
# the real floor is), points at or above the fit get weight 0 (peaks
# can't pull the baseline up). It converges to a curve that hugs the
# lower envelope of the data.
#
# Implementation note: the standard airPLS solves a large banded linear
# system at every point of the spectrum, usually via a sparse-matrix
# package. That package isn't guaranteed to be installed here, so this
# runs the same iteratively-reweighted penalized fit on the existing
# control-point grid instead (same grid already used throughout this
# project for speed) -- at a few hundred control points a plain dense
# solve() is instant, so no extra package is needed at all.
#
# control_point_spacing_ppm: resolution of the control-point grid the
#             baseline is actually solved on (then interpolated back to
#             the full spectrum).
# denoise_window_ppm: each control point's starting value is the local
#             median of the raw spectrum in a small window -- just
#             enough to suppress single-point noise before airPLS does
#             the real work of separating floor from peaks.
# lambda:     smoothness penalty. Higher = stiffer/smoother baseline that
#             can't bend to track every ripple; lower = more flexible,
#             risking tracking into real peaks. 1e5 tracked the drifting
#             baseline in testing to within ~0.002 of the true floor
#             everywhere except directly under a very broad peak (where
#             it correctly refuses to chase the peak's hump at all --
#             exactly the desired behavior, since a broad peak there,
#             e.g. the residual water peak, should be excluded from
#             picking anyway, not used to anchor the baseline).
# max_iter:   CHANGED from 15 to 4 after a real, measured finding. In a
#             crowded region (many overlapping real lines close together,
#             the kind a dense multiplet-rich window actually looks like),
#             there are very few points that are genuinely peak-free --
#             most "floor" samples still carry some tail from a nearby
#             line. The reweighting scheme in the loop below only ever
#             pulls the curve DOWN toward whatever's currently below it,
#             iteration after iteration, and in a crowded region the only
#             points consistently available to pull toward are the rare,
#             narrow true gaps -- which are themselves a biased sample of
#             the lowest local values. Run this to full convergence (the
#             old default, 15, converges by ~iteration 5 and then does
#             nothing further) and the fit keeps chasing those sparse low
#             points past the real floor, landing measurably BELOW it --
#             confirmed directly: on a synthetic crowded-multiplet
#             benchmark built to match a real screenshot, max_iter=15
#             left the baseline sitting 0.015-0.021 below the true floor
#             on average, consistently, across 3 random seeds and across
#             a 300x range of lambda (so lambda was never the lever here).
#             Stopping earlier, before the curve has fully committed to
#             those sparse low anchors, measurably closes that gap
#             (max_iter=4 roughly halves the bias on the same benchmark,
#             every seed tested) WITHOUT hurting a normal, less-crowded
#             region -- re-checked on this file's own gentle synthetic
#             self-test below (excluding the deliberately-untracked
#             broad-water zone for a fair comparison): max_iter=4 scores
#             about the same RMSE there as the old max_iter=15 did, no
#             regression. If your own region is even MORE densely packed
#             than either benchmark, this may still run a little low --
#             try max_iter=3 and compare against your own TopSpin view
#             rather than assume 4 is universally exact.
detect_solvent_exclusion_range <- function(ppm, intensity, seed_center_ppm = 4.7,
                                            seed_half_width_ppm = 0.7, floor_mult = 3) {
  # Finds how far the water (or other dominant solvent) peak's tails
  # actually extend on THIS sample, instead of assuming a fixed ppm
  # window. Why this exists: a fixed window like c(3.9, 5.3) is a guess
  # at how wide the water peak's tails are, and on real data it can be
  # badly wrong in either direction -- too narrow leaves a strong
  # residual tail just outside the window, which (once that ppm range is
  # excluded from the baseline FIT specifically, not just peak
  # detection) no longer gets smoothly bridged, since the fit now has to
  # connect across the gap from a point that's still partway down the
  # tail, not yet back at the real background level, producing a sharp
  # spike in the corrected spectrum right at the edge; too wide throws
  # away real peak territory the baseline fit didn't need to lose.
  # Confirmed on real HAIRD/SWEPF data: c(3.9, 5.3) is too narrow --
  # intensity at both 3.9 ppm and 5.3 ppm is still several times the
  # spectrum's typical background level, nowhere near flat yet.
  #
  # Method: start from the tallest point near seed_center_ppm (the
  # water peak itself), then walk outward on each side until intensity
  # drops below floor_mult times a robust "typical background" estimate.
  # That estimate is the 30th percentile taken ONLY over the
  # conventional 0-10 ppm 1H window, not the full digitized sweep --
  # this spectrometer digitizes well beyond that (roughly -3.75 to
  # 16.3 ppm here), and most of that extra range is near-pure
  # electronic noise close to zero, which otherwise drags a whole-
  # spectrum percentile down to a tiny fraction of the real in-spectrum
  # background level (confirmed: 87 a.u. whole-spectrum vs. ~1000-5000
  # a.u. actually typical between real peaks in 0-10 ppm) and made this
  # function's walk run away across almost the entire spectrum instead
  # of stopping at the water peak's real extent.
  typical_window <- ppm >= 0 & ppm <= 10
  floor_source <- if (any(typical_window)) intensity[typical_window] else intensity
  floor_val <- stats::quantile(floor_source, 0.30, names = FALSE)
  thresh <- floor_mult * floor_val

  seed_mask <- ppm >= (seed_center_ppm - seed_half_width_ppm) & ppm <= (seed_center_ppm + seed_half_width_ppm)
  if (!any(seed_mask)) return(c(seed_center_ppm - seed_half_width_ppm, seed_center_ppm + seed_half_width_ppm))
  center_idx <- which(seed_mask)[which.max(intensity[seed_mask])]

  ord <- order(ppm)
  ppm_s <- ppm[ord]; y_s <- intensity[ord]
  center_pos <- which(ppm_s == ppm[center_idx])[1]
  grid_spacing <- stats::median(diff(ppm_s))

  # A single point below thresh isn't enough to call the tail over --
  # a real, isolated analyte peak taller than thresh (common; that's
  # what makes it a peak) would otherwise fool a point-by-point walk
  # into reading right through the first brief dip after it and
  # continuing outward, well past where the water tail actually ended
  # (this is exactly what happened with a plain point-by-point walk:
  # the detected range came out wider at higher concentration, where
  # real aliphatic peaks are taller -- a tell that real peaks, not the
  # water tail, were extending it). Requiring a SUSTAINED run below
  # thresh (run_ppm wide, entirely below threshold) before declaring
  # "back to baseline" is robust to that: an isolated real peak doesn't
  # have run_ppm of genuinely quiet spectrum sitting right on top of it.
  run_ppm <- 0.15
  run_pts <- max(3, round(run_ppm / grid_spacing))

  find_edge <- function(start_pos, step) {
    pos <- start_pos
    n <- length(y_s)
    repeat {
      pos <- pos + step
      if (pos < 1 || pos > n) return(max(1, min(pos, n)))
      window_idx <- if (step > 0) pos:min(pos + run_pts - 1, n) else max(pos - run_pts + 1, 1):pos
      if (max(y_s[window_idx]) <= thresh) return(pos)
    }
  }
  hi_pos <- find_edge(center_pos, 1L)
  lo_pos <- find_edge(center_pos, -1L)

  c(ppm_s[lo_pos], ppm_s[hi_pos])
}

fit_airpls_baseline <- function(ppm, intensity,
                                 control_point_spacing_ppm = 0.02,
                                 denoise_window_ppm = 0.01,
                                 lambda = 1e5,
                                 max_iter = 4,
                                 ratio = 0.001,
                                 exclude = NULL) {
  # exclude: optional list of c(lo, hi) ppm ranges (e.g. the water/DSS
  # ranges find_shape_based_peaks() already excludes by default) to drop
  # BEFORE the IRLS fit, not just from peak detection afterward.
  #
  # Why this matters: the asymmetric reweighting below (w_new[d>=0] <- 0)
  # assumes peaks are reasonably narrow relative to the baseline's own
  # curvature, so a genuine peak quickly gets zero weight and stops
  # pulling the fit upward. A solvent peak (water, ~4.7-4.9 ppm) breaks
  # that assumption: it can be 100-1000x taller than real signal and
  # wide at its base (its shoulders bleed out to roughly 3.9-5.3 ppm),
  # and with max_iter capped at 4 (deliberately low -- see this
  # function's other callers' notes on full convergence biasing low in
  # crowded regions), there isn't enough reweighting headroom to fully
  # flatten the fit's response to something that extreme. The result,
  # confirmed on real data (HAIRD/SWEPF standards): a broad hump in the
  # fitted baseline centered on the water peak, bleeding into the
  # neighboring signal region, with a small negative dip in the
  # corrected spectrum right at the edge of whatever ppm range peak
  # detection excludes.
  #
  # The direct fix is upstream of the reweighting entirely: never give
  # the solver those control points in the first place. Dropping them
  # here (same NA-then-interpolate path already used below for windows
  # with no data at all) lets stats::approx()'s existing bridging step
  # fill the gap smoothly from the real baseline on either side, instead
  # of the fit ever seeing -- and partially tracking -- the solvent peak.
  # NULL (the default) reproduces the exact old behavior for any
  # existing caller that doesn't pass this.
  cp_ppm <- seq(min(ppm), max(ppm), by = control_point_spacing_ppm)
  half_win <- denoise_window_ppm / 2
  bounds <- .window_bounds(ppm, cp_ppm, half_win)
  excl_mask <- rep(FALSE, length(cp_ppm))
  for (rng in exclude) excl_mask <- excl_mask | (cp_ppm >= rng[1] & cp_ppm <= rng[2])
  cp_y <- numeric(length(cp_ppm))
  for (i in seq_along(cp_ppm)) {
    cp_y[i] <- if (!excl_mask[i] && bounds$hi[i] >= bounds$lo[i]) {
      stats::median(intensity[bounds$lo[i]:bounds$hi[i]])
    } else {
      NA
    }
  }
  ok <- !is.na(cp_y)
  cp_ppm <- cp_ppm[ok]; y <- cp_y[ok]
  m <- length(y)

  # second-difference operator -- the "how smooth" penalty matrix
  D2 <- matrix(0, nrow = m - 2, ncol = m)
  for (i in seq_len(m - 2)) { D2[i, i] <- 1; D2[i, i + 1] <- -2; D2[i, i + 2] <- 1 }
  penalty <- lambda * crossprod(D2)

  w <- rep(1, m)
  z <- y
  for (iter in seq_len(max_iter)) {
    A <- diag(w) + penalty
    z <- as.numeric(solve(A, w * y))
    d <- y - z
    neg <- d[d < 0]
    dssn <- sum(abs(neg))
    if (dssn < ratio * sum(abs(y)) || length(neg) == 0) break
    w_new <- numeric(m)
    w_new[d >= 0] <- 0                              # at/above the fit -- a peak, don't let it pull the baseline up
    w_new[d < 0] <- exp(iter * d[d < 0] / dssn)      # below the fit -- real floor, pull harder each iteration
    w_new[1] <- exp(iter * max(abs(neg)) / dssn)     # keep the two ends anchored the same way
    w_new[m] <- w_new[1]
    w <- w_new
  }
  stats::approx(cp_ppm, z, xout = ppm, rule = 2)$y
}

# ==================================================================
# PART 1b -- fit_binned_baseline()
#
# A different baseline estimator from fit_airpls_baseline() above, built
# specifically to fix a real, confirmed problem with it on this project's
# real HAIRD/SWEPF data: the water/solvent peak's broad hump in the fitted
# baseline doesn't stay local. fit_airpls_baseline() solves ONE penalized
# least-squares system over every control point at once (the `penalty`
# matrix couples each point to its neighbours via a second-difference
# smoothness term), so a huge, broad excursion anywhere in that system
# (the water peak) pulls on the WHOLE fitted curve's shape to keep it
# smooth overall, not just its own immediate neighbourhood. Confirmed
# directly on real HAIRD-60-1H: with the plain airPLS fit, the baseline
# sat measurably ABOVE the true local floor even 4-5 ppm away from the
# water peak (far-field overshoot up to ~20,000 a.u., nowhere near the
# peak), and right on the water peak's own falling shoulder (5.2-5.6 ppm)
# it overshot the real signal there by 60-89%, which is exactly the
# "baseline hallucinating / going beyond where it should" behavior
# reported after the 0-11 ppm reporting-range fix (that fix cleaned up
# digitization padding outside the real spectrum, it never touched this).
#
# This fixes it by doing the opposite of a single global smoothness-
# penalized solve: split the spectrum into BIG ppm bins (bin_width_ppm,
# much wider than any real single peak or multiplet cluster), take a
# robust LOW PERCENTILE of the raw intensity within each bin as that
# bin's background estimate, and interpolate LINEARLY between bin centers
# -- no global penalty matrix, no cross-bin smoothing pass at all. Each
# bin's own value only ever depends on what's actually inside it; nothing
# 4 ppm away can pull it up or down, by construction. This is the direct
# implementation of a two-part design goal:
#   1. Ignore small peaks by using big bins --
#      a real analyte peak (a single line or multiplet cluster, usually
#      well under ~0.1 ppm wide) occupies a small fraction of a
#      bin_width_ppm=0.2 bin, so a LOW percentile (pct=0.10 below) of
#      that bin's ~600-1200 points lands close to the genuine quiet
#      floor around it, the same way a median would, but more robust to
#      a bin that's partly covered by peak -- the narrow peak just can't
#      drag a low-enough percentile up much. Peak PICKING is untouched:
#      find_shape_based_peaks() still runs on the exact same fine-grained
#      raw (ppm, intensity) it always did, this only changes what
#      baseline gets subtracted from it first.
#   2. Let the baseline follow the broad water peak, then apply it to
#      the sample -- there's no separate water-specific model here
#      (every attempt at one, documented above and in fit_airpls_baseline's
#      own notes, either fit the wrong lineshape or needed a hand-tuned
#      window that was wrong on some sample/concentration). Instead, the
#      SAME big-bin-plus-percentile rule is applied everywhere, and it
#      naturally does both jobs on its own: where the water peak genuinely
#      dominates a bin (its base is broad, several ppm wide -- unlike a
#      real small-molecule line), even the 10th percentile of that bin
#      reflects the peak's own real, broad shape, so the baseline still
#      visibly follows it there, as intended.
#      Where it doesn't dominate (everywhere else), the percentile reads
#      the genuine quiet floor instead.
#
# IMPORTANT, found by direct testing: do NOT add a cross-bin smoothing
# pass (a moving average across neighbouring bins' percentile values) on
# top of this, even though the raw bin-to-bin curve looks a little
# stair-steppy without one. Tried first with a 3-bin moving average: it
# reintroduces almost exactly the same failure this function exists to
# avoid, just over a much shorter range -- one bin sitting right on the
# water peak's tall apex gets blended into its lower neighbouring bins,
# so the falling shoulder (5.2-5.6 ppm on real HAIRD-120-1H) still
# overshot the real signal there by 42-174% with smoothing on, vs. 0%
# (no overshoot at all, confirmed numerically and visually against the
# raw trace) with it off. Piecewise-linear interpolation between raw,
# unsmoothed bin percentiles is the one that tested clean; use it as-is.
#
# Confirmed effect across all 6 real samples (HAIRD-30/60/120,
# SWEPF-30/60/120), bin_width_ppm=0.2/pct=0.10 vs. the plain airPLS fit:
# far-field overshoot (0-11 ppm, outside 3.3-6.3) dropped from
# 2,367-20,123 a.u. down to 657-2,835; near-water relative overshoot
# (3.3-6.3 ppm) dropped from +60%..+89% down to -7%..+6% (i.e. essentially
# gone, baseline no longer runs above the real signal). The direct
# consequence, also confirmed and worth knowing before trusting this:
# peak counts (0-11 ppm) roughly doubled on some samples (e.g. HAIRD-60-1H
# 60->200). Checked this wasn't a new artifact two ways -- (a) the extra
# peaks' ppm values are uniformly distributed relative to the 0.2 ppm bin
# edges (mean distance-to-nearest-edge ~0.052 ppm, close to the 0.05
# expected by pure chance), not clustered at edges the way an
# interpolation-kink artifact would be; (b) visually, in the two zones
# that account for most of the increase (3.4-4.1 and 6.4-8.8 ppm), the
# "new" peaks are tall, sharp, cleanly-isolated real multiplets sitting on
# a near-zero corrected floor -- exactly the shape of real signal that the
# old overshooting baseline was suppressing/flattening before, not noise
# crossing a lowered threshold. The far-field 8.5-11 ppm zone (pure
# digitization noise, same finding as the earlier 0-11 ppm reporting-range
# fix) does pick up a handful of extra marginal detections there too; this
# is the known limit of any baseline fit that has no real floor to compare
# against in a region with no real signal at all, and -- same as before --
# nothing in that range is expected to match a real compound downstream.
#
# bin_width_ppm: width of each big bin. 0.2 ppm tested well -- wide
#             enough that a typical single line or multiplet cluster
#             (usually <0.1 ppm) is a minority of the bin's points, narrow
#             enough to still track real, slower background drift.
# pct:        percentile of each bin's raw intensity used as its
#             background estimate. 0.10 tested well; a real bin has
#             hundreds to low thousands of points at typical real spectral
#             density, so this is not a small-sample, noisy estimate.
#             Lower (e.g. 0.05) pulls slightly closer to the true floor
#             but starts to cost a little on far-field accuracy; higher
#             (e.g. 0.25) is more sensitive to a bin where a peak takes up
#             a larger share of it. See this function's own header note
#             above for the numbers each choice produced on real data.
fit_binned_baseline <- function(ppm, intensity, bin_width_ppm = 0.2, pct = 0.10) {
  ord <- order(ppm)
  ppm <- ppm[ord]; intensity <- intensity[ord]
  bin_edges <- seq(min(ppm), max(ppm) + bin_width_ppm, by = bin_width_ppm)
  bin_idx <- findInterval(ppm, bin_edges)
  bin_centers <- as.numeric(tapply(ppm, bin_idx, mean))
  bin_vals <- as.numeric(tapply(intensity, bin_idx, function(v) stats::quantile(v, pct, names = FALSE)))
  o <- order(bin_centers)
  bin_centers <- bin_centers[o]; bin_vals <- bin_vals[o]
  stats::approx(bin_centers, bin_vals, xout = ppm, rule = 2)$y
}

# ==================================================================
# PART 2 -- smoothing for shape detection: savgol_filter()
#
# Replaces the plain moving average with a Savitzky-Golay filter
# (local polynomial fit, evaluated at the center point) -- the standard
# choice for smoothing spectra without flattening real peak shapes the
# way a plain moving average can. Coefficients are computed from
# scratch via a small least-squares fit, so no extra package is needed.
savgol_coeffs <- function(half_win, poly_order) {
  win <- -half_win:half_win
  Van <- outer(win, 0:poly_order, `^`)
  # QR-based pseudo-inverse, not solve(crossprod(Van), t(Van)) -- crossprod
  # SQUARES Van's condition number, and raw powers of a wide index window
  # (win^poly_order) are already badly scaled on their own. Confirmed this
  # matters for real, not just in theory: for a 273-point window (a 0.010
  # ppm smoothing width on a 30000-point, 1.1 ppm-wide run) at poly_order
  # 4, crossprod(Van)'s condition number hit ~4.7e16, right at double-
  # precision's own limit, and solve() failed outright ("computationally
  # singular"). QR avoids squaring that condition number and stayed stable
  # on the exact same case; verified it reproduces the old crossprod-based
  # coefficients to within 1e-16 wherever the old approach still worked,
  # so this is a numerically safer path to the same answer, not a
  # different method.
  qr_van <- qr(Van)
  qr.coef(qr_van, diag(nrow(Van)))[1, ]
}

savgol_filter <- function(y, half_win, poly_order = 3) {
  n <- length(y)
  win_pts <- 2 * half_win + 1
  if (win_pts > n) return(y)
  coeffs <- savgol_coeffs(half_win, poly_order)
  padded <- c(rep(y[1], half_win), y, rep(y[n], half_win))
  out <- numeric(n)
  for (i in seq_len(n)) out[i] <- sum(coeffs * padded[i:(i + win_pts - 1)])
  out
}

# ==================================================================
# Noise estimate shared by the functions below: median absolute value
# of NEGATIVE residuals (a real peak against a correctly fit baseline
# can never produce a negative residual, so these points are always
# pure noise, however much of the window a peak occupies), divided by
# qnorm(0.75) = 0.6745 to correct a real, measured bias: for a
# symmetric, zero-mean noise distribution, |residual| given a negative
# residual has exactly the same distribution as |residual| overall, so
# it should be converted to sigma with the same constant as a two-sided
# MAD. Tested against a two-sided iterative sigma-clipping scheme and a
# low-quantile-of-|residual| scheme first -- both failed: a real peak's
# shape tapers continuously to near zero rather than switching on/off,
# so there's no clean "signal vs. noise" split by magnitude alone, and
# the iterative version was outright unstable (its own threshold grew
# every round instead of shrinking). This estimator is the one that
# held up: stable (single pass, no iteration to diverge) and tested
# close to the true noise level in every window-occupancy scenario
# except the extreme case of one very broad peak filling nearly the
# entire window with no clean edges at all -- see the note on
# find_shape_based_peaks() below for why that extreme case doesn't
# apply to a real excluded-water-region spectrum.
.noise_from_negative_residuals <- function(resid) {
  neg <- resid[resid < 0]
  if (length(neg) < 2) return(NA_real_)
  stats::median(abs(neg)) / 0.6744898
}

# ==================================================================
# PART 3 -- find_shape_based_peaks()  (RECOMMENDED)
#
# Picks peaks by their actual shape -- each peak's own local minima
# (valleys) on either side -- instead of testing each point in
# isolation against a noise floor, now built on the airPLS baseline and
# Savitzky-Golay smoothing above.
#
# ppm, intensity: the full continuous spectrum for ONE sample (from
#                 getPpmGrid(spectra_obj) and one column of
#                 getSpectra(spectra_obj))
# exclude:        ppm ranges to skip entirely (default: the 3.9-5.3
#                 residual water region, and a narrow band around
#                 0.000 ppm where the DSS/TMS calibration standard
#                 sits -- not a metabolite, so it shouldn't be picked
#                 up as a feature).
# airpls_lambda, control_point_spacing_ppm: passed through to
#                 fit_airpls_baseline() -- see that function's notes.
# sg_window_ppm, sg_poly: the Savitzky-Golay smoothing window (in ppm)
#                 and polynomial order, used ONLY to find genuine
#                 peak/valley topology. Real, densely-sampled spectra
#                 have noise on every single point, so every tiny
#                 up-down wiggle between adjacent points creates its
#                 own miniature "local max next to a local min" --
#                 valley-to-valley prominence measured directly on the
#                 raw trace just measures single-point noise amplitude,
#                 not real peak shape, and floods the output with
#                 thousands of spurious noise-level "peaks" (this is
#                 exactly what happened on the first attempt at this
#                 function, before smoothing was added). Peak height
#                 and position are still read back off the UNSMOOTHED
#                 trace (see "refine" in the code), so a real peak's
#                 reported intensity isn't flattened by the smoothing.
#                 Default 0.006 ppm is small relative to a typical
#                 resolvable linewidth (usually >=0.007-0.01 ppm), so
#                 distinct multiplet lines stay distinct; raise it if
#                 your real data is noisier than the synthetic test
#                 below and still over-picks.
# prominence_mult: a peak's summit must rise at least this many times
#                 the local noise level above the higher of its two
#                 flanking valleys to count as real. Main knob for peak
#                 count -- raise for fewer, more conservative peaks.
# min_height_mult: the summit itself must also clear this many times
#                 the local noise level, not just its rise above the
#                 valleys.
# min_width_ppm:  NEW. The peak's width at half its prominence (how far
#                 across, in ppm, between where the smoothed trace
#                 crosses halfway back down on each side) must be at
#                 least this wide. This is the "peak separation" guard
#                 against anything that's tall and locally prominent
#                 but only a point or two wide -- a real resolvable NMR
#                 line always has some minimum physical width; a true
#                 single-point artifact that slipped past smoothing
#                 doesn't.
#
# Known limits, found while testing this against the comparison table
# you put together:
#  - It can only find a peak that HAS its own valley on both sides. A
#    genuine shoulder with no dip at all between it and a much taller
#    neighbour (confirmed with a synthetic test: a peak at 2.50 and a
#    smaller one at 2.47, close enough that there's no valley between
#    them) is reported as ONE peak, not two. refine_peaks_lorentzian()
#    below can refine what's already found, including splitting
#    already-detected close peaks apart more precisely, but it can't by
#    itself invent a component nothing here ever flagged -- an attempt
#    at that (fit one peak, look for a leftover bump in the residual,
#    add it as a second peak, repeat) was tried and scrapped: with
#    thousands of points in a window, ordinary fit imperfection plus
#    noise reliably produces a residual bump somewhere that clears a
#    fixed threshold by chance, and it ended up splitting even a single
#    clean, isolated peak into 2-3 fake components in testing. True
#    hidden-shoulder recovery needs something more principled (e.g. a
#    reference database of known lineshapes/multiplet patterns), not a
#    residual-chasing heuristic.
#  - In the one synthetic case built specifically to stress it (a
#    single very broad peak filling almost the ENTIRE analyzed range,
#    with no normal baseline anywhere nearby to anchor the local noise
#    estimate), this over-detected -- several small spurious peaks
#    trailing off the real peak's flank. That specific shape doesn't
#    occur in your real spectra outside the excluded water region (a
#    real 1D proton spectrum always has normal-width metabolite peaks
#    and genuine flat stretches throughout the rest of the range), and
#    the realistic full-spectrum synthetic test below (55 planted peaks
#    across a real-sized 9 ppm range, including ones right next to a
#    broad excluded water blob) shows no sign of this failure mode --
#    but it's a known boundary condition worth knowing about rather
#    than claiming this handles every possible input.
#  - CONFIRMED, explained finding from testing (not a guess): airPLS
#    deliberately refuses to absorb ANY peak-shaped feature into the
#    baseline, including a broad one -- that's the whole point of
#    switching to it. The old percentile-in-window baseline did, by
#    accident, partly absorb a VERY broad hump like the residual water
#    peak into itself, simply because its window was much narrower than
#    the hump (so a narrow window sitting on the hump's slow slope has
#    its own local 5th percentile close to the hump's local height,
#    wrongly treating the broad trend as "baseline" there). Switching
#    to airPLS loses that accidental side benefit: the water peak's own
#    slowly-decaying tail, past the hard-coded edge of `exclude`, is now
#    real (if small) residual signal, and a noise wiggle riding on that
#    slow decay can clear the detection threshold right next to the
#    boundary, as confirmed in testing. Rather than silently widen
#    `exclude` with a guessed number (your real water peak's true width
#    was established from real TopSpin inspection, not a synthetic
#    test, and I can't verify a better number without your real data),
#    every peak within `boundary_margin_ppm` of an exclude boundary is
#    flagged in the `near_exclusion_boundary` output column instead --
#    check those specifically against TopSpin before trusting them, and
#    widen `exclude` yourself if they turn out to be the water tail
#    rather than a real metabolite peak.
find_shape_based_peaks <- function(ppm, intensity,
                                    exclude = list(c(3.9, 5.3), c(-0.05, 0.05)),
                                    airpls_lambda = 1e5,
                                    control_point_spacing_ppm = 0.02,
                                    sg_window_ppm = 0.006,
                                    sg_poly = 3,
                                    prominence_mult = 5,
                                    min_height_mult = 2,
                                    min_width_ppm = 0.003,
                                    boundary_margin_ppm = 0.3,
                                    precomputed_baseline = NULL) {
  # precomputed_baseline: pass in a baseline already fit for this exact
  # (ppm, intensity) pair (same length/order, raw not corrected) to skip
  # airPLS entirely here -- airPLS is the expensive step in this whole
  # function, and a caller that ALSO needs the baseline for something
  # else (another detector, Lorentzian refinement, a raw-signal fallback
  # check) would otherwise pay for it again from scratch every time. NULL
  # (the default) fits it here exactly as before, nothing changes for any
  # existing caller that doesn't pass this.
  ord <- order(ppm)
  ppm <- ppm[ord]; intensity <- intensity[ord]
  n <- length(ppm)
  grid_spacing <- stats::median(diff(ppm))

  excl_mask <- rep(FALSE, n)
  for (rng in exclude) excl_mask <- excl_mask | (ppm >= rng[1] & ppm <= rng[2])

  if (!is.null(precomputed_baseline)) {
    if (length(precomputed_baseline) != n) {
      stop("precomputed_baseline length (", length(precomputed_baseline),
           ") doesn't match ppm/intensity length (", n, ") -- make sure it ",
           "was fit on this exact (ppm, intensity) pair, same order, before passing it in.")
    }
    baseline <- precomputed_baseline[ord]
  } else {
    baseline <- fit_airpls_baseline(ppm, intensity,
                                     control_point_spacing_ppm = control_point_spacing_ppm,
                                     lambda = airpls_lambda)
  }
  corrected <- intensity - baseline

  # local noise, same control-point + interpolation pattern used
  # throughout this project for speed on a real 30-60k point spectrum --
  # via .window_bounds()'s binary search, not a full rescan per point
  cp_ppm <- seq(min(ppm), max(ppm), by = control_point_spacing_ppm)
  half_win <- 0.04
  bounds <- .window_bounds(ppm, cp_ppm, half_win)
  cp_noise <- numeric(length(cp_ppm))
  for (i in seq_along(cp_ppm)) {
    cp_noise[i] <- if (bounds$hi[i] >= bounds$lo[i]) {
      .noise_from_negative_residuals(corrected[bounds$lo[i]:bounds$hi[i]])
    } else {
      NA
    }
  }
  ok <- !is.na(cp_noise)
  local_noise <- stats::approx(cp_ppm[ok], cp_noise[ok], xout = ppm, rule = 2)$y

  sg_half_win <- max(2, round((sg_window_ppm / grid_spacing) / 2))
  smoothed <- savgol_filter(corrected, sg_half_win, sg_poly)

  # work within each contiguous non-excluded run separately, so a run
  # boundary (next to an excluded zone) is always treated as a valley
  keep_mask <- !excl_mask
  runs <- rle(keep_mask)
  starts <- cumsum(c(1, utils::head(runs$lengths, -1)))
  ends <- cumsum(runs$lengths)

  all_peaks <- data.frame(ppm = numeric(0), intensity = numeric(0),
                           prominence = numeric(0), width_ppm = numeric(0))

  for (r in seq_along(runs$lengths)) {
    if (!runs$values[r]) next
    idx <- starts[r]:ends[r]
    if (length(idx) < 5) next
    ys <- smoothed[idx]   # smoothed trace, for shape/topology only
    yr <- corrected[idx]  # unsmoothed baseline-corrected trace, for reported height
    x <- ppm[idx]
    nn <- length(ys)

    is_max <- rep(FALSE, nn); is_min <- rep(FALSE, nn)
    for (i in 2:(nn - 1)) {
      if (ys[i] > ys[i - 1] && ys[i] >= ys[i + 1]) is_max[i] <- TRUE
      if (ys[i] < ys[i - 1] && ys[i] <= ys[i + 1]) is_min[i] <- TRUE
    }
    max_idx <- which(is_max)
    min_idx <- c(1, which(is_min), nn)  # run boundaries count as valleys too

    for (mi in max_idx) {
      lefts <- min_idx[min_idx < mi]; rights <- min_idx[min_idx > mi]
      if (length(lefts) == 0 || length(rights) == 0) next
      left_valley <- max(lefts); right_valley <- min(rights)
      flank_height <- max(ys[left_valley], ys[right_valley])
      prominence <- ys[mi] - flank_height
      noise_here <- local_noise[idx[mi]]

      # width at half-prominence: walk out from the summit until the
      # smoothed trace drops to halfway between the summit and the
      # higher flanking valley, on each side
      half_h <- flank_height + prominence / 2
      li <- mi; while (li > left_valley && ys[li] > half_h) li <- li - 1
      ri <- mi; while (ri < right_valley && ys[ri] > half_h) ri <- ri + 1
      width_ppm <- abs(x[ri] - x[li])

      if (prominence > prominence_mult * noise_here &&
          ys[mi] > min_height_mult * noise_here &&
          width_ppm >= min_width_ppm) {
        # refine: report this peak's position/height from the RAW
        # (unsmoothed) corrected trace in a small window around the
        # smoothed summit, since smoothing can shave a little off a
        # sharp peak's true height
        refine_win <- max(1, sg_half_win)
        lo <- max(1, mi - refine_win); hi <- min(nn, mi + refine_win)
        best <- which.max(yr[lo:hi]) + lo - 1
        all_peaks <- rbind(all_peaks,
                            data.frame(ppm = x[best], intensity = yr[best],
                                       prominence = prominence, width_ppm = width_ppm))
      }
    }
  }

  all_peaks <- all_peaks[order(-all_peaks$ppm), ]
  rownames(all_peaks) <- NULL
  all_peaks <- all_peaks[!duplicated(all_peaks$ppm), ]  # the refine step can land two summits on the same raw point

  near_boundary <- rep(FALSE, nrow(all_peaks))
  for (rng in exclude) {
    near_boundary <- near_boundary |
      (all_peaks$ppm >= rng[1] - boundary_margin_ppm & all_peaks$ppm <= rng[1]) |
      (all_peaks$ppm >= rng[2] & all_peaks$ppm <= rng[2] + boundary_margin_ppm)
  }
  all_peaks$near_exclusion_boundary <- near_boundary

  all_peaks[, c("ppm", "intensity", "width_ppm", "near_exclusion_boundary")]
}

# ==================================================================
# PART 2b -- second-derivative helpers for find_peaks_second_derivative()
#
# Same local-polynomial machinery as savgol_coeffs()/savgol_filter()
# above, but extracting the QUADRATIC term's coefficient instead of the
# constant term -- for a local fit y = a0 + a1*u + a2*u^2 + ... in index
# units u, d2y/du2 = 2*a2, and converting index units to real ppm units
# needs one more factor of grid_spacing^-2 (chain rule, u = (ppm-center)/
# grid_spacing).
savgol_second_deriv_coeffs <- function(half_win, poly_order) {
  win <- -half_win:half_win
  Van <- outer(win, 0:poly_order, `^`)
  # Same QR-based fix as savgol_coeffs() above, and more important here:
  # this function needs wider windows (curvature estimates need more
  # smoothing) and a higher poly_order (4, not 3) by default, which is
  # exactly the combination that broke the old crossprod-based approach
  # outright on first real use against a fine-grid spectrum (a 273-point
  # window at poly_order 4 hit "computationally singular").
  qr_van <- qr(Van)
  qr.coef(qr_van, diag(nrow(Van)))[3, ]
}

savgol_second_deriv_filter <- function(y, half_win, poly_order, grid_spacing) {
  n <- length(y)
  win_pts <- 2 * half_win + 1
  if (win_pts > n || poly_order < 2) return(rep(0, n))
  coeffs <- savgol_second_deriv_coeffs(half_win, poly_order)
  padded <- c(rep(y[1], half_win), y, rep(y[n], half_win))
  out <- numeric(n)
  for (i in seq_len(n)) out[i] <- sum(coeffs * padded[i:(i + win_pts - 1)])
  out * 2 / (grid_spacing^2)
}

# ==================================================================
# PART 3b -- find_peaks_second_derivative()
#
# A DIFFERENT detection strategy from find_shape_based_peaks() above,
# not a retuned version of it -- built specifically for the case that
# method cannot solve by construction: two or more real lines close
# enough together that their SUM shows no dip between them anywhere in
# the raw trace. find_shape_based_peaks() can only ever report a peak
# that has its own valley on both sides; no threshold on that function
# changes that, it's a structural property of valley-to-valley
# detection, already documented in that function's own notes above.
#
# This looks at the NEGATIVE SECOND DERIVATIVE of the baseline-corrected
# trace instead of the trace itself -- a standard technique in analytical
# spectroscopy (derivative spectroscopy, used for exactly this problem in
# UV-Vis/IR/chromatography). Near the top of ANY individual line, real or
# hidden inside a shoulder, the curve is concave-down, so -y'' has its
# own local maximum there. Curvature is a local, second-order property,
# so two overlapping lines can each still produce their own distinct
# -y'' maximum even when their sum is one smooth hump with no valley at
# all in the raw trace -- that's the actual mechanism that resolves a
# shoulder, not a parameter tweak on the old method.
#
# Honest trade-offs, found in testing, not assumed:
#  - Differentiating amplifies noise, but the fix is NOT simply "use a
#    wider smoothing window on the derivative step" -- that was tried
#    first and failed two different ways depending on which direction it
#    was pushed. Too narrow a derivative window lets ordinary point
#    noise dominate the curvature trace, which plants a spurious
#    noise-level valley right next to a real peak's summit (confirmed:
#    a known isolated real peak came back with its "nearest valley"
#    only ~0.001 ppm away -- a noise wiggle, not its true, much more
#    distant flank) and that same wiggle inflates the local noise
#    estimate enough that nothing clears threshold. Too WIDE a window,
#    or too high a polynomial order (4-5), instead causes genuine
#    numerical ringing: confirmed directly on a clean, noise-free,
#    perfectly isolated synthetic Gaussian, which should only ever have
#    ONE curvature maximum, but came back with 3 once the window got
#    comparable to the peak's own width. Widening or narrowing a
#    SINGLE smoothing step cannot satisfy both constraints: one wants
#    more smoothing, the other wants less.
#  - The fix that actually worked, confirmed by test, is TWO separate
#    smoothing steps with two different jobs: first plain-smooth the
#    baseline-corrected VALUE trace (via the existing savgol_filter(),
#    moderate window, poly=3) to knock down point noise BEFORE any
#    derivative is taken, THEN take the second derivative of that
#    already-smoothed trace with its own, narrower window (poly=2,
#    which rang far less than higher orders in the clean-Gaussian test).
#    On the same known real peak used above, this produced a true,
#    well-separated valley (~0.005 ppm away, consistent with the peak's
#    own linewidth) with strong prominence, where either single-stage
#    attempt gave a spurious valley 5-10x closer. presmooth_window_ppm
#    and sg_window_ppm below are these two SEPARATE windows -- they are
#    not redundant with each other, both are needed.
#  - It will not help two lines so close that even their SEPARATE
#    curvature features merge into one -- it pushes the resolution
#    limit further than the valley-based method, it does not remove the
#    limit. Tested honestly against a synthetic stress benchmark (74
#    real lines in a 1.1 ppm window, 20 of them genuine no-valley
#    shoulder pairs) across 3 random seeds, AT MATCHED thresholds: this
#    detector was NOT a clear, reliable win over find_shape_based_peaks()
#    on that benchmark -- total lines recovered and no-valley pairs
#    recovered came out roughly a wash across seeds (sometimes a bit
#    better, sometimes a bit worse), not the one-sided improvement the
#    underlying theory would suggest. Don't take this function on faith
#    because the mechanism sounds right; it earns its place next to
#    find_shape_based_peaks() as a second, differently-biased detector
#    worth comparing against on real data, not a replacement for it.
#  - What IS real and reproducible: unlike find_shape_based_peaks(),
#    where sweeping prominence_mult/min_height_mult across a wide range
#    changed nothing at all (confirmed on this same benchmark AND
#    separately on real user spectra -- those two knobs are simply not
#    the binding constraint for that function), this detector's
#    thresholds have genuine pull: on the 74-line benchmark, recall rose
#    from 24% (prominence_mult=4) to ~45% (prominence_mult=1.2-1.5) with
#    PERFECT precision throughout (every single detection a real planted
#    line, zero false positives) -- that 1.2-1.5 range is where
#    precision first starts to measurably slip (97%), and it falls off a
#    cliff below prominence_mult=0.4 (75% precision, then 52%, then 43%
#    as it's pushed to 0.1) while recall only crawls up to ~60% even
#    there, with HALF of everything "detected" by that point being noise,
#    not real peaks. Defaults below sit at prominence_mult=1.5,
#    min_height_mult=0.75 -- real, tested headroom over the original
#    4/2, short of where precision starts paying for it.
#  - None of this gets anywhere near large three-digit peak counts in a
#    1 ppm-scale window by legitimate detection. Even deliberately
#    accepting substantial noise contamination (prominence_mult=0.1,
#    where under half of what's reported is real) topped out around 100
#    total detections on this 74-line stress benchmark -- which already
#    has far denser real structure than a typical window. A spectrum
#    window where a person visually counts on the order of 300 distinct
#    features is very likely showing unresolved multiplet fine structure
#    and close higher-order splitting beyond what any single-component
#    threshold-based picker, this one included, can respectably call out
#    as 300 separate "peaks" -- getting a number that large would mean
#    lowering thresholds into the noise floor, not finding more signal.
#  - It reports each resolved component's position and a first-pass
#    height read straight off the raw corrected trace at that point, NOT
#    a proper deconvolved height (the true individual contribution of an
#    overlapping line is biased upward by its neighbour's tail at that
#    same point). Feed this function's output into
#    refine_peaks_lorentzian() with a cluster_tol wide enough to catch
#    the close pairs it found -- that function already fits overlapping
#    peaks TOGETHER as one multi-component model, and this function's
#    real contribution is telling it HOW MANY components are really
#    there and roughly where, which it has no way to know on its own.
#
# exclude, airpls_lambda, control_point_spacing_ppm, boundary_margin_ppm:
#             same meaning as find_shape_based_peaks().
# presmooth_window_ppm: window for the FIRST smoothing pass, applied to
#             the plain value trace before any derivative is taken.
# sg_window_ppm, sg_poly: window/order for the SECOND pass, which
#             estimates the curvature (second derivative) of the
#             already-presmoothed trace -- NOT the same tuning as the
#             shape detector's own sg_window_ppm, don't assume one value
#             works for both functions, and don't assume this one
#             value does both smoothing jobs at once (see trade-offs
#             above for why that doesn't work).
# prominence_mult: a candidate's curvature bump must rise at least this
#             many times the LOCAL noise level of the curvature trace
#             itself above its flanking curvature valleys.
# min_height_mult: the candidate's RAW corrected height at that point
#             (not its curvature) must also clear this many times the
#             raw trace's own local noise -- without this check, a
#             perfectly real-looking curvature wiggle sitting on an
#             already-tiny, noise-level bump would still get reported.
find_peaks_second_derivative <- function(ppm, intensity,
                                          exclude = list(c(3.9, 5.3), c(-0.05, 0.05)),
                                          airpls_lambda = 1e5,
                                          control_point_spacing_ppm = 0.02,
                                          presmooth_window_ppm = 0.003,
                                          sg_window_ppm = 0.003,
                                          sg_poly = 2,
                                          prominence_mult = 1.5,
                                          min_height_mult = 0.75,
                                          boundary_margin_ppm = 0.3,
                                          precomputed_baseline = NULL) {
  # precomputed_baseline: same meaning as find_shape_based_peaks()'s own
  # parameter of the same name -- pass in a baseline already fit for this
  # exact (ppm, intensity) pair to skip airPLS here. See that function's
  # comment for why this exists.
  ord <- order(ppm)
  ppm <- ppm[ord]; intensity <- intensity[ord]
  n <- length(ppm)
  grid_spacing <- stats::median(diff(ppm))

  excl_mask <- rep(FALSE, n)
  for (rng in exclude) excl_mask <- excl_mask | (ppm >= rng[1] & ppm <= rng[2])

  if (!is.null(precomputed_baseline)) {
    if (length(precomputed_baseline) != n) {
      stop("precomputed_baseline length (", length(precomputed_baseline),
           ") doesn't match ppm/intensity length (", n, ") -- make sure it ",
           "was fit on this exact (ppm, intensity) pair, same order, before passing it in.")
    }
    baseline <- precomputed_baseline[ord]
  } else {
    baseline <- fit_airpls_baseline(ppm, intensity,
                                     control_point_spacing_ppm = control_point_spacing_ppm,
                                     lambda = airpls_lambda)
  }
  corrected <- intensity - baseline

  cp_ppm <- seq(min(ppm), max(ppm), by = control_point_spacing_ppm)
  half_win_noise <- 0.04
  bounds <- .window_bounds(ppm, cp_ppm, half_win_noise)
  cp_noise <- numeric(length(cp_ppm))
  for (i in seq_along(cp_ppm)) {
    cp_noise[i] <- if (bounds$hi[i] >= bounds$lo[i]) {
      .noise_from_negative_residuals(corrected[bounds$lo[i]:bounds$hi[i]])
    } else {
      NA
    }
  }
  ok <- !is.na(cp_noise)
  local_noise <- stats::approx(cp_ppm[ok], cp_noise[ok], xout = ppm, rule = 2)$y

  presmooth_half_win <- max(3, round((presmooth_window_ppm / grid_spacing) / 2))
  sg_half_win <- max(3, round((sg_window_ppm / grid_spacing) / 2))

  keep_mask <- !excl_mask
  runs <- rle(keep_mask)
  starts <- cumsum(c(1, utils::head(runs$lengths, -1)))
  ends <- cumsum(runs$lengths)

  all_peaks <- data.frame(ppm = numeric(0), intensity = numeric(0))

  for (r in seq_along(runs$lengths)) {
    if (!runs$values[r]) next
    idx <- starts[r]:ends[r]
    nn <- length(idx)
    if (nn < (2 * sg_half_win + 1) || nn < (2 * presmooth_half_win + 1)) next
    x <- ppm[idx]; yr <- corrected[idx]; noise_idx <- local_noise[idx]

    # TWO stages, not one -- found necessary in testing, not assumed.
    # Differentiating straight off the noisy corrected trace amplifies
    # ordinary point-to-point measurement noise badly enough that the
    # resulting curvature trace is dominated by small-scale wiggle: a
    # real peak's TRUE flanking valley (the one that actually separates
    # it from its neighbours) gets shadowed by a spurious noise-level
    # valley sitting just a few points away, which then gets used as
    # "the" valley and robs the peak of almost all its real prominence --
    # confirmed directly: a known, isolated, obviously-real peak in a
    # real-noise benchmark came back with a "nearest valley" just 0.001
    # ppm away and a measured noise level so inflated by that same
    # wiggle that NOTHING in the whole window could clear threshold.
    # Pre-smoothing the trace (plain value smoothing, not yet a
    # derivative) with its own, separate window first suppresses that
    # point noise before curvature is ever computed, which fixed it in
    # testing (the same known peak's valley moved out to a real ~0.005
    # ppm, consistent with its own linewidth, once this was added).
    # A second, independent risk found in testing, going the OTHER way:
    # if the derivative window itself is too WIDE relative to a real
    # peak's own width, the second derivative starts ringing even on a
    # perfectly clean, noise-free, single isolated Gaussian -- a pure
    # numerical artifact of fitting too high a local polynomial order
    # over too wide a window relative to the feature it's measuring, not
    # overlapping real peaks. Confirmed directly against a clean
    # synthetic single peak: poly_order 4 with a window comparable to or
    # wider than the peak's own width reliably produced 3 spurious local
    # maxima where exactly 1 exists; poly_order 2 with a window narrower
    # than the narrowest real peak width in this project's own
    # benchmarks (0.003 ppm) did not. That's why sg_poly defaults to 2
    # here, not a higher order that would track curvature "better" on
    # paper -- higher order measurably rang on real peak shapes in
    # testing, this didn't.
    presmoothed <- savgol_filter(yr, presmooth_half_win, 3)
    d2 <- savgol_second_deriv_filter(presmoothed, sg_half_win, sg_poly, grid_spacing)
    neg_d2 <- -d2

    cp_x <- seq(min(x), max(x), by = control_point_spacing_ppm)
    half_win_noise2 <- max(control_point_spacing_ppm * 3, 0.04)
    bounds2 <- .window_bounds(x, cp_x, half_win_noise2)
    cp_d2noise <- numeric(length(cp_x))
    for (i in seq_along(cp_x)) {
      if (bounds2$hi[i] < bounds2$lo[i]) { cp_d2noise[i] <- NA; next }
      vals <- neg_d2[bounds2$lo[i]:bounds2$hi[i]]
      if (length(vals) < 5) { cp_d2noise[i] <- NA; next }
      cp_d2noise[i] <- .noise_from_negative_residuals(vals - stats::median(vals))
    }
    ok2 <- !is.na(cp_d2noise)
    if (sum(ok2) < 2) next
    local_d2_noise <- stats::approx(cp_x[ok2], cp_d2noise[ok2], xout = x, rule = 2)$y

    is_max <- rep(FALSE, nn); is_min <- rep(FALSE, nn)
    for (i in 2:(nn - 1)) {
      if (neg_d2[i] > neg_d2[i - 1] && neg_d2[i] >= neg_d2[i + 1]) is_max[i] <- TRUE
      if (neg_d2[i] < neg_d2[i - 1] && neg_d2[i] <= neg_d2[i + 1]) is_min[i] <- TRUE
    }
    max_idx <- which(is_max)
    min_idx <- c(1, which(is_min), nn)

    for (mi in max_idx) {
      lefts <- min_idx[min_idx < mi]; rights <- min_idx[min_idx > mi]
      if (length(lefts) == 0 || length(rights) == 0) next
      left_valley <- max(lefts); right_valley <- min(rights)
      flank <- max(neg_d2[left_valley], neg_d2[right_valley])
      prom <- neg_d2[mi] - flank
      d2n <- local_d2_noise[mi]
      if (is.na(d2n) || d2n <= 0) next
      if (prom <= prominence_mult * d2n) next
      if (is.na(noise_idx[mi]) || yr[mi] <= min_height_mult * noise_idx[mi]) next
      if (yr[mi] <= 0) next

      all_peaks <- rbind(all_peaks, data.frame(ppm = x[mi], intensity = yr[mi]))
    }
  }

  all_peaks <- all_peaks[order(-all_peaks$ppm), ]
  rownames(all_peaks) <- NULL
  all_peaks <- all_peaks[!duplicated(all_peaks$ppm), ]

  near_boundary <- rep(FALSE, nrow(all_peaks))
  for (rng in exclude) {
    near_boundary <- near_boundary |
      (all_peaks$ppm >= rng[1] - boundary_margin_ppm & all_peaks$ppm <= rng[1]) |
      (all_peaks$ppm >= rng[2] & all_peaks$ppm <= rng[2] + boundary_margin_ppm)
  }
  all_peaks$near_exclusion_boundary <- near_boundary

  all_peaks[, c("ppm", "intensity", "near_exclusion_boundary")]
}

# ==================================================================
# PART 3c -- find_shape_based_peaks_fragmented()
#
# Same detector as find_shape_based_peaks() above, but instead of
# fitting the baseline and local noise level ONCE across the full
# spectrum and only filtering excluded-region peaks out of the
# RESULT, this physically cuts the excluded range(s) out of the ppm/
# intensity arrays FIRST, and fits a separate airPLS baseline and a
# separate noise estimate independently on each remaining fragment.
#
# Why this can help here, and why it's safe to try (unlike ASICS):
# ASICS's crash when fragmented came from a per-COMPOUND REFERENCE
# LIBRARY correlation step -- a compound with no real signal inside a
# narrow fragment produces a flat/constant vector, and correlating a
# flat vector against anything is mathematically undefined (0/0),
# which is what broke ASICS's internal matrix math, and that's a
# structural fact about ANY two-way split of its library, confirmed
# numerically against both ASICS's own library and your custom one.
# This detector does no such thing -- it never compares one
# compound's shape to another's, it only looks at ONE sample's OWN
# baseline/noise/peak shape, entirely locally. There's no flat-vector-
# correlation failure mode here to break, at any fragment width.
#
# What fragmenting DOES change, honestly: find_shape_based_peaks()
# already treats the boundary of each excluded run as an implicit
# valley for PEAK-SHAPE detection (see its own code -- run boundaries
# count as valleys too), so for topology/shape purposes the two
# approaches already behave the same. The one real difference is in
# airPLS's own solve: fit over the FULL spectrum, every control point
# is coupled to its neighbours by the smoothness penalty, including
# the control points sitting under the excluded hump itself (down-
# weighted during reweighting, but still part of the same global
# solve) -- so there's a theoretical, if small, path for the excluded
# region's shape to influence the fitted baseline just outside it.
# Fragmenting first removes that path completely: a fragment's airPLS
# solve has no control points from the other side of the cut AT ALL,
# so nothing there can leak in, even indirectly.
#
# What fragmenting does NOT fix: the already-identified cause of
# spurious near-boundary peaks (see find_shape_based_peaks()'s own
# notes above) is a REAL residual tail of the excluded peak (e.g.
# residual water) extending, in the raw data, past wherever `exclude`'s
# hard-coded edge happens to be set -- not a computation artifact.
# Fragmenting at that SAME edge still leaves that real tail sitting
# just outside the cut, same as the exclude-based version does, so by
# itself this will not make a too-narrow exclude/fragment boundary
# stop tripping the near-boundary flag -- confirmed side by side in
# the comparison self-test below. Use this when you want the slightly
# more rigorous no-cross-boundary-coupling-at-all guarantee; still
# check near_exclusion_boundary the same way either method, and widen
# the boundary itself (not switch methods) if those turn out to be
# real excluded-peak tail rather than genuine metabolite signal.
#
# Arguments are identical to find_shape_based_peaks(); exclude's
# complement (everything NOT in one of these ranges, within
# [min(ppm), max(ppm)]) defines the fragment boundaries that get
# processed independently.
find_shape_based_peaks_fragmented <- function(ppm, intensity,
                                               exclude = list(c(3.9, 5.3), c(-0.05, 0.05)),
                                               airpls_lambda = 1e5,
                                               control_point_spacing_ppm = 0.02,
                                               sg_window_ppm = 0.006,
                                               sg_poly = 3,
                                               prominence_mult = 5,
                                               min_height_mult = 2,
                                               min_width_ppm = 0.003,
                                               boundary_margin_ppm = 0.3) {
  ord <- order(ppm)
  ppm <- ppm[ord]; intensity <- intensity[ord]
  grid_spacing <- stats::median(diff(ppm))

  lo <- min(ppm); hi <- max(ppm)
  if (length(exclude) == 0) {
    fragments <- list(c(lo, hi))
  } else {
    starts <- sapply(exclude, `[`, 1); ends <- sapply(exclude, `[`, 2)
    o <- order(starts)
    starts <- starts[o]; ends <- ends[o]
    m_starts <- starts[1]; m_ends <- ends[1]
    for (i in seq_along(starts)[-1]) {
      if (starts[i] <= m_ends[length(m_ends)]) {
        m_ends[length(m_ends)] <- max(m_ends[length(m_ends)], ends[i])
      } else {
        m_starts <- c(m_starts, starts[i]); m_ends <- c(m_ends, ends[i])
      }
    }
    fragments <- list(); cur <- lo
    for (i in seq_along(m_starts)) {
      if (m_starts[i] > cur) fragments[[length(fragments) + 1]] <- c(cur, min(m_starts[i], hi))
      cur <- max(cur, m_ends[i])
      if (cur >= hi) break
    }
    if (cur < hi) fragments[[length(fragments) + 1]] <- c(cur, hi)
  }

  all_peaks <- data.frame(ppm = numeric(0), intensity = numeric(0),
                           prominence = numeric(0), width_ppm = numeric(0))

  for (frag in fragments) {
    in_frag <- ppm >= frag[1] & ppm <= frag[2]
    if (sum(in_frag) < 5) next
    fp <- ppm[in_frag]; fi <- intensity[in_frag]

    baseline <- fit_airpls_baseline(fp, fi,
                                     control_point_spacing_ppm = control_point_spacing_ppm,
                                     lambda = airpls_lambda)
    corrected <- fi - baseline

    cp_ppm <- seq(min(fp), max(fp), by = control_point_spacing_ppm)
    half_win <- 0.04
    cp_noise <- numeric(length(cp_ppm))
    for (i in seq_along(cp_ppm)) {
      in_win <- fp >= (cp_ppm[i] - half_win) & fp <= (cp_ppm[i] + half_win)
      cp_noise[i] <- .noise_from_negative_residuals(corrected[in_win])
    }
    ok <- !is.na(cp_noise)
    if (sum(ok) < 2) next
    local_noise <- stats::approx(cp_ppm[ok], cp_noise[ok], xout = fp, rule = 2)$y

    sg_half_win <- max(2, round((sg_window_ppm / grid_spacing) / 2))
    smoothed <- savgol_filter(corrected, sg_half_win, sg_poly)

    nn <- length(fp)
    is_max <- rep(FALSE, nn); is_min <- rep(FALSE, nn)
    for (i in 2:(nn - 1)) {
      if (smoothed[i] > smoothed[i - 1] && smoothed[i] >= smoothed[i + 1]) is_max[i] <- TRUE
      if (smoothed[i] < smoothed[i - 1] && smoothed[i] <= smoothed[i + 1]) is_min[i] <- TRUE
    }
    max_idx <- which(is_max)
    min_idx <- c(1, which(is_min), nn)  # the fragment's own cut edges count as valleys too

    for (mi in max_idx) {
      lefts <- min_idx[min_idx < mi]; rights <- min_idx[min_idx > mi]
      if (length(lefts) == 0 || length(rights) == 0) next
      left_valley <- max(lefts); right_valley <- min(rights)
      flank_height <- max(smoothed[left_valley], smoothed[right_valley])
      prominence <- smoothed[mi] - flank_height
      noise_here <- local_noise[mi]

      half_h <- flank_height + prominence / 2
      li <- mi; while (li > left_valley && smoothed[li] > half_h) li <- li - 1
      ri <- mi; while (ri < right_valley && smoothed[ri] > half_h) ri <- ri + 1
      width_ppm <- abs(fp[ri] - fp[li])

      if (prominence > prominence_mult * noise_here &&
          smoothed[mi] > min_height_mult * noise_here &&
          width_ppm >= min_width_ppm) {
        refine_win <- max(1, sg_half_win)
        lo2 <- max(1, mi - refine_win); hi2 <- min(nn, mi + refine_win)
        best <- which.max(corrected[lo2:hi2]) + lo2 - 1
        all_peaks <- rbind(all_peaks,
                            data.frame(ppm = fp[best], intensity = corrected[best],
                                       prominence = prominence, width_ppm = width_ppm))
      }
    }
  }

  all_peaks <- all_peaks[order(-all_peaks$ppm), ]
  rownames(all_peaks) <- NULL
  all_peaks <- all_peaks[!duplicated(all_peaks$ppm), ]

  near_boundary <- rep(FALSE, nrow(all_peaks))
  for (rng in exclude) {
    near_boundary <- near_boundary |
      (all_peaks$ppm >= rng[1] - boundary_margin_ppm & all_peaks$ppm <= rng[1]) |
      (all_peaks$ppm >= rng[2] & all_peaks$ppm <= rng[2] + boundary_margin_ppm)
  }
  all_peaks$near_exclusion_boundary <- near_boundary

  all_peaks[, c("ppm", "intensity", "width_ppm", "near_exclusion_boundary")]
}

# ==================================================================
# PART 4 -- refine_peaks_lorentzian()  (OPTIONAL)
#
# Fits a real NMR lineshape (Lorentzian: h / (1 + ((ppm-x0)/w)^2)) to
# each peak find_shape_based_peaks() already found, refining its
# center/height/width and giving an actual integrated area (height x
# width x pi) -- use this when you need real quantification rather
# than just a peak list, or when two detected peaks are close enough
# that their tails bias each other's raw reported height.
#
# peak_tbl: the data.frame returned by find_shape_based_peaks() for
#           ONE sample (needs ppm, intensity columns).
# corrected, ppm: the SAME baseline-corrected trace and ppm grid used to
#           detect peak_tbl (corrected = intensity - baseline from
#           fit_airpls_baseline()).
# cluster_tol: peaks within this many ppm of each other are fit
#           TOGETHER as one multi-component model (so their overlap is
#           properly accounted for) rather than one at a time.
# window_pad: how far past the outermost peak in a cluster to include
#           when fitting, so the model has enough of the flanking
#           baseline-corrected trace to constrain the fit.
#
# Honest limit (see find_shape_based_peaks()'s own notes): this refines
# what's already been found. It does not invent a new component for a
# true shoulder that never had its own valley -- that was tried and
# removed because it overfit noise into fake extra peaks even on
# already-clean, well-separated test cases.
refine_peaks_lorentzian <- function(ppm, corrected, peak_tbl,
                                     cluster_tol = 0.03, window_pad = 0.04,
                                     center_bound_ppm = 0.015,
                                     hwhm_bounds_ppm = c(0.001, 0.05)) {
  # center_bound_ppm, hwhm_bounds_ppm: NEW, exposed rather than hardcoded,
  # for one reason found while building the synthetic multiplet benchmark
  # (accuracy/precision panels against known triplet/quartet ground
  # truth -- see that benchmark's own scripts).
  # These bounds were originally picked for THIS project's real 600 MHz
  # E. coli data, where a real linewidth is essentially always
  # >=0.007-0.01 ppm (several Hz) -- 0.001 ppm as a LOWER bound was never
  # meant to be a meaningful physical floor there, just a "don't let the
  # solver collapse to a spike" guard nowhere near real data. A narrower
  # synthetic line (the benchmark's lw=1 Hz default, matching the
  # MultiNMRFit comparison figure it's modeled on) can sit BELOW that
  # floor at a high field (1 Hz FWHM = 0.5 Hz HWHM = 0.001042 ppm at 800
  # MHz, 0.000625 ppm at 600 MHz) -- fitting with the old hardcoded bound
  # would then measure the bound, not the detector's real resolving
  # power, and silently bias every linewidth estimate high regardless of
  # S/N. Defaults are UNCHANGED from the old hardcoded values, so no
  # existing caller's behavior changes; the benchmark passes its own,
  # tighter bound explicitly and reports results both ways (see that
  # script's own notes).
  if (nrow(peak_tbl) == 0) {
    return(data.frame(ppm = numeric(0), intensity = numeric(0),
                       hwhm = numeric(0), area = numeric(0), converged = logical(0)))
  }
  ord <- order(peak_tbl$ppm)
  peak_ppm <- peak_tbl$ppm[ord]; peak_height <- peak_tbl$intensity[ord]
  grp <- cumsum(c(1, diff(peak_ppm) > cluster_tol))

  fit_cluster <- function(x, y, centers0, heights0, widths0 = rep(0.006, length(centers0))) {
    k <- length(centers0)
    start <- list()
    for (i in seq_len(k)) {
      start[[paste0("x0_", i)]] <- centers0[i]
      start[[paste0("h_", i)]]  <- heights0[i]
      start[[paste0("w_", i)]]  <- widths0[i]
    }
    terms <- sapply(seq_len(k), function(i) sprintf("(h_%d / (1 + ((x - x0_%d)/w_%d)^2))", i, i, i))
    form <- as.formula(paste("y ~", paste(terms, collapse = " + ")))
    lower <- c(); upper <- c()
    for (i in seq_len(k)) {
      lower <- c(lower, centers0[i] - center_bound_ppm, 0, hwhm_bounds_ppm[1])
      upper <- c(upper, centers0[i] + center_bound_ppm, Inf, hwhm_bounds_ppm[2])
    }
    tryCatch(
      stats::nls(form, data = list(x = x, y = y), start = start,
                 algorithm = "port", lower = lower, upper = upper,
                 control = stats::nls.control(maxiter = 200, warnOnly = TRUE)),
      error = function(e) NULL
    )
  }

  out <- list()
  for (g in unique(grp)) {
    idx <- which(grp == g)
    centers0 <- peak_ppm[idx]; heights0 <- peak_height[idx]
    lo <- min(centers0) - window_pad; hi <- max(centers0) + window_pad
    in_win <- ppm >= lo & ppm <= hi
    x <- ppm[in_win]; y <- corrected[in_win]
    k <- length(centers0)

    fit <- fit_cluster(x, y, centers0, heights0)
    if (is.null(fit)) {
      out[[length(out) + 1]] <- data.frame(ppm = centers0, intensity = heights0,
                                            hwhm = NA_real_, area = NA_real_, converged = FALSE)
    } else {
      cf <- coef(fit)
      rows <- data.frame(
        ppm       = sapply(seq_len(k), function(i) cf[paste0("x0_", i)]),
        intensity = sapply(seq_len(k), function(i) cf[paste0("h_", i)]),
        hwhm      = sapply(seq_len(k), function(i) cf[paste0("w_", i)])
      )
      rows$area <- pi * rows$intensity * rows$hwhm
      rows$converged <- isTRUE(fit$convInfo$isConv)
      out[[length(out) + 1]] <- rows
    }
  }
  res <- do.call(rbind, out)
  res <- res[order(-res$ppm), ]
  rownames(res) <- NULL
  res
}

# ==================================================================
# PART 5 -- find_baseline_corrected_peaks()  (legacy, kept for comparison)
#
# The original point-by-point version: flags a point as a peak if it's
# taller than both immediate neighbours AND clears noise_mult times the
# local noise -- no use of a peak's overall shape (its valleys), which
# is exactly what find_shape_based_peaks() above improves on.
find_baseline_corrected_peaks <- function(ppm, intensity,
                                           exclude = list(c(3.9, 5.3), c(-0.05, 0.05)),
                                           baseline_window_ppm = 0.08,
                                           baseline_quantile = 0.05,
                                           noise_mult = 6,
                                           min_separation_ppm = 0.01,
                                           control_point_spacing_ppm = 0.02) {
  ord <- order(ppm)
  ppm <- ppm[ord]; intensity <- intensity[ord]
  n <- length(ppm)

  excl_mask <- rep(FALSE, n)
  for (rng in exclude) excl_mask <- excl_mask | (ppm >= rng[1] & ppm <= rng[2])

  cp_ppm <- seq(min(ppm), max(ppm), by = control_point_spacing_ppm)
  cp_baseline <- numeric(length(cp_ppm))
  cp_noise <- numeric(length(cp_ppm))
  half_win <- baseline_window_ppm / 2
  for (i in seq_along(cp_ppm)) {
    in_win <- ppm >= (cp_ppm[i] - half_win) & ppm <= (cp_ppm[i] + half_win)
    vals <- intensity[in_win]
    if (length(vals) < 3) { cp_baseline[i] <- NA; cp_noise[i] <- NA; next }
    cp_baseline[i] <- stats::quantile(vals, baseline_quantile, na.rm = TRUE)
    resid <- vals - cp_baseline[i]
    cp_noise[i] <- .noise_from_negative_residuals(resid)
  }
  ok <- !is.na(cp_baseline) & !is.na(cp_noise)
  baseline <- stats::approx(cp_ppm[ok], cp_baseline[ok], xout = ppm, rule = 2)$y
  local_noise <- stats::approx(cp_ppm[ok], cp_noise[ok], xout = ppm, rule = 2)$y

  corrected <- intensity - baseline

  is_local_max <- logical(n)
  for (i in 2:(n - 1)) {
    is_local_max[i] <- corrected[i] > corrected[i - 1] && corrected[i] >= corrected[i + 1]
  }

  candidate <- is_local_max & !excl_mask & corrected > (noise_mult * local_noise) & corrected > 0
  cand_df <- data.frame(ppm = ppm[candidate], intensity = intensity[candidate],
                         corrected = corrected[candidate])
  cand_df <- cand_df[order(-cand_df$corrected), ]

  keep <- data.frame(ppm = numeric(0), intensity = numeric(0))
  claimed <- rep(FALSE, nrow(cand_df))
  for (i in seq_len(nrow(cand_df))) {
    if (claimed[i]) next
    row <- cand_df[i, c("ppm", "corrected")]
    colnames(row) <- c("ppm", "intensity")
    keep <- rbind(keep, row)
    claimed[which(abs(cand_df$ppm - cand_df$ppm[i]) <= min_separation_ppm)] <- TRUE
  }
  keep <- keep[order(-keep$ppm), ]
  rownames(keep) <- NULL
  attr(keep, "baseline") <- data.frame(ppm = ppm, baseline = baseline)
  keep
}

# ------------------------------------------------------------------
# SYNTHETIC VALIDATION -- not real data. Same fake dense spectrum used
# throughout this project: known baseline drift, a big broad "water"
# blob, ~55 known narrow peaks of varying height, plus the shoulder
# case and the pathological dominant-peak case discussed above.
if (interactive() || identical(Sys.getenv("RUN_BASELINE_PEAK_PICKING_SELFTEST"), "1")) {
  set.seed(42)
  ppm_grid <- seq(0.3, 9.5, length.out = 40000)
  drifting_baseline <- 0.02 * sin(ppm_grid / 2) + 0.01 * (ppm_grid - 5)^2 / 20
  water_blob <- 5 * exp(-((ppm_grid - 4.6)^2) / (2 * 0.35^2))
  noise <- rnorm(length(ppm_grid), 0, 0.002)
  true_peak_ppm <- sort(runif(55, 0.4, 9.3))
  true_peak_ppm <- true_peak_ppm[!(true_peak_ppm > 3.8 & true_peak_ppm < 5.4)]
  signal <- rep(0, length(ppm_grid))
  for (p in true_peak_ppm) {
    h <- runif(1, 0.05, 1.5); w <- runif(1, 0.003, 0.012)
    signal <- signal + h * exp(-((ppm_grid - p)^2) / (2 * w^2))
  }
  fake_intensity <- drifting_baseline + water_blob + signal + noise

  for (pm in c(3, 4, 5, 6, 8, 10)) {
    picked <- find_shape_based_peaks(ppm_grid, fake_intensity, prominence_mult = pm)
    recovered <- sapply(true_peak_ppm, function(p) any(abs(picked$ppm - p) < 0.02))
    cat(sprintf("prominence_mult=%2d -> picked=%4d  recovered=%d/%d\n",
                pm, nrow(picked), sum(recovered), length(true_peak_ppm)))
  }
  cat("Any picked in water zone:", any(picked$ppm >= 3.9 & picked$ppm <= 5.3), "\n")
  cat("Any negative intensity:", any(picked$intensity <= 0), "\n")
}

# ------------------------------------------------------------------
# Example usage on your real 6-sample set. Run the "import and stack"
# block first in this same session so spectra_obj, ppm, and spec_mat
# already exist -- this does NOT call ASICS() at all.
#
# ppm <- getPpmGrid(spectra_obj)
# spec_mat <- getSpectra(spectra_obj)
# if (nrow(spec_mat) != length(ppm)) spec_mat <- t(spec_mat)
#
# prominence_choice <- 5  # tune this first -- see sweep below
#
# per_sample <- lapply(colnames(spec_mat), function(s) {
#   find_shape_based_peaks(ppm, spec_mat[, s], prominence_mult = prominence_choice)
# })
# names(per_sample) <- colnames(spec_mat)
# for (s in names(per_sample)) cat(s, ":", nrow(per_sample[[s]]), "peaks\n")
#
# -- OPTIONAL: refine for real quantification (center/height/width/area) --
# baselines <- lapply(colnames(spec_mat), function(s) fit_airpls_baseline(ppm, spec_mat[, s]))
# names(baselines) <- colnames(spec_mat)
# refined <- lapply(colnames(spec_mat), function(s) {
#   corrected <- spec_mat[, s] - baselines[[s]]
#   refine_peaks_lorentzian(ppm, corrected, per_sample[[s]])
# })
# names(refined) <- colnames(spec_mat)
#
# -- merge across samples: reference-guided (recommended) --
# spectra_list <- lapply(colnames(spec_mat), function(s) {
#   list(ppm = ppm, corrected = spec_mat[, s] - baselines[[s]])
# })
# names(spectra_list) <- colnames(spec_mat)
# merged <- merge_peaks_reference_guided(per_sample, spectra_list, tol = 0.01)
#
# groups <- c(
#   "sample_control_01" = "control", "sample_control_02" = "control",
#   "sample_treatment_01" = "sample",   "sample_treatment_02" = "sample"
# )
# write_metaboanalyst_csv(merged, groups, "../output/control_vs_sample/peak_list_control_vs_sample.csv")
#
# -- Tuning prominence_mult to land near the peak count you want --
# for (pm in c(3, 4, 5, 6, 8, 10)) {
#   counts <- sapply(colnames(spec_mat), function(s) {
#     nrow(find_shape_based_peaks(ppm, spec_mat[, s], prominence_mult = pm))
#   })
#   cat("prominence_mult =", pm, "-> peaks per sample:", paste(counts, collapse = ", "), "\n")
# }
#
# If your real data is noisier than the synthetic test above and this
# still over-picks even at a high prominence_mult, raise sg_window_ppm a
# bit (e.g. 0.008-0.012) before raising prominence_mult further -- that's
# the knob controlling how much single-point noise wiggle gets smoothed
# away before shape detection, as opposed to how tall a peak has to be
# once shape detection has already run.
