# simulate_multiplet_overlap_benchmark.R
#
# REFERENCE IMPLEMENTATION, not the final benchmark -- a stand-in synthetic
# generator built to validate run_multiplet_overlap_benchmark.R end-to-end
# before a more detailed generator replaces it. Swap this out freely; the only thing run_multiplet_overlap_
# benchmark.R actually requires is the two output files below, in this
# exact shape (see that script's own header for the full contract).
#
# Models the same scenario as the MultiNMRFit comparison figure this
# benchmark is built against: one triplet and one quartet (1:2:1 and
# 1:3:3:1 relative line intensities, both with the same J), placed
# symmetrically around a common center ppm, separated by a controlled
# number of Hz between their own centers (how much the two multiplets'
# NEAREST lines overlap), at a range of S/N ratios, with several
# independent noise draws per (separation, S/N) cell so
# run_multiplet_overlap_benchmark.R can report BIAS (mean error vs. the
# known ground truth, across draws) and PRECISION (SD of that error) the
# same way the paper does, not just a single lucky/unlucky realization.
#
# field_mhz defaults to 600.17 -- this project's own real spectrometer
# field (see read_bruker_1d_direct.R's output on the real HAIRD/SWEPF
# data), not the paper's 800 MHz. The exact field only matters for
# converting the requested Hz values into ppm; keeping it at this
# project's real field keeps every ppm-unit threshold already tuned
# elsewhere in this codebase (cluster_tol, min_width_ppm, sg_window_ppm,
# etc.) in the same regime they were actually tuned for, so this
# benchmark measures the pipeline's real resolving power in the units it
# will actually be used in, not an artifact of a different field.
#
# S/N is defined exactly as the paper states: on the triplet's OWN
# downfield satellite -- i.e. its two equal, smaller outer lines (the
# 1's in the 1:2:1 pattern), whichever sits at higher ppm ("downfield").
# Noise sd is set so that line's peak HEIGHT divided by noise sd equals
# the target ratio.
#
# Output contract (written under ../example_data/multiplet_benchmark/):
#   ground_truth.csv -- one row per TRUE spectral line:
#     draw_id, hz_sep, snr_level, field_mhz, multiplet_type
#     ("triplet"/"quartet"), line_index (0-based, low to high ppm),
#     true_ppm, true_J_hz, true_lw_hz (FWHM), true_height, true_area
#   spectra_long.csv -- one row per spectral point per draw:
#     draw_id, ppm, intensity
# Any generator that writes these two files/columns can replace this one
# without changing run_multiplet_overlap_benchmark.R at all.
#
# Run from this src/ folder:
#   Rscript simulate_multiplet_overlap_benchmark.R

lorentzian <- function(x, x0, height, hwhm) height / (1 + ((x - x0) / hwhm)^2)

simulate_multiplet_benchmark <- function(
    hz_sep_values   = c(15, 20, 30),
    # Densened from the original c(3, 10, 30) to a published-benchmark-style
    # S/N grid (2.5 up to 50 in steps of 5, plus 2.5 itself)
    # so the new heatmap figures (plot_multiplet_benchmark_heatmap.py) have
    # enough S/N rows to look like that figure rather than a 3-row stub.
    snr_values      = c(2.5, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50),
    n_draws         = 10,
    field_mhz       = 600.17,
    center_ppm      = 7.5,
    j_hz            = 7,
    lw_hz_fwhm      = 1,
    total_intensity = 4e5,      # summed AREA across one multiplet's lines
    n_points        = 1024,
    window_ppm      = 1.0,
    seed_base       = 20261019,
    out_dir         = "../example_data/multiplet_benchmark"
) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  hwhm_ppm <- (lw_hz_fwhm / 2) / field_mhz   # HWHM, the fit's own width convention
  j_ppm    <- j_hz / field_mhz

  # 1:2:1 triplet, 1:3:3:1 quartet -- textbook first-order multiplet
  # intensity ratios for coupling to 2 vs. 3 equivalent protons.
  triplet_wt <- c(1, 2, 1) / 4
  quartet_wt <- c(1, 3, 3, 1) / 8

  ppm_grid <- seq(center_ppm - window_ppm / 2, center_ppm + window_ppm / 2, length.out = n_points)

  ground_truth <- list()
  spectra_long <- list()
  draw_id <- 0L

  for (hz_sep in hz_sep_values) {
    sep_ppm <- hz_sep / field_mhz
    triplet_center <- center_ppm - sep_ppm / 2
    quartet_center <- center_ppm + sep_ppm / 2
    triplet_ppm <- triplet_center + (c(-1, 0, 1)) * j_ppm
    quartet_ppm <- quartet_center + (c(-1.5, -0.5, 0.5, 1.5)) * j_ppm

    triplet_area   <- total_intensity * triplet_wt
    triplet_height <- triplet_area / (pi * hwhm_ppm)
    quartet_area   <- total_intensity * quartet_wt
    quartet_height <- quartet_area / (pi * hwhm_ppm)

    # S/N reference line: the triplet's downfield (higher ppm) satellite,
    # i.e. its last (highest-ppm) line -- one of the two smaller "1"
    # lines in the 1:2:1 pattern.
    satellite_height <- triplet_height[3]

    for (snr in snr_values) {
      noise_sd <- satellite_height / snr

      for (draw in seq_len(n_draws)) {
        draw_id <- draw_id + 1L
        set.seed(seed_base + draw_id)

        y <- rep(0, n_points)
        for (i in 1:3) y <- y + lorentzian(ppm_grid, triplet_ppm[i], triplet_height[i], hwhm_ppm)
        for (i in 1:4) y <- y + lorentzian(ppm_grid, quartet_ppm[i], quartet_height[i], hwhm_ppm)
        y <- y + rnorm(n_points, mean = 0, sd = noise_sd)

        spectra_long[[length(spectra_long) + 1]] <- data.frame(
          draw_id = draw_id, ppm = ppm_grid, intensity = y
        )

        ground_truth[[length(ground_truth) + 1]] <- rbind(
          data.frame(draw_id = draw_id, hz_sep = hz_sep, snr_level = snr,
                     field_mhz = field_mhz, multiplet_type = "triplet",
                     line_index = 0:2, true_ppm = triplet_ppm,
                     true_J_hz = j_hz, true_lw_hz = lw_hz_fwhm,
                     true_height = triplet_height, true_area = triplet_area),
          data.frame(draw_id = draw_id, hz_sep = hz_sep, snr_level = snr,
                     field_mhz = field_mhz, multiplet_type = "quartet",
                     line_index = 0:3, true_ppm = quartet_ppm,
                     true_J_hz = j_hz, true_lw_hz = lw_hz_fwhm,
                     true_height = quartet_height, true_area = quartet_area)
        )
      }
    }
  }

  ground_truth <- do.call(rbind, ground_truth)
  spectra_long <- do.call(rbind, spectra_long)

  write.csv(ground_truth, file.path(out_dir, "ground_truth.csv"), row.names = FALSE)
  write.csv(spectra_long, file.path(out_dir, "spectra_long.csv"), row.names = FALSE)

  message("Wrote ", nrow(ground_truth), " ground-truth line(s) across ", draw_id, " draw(s) to ",
          file.path(out_dir, "ground_truth.csv"))
  message("Wrote ", nrow(spectra_long), " spectral point(s) to ",
          file.path(out_dir, "spectra_long.csv"))

  invisible(list(ground_truth = ground_truth, spectra_long = spectra_long))
}

if (sys.nframe() == 0) {
  simulate_multiplet_benchmark()
}
