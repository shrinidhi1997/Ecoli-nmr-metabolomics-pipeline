# read_bruker_1d_direct.R
#
# A minimal, direct reader for a real Bruker 1D processed spectrum
# (pdata/<procno>/1r + the procs text file next to it), used ONLY because
# this sandbox has no network access to CRAN/Bioconductor and so can't
# install ASICS here (confirmed: install.packages() can't even reach
# cloud.r-project.org from this environment). Everywhere else in this
# project, ASICS::importSpectraBruker() is the real import step and
# should stay that way on a machine that has it installed --
# run_peak_picking_control_vs_sample.R is unchanged.
#
# What this reads is 100% real instrument data -- unlike
# simulate_alanine_standard_dilution_1d.R (which fabricates the spectrum
# itself), this bypasses only the Bruker-binary PARSING step, the same
# one narrow thing that script's own header describes bypassing, just
# for a different reason (no package vs. no real file). Every later
# step (airPLS baseline, peak picking, merging, matching, scoring) is
# the same real, unmodified pipeline code either way.
#
# Format (standard Bruker TopSpin processed-data layout, confirmed against
# this sample's own procs/acqus below, not assumed):
#   pdata/<procno>/1r   -- SI real-valued points, stored as 4-byte
#                          integers (DTYPP=0), byte order given by
#                          BYTORDP (0 = little-endian), each point's true
#                          intensity = stored_int * 2^NC_proc.
#   pdata/<procno>/procs -- text file of "##$KEY= value" lines. OFFSET is
#                          the chemical shift (ppm) of the FIRST point;
#                          the ppm axis then decreases linearly across SI
#                          points, spanning SW_p (Hz) / SF (MHz) ppm total.
#
# Verified on this project's own real data (HAIRD-30/60/120-1H, Oct 2026):
# BYTORDP=0, DTYPP=0, SI=65536, 1r file size = 262144 bytes = 65536 x 4,
# SW_p/SF = 12019.2308/600.1699 = 20.0263 ppm, matching acqus's own SW
# (20.0263 ppm) exactly -- so this formula is reading the real sweep
# width this instrument actually recorded, not a guessed constant.

.read_jcampdx_value <- function(lines, key) {
  pat <- paste0("^##\\$", key, "= ?(.*)$")
  hit <- grep(pat, lines, value = TRUE)
  if (length(hit) == 0) stop("Key ", key, " not found in procs/acqus file.")
  val <- sub(pat, "\\1", hit[1])
  val <- gsub("^<|>$", "", val)
  val
}

read_bruker_1r <- function(sample_dir, expno = "1", procno = "1") {
  pdata_dir <- file.path(sample_dir, expno, "pdata", procno)
  procs_path <- file.path(pdata_dir, "procs")
  onerr_path <- file.path(pdata_dir, "1r")
  if (!file.exists(procs_path)) stop("Missing procs file: ", procs_path)
  if (!file.exists(onerr_path)) stop("Missing 1r file: ", onerr_path)

  procs_lines <- readLines(procs_path, warn = FALSE)
  si        <- as.integer(.read_jcampdx_value(procs_lines, "SI"))
  sw_p      <- as.numeric(.read_jcampdx_value(procs_lines, "SW_p"))
  sf        <- as.numeric(.read_jcampdx_value(procs_lines, "SF"))
  offset    <- as.numeric(.read_jcampdx_value(procs_lines, "OFFSET"))
  nc_proc   <- as.integer(.read_jcampdx_value(procs_lines, "NC_proc"))
  bytordp   <- as.integer(.read_jcampdx_value(procs_lines, "BYTORDP"))
  dtypp     <- as.integer(.read_jcampdx_value(procs_lines, "DTYPP"))

  if (dtypp != 0) stop("DTYPP=", dtypp, " (not 4-byte integer) -- this reader only ",
                        "handles the standard DTYPP=0 int32 format; extend it before ",
                        "trusting this on data processed differently.")
  endian <- if (bytordp == 0) "little" else "big"

  n_bytes <- file.info(onerr_path)$size
  expected_n <- si * 4L
  if (n_bytes != expected_n) {
    stop("1r file is ", n_bytes, " bytes, expected ", expected_n,
         " (SI=", si, " x 4 bytes) -- SI in procs doesn't match this file; ",
         "double check procno/expno point at the same processed dataset.")
  }

  con <- file(onerr_path, "rb")
  raw_int <- readBin(con, what = "integer", n = si, size = 4L, endian = endian)
  close(con)

  intensity <- raw_int * (2^nc_proc)
  ppm <- offset - (seq_len(si) - 1) * (sw_p / sf) / si

  # Bruker stores processed spectra highest-ppm-first (ppm descending);
  # every other function in this pipeline (baseline_peak_picking.R,
  # peak_picking_functions.R) assumes ascending ppm, the same convention
  # simulate_alanine_standard_dilution_1d.R's own synthetic seq() grid
  # already uses. find_shape_based_peaks() re-sorts internally, but
  # fit_airpls_baseline() does not (its control-point search assumes a
  # pre-sorted axis) -- sort once here, at the source, instead of
  # relying on every caller downstream to remember to.
  ord <- order(ppm)
  ppm <- ppm[ord]
  intensity <- intensity[ord]

  list(ppm = ppm, intensity = intensity,
       meta = list(si = si, sw_p = sw_p, sf = sf, offset = offset,
                   nc_proc = nc_proc, sample_dir = sample_dir))
}
