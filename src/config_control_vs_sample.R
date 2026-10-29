# config_control_vs_sample.R
#
# The one file to edit for a new dataset. run_peak_picking_control_vs_sample.R
# and run_stats_and_plots_control_vs_sample.R both source this file first and
# use only what is defined here. For a new study, replace the values
# below. You should not need to touch either runner script itself.
#
# See README.md, "Input data layout", for the folder structure
# importSpectraBruker() expects and for the junction trick used to stage
# real Bruker data into that structure without copying it.

# ---------------------------------------------------------------------
# 1. WHERE YOUR STAGED DATA LIVES
# One real folder per sample, directly under this path, each containing a
# numbered experiment folder (the Bruker "expno") with its own
# pdata/<procno> inside.
# ---------------------------------------------------------------------
data_dir <- "/path/to/your/staged/bruker/samples"

# ---------------------------------------------------------------------
# 2. WHICH SAMPLE GOES IN WHICH GROUP
# The name on the left must exactly match that sample's folder name under
# data_dir. Use any two group labels you like on the right (e.g.
# "untreated"/"treated", "wild_type"/"mutant") as long as there are
# exactly two distinct labels across the whole vector.
# ---------------------------------------------------------------------
sample_group_map <- c(
  "sample_control_01" = "control",
  "sample_control_02" = "control",
  "sample_control_03" = "control",
  "sample_treatment_01" = "treatment",
  "sample_treatment_02" = "treatment",
  "sample_treatment_03" = "treatment"
)

group_labels <- unique(sample_group_map)
if (length(group_labels) != 2) {
  stop("sample_group_map must use exactly two distinct group labels, found ",
       length(group_labels), ": ", paste(group_labels, collapse = ", "))
}
control_samples <- names(sample_group_map)[sample_group_map == group_labels[1]]
# test_samples holds the second group's sample names, whatever label you gave
# that group above (named test_samples rather than sample_samples, since
# "sample" already means one NMR specimen everywhere else in this pipeline).
test_samples <- names(sample_group_map)[sample_group_map == group_labels[2]]
all_samples <- c(control_samples, test_samples)

# ---------------------------------------------------------------------
# 3. DETECTION SENSITIVITY
# What these actually change is documented where they're used in
# baseline_peak_picking.R. The values below are the ones this pipeline's
# own tuning run (a dense 74-line synthetic benchmark plus a real sample)
# settled on; a cleaner, less crowded spectrum can usually afford to
# raise prominence_mult/min_height_mult back up for fewer, more
# confident peaks.
# ---------------------------------------------------------------------
prominence_mult <- 3
min_height_mult <- 1.5
deriv_prominence_mult <- 1.5
deriv_min_height_mult <- 0.75
deriv_dedup_tol_ppm <- 0.004
use_lorentzian_refinement <- FALSE
alignment_tol_ppm <- 0.008

# NOTE on the water/solvent peak's effect on the baseline fit: on real
# data, fit_airpls_baseline() can form a broad hump tracking the water
# peak's shoulders (visible in a full-spectrum plot of raw vs. fitted
# baseline). A fix was tried here -- excluding a fixed ppm range around
# the water peak from the fit itself -- and reverted: on real HAIRD/
# SWEPF standards, the water peak's actual tail extent varies by sample
# and concentration (confirmed: a fixed c(3.9, 5.3) window is too narrow
# on some samples, creating a worse spike right at its edge once that
# range is excluded from the fit; a data-driven attempt at detecting the
# tail's real extent came out unstable, expanding to swallow large real
# peak regions on some samples). fit_airpls_baseline() does accept an
# `exclude` argument (see its own comment) if you want to experiment,
# but nothing here passes it by default -- the hump is a known, visible,
# bounded cosmetic artifact; the two things tried to fix it were worse.

# ---------------------------------------------------------------------
# 4. WHERE OUTPUT GOES
# ---------------------------------------------------------------------
out_base <- "../output/control_vs_sample"
