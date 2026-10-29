# run_pipeline_haird_real_dilution.R
#
# Runs the real pipeline (pipeline_standard_dilution_core.R -- see that
# file for the full rationale: direct Bruker reader standing in for
# ASICS, shape-detector-only feature table, water/DSS-excluded baseline
# fit) on three real Bruker 1H spectra: a standard dilution series
# labelled HAIRD-30-1H, HAIRD-60-1H, HAIRD-120-1H, laid out as
# <data_root>/<name>/1/pdata/1/... (set data_root below to wherever
# your copy of these spectra lives).
#
# This used to duplicate pipeline_standard_dilution_core.R's logic
# in-line; pulled out into that shared function (same change
# run_pipeline_swepf_real_dilution.R already used) so a fix made once
# (like the baseline water-exclusion fix) applies to every dataset
# automatically instead of needing to be repeated per copy.
#
# Run from this src/ folder:
#   Rscript run_pipeline_haird_real_dilution.R

source("pipeline_standard_dilution_core.R")

data_root <- "/path/to/your/haird_real_data/extracted"
sample_ids <- c("HAIRD-30-1H", "HAIRD-60-1H", "HAIRD-120-1H")
# label_value: the sample folders are named by uL of stock solution added
# (30/60/120), NOT the final concentration. The actual final
# concentrations are 4, 8, and 16 mM respectively (same
# linear rescaling, x' = x / 7.5, so R^2 on any linearity fit is
# unaffected); see label_units_confirmed below.
label_value <- c("HAIRD-30-1H" = 4, "HAIRD-60-1H" = 8, "HAIRD-120-1H" = 16)
out_dir <- "../example_results/haird_real_dilution"

run_standard_dilution_pipeline("haird", data_root, sample_ids, label_value, out_dir,
                                label_units_confirmed = TRUE)

message("\nDone. Next: nmrmatch-match-feature-table (or python match_r_feature_table.py) ",
        "against the reference library, then the plotting script -- see README.md ",
        "for the exact commands.")
