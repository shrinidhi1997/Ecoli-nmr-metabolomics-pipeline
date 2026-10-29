# run_pipeline_swepf_real_dilution.R
#
# Same real pipeline as run_pipeline_haird_real_dilution.R (see that
# file's header, or pipeline_standard_dilution_core.R, for the full
# rationale), run on a second real standard: a 3-point
# SWEPF-30/60/120-1H dilution series, laid out as
# <data_root>/<sample>/1/pdata/1/... (set data_root below)
# (the zip's own top-level folder is named SWEPF-60-1H and contains all
# three sample subfolders -- same "outer folder named after one of the
# concentrations" layout the HAIRD zip used).
#
# Run from this src/ folder:
#   Rscript run_pipeline_swepf_real_dilution.R

source("pipeline_standard_dilution_core.R")

data_root <- "/path/to/your/swepf_real_data/extracted/SWEPF-60-1H"
sample_ids <- c("SWEPF-30-1H", "SWEPF-60-1H", "SWEPF-120-1H")
# label_value: sample folders are named by uL of stock added (30/60/120),
# not final concentration. Final concentrations are 4, 8,
# and 16 mM respectively (linear rescaling, R^2 unaffected).
label_value <- c("SWEPF-30-1H" = 4, "SWEPF-60-1H" = 8, "SWEPF-120-1H" = 16)
out_dir <- "../example_results/swepf_real_dilution"

run_standard_dilution_pipeline("swepf", data_root, sample_ids, label_value, out_dir,
                                label_units_confirmed = TRUE)

message("\nDone. Next: nmrmatch-match-feature-table (or python match_r_feature_table.py) ",
        "against the reference library, then the plotting script -- see README.md ",
        "for the exact commands.")
