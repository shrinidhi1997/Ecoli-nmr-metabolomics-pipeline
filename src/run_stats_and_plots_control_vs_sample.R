# run_stats_and_plots_control_vs_sample.R
#
# Reads match_r_feature_table.py's output (compound_matches_control_vs_sample.csv,
# one ranked row per sample per candidate compound, with matched_intensity_mean
# -- see utils.py's score_compound() docstring for exactly what that value is
# and isn't) and runs the statistics/figures for the two-group comparison:
# per-compound Welch t-test + Wilcoxon rank-sum + BH-FDR, fold-change, PCA, a
# sample correlation heatmap, a volcano plot, a group-averaged heatmap of the
# most differential compounds, individual-replicate boxplots for the top
# compounds, a confidence distribution histogram, a detected-compound-count
# bar per sample, and a stacked composition bar chart across all samples
# (Sections F, G, H).
#
# score_compound() (in utils.py) rewards both the FRACTION of a compound's
# reference peaks that matched and the ABSOLUTE number of matched peaks and
# how prominent those matched peaks are in the sample's own spectrum. Without
# that second part, a compound with only 2 reference peaks, both matched,
# could automatically outrank a compound matched on far more peaks, since a
# ratio alone can't tell "2 of 2" apart from "27 of 28". If your confidence
# numbers for thin (1-3-reference-peak) compounds run lower than you'd
# naively expect, that's this scoring choice working as intended, not a bug
# (see utils.py's score_compound() docstring for the full rationale).
#
# All dataset-specific settings, including control_samples/test_samples, live
# in config_control_vs_sample.R, which this script sources below. Edit that
# file for a new dataset; this script should not need to change.
#
# NOTE ON STATISTICAL POWER: with n replicates per group, an exact two-sided
# Wilcoxon rank-sum test can reach at best p = 2/choose(2n, n), even for the
# most extreme possible split. At n=3 that floor is 2/20 = 0.1, ABOVE the
# conventional p<0.05 threshold, so no compound can clear raw Wilcoxon
# significance at all with only 3 replicates/group, let alone survive BH-FDR
# correction across the dozens to low hundreds of compounds this pipeline
# typically matches. The Welch t-test doesn't have that exact-test floor and
# can still land under p<0.05, but with few replicates it is far more
# sensitive to one noisy sample. Treat any single-compound p-value from a
# small study as a lead to re-check with more replicates, not a result to
# report on its own.
#
# Run from the project's src/ folder, after match_r_feature_table.py:
#   source("run_stats_and_plots_control_vs_sample.R")

if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(ggplot2)

if (!requireNamespace("nmrpick", quietly = TRUE)) {
  stop("Package 'nmrpick' is not installed. From an R session: ",
       'devtools::install_github("shrinidhi1997/Ecoli-nmr-metabolomics-pipeline", subdir = "nmrpick")')
}
library(nmrpick)

source("config_control_vs_sample.R")
matches_file <- file.path(out_base, "compound_matches_control_vs_sample.csv")
out_dir <- file.path(out_base, "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

while (!is.null(dev.list())) dev.off()

# png() itself (opening the file) is wrapped in tryCatch too, not just the
# drawing call, so a locked/read-only/still-syncing target file skips just
# that one figure instead of killing the whole script.
safe_plot <- function(path, plot_fn, width = 900, height = 800) {
  opened <- tryCatch({
    png(path, width = width, height = height, res = 150)
    TRUE
  }, error = function(e) {
    message("Could not open ", basename(path), " for writing: ", conditionMessage(e),
            " -- likely causes: the file is open in another program, the folder ",
            "is read-only, or a cloud-sync tool hasn't finished syncing it. ",
            "Skipping this one figure; the rest of the run continues.")
    FALSE
  })
  if (!opened) return(invisible(FALSE))
  ok <- tryCatch({
    plot_fn()
    TRUE
  }, error = function(e) {
    message("Could not draw ", basename(path), ": ", conditionMessage(e))
    FALSE
  })
  dev.off()
  invisible(ok)
}

if (!file.exists(matches_file)) {
  stop("Can't find ", matches_file, " -- run match_r_feature_table.py first ",
       "(see its own header for the exact command).")
}
matches <- read.csv(matches_file, stringsAsFactors = FALSE)

missing_samples <- setdiff(all_samples, unique(matches$sample_id))
if (length(missing_samples) > 0) {
  message("Warning: these samples aren't in the match file at all (typo in ",
          "the sample name, or that sample had zero matches): ",
          paste(missing_samples, collapse = ", "))
}
if (!"matched_intensity_mean" %in% colnames(matches)) {
  stop("No matched_intensity_mean column in ", matches_file, " -- re-run ",
       "match_r_feature_table.py with the current utils.py/peak_matcher.py ",
       "(an older match file won't have this column).")
}
matches <- matches[matches$sample_id %in% all_samples, ]
if (nrow(matches) == 0) {
  stop("None of the ", length(all_samples), " sample names in control_samples/",
       "test_samples (config_control_vs_sample.R) appear in ", matches_file,
       " at all. Check the sample names match exactly what ",
       "match_r_feature_table.py was run against.")
}

# -------------------------------------------------------------------
# Build the compound x sample matrix via nmrpick::build_compound_matrix()
# -- one row per compound (keyed on compound_id, the real unique key --
# never display_name/compound_name for lookups), one column per sample,
# value = matched_intensity_mean, NA where that sample had no match.
# compound_keys$display_name disambiguates compound names shared by more
# than one compound_id (a confirmed real case earlier in this project:
# "taurine" = bmse000120 and bmse000863), appending "[compound_id]" ONLY
# to the colliding names so the common case stays readable as-is. This
# used to be inline code in this script; it is now one tested
# package function.
# -------------------------------------------------------------------
built <- build_compound_matrix(matches, all_samples)
wide <- built$wide
compound_keys <- built$compound_keys

message("Loaded ", nrow(matches), " compound-match rows across ", length(all_samples),
        " sample(s), ", nrow(compound_keys), " unique compound(s).")

# ============================================================
# SECTION A -- per-compound Welch t-test + Wilcoxon + BH-FDR + fold-change,
# via nmrpick::compound_stats_table(). Only compounds with >= 2 non-NA
# samples per group get a real p-value -- same floor used every previous
# round of this project (a "p-value" from a single sample isn't one).
# See ?nmrpick::wilcoxon_exact_floor for why the exact Wilcoxon test
# structurally can't clear p<0.05 at small per-group sample sizes; the
# message below reports that floor for this exact dataset's group sizes.
# ============================================================
results <- compound_stats_table(wide, compound_keys, control_samples, test_samples,
                                 min_n_per_group = 2)

out_csv <- file.path(out_base, "stats_control_vs_sample.csv")
write.csv(results, out_csv, row.names = FALSE)
message("Saved per-compound stats to ", normalizePath(out_csv))

# ============================================================
# SECTION B -- volcano plot (Welch t-test p; with a small number of
# replicates per group, the Wilcoxon exact floor from the header note above
# makes most points pile up at a handful of discrete p-values, which is a
# worse y-axis choice than the t-test p here. wilcox_BH_q_value is still in
# the CSV for anyone who wants the distribution-free version).
# ============================================================
top_n_label <- 20
results$significant <- !is.na(results$t_p_value) & results$t_p_value < 0.05 & abs(results$log2fc) > 1
results$neg_log10_p <- -log10(pmax(results$t_p_value, 1e-300))
label_set <- head(results[order(-abs(results$log2fc)), ], top_n_label)

have_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
if (!have_ggrepel) {
  have_ggrepel <- tryCatch({install.packages("ggrepel"); requireNamespace("ggrepel", quietly = TRUE)},
                            error = function(e) FALSE)
}
p_volcano <- ggplot(results, aes(x = log2fc, y = neg_log10_p, color = significant)) +
  geom_point(alpha = 0.75, size = 2) +
  scale_color_manual(values = c("TRUE" = "#d9534f", "FALSE" = "#9aa5b1"), guide = "none") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  labs(title = paste0("Volcano plot, sample vs control (Welch t-test, raw p, ",
                       length(control_samples), " vs ", length(test_samples), ")"),
       subtitle = paste0("Dashed line = uncorrected p = 0.05 (reference only -- with n=",
                          length(control_samples), "/group the exact Wilcoxon floor is ~",
                          wilcoxon_exact_floor(length(control_samples), length(test_samples)),
                          "; see stats_control_vs_sample.csv's t_BH_q_value and ",
                          "wilcox_BH_q_value columns, which are the ones to actually trust)"),
       x = "log2 fold change (sample / control)", y = "-log10(raw t-test p-value)") +
  theme_minimal(base_size = 11)
if (have_ggrepel) {
  p_volcano <- p_volcano + ggrepel::geom_text_repel(data = label_set, aes(label = display_name),
                                                      size = 2.6, max.overlaps = 30, color = "black")
} else {
  p_volcano <- p_volcano + geom_text(data = label_set, aes(label = display_name),
                                       size = 2.6, vjust = -0.7, check_overlap = TRUE, color = "black")
}
ggsave(file.path(out_dir, "volcano_control_vs_sample.png"), p_volcano, width = 9, height = 7, dpi = 150)
message("Section B done: volcano_control_vs_sample.png written to ", normalizePath(out_dir))

# ============================================================
# SECTION C -- group-averaged heatmap of the top differential compounds.
# ============================================================
top_n_heat <- min(40, nrow(results))
top_heat <- head(results, top_n_heat)
heat_long <- data.frame(
  display_name = rep(factor(top_heat$display_name, levels = rev(top_heat$display_name)), 2),
  group = rep(c("control", "sample"), each = nrow(top_heat)),
  mean_intensity = c(top_heat$mean_control, top_heat$mean_sample)
)
p_group_heat <- ggplot(heat_long, aes(x = group, y = display_name, fill = mean_intensity)) +
  geom_tile() +
  scale_fill_gradient(low = "#f7f7f7", high = "#d9534f") +
  theme_minimal(base_size = 8) +
  theme(axis.text.y = element_text(size = 6)) +
  labs(title = paste0("Group-mean matched intensity, top ", top_n_heat, " most-changed compounds"),
       x = NULL, y = NULL, fill = "mean\nintensity")
ggsave(file.path(out_dir, "heatmap_group_averages.png"), p_group_heat,
       width = 6, height = max(5, top_n_heat * 0.18), dpi = 150, limitsize = FALSE)
message("Section C done: heatmap_group_averages.png written to ", normalizePath(out_dir))

# ============================================================
# SECTION D -- individual-replicate boxplots, top 12 most-changed compounds.
# ============================================================
top_n_box <- min(12, nrow(results))
top_box <- head(results, top_n_box)
box_long <- do.call(rbind, lapply(seq_len(nrow(top_box)), function(i) {
  cid <- top_box$compound_id[i]
  dname <- top_box$display_name[i]
  rbind(
    data.frame(display_name = dname, group = "control",
               value = as.numeric(wide[cid, control_samples]),
               sample_id = control_samples),
    data.frame(display_name = dname, group = "sample",
               value = as.numeric(wide[cid, test_samples]),
               sample_id = test_samples)
  )
}))
box_long$display_name <- factor(box_long$display_name, levels = top_box$display_name)
box_long$group <- factor(box_long$group, levels = c("control", "sample"))
p_box <- ggplot(box_long, aes(x = display_name, y = value, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, position = position_dodge(width = 0.75)) +
  geom_jitter(aes(color = group), position = position_dodge(width = 0.75), size = 1.6, na.rm = TRUE) +
  coord_flip() +
  labs(title = paste0("Individual replicates, top ", top_n_box, " most-changed compounds"),
       subtitle = paste0("Only ", length(control_samples), " vs ", length(test_samples),
                          " points per group -- read spread, not just the box, at this sample size"),
       x = NULL, y = "matched_intensity_mean") +
  theme_minimal(base_size = 11)
ggsave(file.path(out_dir, "top12_compound_replicates_boxplot.png"), p_box, width = 8, height = 7, dpi = 150)
message("Section D done: top12_compound_replicates_boxplot.png written to ", normalizePath(out_dir))

# ============================================================
# SECTION E -- sample correlation heatmap + PCA, both restricted to
# compounds with a value in ALL samples (prcomp()/cor() need complete
# cases; with only a few reps/group, imputing the gaps would do more harm
# than restricting to compounds every sample actually shares).
# ============================================================
complete_mask <- stats::complete.cases(wide[, all_samples])
n_complete <- sum(complete_mask)
message(n_complete, " of ", nrow(compound_keys),
        " compound(s) have a value in all ", length(all_samples),
        " samples; PCA and the sample correlation heatmap use only those ",
        "(everything else in this script uses all available data per compound).")

if (n_complete >= 3) {
  complete_wide <- wide[complete_mask, all_samples, drop = FALSE]

  cor_mat <- cor(complete_wide, method = "pearson")
  cor_long <- as.data.frame(as.table(cor_mat))
  names(cor_long) <- c("sample_a", "sample_b", "correlation")
  cor_long$sample_a <- factor(cor_long$sample_a, levels = all_samples)
  cor_long$sample_b <- factor(cor_long$sample_b, levels = all_samples)
  p_cor <- ggplot(cor_long, aes(x = sample_a, y = sample_b, fill = correlation)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", correlation)), size = 3) +
    scale_fill_gradient2(low = "#3b6ea5", mid = "white", high = "#d9534f", midpoint = 0.9,
                          limits = c(min(cor_long$correlation, na.rm = TRUE), 1)) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste0("Sample-sample correlation (", n_complete, " compounds present in all ",
                         length(all_samples), " samples)"),
         x = NULL, y = NULL)
  ggsave(file.path(out_dir, "sample_correlation_heatmap.png"), p_cor, width = 7, height = 6, dpi = 150)
  message("Section E (correlation) done: sample_correlation_heatmap.png written to ", normalizePath(out_dir))

  condition <- ifelse(all_samples %in% control_samples, "control", "sample")
  pca_input <- t(log2(complete_wide + 1e-6))
  resPCA <- tryCatch(prcomp(pca_input, scale. = TRUE), error = function(e) {
    message("PCA failed: ", conditionMessage(e)); NULL
  })
  if (!is.null(resPCA)) {
    var_explained <- round(100 * (resPCA$sdev^2 / sum(resPCA$sdev^2))[1:2], 1)
    pca_df <- data.frame(sample_id = all_samples, condition = condition,
                          PC1 = resPCA$x[, 1], PC2 = resPCA$x[, 2])
    p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, label = sample_id)) +
      geom_point(size = 3) +
      labs(title = paste0("PCA (", n_complete, " compounds present in all ",
                           length(all_samples), " samples)"),
           x = paste0("PC1 (", var_explained[1], "%)"),
           y = paste0("PC2 (", var_explained[2], "%)")) +
      theme_minimal(base_size = 11)
    if (have_ggrepel) {
      p_pca <- p_pca + ggrepel::geom_text_repel(size = 3, show.legend = FALSE)
    } else {
      p_pca <- p_pca + geom_text(size = 3, vjust = -0.9, check_overlap = TRUE, show.legend = FALSE)
    }
    ggsave(file.path(out_dir, "pca_samples.png"), p_pca, width = 7, height = 6, dpi = 150)
    message("Section E (PCA) done: pca_samples.png written to ", normalizePath(out_dir))
  }
} else {
  message("Skipping PCA and the correlation heatmap: fewer than 3 compounds are ",
          "present in all ", length(all_samples), " samples, not enough to say ",
          "anything meaningful.")
}

# ============================================================
# SECTION F -- confidence distribution across every kept candidate match
# (not just the compounds that made it into the per-compound stats table
# above). A quick way to see at a glance whether this run's matches lean
# toward strong evidence or toward the thin, few-peak end -- especially
# relevant now that confidence also reflects evidence_score/peak_prominence,
# not just coverage/shift_score (see this file's header).
# ============================================================
matches$group <- ifelse(matches$sample_id %in% control_samples, "control", "sample")
p_conf <- ggplot(matches, aes(x = confidence, fill = group)) +
  geom_histogram(binwidth = 0.02, alpha = 0.6, position = "identity") +
  scale_fill_manual(values = c("control" = "#3b6ea5", "sample" = "#d9534f")) +
  theme_minimal(base_size = 11) +
  labs(title = "Confidence distribution, every kept candidate match",
       subtitle = paste0(nrow(matches), " match rows across ", length(all_samples), " samples"),
       x = "confidence", y = "count of candidate matches")
ggsave(file.path(out_dir, "confidence_distribution.png"), p_conf, width = 7, height = 5, dpi = 150)
message("Section F done: confidence_distribution.png written to ", normalizePath(out_dir))

# ============================================================
# SECTION G -- how many candidate compounds got kept per sample. A quick
# sanity check, not a result in its own right -- a sample far out of line
# with the rest here usually means a run problem (bad tolerance, a sample
# that didn't import correctly), not a biological finding.
# ============================================================
counts_df <- as.data.frame(table(matches$sample_id), stringsAsFactors = FALSE)
names(counts_df) <- c("sample_id", "n_compounds")
counts_df$sample_id <- factor(counts_df$sample_id, levels = rev(all_samples))
counts_df$group <- ifelse(as.character(counts_df$sample_id) %in% control_samples, "control", "sample")
p_counts <- ggplot(counts_df, aes(x = sample_id, y = n_compounds, fill = group)) +
  geom_col(width = 0.6) +
  coord_flip() +
  scale_fill_manual(values = c("control" = "#3b6ea5", "sample" = "#d9534f")) +
  theme_minimal(base_size = 11) +
  labs(title = "Candidate compounds kept per sample", x = NULL, y = "n_compounds")
ggsave(file.path(out_dir, "detected_compound_counts.png"), p_counts, width = 7, height = 4.5, dpi = 150)
message("Section G done: detected_compound_counts.png written to ", normalizePath(out_dir))

# ============================================================
# SECTION H -- stacked composition bar chart. For the top differential
# compounds (same list as Sections C/D), shows each one's SHARE of that
# group's combined matched_intensity_mean in each sample, stacked to 100%
# per sample -- a relative-composition view, not an absolute-amount one, so
# it answers "does the MIX of these compounds shift between groups" rather
# than "how much of each is there" (the boxplots in Section D already
# answer the absolute-amount question).
# ============================================================
top_n_stack <- min(10, nrow(results))
top_stack <- head(results, top_n_stack)
stack_long <- do.call(rbind, lapply(seq_len(nrow(top_stack)), function(i) {
  cid <- top_stack$compound_id[i]
  vals <- as.numeric(wide[cid, all_samples])
  vals[is.na(vals)] <- 0
  data.frame(sample_id = all_samples, display_name = top_stack$display_name[i], value = vals,
             stringsAsFactors = FALSE)
}))
sample_totals <- tapply(stack_long$value, stack_long$sample_id, sum)
stack_long$fraction <- stack_long$value / sample_totals[stack_long$sample_id]
stack_long$sample_id <- factor(stack_long$sample_id, levels = all_samples)
stack_long$display_name <- factor(stack_long$display_name, levels = rev(top_stack$display_name))
p_stack <- ggplot(stack_long, aes(x = sample_id, y = fraction, fill = display_name)) +
  geom_col(width = 0.7) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = paste0("Relative composition, top ", top_n_stack, " most-changed compounds"),
       subtitle = "Each bar sums to 100% of just these compounds, not the whole spectrum",
       x = NULL, y = "share of these compounds' combined intensity", fill = "compound")
ggsave(file.path(out_dir, "stacked_composition_top_compounds.png"), p_stack, width = 8, height = 6, dpi = 150)
message("Section H done: stacked_composition_top_compounds.png written to ", normalizePath(out_dir))

message("\nAll done. Figures are in ", normalizePath(out_dir),
        "; CSVs are in ", normalizePath(out_base))
