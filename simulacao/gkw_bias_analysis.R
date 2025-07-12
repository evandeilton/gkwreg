# ---------------------------------------------------------------------------- #
# Title: Analysis Script for GKw MLE Bias Study
# Description: This script analyzes the results of the GKw MLE bias simulation
#              study, generates plots and tables for all parameters and methods.
# ---------------------------------------------------------------------------- #

# Load required packages
library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggh4x)    # For facet_nested
library(latex2exp) # For rendering LaTeX in plots
library(knitr)
library(kableExtra)

# Load helper functions
source("gkw_bias_extra_functions.R")

# ---------------------------------------------------------------------------- #
# Load Results
# ---------------------------------------------------------------------------- #

# Create output directories if they don't exist
plots_dir <- "plots"
param_plots_dir <- file.path(plots_dir, "parameter_focus")
csv_dir <- "csv"

dir_list <- c(plots_dir, param_plots_dir, csv_dir)
for (dir_path in dir_list) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
}

# Load the aggregated results
param_results <- readRDS("rds/gkw_bias_param_results.rds")
convergence_df <- readRDS("rds/gkw_bias_convergence_results.rds")
run_config <- readRDS("rds/run_config.rds")

# Extract simulation parameters from config
methods <- run_config$methods
sample_sizes <- run_config$sample_sizes
configurations <- run_config$configurations

# Verify data loaded correctly
cat("Loaded parameter results:", nrow(param_results), "rows\n")
cat("Loaded convergence results:", nrow(convergence_df), "rows\n")

# Extract unique parameters, families, and configurations
unique_params <- unique(param_results$parameter)
unique_families <- unique(param_results$family)
unique_configs <- unique(param_results$config)

cat("Parameters in results:", paste(unique_params, collapse=", "), "\n")
cat("Families in results:", paste(unique_families, collapse=", "), "\n")

# ---------------------------------------------------------------------------- #
# Generate Summary Tables
# ---------------------------------------------------------------------------- #

cat("Generating summary tables...\n")

# 1. Bias summary table
bias_summary <- create_bias_summary(param_results)
write_csv(bias_summary, file.path(csv_dir, "gkw_bias_summary.csv"))

# 2. Coverage summary table
coverage_summary <- create_coverage_summary(param_results)
write_csv(coverage_summary, file.path(csv_dir, "gkw_coverage_summary.csv"))

# 3. Convergence and timing summary
convergence_summary <- create_convergence_summary(convergence_df)
write_csv(convergence_summary, file.path(csv_dir, "gkw_convergence_summary.csv"))

# 4. Method comparison summary
method_comparison <- create_method_comparison(param_results, convergence_df)
write_csv(method_comparison, file.path(csv_dir, "gkw_method_comparison.csv"))

# 5. Method rankings
method_rankings <- rank_methods(param_results, convergence_df,
                                rel_bias_weight = 0.4,
                                coverage_weight = 0.3,
                                convergence_weight = 0.2,
                                time_weight = 0.1)
write_csv(method_rankings, file.path(csv_dir, "gkw_method_rankings.csv"))

# 6. Identify problematic configurations
problematic_configs <- identify_problematic_configs(param_results, convergence_df,
                                                    convergence_threshold = 0.8,
                                                    rel_bias_threshold = 20.0)

write_csv(problematic_configs$low_convergence,
          file.path(csv_dir, "gkw_low_convergence_configs.csv"))
write_csv(problematic_configs$high_bias,
          file.path(csv_dir, "gkw_high_bias_configs.csv"))

# 7. Relative bias comparison table
rel_bias_comparison <- param_results %>%
  group_by(parameter, family, method, sample_size) %>%
  summarize(
    true_value = mean(true_value),
    abs_bias = mean(bias, na.rm = TRUE),
    rel_bias_pct = mean(rel_bias, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(sample_size == max(sample_size)) %>%  # Focus on largest sample size
  select(parameter, family, method, true_value, abs_bias, rel_bias_pct) %>%
  mutate(
    true_value = round(true_value, 4),
    abs_bias = round(abs_bias, 4),
    rel_bias_pct = round(rel_bias_pct, 2)
  ) %>%
  arrange(parameter, family, method)

write_csv(rel_bias_comparison, file.path(csv_dir, "gkw_rel_bias_comparison.csv"))

# Generate a table comparing different bias measures across parameters and methods
bias_measures_comparison <- param_results %>%
  group_by(parameter, method) %>%
  summarize(
    mean_abs_bias = mean(abs(bias), na.rm = TRUE),
    mean_abs_rel_bias = mean(abs(rel_bias), na.rm = TRUE),
    mean_abs_bias_z = mean(abs(bias_z), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    mean_abs_bias = round(mean_abs_bias, 4),
    mean_abs_rel_bias = round(mean_abs_rel_bias, 2),
    mean_abs_bias_z = round(mean_abs_bias_z, 2)
  ) %>%
  arrange(parameter, method)

write_csv(bias_measures_comparison, file.path(csv_dir, "gkw_bias_measures_comparison.csv"))

# ---------------------------------------------------------------------------- #
# Generate Traditional Plots (Original Layout)
# ---------------------------------------------------------------------------- #

cat("Generating traditional plots...\n")

# Create a helper function to create and save plots
create_and_save_plot <- function(plot_func, param, filename_prefix, ...) {
  cat("  Creating", filename_prefix, "plot for parameter", param, "...\n")

  p <- plot_func(param_results, param, ...)

  # Save the plot
  filename <- file.path(plots_dir, paste0(filename_prefix, "_", param, ".pdf"))
  ggsave(filename, p, width = 10, height = 8)

  return(p)
}

# Plot each type of graph for each parameter
for(param in unique_params) {
  # 1. Absolute bias plots
  abs_bias_plot <- create_and_save_plot(
    plot_abs_bias,
    param,
    "abs_bias",
    title = paste("Absolute Bias for Parameter", param)
  )

  # 2. Relative bias plots
  # Set appropriate limits based on parameter
  bias_limits <- NULL
  if(param %in% c("alpha", "beta")) {
    bias_limits <- c(-50, 50)  # Adjust as needed
  } else if(param %in% c("gamma", "delta")) {
    bias_limits <- c(-75, 75)  # Adjust as needed
  } else if(param == "lambda") {
    bias_limits <- c(-60, 60)  # Adjust as needed
  }

  rel_bias_plot <- create_and_save_plot(
    plot_rel_bias,
    param,
    "rel_bias",
    limits = bias_limits,
    title = paste("Relative Bias (%) for Parameter", param)
  )

  # 3. Standardized bias plots
  bias_z_plot <- create_and_save_plot(
    plot_bias_z,
    param,
    "bias_z",
    limits = c(-2, 2),
    title = paste("Standardized Bias for Parameter", param)
  )

  # 4. Coverage plots
  coverage_plot <- create_and_save_plot(
    plot_coverage,
    param,
    "coverage",
    limits = c(0.7, 1.0),
    title = paste("Coverage Rate for Parameter", param)
  )

  # 5. Relative bias trend plots
  rel_bias_trend_plot <- analyze_rel_bias_trend(param_results, param)
  ggsave(file.path(plots_dir, paste0("rel_bias_trend_", param, ".pdf")),
         rel_bias_trend_plot, width = 8, height = 6)
}

# Generate convergence plot
cat("  Creating convergence plot...\n")
conv_plot <- plot_convergence(convergence_df)
ggsave(file.path(plots_dir, "convergence.pdf"), conv_plot, width = 10, height = 8)

# Generate computation time plot
cat("  Creating computation time plot...\n")
time_plot <- plot_computation_time(convergence_df, log_scale = TRUE)
ggsave(file.path(plots_dir, "computation_time.pdf"), time_plot, width = 10, height = 8)

# Generate consolidated plots for manuscript
# Create a plot that shows relative bias for all parameters on the same grid
cat("  Creating consolidated relative bias plot...\n")

# Filter data by parameter and create separate plots
rel_bias_plots <- lapply(unique_params, function(param) {
  plot_rel_bias(
    param_results %>% filter(parameter == param),
    param,
    title = NULL
  ) +
    theme(legend.position = "none") +
    labs(title = param)
})

# Combine all plots using patchwork
combined_rel_bias <- wrap_plots(rel_bias_plots, ncol = 2)

# Add a common legend
combined_rel_bias <- combined_rel_bias +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Save the combined plot
ggsave(file.path(plots_dir, "combined_rel_bias.pdf"),
       combined_rel_bias, width = 16, height = 12)

# Create a consolidated method comparison plot
cat("  Creating method comparison plot...\n")

# Create a long format dataset for the bias measures comparison
method_comparison_long <- bias_measures_comparison %>%
  pivot_longer(
    cols = c(mean_abs_bias, mean_abs_rel_bias, mean_abs_bias_z),
    names_to = "bias_measure",
    values_to = "value"
  ) %>%
  mutate(
    bias_measure = factor(
      bias_measure,
      levels = c("mean_abs_bias", "mean_abs_rel_bias", "mean_abs_bias_z"),
      labels = c("Absolute Bias", "Relative Bias (%)", "Standardized Bias")
    )
  )

# Create the method comparison plot
method_comp_plot <- ggplot(method_comparison_long,
                           aes(x = method, y = value, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(bias_measure ~ parameter, scales = "free_y") +
  labs(
    title = "Comparison of Bias Measures by Method and Parameter",
    y = "Value",
    x = "Method",
    fill = "Method"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    strip.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 12)
  )

# Save the method comparison plot
ggsave(file.path(plots_dir, "method_comparison.pdf"),
       method_comp_plot, width = 12, height = 8)

# Create a heatmap of relative bias by parameter, family, and method
cat("  Creating relative bias heatmap...\n")

# Create a summary for the heatmap
rel_bias_heatmap_data <- param_results %>%
  group_by(parameter, family, method) %>%
  summarize(
    mean_abs_rel_bias = mean(abs(rel_bias), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Cap extreme values for better visualization
  mutate(
    mean_abs_rel_bias_capped = pmin(mean_abs_rel_bias, 50)
  )

# Create the heatmap
rel_bias_heatmap <- ggplot(rel_bias_heatmap_data,
                           aes(x = method, y = family, fill = mean_abs_rel_bias_capped)) +
  geom_tile() +
  facet_wrap(~ parameter, ncol = 3) +
  scale_fill_viridis_c(option = "plasma", name = "Mean Absolute\nRelative Bias (%)") +
  labs(
    title = "Mean Absolute Relative Bias (%) by Parameter, Family, and Method",
    x = "Method",
    y = "Distribution Family"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 12)
  )

# Save the heatmap
ggsave(file.path(plots_dir, "rel_bias_heatmap.pdf"),
       rel_bias_heatmap, width = 12, height = 8)

# ---------------------------------------------------------------------------- #
# Generate Parameter-Focused Plots (Matrix Layout matching reference image)
# ---------------------------------------------------------------------------- #

cat("Generating parameter-focused matrix plots...\n")

# Helper function to create and save parameter-focused matrix plots
create_and_save_param_matrix_plot <- function(plot_func, param, filename_prefix, dir = param_plots_dir, ...) {
  cat("  Creating", filename_prefix, "matrix plot for parameter", param, "...\n")

  p <- plot_func(param_results, param, ...)

  # Save the plot
  filename <- file.path(dir, paste0(filename_prefix, "_", param, ".pdf"))
  ggsave(filename, p, width = 14, height = 12)

  return(p)
}

# Define family ordering
family_order <- c("gkw", "bkw", "kkw", "ekw", "mc", "kw", "beta")

# For each parameter, create the new parameter-focused matrix plots
for(param in unique_params) {
  # Set appropriate limits based on parameter
  abs_limits <- NULL
  rel_limits <- NULL

  if(param %in% c("alpha", "beta")) {
    abs_limits <- c(-0.5, 0.5)
    rel_limits <- c(-50, 50)
  } else if(param %in% c("gamma", "delta")) {
    abs_limits <- c(-0.75, 0.75)
    rel_limits <- c(-75, 75)
  } else if(param == "lambda") {
    abs_limits <- c(-1, 1)
    rel_limits <- c(-60, 60)
  }

  # 1. Absolute bias matrix plot
  abs_bias_plot <- create_and_save_param_matrix_plot(
    plot_parameter_bias_matrix,
    param,
    "abs_bias_matrix",
    limits = abs_limits,
    family_sort = family_order,
    title = paste("Absolute Bias for Parameter", param, "by Family and Method")
  )

  # 2. Relative bias matrix plot
  rel_bias_plot <- create_and_save_param_matrix_plot(
    plot_parameter_rel_bias_matrix,
    param,
    "rel_bias_matrix",
    limits = rel_limits,
    family_sort = family_order,
    title = paste("Relative Bias (%) for Parameter", param, "by Family and Method")
  )
}

# Create parameter-focused heatmaps for different sample sizes
for (sample_size in sort(unique(as.numeric(as.character(unique(param_results$sample_size)))))) {
  cat("  Creating parameter bias heatmaps for sample size", sample_size, "...\n")

  # 1. Relative bias heatmap
  rel_bias_heatmap <- plot_parameter_bias_heatmap(
    param_results,
    sample_size_to_show = sample_size,
    measure = "rel_bias",
    use_abs = TRUE,
    title = paste("Mean Absolute Relative Bias (%) by Parameter, Family, and Method (n =", sample_size, ")")
  )

  ggsave(file.path(param_plots_dir, paste0("rel_bias_heatmap_n", sample_size, ".pdf")),
         rel_bias_heatmap, width = 12, height = 8)

  # 2. Standardized bias heatmap
  std_bias_heatmap <- plot_parameter_bias_heatmap(
    param_results,
    sample_size_to_show = sample_size,
    measure = "bias_z",
    use_abs = TRUE,
    title = paste("Mean Absolute Standardized Bias by Parameter, Family, and Method (n =", sample_size, ")")
  )

  ggsave(file.path(param_plots_dir, paste0("std_bias_heatmap_n", sample_size, ".pdf")),
         std_bias_heatmap, width = 12, height = 8)
}

# Create consolidated parameter comparison plots
cat("  Creating consolidated parameter comparison plots...\n")

# Average absolute relative bias across all parameters
param_rel_bias_comparison <- param_results %>%
  group_by(parameter, method, sample_size) %>%
  summarize(
    mean_abs_rel_bias = mean(abs(rel_bias), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    sample_size_num = as.numeric(as.character(sample_size)),
    method_label = case_when(
      method == "tmb" ~ "TMB",
      method == "nr" ~ "Newton-Raphson",
      method == "nlminb" ~ "nlminb",
      method == "optim" ~ "optim (L-BFGS-B)",
      TRUE ~ as.character(method)
    )
  )

# Create plot showing bias trends across all parameters
param_comparison_plot <- ggplot(param_rel_bias_comparison,
                                aes(x = sample_size_num, y = mean_abs_rel_bias,
                                    color = method_label)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.5) +
  facet_wrap(~ parameter, scales = "free_y") +
  scale_x_log10(breaks = unique(param_rel_bias_comparison$sample_size_num),
                labels = unique(param_rel_bias_comparison$sample_size_num)) +
  labs(
    title = "Mean Absolute Relative Bias (%) by Parameter and Method",
    x = "Sample Size (log scale)",
    y = "Mean Absolute Relative Bias (%)",
    color = "Method"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 12),
    legend.position = "top"
  )

ggsave(file.path(param_plots_dir, "param_comparison_rel_bias.pdf"),
       param_comparison_plot, width = 12, height = 8)

# Create method-focused plots for convergence and computation time
cat("  Creating method performance comparison plots...\n")

# 1. Convergence comparison
for (sample_size in sort(unique(as.numeric(as.character(unique(convergence_df$sample_size)))))) {
  conv_comparison <- plot_convergence_by_method_family(
    convergence_df,
    sample_size_to_show = sample_size,
    title = paste("Convergence Rates by Method and Family (n =", sample_size, ")")
  )

  ggsave(file.path(param_plots_dir, paste0("convergence_comparison_n", sample_size, ".pdf")),
         conv_comparison, width = 10, height = 8)
}

# 2. Computation time comparison
for (sample_size in sort(unique(as.numeric(as.character(unique(convergence_df$sample_size)))))) {
  time_comparison <- plot_computation_time_by_method_family(
    convergence_df,
    sample_size_to_show = sample_size,
    log_scale = TRUE,
    title = paste("Computation Time by Method and Family (n =", sample_size, ")")
  )

  ggsave(file.path(param_plots_dir, paste0("time_comparison_n", sample_size, ".pdf")),
         time_comparison, width = 10, height = 8)
}

# ---------------------------------------------------------------------------- #
# Summarize Overall Findings
# ---------------------------------------------------------------------------- #

cat("Creating final summary tables...\n")

# Create a summary of best methods by parameter
best_methods <- method_rankings %>%
  arrange(family, desc(mean_score)) %>%
  group_by(family) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(family, method, mean_score) %>%
  mutate(mean_score = round(mean_score, 3)) %>%
  rename(best_method = method, score = mean_score)

write_csv(best_methods, file.path(csv_dir, "gkw_best_methods.csv"))

# Create a summary of relative bias by parameter
rel_bias_summary <- param_results %>%
  group_by(parameter, family, method) %>%
  summarize(
    mean_rel_bias = mean(rel_bias, na.rm = TRUE),
    mean_abs_rel_bias = mean(abs(rel_bias), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(parameter, family, mean_abs_rel_bias)

write_csv(rel_bias_summary, file.path(csv_dir, "gkw_rel_bias_summary.csv"))

# Create parameter-specific best method rankings
param_method_rankings <- param_results %>%
  group_by(parameter, method) %>%
  summarize(
    mean_abs_rel_bias = mean(abs(rel_bias), na.rm = TRUE),
    mean_bias_z = mean(abs(bias_z), na.rm = TRUE),
    mean_coverage = mean(coverage, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(parameter) %>%
  mutate(
    rel_bias_rank = rank(mean_abs_rel_bias),
    bias_z_rank = rank(mean_bias_z),
    coverage_rank = rank(abs(mean_coverage - 0.95)),
    overall_rank = (rel_bias_rank + bias_z_rank + coverage_rank) / 3
  ) %>%
  arrange(parameter, overall_rank) %>%
  select(parameter, method, mean_abs_rel_bias, mean_bias_z, mean_coverage, overall_rank)

write_csv(param_method_rankings, file.path(csv_dir, "gkw_param_method_rankings.csv"))

# Create a final report summary
cat("Creating final report summary...\n")

# Calculate best method by parameter
best_method_by_param <- param_method_rankings %>%
  group_by(parameter) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(parameter, method) %>%
  rename(best_method = method)

final_summary <- data.frame(
  Metric = c(
    "Total Parameters Evaluated",
    "Total Distribution Families",
    "Total Methods Compared",
    "Sample Sizes Tested",
    "Simulation Replications",
    "Best Overall Method",
    "Lowest Mean Absolute Relative Bias Method",
    "Highest Coverage Rate Method",
    "Fastest Method",
    "Most Problematic Parameters",
    "Most Challenging Configurations"
  ),
  Value = c(
    length(unique_params),
    length(unique_families),
    length(methods),
    paste(sample_sizes, collapse = ", "),
    run_config$n_sim,
    names(sort(table(best_methods$best_method), decreasing = TRUE)[1]),
    names(sort(tapply(rel_bias_summary$mean_abs_rel_bias, rel_bias_summary$method, mean),
               decreasing = FALSE)[1]),
    names(sort(tapply(coverage_summary$mean_coverage, coverage_summary$method, mean),
               decreasing = TRUE)[1]),
    names(sort(tapply(convergence_summary$mean_time, convergence_summary$method, mean),
               decreasing = FALSE)[1]),
    paste(names(sort(tapply(rel_bias_summary$mean_abs_rel_bias, rel_bias_summary$parameter, mean),
                     decreasing = TRUE)[1:2]), collapse = ", "),
    paste(names(sort(tapply(rel_bias_summary$mean_abs_rel_bias, rel_bias_summary$family, mean),
                     decreasing = TRUE)[1:2]), collapse = ", ")
  )
)

# Add parameter-specific best methods
for (param in unique_params) {
  best_method <- best_method_by_param %>%
    filter(parameter == param) %>%
    pull(best_method)

  final_summary <- rbind(
    final_summary,
    data.frame(
      Metric = paste("Best Method for Parameter", param),
      Value = best_method
    )
  )
}

write_csv(final_summary, file.path(csv_dir, "gkw_final_summary.csv"))

# Print summary to console
cat("\nMLE Bias Study Summary:\n")
knitr::kable(final_summary)

cat("\nCompleted analysis. All results saved to", csv_dir, "and", plots_dir, "directories.\n")


# # ---------------------------------------------------------------------------- #
# # Title: Analysis Script for GKw MLE Bias Study
# # Description: This script analyzes the results of the GKw MLE bias simulation
# #              study, generates plots and tables for all parameters and methods.
# # ---------------------------------------------------------------------------- #
#
# # Load required packages
# library(tidyverse)
# library(ggplot2)
# library(patchwork)
# library(ggh4x)    # For facet_nested
# library(latex2exp) # For rendering LaTeX in plots
# library(knitr)
# library(kableExtra)
#
# # Load helper functions
# source("gkw_bias_extra_functions.R")
#
# # ---------------------------------------------------------------------------- #
# # Load Results
# # ---------------------------------------------------------------------------- #
#
# # Create output directories if they don't exist
# plots_dir <- "plots"
# param_plots_dir <- file.path(plots_dir, "parameter_focus")
# csv_dir <- "csv"
#
# dir_list <- c(plots_dir, param_plots_dir, csv_dir)
# for (dir_path in dir_list) {
#   if (!dir.exists(dir_path)) {
#     dir.create(dir_path, recursive = TRUE)
#   }
# }
#
# # Load the aggregated results
# param_results <- readRDS("rds/gkw_bias_param_results.rds")
# convergence_df <- readRDS("rds/gkw_bias_convergence_results.rds")
# run_config <- readRDS("rds/run_config.rds")
#
# # Extract simulation parameters from config
# methods <- run_config$methods
# sample_sizes <- run_config$sample_sizes
# configurations <- run_config$configurations
#
# # Verify data loaded correctly
# cat("Loaded parameter results:", nrow(param_results), "rows\n")
# cat("Loaded convergence results:", nrow(convergence_df), "rows\n")
#
# # Extract unique parameters, families, and configurations
# unique_params <- unique(param_results$parameter)
# unique_families <- unique(param_results$family)
# unique_configs <- unique(param_results$config)
#
# cat("Parameters in results:", paste(unique_params, collapse=", "), "\n")
# cat("Families in results:", paste(unique_families, collapse=", "), "\n")
#
# # ---------------------------------------------------------------------------- #
# # Generate Summary Tables
# # ---------------------------------------------------------------------------- #
#
# cat("Generating summary tables...\n")
#
# # 1. Bias summary table
# bias_summary <- create_bias_summary(param_results)
# write_csv(bias_summary, file.path(csv_dir, "gkw_bias_summary.csv"))
#
# # 2. Coverage summary table
# coverage_summary <- create_coverage_summary(param_results)
# write_csv(coverage_summary, file.path(csv_dir, "gkw_coverage_summary.csv"))
#
# # 3. Convergence and timing summary
# convergence_summary <- create_convergence_summary(convergence_df)
# write_csv(convergence_summary, file.path(csv_dir, "gkw_convergence_summary.csv"))
#
# # 4. Method comparison summary
# method_comparison <- create_method_comparison(param_results, convergence_df)
# write_csv(method_comparison, file.path(csv_dir, "gkw_method_comparison.csv"))
#
# # 5. Method rankings
# method_rankings <- rank_methods(param_results, convergence_df,
#                                 rel_bias_weight = 0.4,
#                                 coverage_weight = 0.3,
#                                 convergence_weight = 0.2,
#                                 time_weight = 0.1)
# write_csv(method_rankings, file.path(csv_dir, "gkw_method_rankings.csv"))
#
# # 6. Identify problematic configurations
# problematic_configs <- identify_problematic_configs(param_results, convergence_df,
#                                                     convergence_threshold = 0.8,
#                                                     rel_bias_threshold = 20.0)
#
# write_csv(problematic_configs$low_convergence,
#           file.path(csv_dir, "gkw_low_convergence_configs.csv"))
# write_csv(problematic_configs$high_bias,
#           file.path(csv_dir, "gkw_high_bias_configs.csv"))
#
# # 7. Relative bias comparison table
# rel_bias_comparison <- param_results %>%
#   group_by(parameter, family, method, sample_size) %>%
#   summarize(
#     true_value = mean(true_value),
#     abs_bias = mean(bias, na.rm = TRUE),
#     rel_bias_pct = mean(rel_bias, na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   filter(sample_size == max(sample_size)) %>%  # Focus on largest sample size
#   select(parameter, family, method, true_value, abs_bias, rel_bias_pct) %>%
#   mutate(
#     true_value = round(true_value, 4),
#     abs_bias = round(abs_bias, 4),
#     rel_bias_pct = round(rel_bias_pct, 2)
#   ) %>%
#   arrange(parameter, family, method)
#
# write_csv(rel_bias_comparison, file.path(csv_dir, "gkw_rel_bias_comparison.csv"))
#
# # Generate a table comparing different bias measures across parameters and methods
# bias_measures_comparison <- param_results %>%
#   group_by(parameter, method) %>%
#   summarize(
#     mean_abs_bias = mean(abs(bias), na.rm = TRUE),
#     mean_abs_rel_bias = mean(abs(rel_bias), na.rm = TRUE),
#     mean_abs_bias_z = mean(abs(bias_z), na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   mutate(
#     mean_abs_bias = round(mean_abs_bias, 4),
#     mean_abs_rel_bias = round(mean_abs_rel_bias, 2),
#     mean_abs_bias_z = round(mean_abs_bias_z, 2)
#   ) %>%
#   arrange(parameter, method)
#
# write_csv(bias_measures_comparison, file.path(csv_dir, "gkw_bias_measures_comparison.csv"))
#
# # ---------------------------------------------------------------------------- #
# # Generate Traditional Plots (Original Layout)
# # ---------------------------------------------------------------------------- #
#
# cat("Generating traditional plots...\n")
#
# # Create a helper function to create and save plots
# create_and_save_plot <- function(plot_func, param, filename_prefix, ...) {
#   cat("  Creating", filename_prefix, "plot for parameter", param, "...\n")
#
#   p <- plot_func(param_results, param, ...)
#
#   # Save the plot
#   filename <- file.path(plots_dir, paste0(filename_prefix, "_", param, ".pdf"))
#   ggsave(filename, p, width = 10, height = 8)
#
#   return(p)
# }
#
# # Plot each type of graph for each parameter
# for(param in unique_params) {
#   # 1. Absolute bias plots
#   abs_bias_plot <- create_and_save_plot(
#     plot_abs_bias,
#     param,
#     "abs_bias",
#     title = paste("Absolute Bias for Parameter", param)
#   )
#
#   # 2. Relative bias plots
#   # Set appropriate limits based on parameter
#   # bias_limits <- c(-8, 8)
#   # if(param %in% c("alpha", "beta")) {
#   #   bias_limits <- c(-3, 3)  # Adjust as needed
#   # } else if(param %in% c("gamma", "delta")) {
#   #   bias_limits <- c(-3, 3)  # Adjust as needed
#   # } else if(param == "lambda") {
#   #   bias_limits <- c(-3, 3)  # Adjust as needed
#   # }
#
#   rel_bias_plot <- create_and_save_plot(
#     plot_rel_bias,
#     param,
#     "rel_bias",
#     limits = c(-20, 20),
#     title = paste("Relative Bias (%) for Parameter", param)
#   )
#
#   # 3. Standardized bias plots
#   bias_z_plot <- create_and_save_plot(
#     plot_bias_z,
#     param,
#     "bias_z",
#     limits = c(-1.96, 1.96),
#     title = paste("Standardized Bias for Parameter", param)
#   )
#
#   # 4. Coverage plots
#   coverage_plot <- create_and_save_plot(
#     plot_coverage,
#     param,
#     "coverage",
#     limits = c(0.7, 1.0),
#     title = paste("Coverage Rate for Parameter", param)
#   )
#
#   # 5. Relative bias trend plots
#   rel_bias_trend_plot <- analyze_rel_bias_trend(param_results, param)
#   ggsave(file.path(plots_dir, paste0("rel_bias_trend_", param, ".pdf")),
#          rel_bias_trend_plot, width = 8, height = 6)
# }
#
# # Generate convergence plot
# cat("  Creating convergence plot...\n")
# conv_plot <- plot_convergence(convergence_df)
# ggsave(file.path(plots_dir, "convergence.pdf"), conv_plot, width = 10, height = 8)
#
# # Generate computation time plot
# cat("  Creating computation time plot...\n")
# time_plot <- plot_computation_time(convergence_df, log_scale = TRUE)
# ggsave(file.path(plots_dir, "computation_time.pdf"), time_plot, width = 10, height = 8)
#
# # Generate consolidated plots for manuscript
# # Create a plot that shows relative bias for all parameters on the same grid
# cat("  Creating consolidated relative bias plot...\n")
#
# # Filter data by parameter and create separate plots
# rel_bias_plots <- lapply(unique_params, function(param) {
#   plot_rel_bias(
#     param_results %>% filter(parameter == param),
#     param,
#     title = NULL
#   ) +
#     theme(legend.position = "none") +
#     labs(title = param)
# })
#
# # Combine all plots using patchwork
# combined_rel_bias <- wrap_plots(rel_bias_plots, ncol = 2)
#
# # Add a common legend
# combined_rel_bias <- combined_rel_bias +
#   plot_layout(guides = "collect") &
#   theme(legend.position = "bottom")
#
# # Save the combined plot
# ggsave(file.path(plots_dir, "combined_rel_bias.pdf"),
#        combined_rel_bias, width = 16, height = 12)
#
# # Create a consolidated method comparison plot
# cat("  Creating method comparison plot...\n")
#
# # Create a long format dataset for the bias measures comparison
# method_comparison_long <- bias_measures_comparison %>%
#   pivot_longer(
#     cols = c(mean_abs_bias, mean_abs_rel_bias, mean_abs_bias_z),
#     names_to = "bias_measure",
#     values_to = "value"
#   ) %>%
#   mutate(
#     bias_measure = factor(
#       bias_measure,
#       levels = c("mean_abs_bias", "mean_abs_rel_bias", "mean_abs_bias_z"),
#       labels = c("Absolute Bias", "Relative Bias (%)", "Standardized Bias")
#     )
#   )
#
# # Create the method comparison plot
# method_comp_plot <- ggplot(method_comparison_long,
#                            aes(x = method, y = value, fill = method)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   facet_grid(bias_measure ~ parameter, scales = "free_y") +
#   labs(
#     title = "Comparison of Bias Measures by Method and Parameter",
#     y = "Value",
#     x = "Method",
#     fill = "Method"
#   ) +
#   theme_bw() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
#     strip.text = element_text(size = 10),
#     plot.title = element_text(hjust = 0.5, size = 12)
#   )
#
# # Save the method comparison plot
# ggsave(file.path(plots_dir, "method_comparison.pdf"),
#        method_comp_plot, width = 12, height = 8)
#
# # Create a heatmap of relative bias by parameter, family, and method
# cat("  Creating relative bias heatmap...\n")
#
# # Create a summary for the heatmap
# rel_bias_heatmap_data <- param_results %>%
#   group_by(parameter, family, method) %>%
#   summarize(
#     mean_abs_rel_bias = mean(abs(rel_bias), na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   # Cap extreme values for better visualization
#   mutate(
#     mean_abs_rel_bias_capped = pmin(mean_abs_rel_bias, 50)
#   )
#
# # Create the heatmap
# rel_bias_heatmap <- ggplot(rel_bias_heatmap_data,
#                            aes(x = method, y = family, fill = mean_abs_rel_bias_capped)) +
#   geom_tile() +
#   facet_wrap(~ parameter, ncol = 3) +
#   scale_fill_viridis_c(option = "plasma", name = "Mean Absolute\nRelative Bias (%)") +
#   labs(
#     title = "Mean Absolute Relative Bias (%) by Parameter, Family, and Method",
#     x = "Method",
#     y = "Distribution Family"
#   ) +
#   theme_bw() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     strip.text = element_text(size = 10),
#     plot.title = element_text(hjust = 0.5, size = 12)
#   )
#
# # Save the heatmap
# ggsave(file.path(plots_dir, "rel_bias_heatmap.pdf"),
#        rel_bias_heatmap, width = 12, height = 8)
#
# # ---------------------------------------------------------------------------- #
# # Generate Parameter-Focused Plots (Matrix Layout matching reference image)
# # ---------------------------------------------------------------------------- #
#
# cat("Generating parameter-focused matrix plots...\n")
#
# # Helper function to create and save parameter-focused matrix plots
# create_and_save_param_matrix_plot <- function(plot_func, param, filename_prefix, dir = param_plots_dir, ...) {
#   cat("  Creating", filename_prefix, "matrix plot for parameter", param, "...\n")
#
#   p <- plot_func(param_results, param, ...)
#
#   # Save the plot
#   filename <- file.path(dir, paste0(filename_prefix, "_", param, ".pdf"))
#   ggsave(filename, p, width = 14, height = 12)
#
#   return(p)
# }
#
# # For each parameter, create the new parameter-focused matrix plots
# for(param in unique_params) {
#   # Set appropriate limits based on parameter
#   abs_limits <- NULL
#   rel_limits <- NULL
#
#   if(param %in% c("alpha", "beta")) {
#     abs_limits <- c(-0.5, 0.5)
#     rel_limits <- c(-50, 50)
#   } else if(param %in% c("gamma", "delta")) {
#     abs_limits <- c(-0.75, 0.75)
#     rel_limits <- c(-75, 75)
#   } else if(param == "lambda") {
#     abs_limits <- c(-1, 1)
#     rel_limits <- c(-60, 60)
#   }
#
#   # 1. Absolute bias matrix plot
#   abs_bias_plot <- create_and_save_param_matrix_plot(
#     plot_parameter_bias_matrix,
#     param,
#     "abs_bias_matrix",
#     limits = abs_limits,
#     title = paste("Absolute Bias for Parameter", param, "by Family and Method")
#   )
#
#   # 2. Relative bias matrix plot
#   rel_bias_plot <- create_and_save_param_matrix_plot(
#     plot_parameter_rel_bias_matrix,
#     param,
#     "rel_bias_matrix",
#     limits = rel_limits,
#     title = paste("Relative Bias (%) for Parameter", param, "by Family and Method")
#   )
# }
#
# # Create parameter-focused heatmaps for different sample sizes
# for (sample_size in sort(unique(as.numeric(as.character(unique(param_results$sample_size)))))) {
#   cat("  Creating parameter bias heatmaps for sample size", sample_size, "...\n")
#
#   # 1. Relative bias heatmap
#   rel_bias_heatmap <- plot_parameter_bias_heatmap(
#     param_results,
#     sample_size_to_show = sample_size,
#     measure = "rel_bias",
#     use_abs = TRUE,
#     title = paste("Mean Absolute Relative Bias (%) by Parameter, Family, and Method (n =", sample_size, ")")
#   )
#
#   ggsave(file.path(param_plots_dir, paste0("rel_bias_heatmap_n", sample_size, ".pdf")),
#          rel_bias_heatmap, width = 12, height = 8)
#
#   # 2. Standardized bias heatmap
#   std_bias_heatmap <- plot_parameter_bias_heatmap(
#     param_results,
#     sample_size_to_show = sample_size,
#     measure = "bias_z",
#     use_abs = TRUE,
#     title = paste("Mean Absolute Standardized Bias by Parameter, Family, and Method (n =", sample_size, ")")
#   )
#
#   ggsave(file.path(param_plots_dir, paste0("std_bias_heatmap_n", sample_size, ".pdf")),
#          std_bias_heatmap, width = 12, height = 8)
# }
#
# # Create consolidated parameter comparison plots
# cat("  Creating consolidated parameter comparison plots...\n")
#
# # Average absolute relative bias across all parameters
# param_rel_bias_comparison <- param_results %>%
#   group_by(parameter, method, sample_size) %>%
#   summarize(
#     mean_abs_rel_bias = mean(abs(rel_bias), na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   mutate(
#     sample_size_num = as.numeric(as.character(sample_size)),
#     method_label = case_when(
#       method == "tmb" ~ "TMB",
#       method == "nr" ~ "Newton-Raphson",
#       method == "nlminb" ~ "nlminb",
#       method == "optim" ~ "optim (L-BFGS-B)",
#       TRUE ~ as.character(method)
#     )
#   )
#
# # Create plot showing bias trends across all parameters
# param_comparison_plot <- ggplot(param_rel_bias_comparison,
#                                 aes(x = sample_size_num, y = mean_abs_rel_bias,
#                                     color = method_label)) +
#   geom_line(linewidth = 0.6) +
#   geom_point(size = 1.5) +
#   facet_wrap(~ parameter, scales = "free_y") +
#   scale_x_log10(breaks = unique(param_rel_bias_comparison$sample_size_num),
#                 labels = unique(param_rel_bias_comparison$sample_size_num)) +
#   labs(
#     title = "Mean Absolute Relative Bias (%) by Parameter and Method",
#     x = "Sample Size (log scale)",
#     y = "Mean Absolute Relative Bias (%)",
#     color = "Method"
#   ) +
#   theme_bw() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     strip.text = element_text(size = 10),
#     plot.title = element_text(hjust = 0.5, size = 12),
#     legend.position = "top"
#   )
#
# ggsave(file.path(param_plots_dir, "param_comparison_rel_bias.pdf"),
#        param_comparison_plot, width = 12, height = 8)
#
# # Create method-focused plots for convergence and computation time
# cat("  Creating method performance comparison plots...\n")
#
# # 1. Convergence comparison
# for (sample_size in sort(unique(as.numeric(as.character(unique(convergence_df$sample_size)))))) {
#   conv_comparison <- plot_convergence_by_method_family(
#     convergence_df,
#     sample_size_to_show = sample_size,
#     title = paste("Convergence Rates by Method and Family (n =", sample_size, ")")
#   )
#
#   ggsave(file.path(param_plots_dir, paste0("convergence_comparison_n", sample_size, ".pdf")),
#          conv_comparison, width = 10, height = 8)
# }
#
# # 2. Computation time comparison
# for (sample_size in sort(unique(as.numeric(as.character(unique(convergence_df$sample_size)))))) {
#   time_comparison <- plot_computation_time_by_method_family(
#     convergence_df,
#     sample_size_to_show = sample_size,
#     log_scale = TRUE,
#     title = paste("Computation Time by Method and Family (n =", sample_size, ")")
#   )
#
#   ggsave(file.path(param_plots_dir, paste0("time_comparison_n", sample_size, ".pdf")),
#          time_comparison, width = 10, height = 8)
# }
#
# # ---------------------------------------------------------------------------- #
# # Summarize Overall Findings
# # ---------------------------------------------------------------------------- #
#
# cat("Creating final summary tables...\n")
#
# # Create a summary of best methods by parameter
# best_methods <- method_rankings %>%
#   arrange(family, desc(mean_score)) %>%
#   group_by(family) %>%
#   slice_head(n = 1) %>%
#   ungroup() %>%
#   select(family, method, mean_score) %>%
#   mutate(mean_score = round(mean_score, 3)) %>%
#   rename(best_method = method, score = mean_score)
#
# write_csv(best_methods, file.path(csv_dir, "gkw_best_methods.csv"))
#
# # Create a summary of relative bias by parameter
# rel_bias_summary <- param_results %>%
#   group_by(parameter, family, method) %>%
#   summarize(
#     mean_rel_bias = mean(rel_bias, na.rm = TRUE),
#     mean_abs_rel_bias = mean(abs(rel_bias), na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   arrange(parameter, family, mean_abs_rel_bias)
#
# write_csv(rel_bias_summary, file.path(csv_dir, "gkw_rel_bias_summary.csv"))
#
# # Create parameter-specific best method rankings
# param_method_rankings <- param_results %>%
#   group_by(parameter, method) %>%
#   summarize(
#     mean_abs_rel_bias = mean(abs(rel_bias), na.rm = TRUE),
#     mean_bias_z = mean(abs(bias_z), na.rm = TRUE),
#     mean_coverage = mean(coverage, na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   group_by(parameter) %>%
#   mutate(
#     rel_bias_rank = rank(mean_abs_rel_bias),
#     bias_z_rank = rank(mean_bias_z),
#     coverage_rank = rank(abs(mean_coverage - 0.95)),
#     overall_rank = (rel_bias_rank + bias_z_rank + coverage_rank) / 3
#   ) %>%
#   arrange(parameter, overall_rank) %>%
#   select(parameter, method, mean_abs_rel_bias, mean_bias_z, mean_coverage, overall_rank)
#
# write_csv(param_method_rankings, file.path(csv_dir, "gkw_param_method_rankings.csv"))
#
# # Create a final report summary
# cat("Creating final report summary...\n")
#
# # Calculate best method by parameter
# best_method_by_param <- param_method_rankings %>%
#   group_by(parameter) %>%
#   slice_head(n = 1) %>%
#   ungroup() %>%
#   select(parameter, method) %>%
#   rename(best_method = method)
#
# final_summary <- data.frame(
#   Metric = c(
#     "Total Parameters Evaluated",
#     "Total Distribution Families",
#     "Total Methods Compared",
#     "Sample Sizes Tested",
#     "Simulation Replications",
#     "Best Overall Method",
#     "Lowest Mean Absolute Relative Bias Method",
#     "Highest Coverage Rate Method",
#     "Fastest Method",
#     "Most Problematic Parameters",
#     "Most Challenging Configurations"
#   ),
#   Value = c(
#     length(unique_params),
#     length(unique_families),
#     length(methods),
#     paste(sample_sizes, collapse = ", "),
#     run_config$n_sim,
#     names(sort(table(best_methods$best_method), decreasing = TRUE)[1]),
#     names(sort(tapply(rel_bias_summary$mean_abs_rel_bias, rel_bias_summary$method, mean),
#                decreasing = FALSE)[1]),
#     names(sort(tapply(coverage_summary$mean_coverage, coverage_summary$method, mean),
#                decreasing = TRUE)[1]),
#     names(sort(tapply(convergence_summary$mean_time, convergence_summary$method, mean),
#                decreasing = FALSE)[1]),
#     paste(names(sort(tapply(rel_bias_summary$mean_abs_rel_bias, rel_bias_summary$parameter, mean),
#                      decreasing = TRUE)[1:2]), collapse = ", "),
#     paste(names(sort(tapply(rel_bias_summary$mean_abs_rel_bias, rel_bias_summary$family, mean),
#                      decreasing = TRUE)[1:2]), collapse = ", ")
#   )
# )
#
# # Add parameter-specific best methods
# for (param in unique_params) {
#   best_method <- best_method_by_param %>%
#     filter(parameter == param) %>%
#     pull(best_method)
#
#   final_summary <- rbind(
#     final_summary,
#     data.frame(
#       Metric = paste("Best Method for Parameter", param),
#       Value = best_method
#     )
#   )
# }
#
# write_csv(final_summary, file.path(csv_dir, "gkw_final_summary.csv"))
#
# # Print summary to console
# cat("\nMLE Bias Study Summary:\n")
# knitr::kable(final_summary)
#
# cat("\nCompleted analysis. All results saved to", csv_dir, "and", plots_dir, "directories.\n")
