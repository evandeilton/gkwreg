# ---------------------------------------------------------------------------- #
# Title: Main Execution Script for GKw MLE Bias Study
# Description: This script orchestrates the complete GKw MLE bias study workflow
#              with enhanced parameter-focused visualization analysis
# ---------------------------------------------------------------------------- #

# Record start time for the entire study
total_start_time <- Sys.time()

# Create timestamp for this run
run_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# Create master output directory
master_dir <- paste0("run_", run_timestamp)
dir.create(master_dir, showWarnings = FALSE)

# Set working directory to the master directory
original_dir <- getwd()
setwd(master_dir)

# Create log file
log_file <- file("gkw_bias_study_master_log.txt", "w")
writeLines(paste("GKw Bias Study Master Log -", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), log_file)
writeLines("=================================================================", log_file)
close(log_file)

# Function to log messages
log_message <- function(message) {
  timestamp <- format(Sys.time(), "%H:%M:%S")
  message_with_time <- paste0("[", timestamp, "] ", message)

  # Print to console
  cat(message_with_time, "\n")

  # Append to log file
  log_file <- file("gkw_bias_study_master_log.txt", "a")
  writeLines(message_with_time, log_file)
  close(log_file)
}

# Copy all required script files to the working directory
log_message("Setting up environment...")

source_files <- c(
  "gkw_bias_extra_functions.R",
  "gkw_bias_simulation.R",
  "gkw_bias_analysis.R",
  "gkw_bias_report.Rmd"
)

for (file in source_files) {
  if (file.exists(file.path(original_dir, file))) {
    file.copy(file.path(original_dir, file), ".", overwrite = TRUE)
    log_message(paste("Copied", file, "to working directory"))
  } else {
    log_message(paste("WARNING: Could not find source file:", file))
  }
}

# Step 1: Check for required packages and install if missing
log_message("Checking for required packages...")

required_packages <- c(
  "tidyverse", "ggplot2", "patchwork", "parallel", "foreach", "doParallel",
  "ggh4x", "latex2exp", "knitr", "kableExtra", "rmarkdown", "gkwreg"
)

missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if (length(missing_packages) > 0) {
  log_message(paste("Installing missing packages:", paste(missing_packages, collapse = ", ")))
  install.packages(missing_packages)
}

# Load utility functions first
source("gkw_bias_extra_functions.R")

# Step 2: Run the simulation
log_message("Starting simulation phase...")
simulation_start <- Sys.time()

# Run the full simulation
source("gkw_bias_simulation.R")

simulation_end <- Sys.time()
simulation_time <- difftime(simulation_end, simulation_start, units = "mins")
log_message(paste("Simulation phase completed in", round(simulation_time, 2), "minutes"))

# Step 3: Run the analysis
log_message("Starting analysis phase...")
analysis_start <- Sys.time()

source("gkw_bias_analysis.R")

analysis_end <- Sys.time()
analysis_time <- difftime(analysis_end, analysis_start, units = "mins")
log_message(paste("Analysis phase completed in", round(analysis_time, 2), "minutes"))

# Step 4: Generate the report
log_message("Generating final report...")
report_start <- Sys.time()

rmarkdown::render("gkw_bias_report.Rmd",
                  output_file = "gkw_bias_report.html",
                  quiet = TRUE)

report_end <- Sys.time()
report_time <- difftime(report_end, report_start, units = "mins")
log_message(paste("Report generation completed in", round(report_time, 2), "minutes"))

# Calculate total execution time
total_end_time <- Sys.time()
total_time <- difftime(total_end_time, total_start_time, units = "mins")

log_message("=================================================================")
log_message(paste("GKw Bias Study completed successfully in", round(total_time, 2), "minutes"))
log_message(paste("Results available in:", getwd()))
log_message("=================================================================")

# Create a summary file with execution statistics
summary_data <- data.frame(
  Phase = c("Simulation", "Analysis", "Report Generation", "Total"),
  Duration_Minutes = c(
    round(as.numeric(simulation_time), 2),
    round(as.numeric(analysis_time), 2),
    round(as.numeric(report_time), 2),
    round(as.numeric(total_time), 2)
  )
)

write.csv(summary_data, "execution_summary.csv", row.names = FALSE)

# Return to the original directory
setwd(original_dir)

# Print final message
cat("\n=================================================================\n")
cat("GKw Bias Study completed successfully!\n")
cat("Results available in:", master_dir, "\n")
cat("=================================================================\n")


# # ---------------------------------------------------------------------------- #
# # Title: Main Execution Script for GKw MLE Bias Study
# # Description: This script orchestrates the complete GKw MLE bias study workflow
# # ---------------------------------------------------------------------------- #
#
# # Record start time for the entire study
# total_start_time <- Sys.time()
#
# # Create timestamp for this run
# run_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
#
# # Create master output directory
# master_dir <- paste0("run_", run_timestamp)
# dir.create(master_dir, showWarnings = FALSE)
#
# # Set working directory to the master directory
# original_dir <- getwd()
# setwd(master_dir)
#
# # Create log file
# log_file <- file("gkw_bias_study_master_log.txt", "w")
# writeLines(paste("GKw Bias Study Master Log -", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), log_file)
# writeLines("=================================================================", log_file)
# close(log_file)
#
# # Function to log messages
# log_message <- function(message) {
#   timestamp <- format(Sys.time(), "%H:%M:%S")
#   message_with_time <- paste0("[", timestamp, "] ", message)
#
#   # Print to console
#   cat(message_with_time, "\n")
#
#   # Append to log file
#   log_file <- file("gkw_bias_study_master_log.txt", "a")
#   writeLines(message_with_time, log_file)
#   close(log_file)
# }
#
# # Copy all required script files to the working directory
# log_message("Setting up environment...")
#
# source_files <- c(
#   "gkw_bias_extra_functions.R",
#   "gkw_bias_simulation.R",
#   "gkw_bias_analysis.R",
#   "gkw_bias_report.Rmd"
# )
#
# for (file in source_files) {
#   if (file.exists(file.path(original_dir, file))) {
#     file.copy(file.path(original_dir, file), ".", overwrite = TRUE)
#     log_message(paste("Copied", file, "to working directory"))
#   } else {
#     log_message(paste("WARNING: Could not find source file:", file))
#   }
# }
#
# # Step 1: Check for required packages and install if missing
# log_message("Checking for required packages...")
#
# required_packages <- c(
#   "tidyverse", "ggplot2", "patchwork", "parallel", "foreach", "doParallel",
#   "ggh4x", "latex2exp", "knitr", "kableExtra", "rmarkdown", "gkwreg"
# )
#
# missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
#
# if (length(missing_packages) > 0) {
#   log_message(paste("Installing missing packages:", paste(missing_packages, collapse = ", ")))
#   install.packages(missing_packages)
# }
#
# # Load utility functions first
# source("gkw_bias_extra_functions.R")
#
#
# # Step 2: Run the simulation
# log_message("Starting simulation phase...")
# simulation_start <- Sys.time()
#
# # Run the full simulation
# source("gkw_bias_simulation.R")
#
#
# simulation_end <- Sys.time()
# simulation_time <- difftime(simulation_end, simulation_start, units = "mins")
# log_message(paste("Simulation phase completed in", round(simulation_time, 2), "minutes"))
#
# # Step 3: Run the analysis
# log_message("Starting analysis phase...")
# analysis_start <- Sys.time()
#
# source("gkw_bias_analysis.R")
#
# analysis_end <- Sys.time()
# analysis_time <- difftime(analysis_end, analysis_start, units = "mins")
# log_message(paste("Analysis phase completed in", round(analysis_time, 2), "minutes"))
#
# # Step 4: Generate the report
# log_message("Generating final report...")
# report_start <- Sys.time()
#
# rmarkdown::render("gkw_bias_report.Rmd",
#                   output_file = "gkw_bias_report.html",
#                   quiet = TRUE)
#
# report_end <- Sys.time()
# report_time <- difftime(report_end, report_start, units = "mins")
# log_message(paste("Report generation completed in", round(report_time, 2), "minutes"))
#
# # Calculate total execution time
# total_end_time <- Sys.time()
# total_time <- difftime(total_end_time, total_start_time, units = "mins")
#
# log_message("=================================================================")
# log_message(paste("GKw Bias Study completed successfully in", round(total_time, 2), "minutes"))
# log_message(paste("Results available in:", getwd()))
# log_message("=================================================================")
#
# # Create a summary file with execution statistics
# summary_data <- data.frame(
#   Phase = c("Simulation", "Analysis", "Report Generation", "Total"),
#   Duration_Minutes = c(
#     round(as.numeric(simulation_time), 2),
#     round(as.numeric(analysis_time), 2),
#     round(as.numeric(report_time), 2),
#     round(as.numeric(total_time), 2)
#   )
# )
#
# write.csv(summary_data, "execution_summary.csv", row.names = FALSE)
#
# # Return to the original directory
# setwd(original_dir)
#
# # Print final message
# cat("\n=================================================================\n")
# cat("GKw Bias Study completed successfully!\n")
# cat("Results available in:", master_dir, "\n")
# cat("=================================================================\n")
