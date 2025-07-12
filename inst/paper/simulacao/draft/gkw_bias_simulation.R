# ---------------------------------------------------------------------------- #
# Title: Simulation Script for GKw MLE Bias Study
# Description: This script runs the simulation study to evaluate bias in MLEs
#              for the Generalized Kumaraswamy distribution family using four
#              estimation methods: TMB, Newton-Raphson, nlminb, and optim.
# ---------------------------------------------------------------------------- #

# Load required packages
library(tidyverse)
library(parallel)
library(foreach)
library(doParallel)
library(ggh4x)    # For facet_nested
library(latex2exp) # For rendering LaTeX in plots
library(gkwreg)   # Assumed to contain all GKw distribution functions and fit methods

# Load helper functions
source("gkw_bias_extra_functions.R")

# ---------------------------------------------------------------------------- #
# Setup Configuration
# ---------------------------------------------------------------------------- #

# Define the configurations to test
configurations <- list(
  # GKw (5 parameters)
  gkw_bimodal = list(family = "gkw", params = list(alpha = 0.3, beta = 0.5, gamma = 1, delta = 1, lambda = 3)),
  gkw_u_shape = list(family = "gkw", params = list(alpha = 0.5, beta = 0.5, gamma = 0.5, delta = 0.5, lambda = 1)),
  gkw_j_shape = list(family = "gkw", params = list(alpha = 2, beta = 0.5, gamma = 1, delta = 1, lambda = 1)),

  # BKw (4 parameters)
  bkw_l_shape = list(family = "bkw", params = list(alpha = 1, beta = 3, gamma = 0.3, delta = 3)),
  bkw_unimodal = list(family = "bkw", params = list(alpha = 2, beta = 2, gamma = 2, delta = 2)),

  # KKw (4 parameters)
  kkw_j_shape = list(family = "kkw", params = list(alpha = 2, beta = 0.5, delta = 1, lambda = 1)),
  kkw_bimodal = list(family = "kkw", params = list(alpha = 0.5, beta = 5, delta = 30, lambda = 1)),

  # EKw (3 parameters)
  ekw_asym = list(family = "ekw", params = list(alpha = 2, beta = 3, lambda = 1.5)),
  ekw_extreme = list(family = "ekw", params = list(alpha = 0.5, beta = 0.5, lambda = 10)),

  # Mc (3 parameters)
  mc_u_shape = list(family = "mc", params = list(gamma = 0.5, delta = 0.5, lambda = 1)),
  mc_unimodal = list(family = "mc", params = list(gamma = 2, delta = 2, lambda = 1)),

  # Kw (2 parameters)
  kw_j_shape = list(family = "kw", params = list(alpha = 0.5, beta = 3)),
  kw_unimodal = list(family = "kw", params = list(alpha = 2, beta = 2)),

  # Beta (2 parameters)
  beta_u_shape = list(family = "beta", params = list(gamma = 0.5, delta = 0.5)),
  beta_unimodal = list(family = "beta", params = list(gamma = 2, delta = 2))
)

# Estimation methods to compare
methods <- c("tmb", "nr", "nlminb", "optim")

# Sample sizes to test
# sample_sizes <- c(50, 100, 250, 500, 1000)
sample_sizes <- c(250, 500)

# Number of simulation replications
# n_sim <- 1000
n_sim <- 200

# Create timestamp for this run
run_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# Create output directories
rds_dir <- "rds"
if (!dir.exists(rds_dir)) {
  dir.create(rds_dir, recursive = TRUE)
}

# Create a log file
log_file <- file(paste0("gkw_bias_study_log_", run_timestamp, ".txt"), "w")
writeLines(paste("GKw Bias Study Log -", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), log_file)
writeLines(paste("Number of simulations:", n_sim), log_file)
writeLines(paste("Sample sizes:", paste(sample_sizes, collapse=", ")), log_file)
writeLines(paste("Methods:", paste(methods, collapse=", ")), log_file)
writeLines(paste("Number of configurations:", length(configurations)), log_file)
close(log_file)

# Save run configuration
run_config <- list(
  timestamp = run_timestamp,
  n_sim = n_sim,
  sample_sizes = sample_sizes,
  methods = methods,
  configurations = configurations
)
saveRDS(run_config, file.path(rds_dir, "run_config.rds"))

# ---------------------------------------------------------------------------- #
# Execute Simulation
# ---------------------------------------------------------------------------- #

# Create all combinations of configurations
simulation_grid <- expand.grid(
  config = 1:length(configurations),
  method = methods,
  sample_size = sample_sizes,
  stringsAsFactors = FALSE
)

# Add reference columns for easier tracking
simulation_grid$family <- sapply(configurations[simulation_grid$config],
                                 function(x) x$family)
simulation_grid$config_name <- names(configurations)[simulation_grid$config]

# Print simulation overview
cat("Starting simulation with", n_sim, "iterations per configuration...\n")
cat("Total configurations:", nrow(simulation_grid), "\n")
cat("Output will be saved to", rds_dir, "directory\n")

# Set up parallel backend
n_cores <- max(1, detectCores() - 1)
cat("Using", n_cores, "cores for parallel processing...\n")
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Export necessary functions and data to the cluster
clusterExport(cl, c("configurations", "methods", "sample_sizes", "n_sim",
                    "rds_dir", "create_param_df", "compute_coverage",
                    "run_simulation", "generate_seed"))

# Load required libraries on the cluster
clusterEvalQ(cl, {
  library(tidyverse)
  library(gkwreg)  # Assumed to contain all GKw distribution functions and fit methods
})

# Run simulations in parallel
start_time_all <- Sys.time()

results <- foreach(i = 1:nrow(simulation_grid),
                   .packages = c("tidyverse", "gkwreg")) %dopar% {
                     config_idx <- simulation_grid$config[i]
                     method <- simulation_grid$method[i]
                     sample_size <- simulation_grid$sample_size[i]

                     # Print progress update (will be seen in log file)
                     cat("Running simulation for", names(configurations)[config_idx],
                         "with method", method, "and sample size", sample_size, "...\n")

                     # Run simulation for this configuration
                     sim_result <- run_simulation(configurations[[config_idx]], method, sample_size, n_sim)

                     # Save individual result to disk
                     result_filename <- paste0(
                       "gkw_bias_", simulation_grid$family[i], "_",
                       simulation_grid$config_name[i], "_",
                       method, "_n", sample_size, ".rds"
                     )
                     saveRDS(sim_result, file.path(rds_dir, result_filename))

                     # Return the result for aggregation
                     sim_result
                   }

end_time_all <- Sys.time()
total_time_all <- difftime(end_time_all, start_time_all, units = "mins")

# Clean up parallel cluster
stopCluster(cl)

# ---------------------------------------------------------------------------- #
# Process and Save Aggregated Results
# ---------------------------------------------------------------------------- #

# Extract parameter results from all simulations
param_results <- do.call(rbind, lapply(results, function(x) x$param_results))

# Extract convergence rates
convergence_df <- data.frame(
  config = simulation_grid$config,
  method = simulation_grid$method,
  sample_size = simulation_grid$sample_size,
  family = simulation_grid$family,
  config_name = simulation_grid$config_name,
  convergence_rate = sapply(results, function(x) x$convergence_rate),
  avg_time = sapply(results, function(x) x$avg_time)
)

# Save aggregated results to disk
saveRDS(param_results, file.path(rds_dir, "gkw_bias_param_results.rds"))
saveRDS(convergence_df, file.path(rds_dir, "gkw_bias_convergence_results.rds"))

# Log completion
log_file <- file(paste0("gkw_bias_study_log_", run_timestamp, ".txt"), "a")
writeLines(paste("Simulation completed in", round(total_time_all, 2), "minutes."), log_file)
writeLines(paste("Results saved to", rds_dir, "directory."), log_file)
writeLines(paste("Total parameter results:", nrow(param_results)), log_file)
writeLines(paste("Total convergence results:", nrow(convergence_df)), log_file)
close(log_file)

# Print summary
cat("\nSimulation completed in", round(total_time_all, 2), "minutes.\n")
cat("Results saved to", rds_dir, "directory.\n")
cat("Total parameter results:", nrow(param_results), "\n")
cat("Total convergence results:", nrow(convergence_df), "\n")





# # ---------------------------------------------------------------------------- #
# # Title: Simulation Script for GKw MLE Bias Study
# # Description: This script runs the simulation study to evaluate bias in MLEs
# #              for the Generalized Kumaraswamy distribution family using four
# #              estimation methods: TMB, Newton-Raphson, nlminb, and optim.
# # ---------------------------------------------------------------------------- #
#
# # Load required packages
# library(tidyverse)
# library(parallel)
# library(foreach)
# library(doParallel)
# library(ggh4x)    # For facet_nested
# library(latex2exp) # For rendering LaTeX in plots
# library(gkwreg)   # Assumed to contain all GKw distribution functions and fit methods
#
# # Load helper functions
# source("gkw_bias_extra_functions.R")
#
# # ---------------------------------------------------------------------------- #
# # Setup Configuration
# # ---------------------------------------------------------------------------- #
#
# # Define the configurations to test
# configurations <- list(
#   # GKw (5 parameters)
#   gkw_bimodal = list(family = "gkw", params = list(alpha = 0.3, beta = 0.5, gamma = 1, delta = 1, lambda = 3)),
#   gkw_u_shape = list(family = "gkw", params = list(alpha = 0.5, beta = 0.5, gamma = 0.5, delta = 0.5, lambda = 1)),
#   gkw_j_shape = list(family = "gkw", params = list(alpha = 2, beta = 0.5, gamma = 1, delta = 1, lambda = 1)),
#
#   # BKw (4 parameters)
#   bkw_l_shape = list(family = "bkw", params = list(alpha = 1, beta = 3, gamma = 0.3, delta = 3)),
#   bkw_unimodal = list(family = "bkw", params = list(alpha = 2, beta = 2, gamma = 2, delta = 2)),
#
#   # KKw (4 parameters)
#   kkw_j_shape = list(family = "kkw", params = list(alpha = 2, beta = 0.5, delta = 1, lambda = 1)),
#   kkw_bimodal = list(family = "kkw", params = list(alpha = 0.5, beta = 5, delta = 30, lambda = 1)),
#
#   # EKw (3 parameters)
#   ekw_asym = list(family = "ekw", params = list(alpha = 2, beta = 3, lambda = 1.5)),
#   ekw_extreme = list(family = "ekw", params = list(alpha = 0.5, beta = 0.5, lambda = 10)),
#
#   # Mc (3 parameters)
#   mc_u_shape = list(family = "mc", params = list(gamma = 0.5, delta = 0.5, lambda = 1)),
#   mc_unimodal = list(family = "mc", params = list(gamma = 2, delta = 2, lambda = 1)),
#
#   # Kw (2 parameters)
#   kw_j_shape = list(family = "kw", params = list(alpha = 0.5, beta = 3)),
#   kw_unimodal = list(family = "kw", params = list(alpha = 2, beta = 2)),
#
#   # Beta (2 parameters)
#   beta_u_shape = list(family = "beta", params = list(gamma = 0.5, delta = 0.5)),
#   beta_unimodal = list(family = "beta", params = list(gamma = 2, delta = 2))
# )
#
# # Estimation methods to compare
# methods <- c("tmb", "nr", "nlminb", "optim")
#
# # Sample sizes to test
# # sample_sizes <- c(50, 100, 250, 500, 1000)
# sample_sizes <- c(250, 500)
#
# # Number of simulation replications
# # n_sim <- 1000
# n_sim <- 200
#
# # Create timestamp for this run
# run_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
#
# # Create output directories
# rds_dir <- "rds"
# if (!dir.exists(rds_dir)) {
#   dir.create(rds_dir, recursive = TRUE)
# }
#
# # Create a log file
# log_file <- file(paste0("gkw_bias_study_log_", run_timestamp, ".txt"), "w")
# writeLines(paste("GKw Bias Study Log -", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), log_file)
# writeLines(paste("Number of simulations:", n_sim), log_file)
# writeLines(paste("Sample sizes:", paste(sample_sizes, collapse=", ")), log_file)
# writeLines(paste("Methods:", paste(methods, collapse=", ")), log_file)
# writeLines(paste("Number of configurations:", length(configurations)), log_file)
# close(log_file)
#
# # Save run configuration
# run_config <- list(
#   timestamp = run_timestamp,
#   n_sim = n_sim,
#   sample_sizes = sample_sizes,
#   methods = methods,
#   configurations = configurations
# )
# saveRDS(run_config, file.path(rds_dir, "run_config.rds"))
#
# # ---------------------------------------------------------------------------- #
# # Execute Simulation
# # ---------------------------------------------------------------------------- #
#
# # Create all combinations of configurations
# simulation_grid <- expand.grid(
#   config = 1:length(configurations),
#   method = methods,
#   sample_size = sample_sizes,
#   stringsAsFactors = FALSE
# )
#
# # Add reference columns for easier tracking
# simulation_grid$family <- sapply(configurations[simulation_grid$config],
#                                  function(x) x$family)
# simulation_grid$config_name <- names(configurations)[simulation_grid$config]
#
# # Print simulation overview
# cat("Starting simulation with", n_sim, "iterations per configuration...\n")
# cat("Total configurations:", nrow(simulation_grid), "\n")
# cat("Output will be saved to", rds_dir, "directory\n")
#
# # Set up parallel backend
# n_cores <- max(1, detectCores() - 1)
# cat("Using", n_cores, "cores for parallel processing...\n")
# cl <- makeCluster(n_cores)
# registerDoParallel(cl)
#
# # Export necessary functions and data to the cluster
# clusterExport(cl, c("configurations", "methods", "sample_sizes", "n_sim",
#                     "rds_dir", "create_param_df", "compute_coverage",
#                     "run_simulation", "generate_seed"))
#
# # Load required libraries on the cluster
# clusterEvalQ(cl, {
#   library(tidyverse)
#   library(gkwreg)  # Assumed to contain all GKw distribution functions and fit methods
# })
#
# # Run simulations in parallel
# start_time_all <- Sys.time()
#
# results <- foreach(i = 1:nrow(simulation_grid),
#                    .packages = c("tidyverse", "gkwreg")) %dopar% {
#                      config_idx <- simulation_grid$config[i]
#                      method <- simulation_grid$method[i]
#                      sample_size <- simulation_grid$sample_size[i]
#
#                      # Print progress update (will be seen in log file)
#                      cat("Running simulation for", names(configurations)[config_idx],
#                          "with method", method, "and sample size", sample_size, "...\n")
#
#                      # Run simulation for this configuration
#                      sim_result <- run_simulation(configurations[[config_idx]], method, sample_size, n_sim)
#
#                      # Save individual result to disk
#                      result_filename <- paste0(
#                        "gkw_bias_", simulation_grid$family[i], "_",
#                        simulation_grid$config_name[i], "_",
#                        method, "_n", sample_size, ".rds"
#                      )
#                      saveRDS(sim_result, file.path(rds_dir, result_filename))
#
#                      # Return the result for aggregation
#                      sim_result
#                    }
#
# end_time_all <- Sys.time()
# total_time_all <- difftime(end_time_all, start_time_all, units = "mins")
#
# # Clean up parallel cluster
# stopCluster(cl)
#
# # ---------------------------------------------------------------------------- #
# # Process and Save Aggregated Results
# # ---------------------------------------------------------------------------- #
#
# # Extract parameter results from all simulations
# param_results <- do.call(rbind, lapply(results, function(x) x$param_results))
#
# # Extract convergence rates
# convergence_df <- data.frame(
#   config = simulation_grid$config,
#   method = simulation_grid$method,
#   sample_size = simulation_grid$sample_size,
#   family = simulation_grid$family,
#   config_name = simulation_grid$config_name,
#   convergence_rate = sapply(results, function(x) x$convergence_rate),
#   avg_time = sapply(results, function(x) x$avg_time)
# )
#
# # Save aggregated results to disk
# saveRDS(param_results, file.path(rds_dir, "gkw_bias_param_results.rds"))
# saveRDS(convergence_df, file.path(rds_dir, "gkw_bias_convergence_results.rds"))
#
# # Log completion
# log_file <- file(paste0("gkw_bias_study_log_", run_timestamp, ".txt"), "a")
# writeLines(paste("Simulation completed in", round(total_time_all, 2), "minutes."), log_file)
# writeLines(paste("Results saved to", rds_dir, "directory."), log_file)
# writeLines(paste("Total parameter results:", nrow(param_results)), log_file)
# writeLines(paste("Total convergence results:", nrow(convergence_df)), log_file)
# close(log_file)
#
# # Print summary
# cat("\nSimulation completed in", round(total_time_all, 2), "minutes.\n")
# cat("Results saved to", rds_dir, "directory.\n")
# cat("Total parameter results:", nrow(param_results), "\n")
# cat("Total convergence results:", nrow(convergence_df), "\n")
