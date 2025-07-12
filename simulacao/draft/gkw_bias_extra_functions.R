# ---------------------------------------------------------------------------- #
# Title: Utility Functions for GKw MLE Bias Study
# Description: Contains all functions needed for analysis, visualization, and
#              reporting in the GKw MLE bias simulation study
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Utility Functions
# ---------------------------------------------------------------------------- #

#' Generate a deterministic seed based on configuration
#'
#' @param i Iteration number
#' @param config_name Configuration name
#' @param method Estimation method
#' @param sample_size Sample size
#' @param base_seed Base seed to start from
#' @return Integer seed value
generate_seed <- function(i, config_name, method, sample_size, base_seed = 12345) {
  # Encode method as a number (1-4)
  method_num <- match(method, c("tmb", "nr", "nlminb", "optim"))

  # Encode config_name as a number based on its position
  config_num <- which(names(configurations) == config_name)

  # Create a deterministic but unique seed
  seed <- base_seed + i * 10000 + config_num * 100 + method_num * 10 +
    as.numeric(as.factor(sample_size))

  # Ensure it's within integer limits
  seed <- seed %% 2147483647  # Max integer in R

  return(seed)
}

#' Create a data frame from the parameter results with relative bias
#'
#' @param param_name Parameter name
#' @param param_value True parameter value
#' @param method Estimation method
#' @param family Distribution family
#' @param config Configuration name
#' @param sample_size Sample size
#' @param bias Absolute bias
#' @param se Standard error
#' @param emp_se Empirical standard error
#' @param coverage Coverage probability
#' @return Data frame with all bias metrics
create_param_df <- function(param_name, param_value, method, family, config,
                            sample_size, bias, se, emp_se, coverage) {
  # Calculate bias_z outside the data.frame call
  bias_z_val <- bias / se

  # Calculate relative bias (as percentage of true value)
  rel_bias <- ifelse(param_value != 0,
                     (bias / param_value) * 100,
                     NA)  # Avoid division by zero

  # Calculate relative standard error as percentage of true value
  rel_se <- ifelse(param_value != 0,
                   (se / param_value) * 100,
                   NA)

  # Calculate relative empirical standard error
  rel_emp_se <- ifelse(param_value != 0,
                       (emp_se / param_value) * 100,
                       NA)

  # Create data frame with all metrics
  data.frame(
    parameter = param_name,
    true_value = param_value,
    method = method,
    family = family,
    config = config,
    sample_size = sample_size,
    bias = bias,                    # Absolute bias
    rel_bias = rel_bias,            # Relative bias (%)
    se = se,                        # Standard error
    emp_se = emp_se,                # Empirical standard error
    rel_se = rel_se,                # Relative standard error (%)
    rel_emp_se = rel_emp_se,        # Relative empirical standard error (%)
    coverage = coverage,            # Coverage probability
    bias_z = bias_z_val,            # Standardized bias
    bias_lower_z = bias_z_val - 1,  # Lower bound for standardized bias
    bias_upper_z = bias_z_val + 1,  # Upper bound for standardized bias
    stringsAsFactors = FALSE
  )
}

#' Compute coverage rate
#'
#' @param estimates Vector of parameter estimates
#' @param std_errors Vector of standard errors
#' @param true_value True parameter value
#' @return Coverage rate (proportion of intervals containing the true value)
compute_coverage <- function(estimates, std_errors, true_value) {
  lower_ci <- estimates - qnorm(0.975) * std_errors
  upper_ci <- estimates + qnorm(0.975) * std_errors
  fora <- sum(lower_ci > true_value) + sum(upper_ci < true_value)
  coverage <- 1 - fora / length(estimates)
  return(coverage)
}

# ---------------------------------------------------------------------------- #
# Simulation Function
# ---------------------------------------------------------------------------- #

#' Run a single simulation configuration
#'
#' @param config Configuration list with family and parameters
#' @param method Estimation method
#' @param sample_size Sample size
#' @param n_sim Number of simulation replications
#' @return List with simulation results
run_simulation <- function(config, method, sample_size, n_sim) {

  # Extract configuration details
  family <- config$family
  params <- config$params
  config_name <- names(configurations)[which(sapply(configurations, function(x)
    identical(x$family, family) && identical(x$params, params)))]

  # Function to generate samples for the appropriate family
  generate_sample <- function(n, family, params) {
    # Switch to the right random generation function
    rand_func <- switch(family,
                        "gkw" = rgkw,
                        "bkw" = rbkw,
                        "kkw" = rkkw,
                        "ekw" = rekw,
                        "mc" = rmc,
                        "kw" = rkw,
                        "beta" = rbeta_)

    # Generate random sample
    do.call(rand_func, c(list(n = n), params))
  }

  # Initialize results storage
  n_params <- length(params)
  param_names <- names(params)

  # Create matrices to store results
  estimates <- matrix(NA, nrow = n_sim, ncol = n_params)
  std_errors <- matrix(NA, nrow = n_sim, ncol = n_params)
  converged <- rep(FALSE, n_sim)
  execution_times <- rep(NA, n_sim)

  # Run simulations
  for (i in 1:n_sim) {
    # Set seed for reproducibility (unique to this configuration)
    seed <- generate_seed(i, config_name, method, sample_size)
    set.seed(seed)

    # Generate sample
    data <- generate_sample(sample_size, family, params)

    # Measure execution time
    start_time <- Sys.time()

    # Fit model using specified method
    fit_result <- tryCatch({
      fit <- gkwfit(data = data,
                    family = family,
                    fit = method,
                    hessian = TRUE,
                    silent = TRUE)

      list(success = TRUE, fit = fit)
    },
    error = function(e) {
      return(list(success = FALSE, error = e))
    },
    warning = function(w) {
      return(list(success = FALSE, warning = w))
    })

    # Record execution time
    end_time <- Sys.time()
    execution_times[i] <- as.numeric(difftime(end_time, start_time, units = "secs"))

    # Process results if successful
    if (fit_result$success) {
      fit <- fit_result$fit

      # Check if convergence was achieved
      if (fit$convergence) {
        converged[i] <- TRUE

        # Extract parameter estimates and standard errors
        for (j in 1:n_params) {
          param <- param_names[j]
          if (param %in% names(fit$coefficients)) {
            estimates[i, j] <- fit$coefficients[param]
            std_errors[i, j] <- fit$std.errors[param]
          }
        }
      }
    }
  }

  # Compute summary statistics for each parameter
  results <- list()

  for (j in 1:n_params) {
    param <- param_names[j]
    true_value <- params[[param]]

    # Filter for converged runs
    valid_estimates <- estimates[converged, j]
    valid_std_errors <- std_errors[converged, j]

    # Compute bias and other metrics (if we have enough converged runs)
    if (sum(converged) > 0) {
      bias <- mean(valid_estimates, na.rm = TRUE) - true_value
      se <- mean(valid_std_errors, na.rm = TRUE)
      emp_se <- sd(valid_estimates, na.rm = TRUE)
      coverage <- compute_coverage(valid_estimates, valid_std_errors, true_value)

      # Add results to the list
      results[[param]] <- create_param_df(
        param_name = param,
        param_value = true_value,
        method = method,
        family = family,
        config = config_name,
        sample_size = sample_size,
        bias = bias,
        se = se,
        emp_se = emp_se,
        coverage = coverage
      )
    }
  }

  # Add convergence rate and timing information
  convergence_rate <- mean(converged)
  avg_time <- mean(execution_times[converged], na.rm = TRUE)

  return(list(
    param_results = do.call(rbind, results),
    convergence_rate = convergence_rate,
    avg_time = avg_time,
    raw_estimates = estimates[converged, ],
    raw_std_errors = std_errors[converged, ],
    config_name = config_name,
    family = family,
    method = method,
    sample_size = sample_size
  ))
}

# ---------------------------------------------------------------------------- #
# Visualization Functions
# ---------------------------------------------------------------------------- #

#' Plot absolute bias using forest plots
#'
#' @param data Data frame with simulation results
#' @param param_name Parameter name to plot
#' @param methods_to_show Methods to include in the plot
#' @param families_to_show Families to include in the plot
#' @param configs_to_show Configurations to include in the plot
#' @param limits Y-axis limits
#' @param title Plot title
#' @param xlab X-axis label
#' @return ggplot object
plot_abs_bias <- function(data, param_name, methods_to_show = methods,
                          families_to_show = unique(data$family),
                          configs_to_show = unique(data$config),
                          limits = NULL,
                          title = NULL,
                          xlab = "Sample Size") {

  # Filter data for the specific parameter
  filtered_data <- data %>%
    filter(parameter == param_name,
           method %in% methods_to_show,
           family %in% families_to_show,
           config %in% configs_to_show) %>%
    mutate(
      # Convert to factors for proper ordering
      sample_size = factor(sample_size, levels = sort(unique(sample_size))),
      method = factor(method, levels = methods_to_show),

      # Create nice labels for the plots
      method_label = case_when(
        method == "tmb" ~ "TMB",
        method == "nr" ~ "Newton-Raphson",
        method == "nlminb" ~ "nlminb",
        method == "optim" ~ "optim (L-BFGS-B)",
        TRUE ~ as.character(method)
      ),

      config_label = gsub("_", "-", config),

      # Create LaTeX-formatted parameter name
      param_latex = TeX(paste0("$\\", param_name, "$"), output = "character")
    )

  # Create standard error range for bias
  filtered_data <- filtered_data %>%
    mutate(
      bias_lower = bias - se,
      bias_upper = bias + se
    )

  # Create the ggplot
  g <- ggplot(filtered_data, aes(x = sample_size, y = bias, color = method_label)) +
    geom_pointrange(
      aes(ymin = bias_lower, ymax = bias_upper),
      position = position_dodge(width = 0.5),
      size = 0.2,
      linewidth = 0.3
    ) +
    facet_nested(family + config_label ~ .,
                 scales = "free_y",
                 switch = "y",
                 labeller = labeller(config_label = label_value, family = label_value)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    labs(
      title = ifelse(is.null(title),
                     paste("Absolute bias for parameter", param_name),
                     title),
      y = expression(hat(bias) %+-% SE),
      x = xlab,
      color = "Method"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 7),
      strip.text.x = element_text(size = 10, margin = margin(0.1, 0, 0.1, 0, "cm")),
      strip.text.y = element_text(size = 10, margin = margin(0.1, 0, 0.1, 0, "cm")),
      panel.spacing.x = unit(0.12, "lines"),
      panel.spacing.y = unit(0.12, "lines"),
      plot.title = element_text(hjust = 0.5),
      legend.position = "top"
    ) +
    coord_flip()

  # Add limits if provided
  if (!is.null(limits)) {
    g <- g + scale_y_continuous(limits = limits)
  }

  return(g)
}

#' Plot relative bias using forest plots
#'
#' @param data Data frame with simulation results
#' @param param_name Parameter name to plot
#' @param methods_to_show Methods to include in the plot
#' @param families_to_show Families to include in the plot
#' @param configs_to_show Configurations to include in the plot
#' @param limits Y-axis limits
#' @param title Plot title
#' @param xlab X-axis label
#' @return ggplot object
plot_rel_bias <- function(data, param_name, methods_to_show = methods,
                          families_to_show = unique(data$family),
                          configs_to_show = unique(data$config),
                          limits = NULL,
                          title = NULL,
                          xlab = "Sample Size") {

  # Filter data for the specific parameter
  filtered_data <- data %>%
    filter(parameter == param_name,
           method %in% methods_to_show,
           family %in% families_to_show,
           config %in% configs_to_show) %>%
    mutate(
      # Convert to factors for proper ordering
      sample_size = factor(sample_size, levels = sort(unique(sample_size))),
      method = factor(method, levels = methods_to_show),

      # Create nice labels for the plots
      method_label = case_when(
        method == "tmb" ~ "TMB",
        method == "nr" ~ "Newton-Raphson",
        method == "nlminb" ~ "nlminb",
        method == "optim" ~ "optim (L-BFGS-B)",
        TRUE ~ as.character(method)
      ),

      config_label = gsub("_", "-", config),

      # Create LaTeX-formatted parameter name
      param_latex = TeX(paste0("$\\", param_name, "$"), output = "character")
    )

  # Calculate relative bias confidence intervals
  # using error propagation for the ratio: Var(X/Y) ≈ (X/Y)^2 * [Var(X)/X^2 + Var(Y)/Y^2]
  filtered_data <- filtered_data %>%
    mutate(
      # For relative bias, we need the relative standard error
      rel_bias_lower = rel_bias - rel_se,
      rel_bias_upper = rel_bias + rel_se
    )

  # Create the ggplot for relative bias
  g <- ggplot(filtered_data, aes(x = sample_size, y = rel_bias, color = method_label)) +
    geom_pointrange(
      aes(ymin = rel_bias_lower, ymax = rel_bias_upper),
      position = position_dodge(width = 0.5),
      size = 0.2,
      linewidth = 0.3
    ) +
    facet_nested(family + config_label ~ .,
                 scales = "free_y",
                 switch = "y",
                 labeller = labeller(config_label = label_value, family = label_value)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    labs(
      title = ifelse(is.null(title),
                     paste("Relative bias (%) for parameter", param_name),
                     title),
      y = "Relative Bias (%) ± SE",
      x = xlab,
      color = "Method"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 7),
      strip.text.x = element_text(size = 10, margin = margin(0.1, 0, 0.1, 0, "cm")),
      strip.text.y = element_text(size = 10, margin = margin(0.1, 0, 0.1, 0, "cm")),
      panel.spacing.x = unit(0.12, "lines"),
      panel.spacing.y = unit(0.12, "lines"),
      plot.title = element_text(hjust = 0.5),
      legend.position = "top"
    ) +
    coord_flip()

  # Add limits if provided
  if (!is.null(limits)) {
    g <- g + scale_y_continuous(limits = limits)
  }

  return(g)
}

# ---------------------------------------------------------------------------- #
# New Parameter-Focused Visualization Functions
# ---------------------------------------------------------------------------- #

#' Plot bias with format matching the reference image
#'
#' @param data Data frame with simulation results
#' @param bias_type Type of bias to plot ("absolute" or "relative")
#' @param methods_to_show Methods to include in the plot (columns)
#' @param families_to_show Families to include in the plot (rows)
#' @param params_to_show Parameters to include in the plot (rows, combined with family)
#' @param limits X-axis limits
#' @param title Plot title
#' @return ggplot object
plot_bias_matrix <- function(data,
                             bias_type = "absolute",
                             methods_to_show = methods,
                             families_to_show = unique(data$family),
                             params_to_show = unique(data$parameter),
                             limits = NULL,
                             title = NULL) {

  # Filter data for the specified parameters and families
  filtered_data <- data %>%
    filter(parameter %in% params_to_show,
           method %in% methods_to_show,
           family %in% families_to_show) %>%
    mutate(
      # Convert sample size to factor for ordering on y-axis
      sample_size = factor(sample_size, levels = sort(unique(as.numeric(as.character(sample_size))))),

      # Create nice labels for the plots
      method_label = case_when(
        method == "tmb" ~ "TMB",
        method == "nr" ~ "Newton-Raphson",
        method == "nlminb" ~ "nlminb",
        method == "optim" ~ "optim (L-BFGS-B)",
        TRUE ~ as.character(method)
      ),

      # Create combined family-parameter label
      family_param = paste(family, parameter, sep = " - ")
    )

  # Set up bias variables based on bias type
  if (bias_type == "relative") {
    filtered_data <- filtered_data %>%
      mutate(
        bias_value = rel_bias,
        bias_lower = rel_bias - rel_se,
        bias_upper = rel_bias + rel_se
      )
    x_label <- "Relative Bias (%) ± SE"
  } else {
    filtered_data <- filtered_data %>%
      mutate(
        bias_value = bias,
        bias_lower = bias - se,
        bias_upper = bias + se
      )
    x_label <- expression(hat(bias) %+-% SE)
  }

  # Create the ggplot with family-parameter combinations in rows and methods in columns
  g <- ggplot(filtered_data,
              aes(x = bias_value, y = sample_size)) +
    geom_vline(xintercept = 0, linetype = "dotted") +
    geom_pointrange(
      aes(xmin = bias_lower, xmax = bias_upper),
      size = 0.3,
      linewidth = 0.3
    ) +
    facet_grid(family_param ~ method_label,
               scales = "free_y",
               switch = "y") +
    labs(
      title = title,
      x = x_label,
      y = "Sample Size"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 7),
      strip.text.x = element_text(size = 9),
      strip.text.y = element_text(size = 9),
      panel.spacing = unit(0.2, "lines"),
      plot.title = element_text(hjust = 0.5),
      legend.position = "none"
    )

  # Add limits if provided
  if (!is.null(limits)) {
    g <- g + scale_x_continuous(limits = limits)
  }

  return(g)
}

#' Plot parameter absolute bias by family and method using the new matrix format
#'
#' @param data Data frame with simulation results
#' @param param_name Parameter name to plot
#' @param methods_to_show Methods to include in the plot (columns)
#' @param families_to_show Families to include in the plot (rows)
#' @param limits X-axis limits
#' @param title Plot title
#' @return ggplot object
plot_parameter_bias_matrix <- function(data,
                                       param_name,
                                       methods_to_show = methods,
                                       families_to_show = unique(data$family),
                                       limits = NULL,
                                       title = NULL) {

  # Filter data for the specific parameter
  filtered_data <- data %>%
    filter(parameter == param_name)

  # Create the matrix plot
  plot_bias_matrix(
    filtered_data,
    bias_type = "absolute",
    methods_to_show = methods_to_show,
    families_to_show = families_to_show,
    params_to_show = param_name,
    limits = limits,
    title = ifelse(is.null(title),
                   paste("Absolute bias for parameter", param_name),
                   title)
  )
}

#' Plot parameter relative bias by family and method using the new matrix format
#'
#' @param data Data frame with simulation results
#' @param param_name Parameter name to plot
#' @param methods_to_show Methods to include in the plot (columns)
#' @param families_to_show Families to include in the plot (rows)
#' @param limits X-axis limits
#' @param title Plot title
#' @return ggplot object
plot_parameter_rel_bias_matrix <- function(data,
                                           param_name,
                                           methods_to_show = methods,
                                           families_to_show = unique(data$family),
                                           limits = NULL,
                                           title = NULL) {

  # Filter data for the specific parameter
  filtered_data <- data %>%
    filter(parameter == param_name)

  # Create the matrix plot
  plot_bias_matrix(
    filtered_data,
    bias_type = "relative",
    methods_to_show = methods_to_show,
    families_to_show = families_to_show,
    params_to_show = param_name,
    limits = limits,
    title = ifelse(is.null(title),
                   paste("Relative bias (%) for parameter", param_name),
                   title)
  )
}

#' Plot parameter bias across sample sizes by method
#'
#' @param data Data frame with simulation results
#' @param param_name Parameter name to plot
#' @param methods_to_show Methods to include in the plot
#' @param families_to_show Families to include in the plot
#' @param bias_type Type of bias to plot ('absolute' or 'relative')
#' @param limits Y-axis limits
#' @param title Plot title
#' @return ggplot object
plot_parameter_bias_by_sample_size <- function(data, param_name,
                                               methods_to_show = methods,
                                               families_to_show = unique(data$family),
                                               bias_type = "relative",
                                               limits = NULL,
                                               title = NULL) {

  # Filter data for the specific parameter
  filtered_data <- data %>%
    filter(parameter == param_name,
           method %in% methods_to_show,
           family %in% families_to_show) %>%
    mutate(
      # Convert sample size to numeric for proper plotting
      sample_size_num = as.numeric(as.character(sample_size)),

      # Convert to factors for proper ordering
      method = factor(method, levels = methods_to_show),
      family = factor(family, levels = families_to_show),

      # Create nice labels for the plots
      method_label = case_when(
        method == "tmb" ~ "TMB",
        method == "nr" ~ "Newton-Raphson",
        method == "nlminb" ~ "nlminb",
        method == "optim" ~ "optim (L-BFGS-B)",
        TRUE ~ as.character(method)
      ),

      # Create configuration label
      config_label = gsub("_", "-", config),

      # Create LaTeX-formatted parameter name for title
      param_latex = TeX(paste0("$\\", param_name, "$"), output = "character")
    )

  # Set up the plot based on bias type
  if (bias_type == "relative") {
    y_var <- "rel_bias"
    y_lab <- "Relative Bias (%)"
    default_title <- paste("Relative bias trends for parameter", param_name)
  } else {
    y_var <- "bias"
    y_lab <- "Absolute Bias"
    default_title <- paste("Absolute bias trends for parameter", param_name)
  }

  # Create the ggplot showing bias trends by sample size
  g <- ggplot(filtered_data, aes_string(x = "sample_size_num", y = y_var, color = "method_label")) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_line(linewidth = 0.5) +
    geom_point(size = 1.5) +
    facet_grid(family ~ config_label, scales = "free_y") +
    scale_x_log10(breaks = unique(filtered_data$sample_size_num),
                  labels = unique(filtered_data$sample_size_num)) +
    labs(
      title = ifelse(is.null(title), default_title, title),
      x = "Sample Size (log scale)",
      y = y_lab,
      color = "Method"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 7),
      strip.text = element_text(size = 9),
      panel.spacing = unit(0.2, "lines"),
      plot.title = element_text(hjust = 0.5),
      legend.position = "top"
    )

  # Add limits if provided
  if (!is.null(limits)) {
    g <- g + scale_y_continuous(limits = limits)
  }

  return(g)
}

#' Create a comparison heatmap of bias measures across parameters and methods
#'
#' @param data Data frame with parameter results
#' @param families_to_show Families to include
#' @param methods_to_show Methods to include
#' @param sample_size_to_show Sample size to show (default: largest)
#' @param measure Bias measure to use ('rel_bias', 'bias_z', or 'bias')
#' @param use_abs Whether to use absolute values of the measure
#' @param title Plot title
#' @return ggplot object
plot_parameter_bias_heatmap <- function(data,
                                        families_to_show = unique(data$family),
                                        methods_to_show = methods,
                                        sample_size_to_show = NULL,
                                        measure = "rel_bias",
                                        use_abs = TRUE,
                                        title = NULL) {

  # If no sample size provided, use largest available
  if (is.null(sample_size_to_show)) {
    sample_size_to_show <- max(as.numeric(as.character(unique(data$sample_size))))
  }

  # Filter data
  filtered_data <- data %>%
    filter(method %in% methods_to_show,
           family %in% families_to_show,
           sample_size == sample_size_to_show)

  # Apply abs function if requested
  if (use_abs) {
    filtered_data[[measure]] <- abs(filtered_data[[measure]])
  }

  # Set up measure label
  measure_label <- case_when(
    measure == "rel_bias" & use_abs ~ "Mean Absolute Relative Bias (%)",
    measure == "rel_bias" & !use_abs ~ "Mean Relative Bias (%)",
    measure == "bias_z" & use_abs ~ "Mean Absolute Standardized Bias",
    measure == "bias_z" & !use_abs ~ "Mean Standardized Bias",
    measure == "bias" & use_abs ~ "Mean Absolute Bias",
    measure == "bias" & !use_abs ~ "Mean Bias",
    TRUE ~ "Bias Measure"
  )

  # Group and summarize data
  heatmap_data <- filtered_data %>%
    group_by(parameter, family, method) %>%
    summarize(
      mean_measure = mean(!!sym(measure), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # Create method labels
    mutate(
      method_label = case_when(
        method == "tmb" ~ "TMB",
        method == "nr" ~ "Newton-Raphson",
        method == "nlminb" ~ "nlminb",
        method == "optim" ~ "optim (L-BFGS-B)",
        TRUE ~ as.character(method)
      )
    )

  # Create the heatmap
  g <- ggplot(heatmap_data, aes(x = method_label, y = family, fill = mean_measure)) +
    geom_tile() +
    facet_wrap(~ parameter, scales = "free") +
    scale_fill_viridis_c(option = "plasma", name = measure_label) +
    labs(
      title = ifelse(is.null(title),
                     paste(measure_label, "by Parameter, Family, and Method (n =", sample_size_to_show, ")"),
                     title),
      x = "Method",
      y = "Distribution Family"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, size = 12)
    )

  return(g)
}

#' Plot standardized bias (bias/se) using forest plots
#'
#' @param data Data frame with simulation results
#' @param param_name Parameter name to plot
#' @param methods_to_show Methods to include in the plot
#' @param families_to_show Families to include in the plot
#' @param configs_to_show Configurations to include in the plot
#' @param limits Y-axis limits
#' @param title Plot title
#' @param xlab X-axis label
#' @return ggplot object
plot_bias_z <- function(data, param_name, methods_to_show = methods,
                        families_to_show = unique(data$family),
                        configs_to_show = unique(data$config),
                        limits = c(-2, 2),
                        title = NULL,
                        xlab = "Sample Size") {

  # Filter data for the specific parameter
  filtered_data <- data %>%
    filter(parameter == param_name,
           method %in% methods_to_show,
           family %in% families_to_show,
           config %in% configs_to_show) %>%
    mutate(
      # Convert to factors for proper ordering
      sample_size = factor(sample_size, levels = sort(unique(sample_size))),
      method = factor(method, levels = methods_to_show),

      # Create nice labels for the plots
      method_label = case_when(
        method == "tmb" ~ "TMB",
        method == "nr" ~ "Newton-Raphson",
        method == "nlminb" ~ "nlminb",
        method == "optim" ~ "optim (L-BFGS-B)",
        TRUE ~ as.character(method)
      ),

      config_label = gsub("_", "-", config),

      # Create LaTeX-formatted parameter name
      param_latex = TeX(paste0("$\\", param_name, "$"), output = "character")
    )

  # Create the ggplot
  g <- ggplot(filtered_data, aes(x = sample_size, y = bias_z, color = method_label)) +
    geom_pointrange(
      aes(ymin = bias_lower_z, ymax = bias_upper_z),
      position = position_dodge(width = 0.5),
      size = 0.2,
      linewidth = 0.3
    ) +
    facet_nested(family + config_label ~ .,
                 scales = "free_y",
                 switch = "y",
                 labeller = labeller(config_label = label_value, family = label_value)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    labs(
      title = ifelse(is.null(title),
                     paste("Standardized bias for parameter", param_name),
                     title),
      y = expression(hat(bias) / SE %+-% 1),
      x = xlab,
      color = "Method"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 7),
      strip.text.x = element_text(size = 10, margin = margin(0.1, 0, 0.1, 0, "cm")),
      strip.text.y = element_text(size = 10, margin = margin(0.1, 0, 0.1, 0, "cm")),
      panel.spacing.x = unit(0.12, "lines"),
      panel.spacing.y = unit(0.12, "lines"),
      plot.title = element_text(hjust = 0.5),
      legend.position = "top"
    ) +
    coord_flip()

  # Add limits if provided
  if (!is.null(limits)) {
    g <- g + scale_y_continuous(limits = limits)
  }

  return(g)
}

#' Plot parameter standardized bias by family and method
#'
#' @param data Data frame with simulation results
#' @param param_name Parameter name to plot
#' @param methods_to_show Methods to include in the plot
#' @param families_to_show Families to include in the plot
#' @param sample_size_to_show Sample size to show (default: largest)
#' @param limits X-axis limits
#' @param title Plot title
#' @return ggplot object
plot_parameter_bias_z_by_family <- function(data, param_name,
                                            methods_to_show = methods,
                                            families_to_show = unique(data$family),
                                            sample_size_to_show = NULL,
                                            limits = c(-2, 2),
                                            title = NULL) {

  # If no sample size provided, use largest available
  if (is.null(sample_size_to_show)) {
    sample_size_to_show <- max(as.numeric(as.character(unique(data$sample_size))))
  }

  # Filter data for the specific parameter and sample size
  filtered_data <- data %>%
    filter(parameter == param_name,
           method %in% methods_to_show,
           family %in% families_to_show,
           sample_size == sample_size_to_show) %>%
    mutate(
      # Convert to factors for proper ordering
      method = factor(method, levels = methods_to_show),
      family = factor(family, levels = families_to_show),

      # Create nice labels for the plots
      method_label = case_when(
        method == "tmb" ~ "TMB",
        method == "nr" ~ "Newton-Raphson",
        method == "nlminb" ~ "nlminb",
        method == "optim" ~ "optim (L-BFGS-B)",
        TRUE ~ as.character(method)
      ),

      # Create configuration label
      config_label = gsub("_", "-", config),

      # Create LaTeX-formatted parameter name for title
      param_latex = TeX(paste0("$\\", param_name, "$"), output = "character")
    )

  # Create the ggplot with parameters in rows
  g <- ggplot(filtered_data, aes(x = bias_z, y = config_label, color = method_label)) +
    geom_vline(xintercept = 0, linetype = "dotted") +
    geom_pointrange(
      aes(xmin = bias_lower_z, xmax = bias_upper_z),
      position = position_dodge(width = 0.5),
      size = 0.3,
      linewidth = 0.4
    ) +
    facet_grid(family ~ method_label, scales = "free_y", space = "free_y") +
    labs(
      title = ifelse(is.null(title),
                     paste("Standardized bias for parameter", param_name, "(n =", sample_size_to_show, ")"),
                     title),
      x = expression(hat(bias) / SE %+-% 1),
      y = "Configuration",
      color = "Method"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 7),
      strip.text = element_text(size = 10),
      panel.spacing = unit(0.2, "lines"),
      plot.title = element_text(hjust = 0.5),
      legend.position = "top"
    )

  # Add limits if provided
  if (!is.null(limits)) {
    g <- g + scale_x_continuous(limits = limits)
  }

  return(g)
}

#' Plot coverage rates
#'
#' @param data Data frame with simulation results
#' @param param_name Parameter name to plot
#' @param methods_to_show Methods to include in the plot
#' @param families_to_show Families to include in the plot
#' @param configs_to_show Configurations to include in the plot
#' @param target_coverage Target coverage rate (usually 0.95)
#' @param limits Y-axis limits
#' @param title Plot title
#' @param xlab X-axis label
#' @return ggplot object
plot_coverage <- function(data, param_name, methods_to_show = methods,
                          families_to_show = unique(data$family),
                          configs_to_show = unique(data$config),
                          target_coverage = 0.95,
                          limits = c(0.7, 1.0),
                          title = NULL,
                          xlab = "Sample Size") {

  # Filter data for the specific parameter
  filtered_data <- data %>%
    filter(parameter == param_name,
           method %in% methods_to_show,
           family %in% families_to_show,
           config %in% configs_to_show) %>%
    mutate(
      # Convert to factors for proper ordering
      sample_size = factor(sample_size, levels = sort(unique(sample_size))),
      method = factor(method, levels = methods_to_show),

      # Create nice labels for the plots
      method_label = case_when(
        method == "tmb" ~ "TMB",
        method == "nr" ~ "Newton-Raphson",
        method == "nlminb" ~ "nlminb",
        method == "optim" ~ "optim (L-BFGS-B)",
        TRUE ~ as.character(method)
      ),

      config_label = gsub("_", "-", config),

      # Create LaTeX-formatted parameter name
      param_latex = TeX(paste0("$\\", param_name, "$"), output = "character")
    )

  # Create the coverage plot
  g <- ggplot(filtered_data,
              aes(x = sample_size, y = coverage,
                  group = interaction(method_label, family),
                  color = method_label)) +
    geom_line(linewidth = 0.3) +
    geom_point(aes(shape = method_label)) +
    facet_nested(family + config_label ~ .,
                 switch = "y",
                 labeller = labeller(config_label = label_value, family = label_value)) +
    geom_hline(yintercept = target_coverage, linewidth = 0.2, linetype = "dashed") +
    labs(
      title = ifelse(is.null(title),
                     paste("Coverage rate for parameter", param_name),
                     title),
      y = "Coverage Rate",
      x = xlab,
      color = "Method",
      shape = "Method"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 7),
      strip.text.x = element_text(size = 10, margin = margin(0.1, 0, 0.1, 0, "cm")),
      strip.text.y = element_text(size = 10, margin = margin(0.1, 0, 0.1, 0, "cm")),
      panel.spacing.x = unit(0.12, "lines"),
      panel.spacing.y = unit(0.12, "lines"),
      plot.title = element_text(hjust = 0.5),
      legend.position = "top"
    ) +
    scale_y_continuous(limits = limits) +
    coord_flip()

  return(g)
}

#' Plot parameter coverage by family and method
#'
#' @param data Data frame with simulation results
#' @param param_name Parameter name to plot
#' @param methods_to_show Methods to include in the plot
#' @param families_to_show Families to include in the plot
#' @param sample_size_to_show Sample size to show (default: largest)
#' @param target_coverage Target coverage rate (usually 0.95)
#' @param limits X-axis limits
#' @param title Plot title
#' @return ggplot object
plot_parameter_coverage_by_family <- function(data, param_name,
                                              methods_to_show = methods,
                                              families_to_show = unique(data$family),
                                              sample_size_to_show = NULL,
                                              target_coverage = 0.95,
                                              limits = c(0.7, 1.0),
                                              title = NULL) {

  # If no sample size provided, use largest available
  if (is.null(sample_size_to_show)) {
    sample_size_to_show <- max(as.numeric(as.character(unique(data$sample_size))))
  }

  # Filter data for the specific parameter and sample size
  filtered_data <- data %>%
    filter(parameter == param_name,
           method %in% methods_to_show,
           family %in% families_to_show,
           sample_size == sample_size_to_show) %>%
    mutate(
      # Convert to factors for proper ordering
      method = factor(method, levels = methods_to_show),
      family = factor(family, levels = families_to_show),

      # Create nice labels for the plots
      method_label = case_when(
        method == "tmb" ~ "TMB",
        method == "nr" ~ "Newton-Raphson",
        method == "nlminb" ~ "nlminb",
        method == "optim" ~ "optim (L-BFGS-B)",
        TRUE ~ as.character(method)
      ),

      # Create configuration label
      config_label = gsub("_", "-", config),

      # Create LaTeX-formatted parameter name for title
      param_latex = TeX(paste0("$\\", param_name, "$"), output = "character")
    )

  # Create the ggplot with parameters in rows
  g <- ggplot(filtered_data, aes(x = coverage, y = config_label, color = method_label)) +
    geom_vline(xintercept = target_coverage, linetype = "dashed") +
    geom_point(size = 3) +
    facet_grid(family ~ method_label, scales = "free_y", space = "free_y") +
    labs(
      title = ifelse(is.null(title),
                     paste("Coverage rate for parameter", param_name, "(n =", sample_size_to_show, ")"),
                     title),
      x = "Coverage Rate",
      y = "Configuration",
      color = "Method"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 7),
      strip.text = element_text(size = 10),
      panel.spacing = unit(0.2, "lines"),
      plot.title = element_text(hjust = 0.5),
      legend.position = "top"
    )

  # Add limits if provided
  if (!is.null(limits)) {
    g <- g + scale_x_continuous(limits = limits)
  }

  return(g)
}

#' Plot convergence rates
#'
#' @param convergence_data Data frame with convergence results
#' @param methods_to_show Methods to include in the plot
#' @param families_to_show Families to include in the plot
#' @param configs_to_show Configurations to include in the plot
#' @param limits Y-axis limits
#' @param title Plot title
#' @param xlab X-axis label
#' @return ggplot object
plot_convergence <- function(convergence_data,
                             methods_to_show = methods,
                             families_to_show = unique(convergence_data$family),
                             configs_to_show = unique(convergence_data$config_name),
                             limits = c(0, 1),
                             title = "Convergence Rates",
                             xlab = "Sample Size") {

  # Filter data
  filtered_data <- convergence_data %>%
    filter(method %in% methods_to_show,
           family %in% families_to_show,
           config_name %in% configs_to_show) %>%
    mutate(
      # Convert to factors for proper ordering
      sample_size = factor(sample_size, levels = sort(unique(sample_size))),
      method = factor(method, levels = methods_to_show),

      # Create nice labels for the plots
      method_label = case_when(
        method == "tmb" ~ "TMB",
        method == "nr" ~ "Newton-Raphson",
        method == "nlminb" ~ "nlminb",
        method == "optim" ~ "optim (L-BFGS-B)",
        TRUE ~ as.character(method)
      ),

      config_label = gsub("_", " ", config_name)
    )

  # Create the convergence plot
  g <- ggplot(filtered_data,
              aes(x = sample_size, y = convergence_rate,
                  group = interaction(method_label, family),
                  color = method_label)) +
    geom_line(linewidth = 0.3) +
    geom_point(aes(shape = method_label)) +
    facet_nested(family + config_label ~ .,
                 switch = "y",
                 labeller = labeller(config_label = label_value, family = label_value)) +
    labs(
      title = title,
      y = "Convergence Rate",
      x = xlab,
      color = "Method",
      shape = "Method"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 7),
      strip.text.x = element_text(size = 10, margin = margin(0.1, 0, 0.1, 0, "cm")),
      strip.text.y = element_text(size = 10, margin = margin(0.1, 0, 0.1, 0, "cm")),
      panel.spacing.x = unit(0.12, "lines"),
      panel.spacing.y = unit(0.12, "lines"),
      plot.title = element_text(hjust = 0.5),
      legend.position = "top"
    ) +
    scale_y_continuous(limits = limits) +
    coord_flip()

  return(g)
}

#' Plot convergence by method and family
#'
#' @param convergence_data Data frame with convergence results
#' @param methods_to_show Methods to include in the plot
#' @param families_to_show Families to include in the plot
#' @param sample_size_to_show Sample size to show (default: smallest)
#' @param limits X-axis limits
#' @param title Plot title
#' @return ggplot object
plot_convergence_by_method_family <- function(convergence_data,
                                              methods_to_show = methods,
                                              families_to_show = unique(convergence_data$family),
                                              sample_size_to_show = NULL,
                                              limits = c(0, 1),
                                              title = NULL) {

  # If no sample size provided, use smallest available (typically challenging convergence)
  if (is.null(sample_size_to_show)) {
    sample_size_to_show <- min(as.numeric(as.character(unique(convergence_data$sample_size))))
  }

  # Filter data
  filtered_data <- convergence_data %>%
    filter(method %in% methods_to_show,
           family %in% families_to_show,
           sample_size == sample_size_to_show) %>%
    mutate(
      # Convert to factors for proper ordering
      method = factor(method, levels = methods_to_show),

      # Create nice labels for the plots
      method_label = case_when(
        method == "tmb" ~ "TMB",
        method == "nr" ~ "Newton-Raphson",
        method == "nlminb" ~ "nlminb",
        method == "optim" ~ "optim (L-BFGS-B)",
        TRUE ~ as.character(method)
      ),

      config_label = gsub("_", " ", config_name)
    )

  # Create the convergence plot
  g <- ggplot(filtered_data, aes(x = convergence_rate, y = config_label, fill = method_label)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_grid(family ~ ., scales = "free_y", space = "free_y") +
    labs(
      title = ifelse(is.null(title),
                     paste("Convergence rates by method and family (n =", sample_size_to_show, ")"),
                     title),
      x = "Convergence Rate",
      y = "Configuration",
      fill = "Method"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 7),
      strip.text = element_text(size = 10),
      panel.spacing = unit(0.2, "lines"),
      plot.title = element_text(hjust = 0.5),
      legend.position = "top"
    )

  # Add limits if provided
  if (!is.null(limits)) {
    g <- g + scale_x_continuous(limits = limits)
  }

  return(g)
}

#' Plot computation time
#'
#' @param convergence_data Data frame with convergence results
#' @param methods_to_show Methods to include in the plot
#' @param families_to_show Families to include in the plot
#' @param configs_to_show Configurations to include in the plot
#' @param log_scale Whether to use log scale for y-axis
#' @param title Plot title
#' @param xlab X-axis label
#' @return ggplot object
plot_computation_time <- function(convergence_data,
                                  methods_to_show = methods,
                                  families_to_show = unique(convergence_data$family),
                                  configs_to_show = unique(convergence_data$config_name),
                                  log_scale = TRUE,
                                  title = "Computation Time",
                                  xlab = "Sample Size") {

  # Filter data
  filtered_data <- convergence_data %>%
    filter(method %in% methods_to_show,
           family %in% families_to_show,
           config_name %in% configs_to_show) %>%
    mutate(
      # Convert to factors for proper ordering
      sample_size = factor(sample_size, levels = sort(unique(sample_size))),
      method = factor(method, levels = methods_to_show),

      # Create nice labels for the plots
      method_label = case_when(
        method == "tmb" ~ "TMB",
        method == "nr" ~ "Newton-Raphson",
        method == "nlminb" ~ "nlminb",
        method == "optim" ~ "optim (L-BFGS-B)",
        TRUE ~ as.character(method)
      ),

      config_label = gsub("_", " ", config_name)
    )

  # Create the computation time plot
  g <- ggplot(filtered_data,
              aes(x = sample_size, y = avg_time,
                  group = interaction(method_label, family),
                  color = method_label)) +
    geom_line(linewidth = 0.3) +
    geom_point(aes(shape = method_label)) +
    facet_nested(family + config_label ~ .,
                 switch = "y",
                 labeller = labeller(config_label = label_value, family = label_value)) +
    labs(
      title = title,
      y = ifelse(log_scale, "Computation Time (log seconds)", "Computation Time (seconds)"),
      x = xlab,
      color = "Method",
      shape = "Method"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 7),
      strip.text.x = element_text(size = 10, margin = margin(0.1, 0, 0.1, 0, "cm")),
      strip.text.y = element_text(size = 10, margin = margin(0.1, 0, 0.1, 0, "cm")),
      panel.spacing.x = unit(0.12, "lines"),
      panel.spacing.y = unit(0.12, "lines"),
      plot.title = element_text(hjust = 0.5),
      legend.position = "top"
    ) +
    coord_flip()

  # Apply log scale if requested
  if (log_scale) {
    g <- g + scale_y_log10()
  }

  return(g)
}

#' Plot comparison of computation time by method and family
#'
#' @param convergence_data Data frame with convergence results
#' @param methods_to_show Methods to include in the plot
#' @param families_to_show Families to include in the plot
#' @param sample_size_to_show Sample size to show (default: largest)
#' @param log_scale Whether to use log scale for x-axis
#' @param title Plot title
#' @return ggplot object
plot_computation_time_by_method_family <- function(convergence_data,
                                                   methods_to_show = methods,
                                                   families_to_show = unique(convergence_data$family),
                                                   sample_size_to_show = NULL,
                                                   log_scale = TRUE,
                                                   title = NULL) {

  # If no sample size provided, use largest available
  if (is.null(sample_size_to_show)) {
    sample_size_to_show <- max(as.numeric(as.character(unique(convergence_data$sample_size))))
  }

  # Filter data
  filtered_data <- convergence_data %>%
    filter(method %in% methods_to_show,
           family %in% families_to_show,
           sample_size == sample_size_to_show) %>%
    mutate(
      # Convert to factors for proper ordering
      method = factor(method, levels = methods_to_show),

      # Create nice labels for the plots
      method_label = case_when(
        method == "tmb" ~ "TMB",
        method == "nr" ~ "Newton-Raphson",
        method == "nlminb" ~ "nlminb",
        method == "optim" ~ "optim (L-BFGS-B)",
        TRUE ~ as.character(method)
      ),

      config_label = gsub("_", " ", config_name)
    )

  # Create the computation time plot
  g <- ggplot(filtered_data, aes(x = avg_time, y = config_label, fill = method_label)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_grid(family ~ ., scales = "free_y", space = "free_y") +
    labs(
      title = ifelse(is.null(title),
                     paste("Computation time by method and family (n =", sample_size_to_show, ")"),
                     title),
      x = ifelse(log_scale, "Computation Time (log seconds)", "Computation Time (seconds)"),
      y = "Configuration",
      fill = "Method"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 7),
      strip.text = element_text(size = 10),
      panel.spacing = unit(0.2, "lines"),
      plot.title = element_text(hjust = 0.5),
      legend.position = "top"
    )

  # Apply log scale if requested
  if (log_scale) {
    g <- g + scale_x_log10()
  }

  return(g)
}

# ---------------------------------------------------------------------------- #
# Analysis Functions
# ---------------------------------------------------------------------------- #

#' Create bias summary tables
#'
#' @param param_results Data frame with parameter results
#' @return Data frame with bias summary
create_bias_summary <- function(param_results) {
  bias_summary <- param_results %>%
    group_by(parameter, family, method, sample_size) %>%
    summarize(
      mean_bias = mean(bias, na.rm = TRUE),
      mean_rel_bias = mean(rel_bias, na.rm = TRUE),
      mean_bias_z = mean(bias_z, na.rm = TRUE),
      mean_true_value = mean(true_value, na.rm = TRUE),
      min_rel_bias = min(rel_bias, na.rm = TRUE),
      max_rel_bias = max(rel_bias, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(parameter, family, method, sample_size)

  return(bias_summary)
}

#' Create coverage summary tables
#'
#' @param param_results Data frame with parameter results
#' @return Data frame with coverage summary
create_coverage_summary <- function(param_results) {
  coverage_summary <- param_results %>%
    group_by(parameter, family, method, sample_size) %>%
    summarize(
      mean_coverage = mean(coverage, na.rm = TRUE),
      min_coverage = min(coverage, na.rm = TRUE),
      max_coverage = max(coverage, na.rm = TRUE),
      target_coverage = 0.95,
      coverage_diff = mean_coverage - 0.95,
      .groups = "drop"
    ) %>%
    arrange(parameter, family, method, sample_size)

  return(coverage_summary)
}

#' Create convergence and timing summary
#'
#' @param convergence_df Data frame with convergence results
#' @return Data frame with convergence and timing summary
create_convergence_summary <- function(convergence_df) {
  convergence_summary <- convergence_df %>%
    group_by(family, method, sample_size) %>%
    summarize(
      mean_convergence = mean(convergence_rate, na.rm = TRUE),
      min_convergence = min(convergence_rate, na.rm = TRUE),
      max_convergence = max(convergence_rate, na.rm = TRUE),
      mean_time = mean(avg_time, na.rm = TRUE),
      min_time = min(avg_time, na.rm = TRUE),
      max_time = max(avg_time, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(family, method, sample_size)

  return(convergence_summary)
}

#' Compare methods
#'
#' @param param_results Data frame with parameter results
#' @param convergence_df Data frame with convergence results
#' @return Data frame with method comparison
create_method_comparison <- function(param_results, convergence_df) {
  method_comparison <- param_results %>%
    group_by(method, parameter, sample_size) %>%
    summarize(
      mean_abs_bias = mean(abs(bias), na.rm = TRUE),
      mean_abs_rel_bias = mean(abs(rel_bias), na.rm = TRUE),
      mean_abs_bias_z = mean(abs(bias_z), na.rm = TRUE),
      mean_coverage = mean(coverage, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(
      convergence_df %>%
        group_by(method, sample_size) %>%
        summarize(
          mean_convergence = mean(convergence_rate, na.rm = TRUE),
          mean_time = mean(avg_time, na.rm = TRUE),
          .groups = "drop"
        ),
      by = c("method", "sample_size")
    ) %>%
    arrange(parameter, method, sample_size)

  return(method_comparison)
}

#' Rank methods by performance
#'
#' @param param_data Data frame with parameter results
#' @param convergence_data Data frame with convergence results
#' @param rel_bias_weight Weight for relative bias in overall score
#' @param coverage_weight Weight for coverage in overall score
#' @param convergence_weight Weight for convergence in overall score
#' @param time_weight Weight for time in overall score
#' @return Data frame with method rankings
rank_methods <- function(param_data, convergence_data,
                         rel_bias_weight = 0.4,
                         coverage_weight = 0.3,
                         convergence_weight = 0.2,
                         time_weight = 0.1) {

  # Normalize and score relative bias (lower absolute value is better)
  bias_scores <- param_data %>%
    group_by(method, parameter, family, config) %>%
    summarize(
      mean_abs_rel_bias = mean(abs(rel_bias), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(parameter, family, config) %>%
    mutate(
      normalized_bias_score = 1 - (mean_abs_rel_bias / max(mean_abs_rel_bias, na.rm = TRUE))
    ) %>%
    ungroup()

  # Normalize and score coverage (closer to 0.95 is better)
  coverage_scores <- param_data %>%
    group_by(method, parameter, family, config) %>%
    summarize(
      mean_coverage = mean(coverage, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(parameter, family, config) %>%
    mutate(
      normalized_coverage_score = 1 - abs(mean_coverage - 0.95) / max(abs(mean_coverage - 0.95), na.rm = TRUE)
    ) %>%
    ungroup()

  # Normalize and score convergence (higher is better)
  convergence_scores <- convergence_data %>%
    group_by(method, family, config_name) %>%
    summarize(
      mean_convergence = mean(convergence_rate, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(family, config_name) %>%
    mutate(
      normalized_convergence_score = mean_convergence / max(mean_convergence, na.rm = TRUE)
    ) %>%
    ungroup()

  # Normalize and score computation time (lower is better)
  time_scores <- convergence_data %>%
    group_by(method, family, config_name) %>%
    summarize(
      mean_time = mean(avg_time, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(family, config_name) %>%
    mutate(
      normalized_time_score = 1 - (mean_time / max(mean_time, na.rm = TRUE))
    ) %>%
    ungroup()

  # Combine scores
  result <- bias_scores %>%
    rename(config_name = config) %>%
    left_join(coverage_scores %>% rename(config_name = config),
              by = c("method", "parameter", "family", "config_name")) %>%
    left_join(convergence_scores,
              by = c("method", "family", "config_name")) %>%
    left_join(time_scores,
              by = c("method", "family", "config_name")) %>%
    mutate(
      weighted_score = rel_bias_weight * normalized_bias_score +
        coverage_weight * normalized_coverage_score +
        convergence_weight * normalized_convergence_score +
        time_weight * normalized_time_score
    ) %>%
    group_by(method, family) %>%
    summarize(
      mean_score = mean(weighted_score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(family, desc(mean_score))

  return(result)
}

#' Analyze relative bias trends by sample size
#'
#' @param data Data frame with parameter results
#' @param param Parameter name to analyze
#' @return ggplot object with bias trend
analyze_rel_bias_trend <- function(data, param) {
  data %>%
    filter(parameter == param) %>%
    group_by(method, sample_size) %>%
    summarize(
      mean_abs_rel_bias = mean(abs(rel_bias), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    ggplot(aes(x = as.numeric(as.character(sample_size)),
               y = mean_abs_rel_bias, color = method)) +
    geom_line() +
    geom_point() +
    scale_x_log10() +
    labs(
      title = paste("Relative Bias Convergence Rate for", param),
      x = "Sample Size (log scale)",
      y = "Mean Absolute Relative Bias (%)",
      color = "Method"
    ) +
    theme_bw()
}

#' Identify problematic configurations
#'
#' @param param_data Data frame with parameter results
#' @param convergence_data Data frame with convergence results
#' @param convergence_threshold Threshold for low convergence
#' @param rel_bias_threshold Threshold for high relative bias
#' @param time_threshold Threshold for high computation time
#' @return List with problematic configurations
identify_problematic_configs <- function(param_data, convergence_data,
                                         convergence_threshold = 0.8,
                                         rel_bias_threshold = 20.0,
                                         time_threshold = NULL) {

  # Identify configurations with low convergence rates
  low_convergence <- convergence_data %>%
    filter(convergence_rate < convergence_threshold) %>%
    select(family, config_name, method, sample_size, convergence_rate) %>%
    arrange(convergence_rate)

  # Identify configurations with high relative bias
  high_bias <- param_data %>%
    filter(abs(rel_bias) > rel_bias_threshold) %>%
    select(family, config, method, sample_size, parameter, rel_bias) %>%
    arrange(desc(abs(rel_bias)))

  # Identify configurations with excessive computation time (if threshold provided)
  high_time <- NULL
  if (!is.null(time_threshold)) {
    high_time <- convergence_data %>%
      filter(avg_time > time_threshold) %>%
      select(family, config_name, method, sample_size, avg_time) %>%
      arrange(desc(avg_time))
  }

  # Combine results
  list(
    low_convergence = low_convergence,
    high_bias = high_bias,
    high_time = high_time
  )
}



















#' # ---------------------------------------------------------------------------- #
#' # Title: Utility Functions for GKw MLE Bias Study
#' # Description: Contains all functions needed for analysis, visualization, and
#' #              reporting in the GKw MLE bias simulation study
#' # ---------------------------------------------------------------------------- #
#'
#' # ---------------------------------------------------------------------------- #
#' # Utility Functions
#' # ---------------------------------------------------------------------------- #
#'
#' #' Generate a deterministic seed based on configuration
#' #'
#' #' @param i Iteration number
#' #' @param config_name Configuration name
#' #' @param method Estimation method
#' #' @param sample_size Sample size
#' #' @param base_seed Base seed to start from
#' #' @return Integer seed value
#' generate_seed <- function(i, config_name, method, sample_size, base_seed = 12345) {
#'   # Encode method as a number (1-4)
#'   method_num <- match(method, c("tmb", "nr", "nlminb", "optim"))
#'
#'   # Encode config_name as a number based on its position
#'   config_num <- which(names(configurations) == config_name)
#'
#'   # Create a deterministic but unique seed
#'   seed <- base_seed + i * 10000 + config_num * 100 + method_num * 10 +
#'     as.numeric(as.factor(sample_size))
#'
#'   # Ensure it's within integer limits
#'   seed <- seed %% 2147483647  # Max integer in R
#'
#'   return(seed)
#' }
#'
#' #' Create a data frame from the parameter results with relative bias
#' #'
#' #' @param param_name Parameter name
#' #' @param param_value True parameter value
#' #' @param method Estimation method
#' #' @param family Distribution family
#' #' @param config Configuration name
#' #' @param sample_size Sample size
#' #' @param bias Absolute bias
#' #' @param se Standard error
#' #' @param emp_se Empirical standard error
#' #' @param coverage Coverage probability
#' #' @return Data frame with all bias metrics
#' create_param_df <- function(param_name, param_value, method, family, config,
#'                             sample_size, bias, se, emp_se, coverage) {
#'   # Calculate bias_z outside the data.frame call
#'   bias_z_val <- bias / se
#'
#'   # Calculate relative bias (as percentage of true value)
#'   rel_bias <- ifelse(param_value != 0,
#'                      (bias / param_value) * 100,
#'                      NA)  # Avoid division by zero
#'
#'   # Calculate relative standard error as percentage of true value
#'   rel_se <- ifelse(param_value != 0,
#'                    (se / param_value) * 100,
#'                    NA)
#'
#'   # Calculate relative empirical standard error
#'   rel_emp_se <- ifelse(param_value != 0,
#'                        (emp_se / param_value) * 100,
#'                        NA)
#'
#'   # Create data frame with all metrics
#'   data.frame(
#'     parameter = param_name,
#'     true_value = param_value,
#'     method = method,
#'     family = family,
#'     config = config,
#'     sample_size = sample_size,
#'     bias = bias,                    # Absolute bias
#'     rel_bias = rel_bias,            # Relative bias (%)
#'     se = se,                        # Standard error
#'     emp_se = emp_se,                # Empirical standard error
#'     rel_se = rel_se,                # Relative standard error (%)
#'     rel_emp_se = rel_emp_se,        # Relative empirical standard error (%)
#'     coverage = coverage,            # Coverage probability
#'     bias_z = bias_z_val,            # Standardized bias
#'     bias_lower_z = bias_z_val - 1,  # Lower bound for standardized bias
#'     bias_upper_z = bias_z_val + 1,  # Upper bound for standardized bias
#'     stringsAsFactors = FALSE
#'   )
#' }
#'
#' #' Compute coverage rate
#' #'
#' #' @param estimates Vector of parameter estimates
#' #' @param std_errors Vector of standard errors
#' #' @param true_value True parameter value
#' #' @return Coverage rate (proportion of intervals containing the true value)
#' compute_coverage <- function(estimates, std_errors, true_value) {
#'   lower_ci <- estimates - qnorm(0.975) * std_errors
#'   upper_ci <- estimates + qnorm(0.975) * std_errors
#'   fora <- sum(lower_ci > true_value) + sum(upper_ci < true_value)
#'   coverage <- 1 - fora / length(estimates)
#'   return(coverage)
#' }
#'
#' # ---------------------------------------------------------------------------- #
#' # Simulation Function
#' # ---------------------------------------------------------------------------- #
#'
#' #' Run a single simulation configuration
#' #'
#' #' @param config Configuration list with family and parameters
#' #' @param method Estimation method
#' #' @param sample_size Sample size
#' #' @param n_sim Number of simulation replications
#' #' @return List with simulation results
#' run_simulation <- function(config, method, sample_size, n_sim) {
#'
#'   # Extract configuration details
#'   family <- config$family
#'   params <- config$params
#'   config_name <- names(configurations)[which(sapply(configurations, function(x)
#'     identical(x$family, family) && identical(x$params, params)))]
#'
#'   # Function to generate samples for the appropriate family
#'   generate_sample <- function(n, family, params) {
#'     # Switch to the right random generation function
#'     rand_func <- switch(family,
#'                         "gkw" = rgkw,
#'                         "bkw" = rbkw,
#'                         "kkw" = rkkw,
#'                         "ekw" = rekw,
#'                         "mc" = rmc,
#'                         "kw" = rkw,
#'                         "beta" = rbeta_)
#'
#'     # Generate random sample
#'     do.call(rand_func, c(list(n = n), params))
#'   }
#'
#'   # Initialize results storage
#'   n_params <- length(params)
#'   param_names <- names(params)
#'
#'   # Create matrices to store results
#'   estimates <- matrix(NA, nrow = n_sim, ncol = n_params)
#'   std_errors <- matrix(NA, nrow = n_sim, ncol = n_params)
#'   converged <- rep(FALSE, n_sim)
#'   execution_times <- rep(NA, n_sim)
#'
#'   # Run simulations
#'   for (i in 1:n_sim) {
#'     # Set seed for reproducibility (unique to this configuration)
#'     seed <- generate_seed(i, config_name, method, sample_size)
#'     set.seed(seed)
#'
#'     # Generate sample
#'     data <- generate_sample(sample_size, family, params)
#'
#'     # Measure execution time
#'     start_time <- Sys.time()
#'
#'     # Fit model using specified method
#'     fit_result <- tryCatch({
#'       fit <- gkwfit(data = data,
#'                     family = family,
#'                     fit = method,
#'                     hessian = TRUE,
#'                     silent = TRUE)
#'
#'       list(success = TRUE, fit = fit)
#'     },
#'     error = function(e) {
#'       return(list(success = FALSE, error = e))
#'     },
#'     warning = function(w) {
#'       return(list(success = FALSE, warning = w))
#'     })
#'
#'     # Record execution time
#'     end_time <- Sys.time()
#'     execution_times[i] <- as.numeric(difftime(end_time, start_time, units = "secs"))
#'
#'     # Process results if successful
#'     if (fit_result$success) {
#'       fit <- fit_result$fit
#'
#'       # Check if convergence was achieved
#'       if (fit$convergence) {
#'         converged[i] <- TRUE
#'
#'         # Extract parameter estimates and standard errors
#'         for (j in 1:n_params) {
#'           param <- param_names[j]
#'           if (param %in% names(fit$coefficients)) {
#'             estimates[i, j] <- fit$coefficients[param]
#'             std_errors[i, j] <- fit$std.errors[param]
#'           }
#'         }
#'       }
#'     }
#'   }
#'
#'   # Compute summary statistics for each parameter
#'   results <- list()
#'
#'   for (j in 1:n_params) {
#'     param <- param_names[j]
#'     true_value <- params[[param]]
#'
#'     # Filter for converged runs
#'     valid_estimates <- estimates[converged, j]
#'     valid_std_errors <- std_errors[converged, j]
#'
#'     # Compute bias and other metrics (if we have enough converged runs)
#'     if (sum(converged) > 0) {
#'       bias <- mean(valid_estimates, na.rm = TRUE) - true_value
#'       se <- mean(valid_std_errors, na.rm = TRUE)
#'       emp_se <- sd(valid_estimates, na.rm = TRUE)
#'       coverage <- compute_coverage(valid_estimates, valid_std_errors, true_value)
#'
#'       # Add results to the list
#'       results[[param]] <- create_param_df(
#'         param_name = param,
#'         param_value = true_value,
#'         method = method,
#'         family = family,
#'         config = config_name,
#'         sample_size = sample_size,
#'         bias = bias,
#'         se = se,
#'         emp_se = emp_se,
#'         coverage = coverage
#'       )
#'     }
#'   }
#'
#'   # Add convergence rate and timing information
#'   convergence_rate <- mean(converged)
#'   avg_time <- mean(execution_times[converged], na.rm = TRUE)
#'
#'   return(list(
#'     param_results = do.call(rbind, results),
#'     convergence_rate = convergence_rate,
#'     avg_time = avg_time,
#'     raw_estimates = estimates[converged, ],
#'     raw_std_errors = std_errors[converged, ],
#'     config_name = config_name,
#'     family = family,
#'     method = method,
#'     sample_size = sample_size
#'   ))
#' }
#'
#' # ---------------------------------------------------------------------------- #
#' # Visualization Functions
#' # ---------------------------------------------------------------------------- #
#'
#' #' Plot absolute bias using forest plots
#' #'
#' #' @param data Data frame with simulation results
#' #' @param param_name Parameter name to plot
#' #' @param methods_to_show Methods to include in the plot
#' #' @param families_to_show Families to include in the plot
#' #' @param configs_to_show Configurations to include in the plot
#' #' @param limits Y-axis limits
#' #' @param title Plot title
#' #' @param xlab X-axis label
#' #' @return ggplot object
#' plot_abs_bias <- function(data, param_name, methods_to_show = methods,
#'                           families_to_show = unique(data$family),
#'                           configs_to_show = unique(data$config),
#'                           limits = NULL,
#'                           title = NULL,
#'                           xlab = "Sample Size") {
#'
#'   # Filter data for the specific parameter
#'   filtered_data <- data %>%
#'     filter(parameter == param_name,
#'            method %in% methods_to_show,
#'            family %in% families_to_show,
#'            config %in% configs_to_show) %>%
#'     mutate(
#'       # Convert to factors for proper ordering
#'       sample_size = factor(sample_size, levels = sort(unique(sample_size))),
#'       method = factor(method, levels = methods_to_show),
#'
#'       # Create nice labels for the plots
#'       method_label = case_when(
#'         method == "tmb" ~ "TMB",
#'         method == "nr" ~ "Newton-Raphson",
#'         method == "nlminb" ~ "nlminb",
#'         method == "optim" ~ "optim (L-BFGS-B)",
#'         TRUE ~ as.character(method)
#'       ),
#'
#'       config_label = gsub("_", "-", config),
#'
#'       # Create LaTeX-formatted parameter name
#'       param_latex = TeX(paste0("$\\", param_name, "$"), output = "character")
#'     )
#'
#'   # Create standard error range for bias
#'   filtered_data <- filtered_data %>%
#'     mutate(
#'       bias_lower = bias - se,
#'       bias_upper = bias + se
#'     )
#'
#'   # Create the ggplot
#'   g <- ggplot(filtered_data, aes(x = sample_size, y = bias, color = method_label)) +
#'     geom_pointrange(
#'       aes(ymin = bias_lower, ymax = bias_upper),
#'       position = position_dodge(width = 0.5),
#'       size = 0.2,
#'       linewidth = 0.3
#'     ) +
#'     facet_nested(family + config_label ~ .,
#'                  scales = "free_y",
#'                  switch = "y",
#'                  labeller = labeller(config_label = label_value, family = label_value)) +
#'     geom_hline(yintercept = 0, linetype = "dotted") +
#'     labs(
#'       title = ifelse(is.null(title),
#'                      paste("Absolute bias for parameter", param_name),
#'                      title),
#'       y = expression(hat(bias) %+-% SE),
#'       x = xlab,
#'       color = "Method"
#'     ) +
#'     theme_bw() +
#'     theme(
#'       axis.text.x = element_text(size = 8),
#'       axis.text.y = element_text(size = 7),
#'       strip.text.x = element_text(size = 10, margin = margin(0.1, 0, 0.1, 0, "cm")),
#'       strip.text.y = element_text(size = 10, margin = margin(0.1, 0, 0.1, 0, "cm")),
#'       panel.spacing.x = unit(0.12, "lines"),
#'       panel.spacing.y = unit(0.12, "lines"),
#'       plot.title = element_text(hjust = 0.5),
#'       legend.position = "top"
#'     ) +
#'     coord_flip()
#'
#'   # Add limits if provided
#'   if (!is.null(limits)) {
#'     g <- g + scale_y_continuous(limits = limits)
#'   }
#'
#'   return(g)
#' }
#'
#' #' Plot relative bias using forest plots
#' #'
#' #' @param data Data frame with simulation results
#' #' @param param_name Parameter name to plot
#' #' @param methods_to_show Methods to include in the plot
#' #' @param families_to_show Families to include in the plot
#' #' @param configs_to_show Configurations to include in the plot
#' #' @param limits Y-axis limits
#' #' @param title Plot title
#' #' @param xlab X-axis label
#' #' @return ggplot object
#' plot_rel_bias <- function(data, param_name, methods_to_show = methods,
#'                           families_to_show = unique(data$family),
#'                           configs_to_show = unique(data$config),
#'                           limits = NULL,
#'                           title = NULL,
#'                           xlab = "Sample Size") {
#'
#'   # Filter data for the specific parameter
#'   filtered_data <- data %>%
#'     filter(parameter == param_name,
#'            method %in% methods_to_show,
#'            family %in% families_to_show,
#'            config %in% configs_to_show) %>%
#'     mutate(
#'       # Convert to factors for proper ordering
#'       sample_size = factor(sample_size, levels = sort(unique(sample_size))),
#'       method = factor(method, levels = methods_to_show),
#'
#'       # Create nice labels for the plots
#'       method_label = case_when(
#'         method == "tmb" ~ "TMB",
#'         method == "nr" ~ "Newton-Raphson",
#'         method == "nlminb" ~ "nlminb",
#'         method == "optim" ~ "optim (L-BFGS-B)",
#'         TRUE ~ as.character(method)
#'       ),
#'
#'       config_label = gsub("_", "-", config),
#'
#'       # Create LaTeX-formatted parameter name
#'       param_latex = TeX(paste0("$\\", param_name, "$"), output = "character")
#'     )
#'
#'   # Calculate relative bias confidence intervals
#'   # using error propagation for the ratio: Var(X/Y) ≈ (X/Y)^2 * [Var(X)/X^2 + Var(Y)/Y^2]
#'   filtered_data <- filtered_data %>%
#'     mutate(
#'       # For relative bias, we need the relative standard error
#'       rel_bias_lower = rel_bias - rel_se,
#'       rel_bias_upper = rel_bias + rel_se
#'     )
#'
#'   # Create the ggplot for relative bias
#'   g <- ggplot(filtered_data, aes(x = sample_size, y = rel_bias, color = method_label)) +
#'     geom_pointrange(
#'       aes(ymin = rel_bias_lower, ymax = rel_bias_upper),
#'       position = position_dodge(width = 0.5),
#'       size = 0.2,
#'       linewidth = 0.3
#'     ) +
#'     facet_nested(family + config_label ~ .,
#'                  scales = "free_y",
#'                  switch = "y",
#'                  labeller = labeller(config_label = label_value, family = label_value)) +
#'     geom_hline(yintercept = 0, linetype = "dotted") +
#'     labs(
#'       title = ifelse(is.null(title),
#'                      paste("Relative bias (%) for parameter", param_name),
#'                      title),
#'       y = "Relative Bias (%) ± SE",
#'       x = xlab,
#'       color = "Method"
#'     ) +
#'     theme_bw() +
#'     theme(
#'       axis.text.x = element_text(size = 8),
#'       axis.text.y = element_text(size = 7),
#'       strip.text.x = element_text(size = 10, margin = margin(0.1, 0, 0.1, 0, "cm")),
#'       strip.text.y = element_text(size = 10, margin = margin(0.1, 0, 0.1, 0, "cm")),
#'       panel.spacing.x = unit(0.12, "lines"),
#'       panel.spacing.y = unit(0.12, "lines"),
#'       plot.title = element_text(hjust = 0.5),
#'       legend.position = "top"
#'     ) +
#'     coord_flip()
#'
#'   # Add limits if provided
#'   if (!is.null(limits)) {
#'     g <- g + scale_y_continuous(limits = limits)
#'   }
#'
#'   return(g)
#' }
#'
#' #' Plot standardized bias (bias/se) using forest plots
#' #'
#' #' @param data Data frame with simulation results
#' #' @param param_name Parameter name to plot
#' #' @param methods_to_show Methods to include in the plot
#' #' @param families_to_show Families to include in the plot
#' #' @param configs_to_show Configurations to include in the plot
#' #' @param limits Y-axis limits
#' #' @param title Plot title
#' #' @param xlab X-axis label
#' #' @return ggplot object
#' plot_bias_z <- function(data, param_name, methods_to_show = methods,
#'                         families_to_show = unique(data$family),
#'                         configs_to_show = unique(data$config),
#'                         limits = c(-2, 2),
#'                         title = NULL,
#'                         xlab = "Sample Size") {
#'
#'   # Filter data for the specific parameter
#'   filtered_data <- data %>%
#'     filter(parameter == param_name,
#'            method %in% methods_to_show,
#'            family %in% families_to_show,
#'            config %in% configs_to_show) %>%
#'     mutate(
#'       # Convert to factors for proper ordering
#'       sample_size = factor(sample_size, levels = sort(unique(sample_size))),
#'       method = factor(method, levels = methods_to_show),
#'
#'       # Create nice labels for the plots
#'       method_label = case_when(
#'         method == "tmb" ~ "TMB",
#'         method == "nr" ~ "Newton-Raphson",
#'         method == "nlminb" ~ "nlminb",
#'         method == "optim" ~ "optim (L-BFGS-B)",
#'         TRUE ~ as.character(method)
#'       ),
#'
#'       config_label = gsub("_", "-", config),
#'
#'       # Create LaTeX-formatted parameter name
#'       param_latex = TeX(paste0("$\\", param_name, "$"), output = "character")
#'     )
#'
#'   # Create the ggplot
#'   g <- ggplot(filtered_data, aes(x = sample_size, y = bias_z, color = method_label)) +
#'     geom_pointrange(
#'       aes(ymin = bias_lower_z, ymax = bias_upper_z),
#'       position = position_dodge(width = 0.5),
#'       size = 0.2,
#'       linewidth = 0.3
#'     ) +
#'     facet_nested(family + config_label ~ .,
#'                  scales = "free_y",
#'                  switch = "y",
#'                  labeller = labeller(config_label = label_value, family = label_value)) +
#'     geom_hline(yintercept = 0, linetype = "dotted") +
#'     labs(
#'       title = ifelse(is.null(title),
#'                      paste("Standardized bias for parameter", param_name),
#'                      title),
#'       y = expression(hat(bias) / SE %+-% 1),
#'       x = xlab,
#'       color = "Method"
#'     ) +
#'     theme_bw() +
#'     theme(
#'       axis.text.x = element_text(size = 8),
#'       axis.text.y = element_text(size = 7),
#'       strip.text.x = element_text(size = 10, margin = margin(0.1, 0, 0.1, 0, "cm")),
#'       strip.text.y = element_text(size = 10, margin = margin(0.1, 0, 0.1, 0, "cm")),
#'       panel.spacing.x = unit(0.12, "lines"),
#'       panel.spacing.y = unit(0.12, "lines"),
#'       plot.title = element_text(hjust = 0.5),
#'       legend.position = "top"
#'     ) +
#'     coord_flip()
#'
#'   # Add limits if provided
#'   if (!is.null(limits)) {
#'     g <- g + scale_y_continuous(limits = limits)
#'   }
#'
#'   return(g)
#' }
#'
#' #' Plot coverage rates
#' #'
#' #' @param data Data frame with simulation results
#' #' @param param_name Parameter name to plot
#' #' @param methods_to_show Methods to include in the plot
#' #' @param families_to_show Families to include in the plot
#' #' @param configs_to_show Configurations to include in the plot
#' #' @param target_coverage Target coverage rate (usually 0.95)
#' #' @param limits Y-axis limits
#' #' @param title Plot title
#' #' @param xlab X-axis label
#' #' @return ggplot object
#' plot_coverage <- function(data, param_name, methods_to_show = methods,
#'                           families_to_show = unique(data$family),
#'                           configs_to_show = unique(data$config),
#'                           target_coverage = 0.95,
#'                           limits = c(0.7, 1.0),
#'                           title = NULL,
#'                           xlab = "Sample Size") {
#'
#'   # Filter data for the specific parameter
#'   filtered_data <- data %>%
#'     filter(parameter == param_name,
#'            method %in% methods_to_show,
#'            family %in% families_to_show,
#'            config %in% configs_to_show) %>%
#'     mutate(
#'       # Convert to factors for proper ordering
#'       sample_size = factor(sample_size, levels = sort(unique(sample_size))),
#'       method = factor(method, levels = methods_to_show),
#'
#'       # Create nice labels for the plots
#'       method_label = case_when(
#'         method == "tmb" ~ "TMB",
#'         method == "nr" ~ "Newton-Raphson",
#'         method == "nlminb" ~ "nlminb",
#'         method == "optim" ~ "optim (L-BFGS-B)",
#'         TRUE ~ as.character(method)
#'       ),
#'
#'       config_label = gsub("_", "-", config),
#'
#'       # Create LaTeX-formatted parameter name
#'       param_latex = TeX(paste0("$\\", param_name, "$"), output = "character")
#'     )
#'
#'   # Create the coverage plot
#'   g <- ggplot(filtered_data,
#'               aes(x = sample_size, y = coverage,
#'                   group = interaction(method_label, family),
#'                   color = method_label)) +
#'     geom_line(linewidth = 0.3) +
#'     geom_point(aes(shape = method_label)) +
#'     facet_nested(family + config_label ~ .,
#'                  switch = "y",
#'                  labeller = labeller(config_label = label_value, family = label_value)) +
#'     geom_hline(yintercept = target_coverage, linewidth = 0.2, linetype = "dashed") +
#'     labs(
#'       title = ifelse(is.null(title),
#'                      paste("Coverage rate for parameter", param_name),
#'                      title),
#'       y = "Coverage Rate",
#'       x = xlab,
#'       color = "Method",
#'       shape = "Method"
#'     ) +
#'     theme_bw() +
#'     theme(
#'       axis.text.x = element_text(size = 8),
#'       axis.text.y = element_text(size = 7),
#'       strip.text.x = element_text(size = 10, margin = margin(0.1, 0, 0.1, 0, "cm")),
#'       strip.text.y = element_text(size = 10, margin = margin(0.1, 0, 0.1, 0, "cm")),
#'       panel.spacing.x = unit(0.12, "lines"),
#'       panel.spacing.y = unit(0.12, "lines"),
#'       plot.title = element_text(hjust = 0.5),
#'       legend.position = "top"
#'     ) +
#'     scale_y_continuous(limits = limits) +
#'     coord_flip()
#'
#'   return(g)
#' }
#'
#' #' Plot convergence rates
#' #'
#' #' @param convergence_data Data frame with convergence results
#' #' @param methods_to_show Methods to include in the plot
#' #' @param families_to_show Families to include in the plot
#' #' @param configs_to_show Configurations to include in the plot
#' #' @param limits Y-axis limits
#' #' @param title Plot title
#' #' @param xlab X-axis label
#' #' @return ggplot object
#' plot_convergence <- function(convergence_data,
#'                              methods_to_show = methods,
#'                              families_to_show = unique(convergence_data$family),
#'                              configs_to_show = unique(convergence_data$config_name),
#'                              limits = c(0, 1),
#'                              title = "Convergence Rates",
#'                              xlab = "Sample Size") {
#'
#'   # Filter data
#'   filtered_data <- convergence_data %>%
#'     filter(method %in% methods_to_show,
#'            family %in% families_to_show,
#'            config_name %in% configs_to_show) %>%
#'     mutate(
#'       # Convert to factors for proper ordering
#'       sample_size = factor(sample_size, levels = sort(unique(sample_size))),
#'       method = factor(method, levels = methods_to_show),
#'
#'       # Create nice labels for the plots
#'       method_label = case_when(
#'         method == "tmb" ~ "TMB",
#'         method == "nr" ~ "Newton-Raphson",
#'         method == "nlminb" ~ "nlminb",
#'         method == "optim" ~ "optim (L-BFGS-B)",
#'         TRUE ~ as.character(method)
#'       ),
#'
#'       config_label = gsub("_", " ", config_name)
#'     )
#'
#'   # Create the convergence plot
#'   g <- ggplot(filtered_data,
#'               aes(x = sample_size, y = convergence_rate,
#'                   group = interaction(method_label, family),
#'                   color = method_label)) +
#'     geom_line(linewidth = 0.3) +
#'     geom_point(aes(shape = method_label)) +
#'     facet_nested(family + config_label ~ .,
#'                  switch = "y",
#'                  labeller = labeller(config_label = label_value, family = label_value)) +
#'     labs(
#'       title = title,
#'       y = "Convergence Rate",
#'       x = xlab,
#'       color = "Method",
#'       shape = "Method"
#'     ) +
#'     theme_bw() +
#'     theme(
#'       axis.text.x = element_text(size = 8),
#'       axis.text.y = element_text(size = 7),
#'       strip.text.x = element_text(size = 10, margin = margin(0.1, 0, 0.1, 0, "cm")),
#'       strip.text.y = element_text(size = 10, margin = margin(0.1, 0, 0.1, 0, "cm")),
#'       panel.spacing.x = unit(0.12, "lines"),
#'       panel.spacing.y = unit(0.12, "lines"),
#'       plot.title = element_text(hjust = 0.5),
#'       legend.position = "top"
#'     ) +
#'     scale_y_continuous(limits = limits) +
#'     coord_flip()
#'
#'   return(g)
#' }
#'
#' #' Plot computation time
#' #'
#' #' @param convergence_data Data frame with convergence results
#' #' @param methods_to_show Methods to include in the plot
#' #' @param families_to_show Families to include in the plot
#' #' @param configs_to_show Configurations to include in the plot
#' #' @param log_scale Whether to use log scale for y-axis
#' #' @param title Plot title
#' #' @param xlab X-axis label
#' #' @return ggplot object
#' plot_computation_time <- function(convergence_data,
#'                                   methods_to_show = methods,
#'                                   families_to_show = unique(convergence_data$family),
#'                                   configs_to_show = unique(convergence_data$config_name),
#'                                   log_scale = TRUE,
#'                                   title = "Computation Time",
#'                                   xlab = "Sample Size") {
#'
#'   # Filter data
#'   filtered_data <- convergence_data %>%
#'     filter(method %in% methods_to_show,
#'            family %in% families_to_show,
#'            config_name %in% configs_to_show) %>%
#'     mutate(
#'       # Convert to factors for proper ordering
#'       sample_size = factor(sample_size, levels = sort(unique(sample_size))),
#'       method = factor(method, levels = methods_to_show),
#'
#'       # Create nice labels for the plots
#'       method_label = case_when(
#'         method == "tmb" ~ "TMB",
#'         method == "nr" ~ "Newton-Raphson",
#'         method == "nlminb" ~ "nlminb",
#'         method == "optim" ~ "optim (L-BFGS-B)",
#'         TRUE ~ as.character(method)
#'       ),
#'
#'       config_label = gsub("_", " ", config_name)
#'     )
#'
#'   # Create the computation time plot
#'   g <- ggplot(filtered_data,
#'               aes(x = sample_size, y = avg_time,
#'                   group = interaction(method_label, family),
#'                   color = method_label)) +
#'     geom_line(linewidth = 0.3) +
#'     geom_point(aes(shape = method_label)) +
#'     facet_nested(family + config_label ~ .,
#'                  switch = "y",
#'                  labeller = labeller(config_label = label_value, family = label_value)) +
#'     labs(
#'       title = title,
#'       y = ifelse(log_scale, "Computation Time (log seconds)", "Computation Time (seconds)"),
#'       x = xlab,
#'       color = "Method",
#'       shape = "Method"
#'     ) +
#'     theme_bw() +
#'     theme(
#'       axis.text.x = element_text(size = 8),
#'       axis.text.y = element_text(size = 7),
#'       strip.text.x = element_text(size = 10, margin = margin(0.1, 0, 0.1, 0, "cm")),
#'       strip.text.y = element_text(size = 10, margin = margin(0.1, 0, 0.1, 0, "cm")),
#'       panel.spacing.x = unit(0.12, "lines"),
#'       panel.spacing.y = unit(0.12, "lines"),
#'       plot.title = element_text(hjust = 0.5),
#'       legend.position = "top"
#'     ) +
#'     coord_flip()
#'
#'   # Apply log scale if requested
#'   if (log_scale) {
#'     g <- g + scale_y_log10()
#'   }
#'
#'   return(g)
#' }
#'
#' # ---------------------------------------------------------------------------- #
#' # Analysis Functions
#' # ---------------------------------------------------------------------------- #
#'
#' #' Create bias summary tables
#' #'
#' #' @param param_results Data frame with parameter results
#' #' @return Data frame with bias summary
#' create_bias_summary <- function(param_results) {
#'   bias_summary <- param_results %>%
#'     group_by(parameter, family, method, sample_size) %>%
#'     summarize(
#'       mean_bias = mean(bias, na.rm = TRUE),
#'       mean_rel_bias = mean(rel_bias, na.rm = TRUE),
#'       mean_bias_z = mean(bias_z, na.rm = TRUE),
#'       mean_true_value = mean(true_value, na.rm = TRUE),
#'       min_rel_bias = min(rel_bias, na.rm = TRUE),
#'       max_rel_bias = max(rel_bias, na.rm = TRUE),
#'       .groups = "drop"
#'     ) %>%
#'     arrange(parameter, family, method, sample_size)
#'
#'   return(bias_summary)
#' }
#'
#' #' Create coverage summary tables
#' #'
#' #' @param param_results Data frame with parameter results
#' #' @return Data frame with coverage summary
#' create_coverage_summary <- function(param_results) {
#'   coverage_summary <- param_results %>%
#'     group_by(parameter, family, method, sample_size) %>%
#'     summarize(
#'       mean_coverage = mean(coverage, na.rm = TRUE),
#'       min_coverage = min(coverage, na.rm = TRUE),
#'       max_coverage = max(coverage, na.rm = TRUE),
#'       target_coverage = 0.95,
#'       coverage_diff = mean_coverage - 0.95,
#'       .groups = "drop"
#'     ) %>%
#'     arrange(parameter, family, method, sample_size)
#'
#'   return(coverage_summary)
#' }
#'
#' #' Create convergence and timing summary
#' #'
#' #' @param convergence_df Data frame with convergence results
#' #' @return Data frame with convergence and timing summary
#' create_convergence_summary <- function(convergence_df) {
#'   convergence_summary <- convergence_df %>%
#'     group_by(family, method, sample_size) %>%
#'     summarize(
#'       mean_convergence = mean(convergence_rate, na.rm = TRUE),
#'       min_convergence = min(convergence_rate, na.rm = TRUE),
#'       max_convergence = max(convergence_rate, na.rm = TRUE),
#'       mean_time = mean(avg_time, na.rm = TRUE),
#'       min_time = min(avg_time, na.rm = TRUE),
#'       max_time = max(avg_time, na.rm = TRUE),
#'       .groups = "drop"
#'     ) %>%
#'     arrange(family, method, sample_size)
#'
#'   return(convergence_summary)
#' }
#'
#' #' Compare methods
#' #'
#' #' @param param_results Data frame with parameter results
#' #' @param convergence_df Data frame with convergence results
#' #' @return Data frame with method comparison
#' create_method_comparison <- function(param_results, convergence_df) {
#'   method_comparison <- param_results %>%
#'     group_by(method, parameter, sample_size) %>%
#'     summarize(
#'       mean_abs_bias = mean(abs(bias), na.rm = TRUE),
#'       mean_abs_rel_bias = mean(abs(rel_bias), na.rm = TRUE),
#'       mean_abs_bias_z = mean(abs(bias_z), na.rm = TRUE),
#'       mean_coverage = mean(coverage, na.rm = TRUE),
#'       .groups = "drop"
#'     ) %>%
#'     left_join(
#'       convergence_df %>%
#'         group_by(method, sample_size) %>%
#'         summarize(
#'           mean_convergence = mean(convergence_rate, na.rm = TRUE),
#'           mean_time = mean(avg_time, na.rm = TRUE),
#'           .groups = "drop"
#'         ),
#'       by = c("method", "sample_size")
#'     ) %>%
#'     arrange(parameter, method, sample_size)
#'
#'   return(method_comparison)
#' }
#'
#' #' Rank methods by performance
#' #'
#' #' @param param_data Data frame with parameter results
#' #' @param convergence_data Data frame with convergence results
#' #' @param rel_bias_weight Weight for relative bias in overall score
#' #' @param coverage_weight Weight for coverage in overall score
#' #' @param convergence_weight Weight for convergence in overall score
#' #' @param time_weight Weight for time in overall score
#' #' @return Data frame with method rankings
#' rank_methods <- function(param_data, convergence_data,
#'                          rel_bias_weight = 0.4,
#'                          coverage_weight = 0.3,
#'                          convergence_weight = 0.2,
#'                          time_weight = 0.1) {
#'
#'   # Normalize and score relative bias (lower absolute value is better)
#'   bias_scores <- param_data %>%
#'     group_by(method, parameter, family, config) %>%
#'     summarize(
#'       mean_abs_rel_bias = mean(abs(rel_bias), na.rm = TRUE),
#'       .groups = "drop"
#'     ) %>%
#'     group_by(parameter, family, config) %>%
#'     mutate(
#'       normalized_bias_score = 1 - (mean_abs_rel_bias / max(mean_abs_rel_bias, na.rm = TRUE))
#'     ) %>%
#'     ungroup()
#'
#'   # Normalize and score coverage (closer to 0.95 is better)
#'   coverage_scores <- param_data %>%
#'     group_by(method, parameter, family, config) %>%
#'     summarize(
#'       mean_coverage = mean(coverage, na.rm = TRUE),
#'       .groups = "drop"
#'     ) %>%
#'     group_by(parameter, family, config) %>%
#'     mutate(
#'       normalized_coverage_score = 1 - abs(mean_coverage - 0.95) / max(abs(mean_coverage - 0.95), na.rm = TRUE)
#'     ) %>%
#'     ungroup()
#'
#'   # Normalize and score convergence (higher is better)
#'   convergence_scores <- convergence_data %>%
#'     group_by(method, family, config_name) %>%
#'     summarize(
#'       mean_convergence = mean(convergence_rate, na.rm = TRUE),
#'       .groups = "drop"
#'     ) %>%
#'     group_by(family, config_name) %>%
#'     mutate(
#'       normalized_convergence_score = mean_convergence / max(mean_convergence, na.rm = TRUE)
#'     ) %>%
#'     ungroup()
#'
#'   # Normalize and score computation time (lower is better)
#'   time_scores <- convergence_data %>%
#'     group_by(method, family, config_name) %>%
#'     summarize(
#'       mean_time = mean(avg_time, na.rm = TRUE),
#'       .groups = "drop"
#'     ) %>%
#'     group_by(family, config_name) %>%
#'     mutate(
#'       normalized_time_score = 1 - (mean_time / max(mean_time, na.rm = TRUE))
#'     ) %>%
#'     ungroup()
#'
#'   # Combine scores
#'   result <- bias_scores %>%
#'     rename(config_name = config) %>%
#'     left_join(coverage_scores %>% rename(config_name = config),
#'               by = c("method", "parameter", "family", "config_name")) %>%
#'     left_join(convergence_scores,
#'               by = c("method", "family", "config_name")) %>%
#'     left_join(time_scores,
#'               by = c("method", "family", "config_name")) %>%
#'     mutate(
#'       weighted_score = rel_bias_weight * normalized_bias_score +
#'         coverage_weight * normalized_coverage_score +
#'         convergence_weight * normalized_convergence_score +
#'         time_weight * normalized_time_score
#'     ) %>%
#'     group_by(method, family) %>%
#'     summarize(
#'       mean_score = mean(weighted_score, na.rm = TRUE),
#'       .groups = "drop"
#'     ) %>%
#'     arrange(family, desc(mean_score))
#'
#'   return(result)
#' }
#'
#' #' Analyze relative bias trends by sample size
#' #'
#' #' @param data Data frame with parameter results
#' #' @param param Parameter name to analyze
#' #' @return ggplot object with bias trend
#' analyze_rel_bias_trend <- function(data, param) {
#'   data %>%
#'     filter(parameter == param) %>%
#'     group_by(method, sample_size) %>%
#'     summarize(
#'       mean_abs_rel_bias = mean(abs(rel_bias), na.rm = TRUE),
#'       .groups = "drop"
#'     ) %>%
#'     ggplot(aes(x = as.numeric(as.character(sample_size)),
#'                y = mean_abs_rel_bias, color = method)) +
#'     geom_line() +
#'     geom_point() +
#'     scale_x_log10() +
#'     labs(
#'       title = paste("Relative Bias Convergence Rate for", param),
#'       x = "Sample Size (log scale)",
#'       y = "Mean Absolute Relative Bias (%)",
#'       color = "Method"
#'     ) +
#'     theme_bw()
#' }
#'
#' #' Identify problematic configurations
#' #'
#' #' @param param_data Data frame with parameter results
#' #' @param convergence_data Data frame with convergence results
#' #' @param convergence_threshold Threshold for low convergence
#' #' @param rel_bias_threshold Threshold for high relative bias
#' #' @param time_threshold Threshold for high computation time
#' #' @return List with problematic configurations
#' identify_problematic_configs <- function(param_data, convergence_data,
#'                                          convergence_threshold = 0.8,
#'                                          rel_bias_threshold = 20.0,
#'                                          time_threshold = NULL) {
#'
#'   # Identify configurations with low convergence rates
#'   low_convergence <- convergence_data %>%
#'     filter(convergence_rate < convergence_threshold) %>%
#'     select(family, config_name, method, sample_size, convergence_rate) %>%
#'     arrange(convergence_rate)
#'
#'   # Identify configurations with high relative bias
#'   high_bias <- param_data %>%
#'     filter(abs(rel_bias) > rel_bias_threshold) %>%
#'     select(family, config, method, sample_size, parameter, rel_bias) %>%
#'     arrange(desc(abs(rel_bias)))
#'
#'   # Identify configurations with excessive computation time (if threshold provided)
#'   high_time <- NULL
#'   if (!is.null(time_threshold)) {
#'     high_time <- convergence_data %>%
#'       filter(avg_time > time_threshold) %>%
#'       select(family, config_name, method, sample_size, avg_time) %>%
#'       arrange(desc(avg_time))
#'   }
#'
#'   # Combine results
#'   list(
#'     low_convergence = low_convergence,
#'     high_bias = high_bias,
#'     high_time = high_time
#'   )
#' }
