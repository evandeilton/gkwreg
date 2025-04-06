#' @title Fit All Generalized Kumaraswamy Family Distributions and Compare Them
#'
#' @description
#' Fits all seven distributions from the Generalized Kumaraswamy (GKw) family to data
#' using maximum likelihood estimation through Template Model Builder (TMB).
#' It provides a comprehensive comparison of fit quality across all families
#' through statistics and visualization.
#'
#' @param data A numeric vector with values strictly between 0 and 1.
#' @param method Optimization method to use. One of: \code{"nlminb"} (default), \code{"Nelder-Mead"},
#'   \code{"BFGS"}, \code{"CG"}, \code{"L-BFGS-B"} or \code{"SANN"}.
#' @param use_moments Logical; if \code{TRUE}, attempts to use method of moments estimates
#'   as initial values. Default: \code{FALSE}.
#' @param profile Logical; if \code{TRUE}, computes likelihood profiles for parameters. Default: \code{TRUE}.
#' @param npoints Integer; number of points to use in profile likelihood calculations. Default: 20.
#' @param plot Logical; if \code{TRUE}, generates comparison plots. Default: \code{TRUE}.
#'
#' @details
#' This function fits all seven distributions in the GKw family:
#' \itemize{
#'   \item \strong{GKw}: 5 parameters (\eqn{\alpha, \beta, \gamma, \delta, \lambda}) - All positive.
#'   \item \strong{BKw}: 4 parameters (\eqn{\alpha, \beta, \gamma, \delta}), \eqn{\lambda = 1} fixed - All positive.
#'   \item \strong{KKw}: 4 parameters (\eqn{\alpha, \beta, \delta, \lambda}), \eqn{\gamma = 1} fixed - All positive.
#'   \item \strong{EKw}: 3 parameters (\eqn{\alpha, \beta, \lambda}), \eqn{\gamma = 1, \delta = 0} fixed - All positive.
#'   \item \strong{Mc} (McDonald / Beta Power): 3 parameters (\eqn{\gamma, \delta, \lambda}), \eqn{\alpha = 1, \beta = 1} fixed - All positive.
#'   \item \strong{Kw} (Kumaraswamy): 2 parameters (\eqn{\alpha, \beta}), \eqn{\gamma = 1, \delta = 0, \lambda = 1} fixed - All positive.
#'   \item \strong{Beta}: 2 parameters (\eqn{\gamma, \delta}), \eqn{\alpha = 1, \beta = 1, \lambda = 1} fixed - All positive.
#' }
#'
#' The function generates comparison statistics including AIC, BIC, AICc, log-likelihood values,
#' and various goodness-of-fit measures. It also produces visualizations with all fitted
#' distributions overlaid on diagnostic plots.
#'
#' @return A list containing:
#' \item{fits}{A list of \code{gkwfit} objects for all seven distribution families.}
#' \item{comparison}{A data frame with comparison statistics (AIC, BIC, log-likelihood, etc.) for all models.}
#' \item{plots}{A ggplot2 object with diagnostic plots for all models if \code{plot = TRUE}.}
#'
#' @examples
#' \dontrun{
#' # Create sample data
#' set.seed(123)
#' sample_data <- rbeta(200, 2, 3)
#'
#' # Fit all GKw family distributions and compare
#' results <- gkwfitall(sample_data)
#'
#' # View comparison statistics
#' results$comparison
#'
#' # View comparison plots
#' print(results$plots)
#'
#' # Examine the best-fitting model
#' best_model <- results$fits[[which.min(results$comparison$AIC)]]
#' summary(best_model)
#' }
#'
#' @importFrom stats AIC BIC logLik
#' @importFrom ggplot2 ggplot aes geom_histogram geom_line geom_point geom_abline
#'   geom_text labs theme_minimal facet_wrap scale_color_brewer
#' @importFrom patchwork wrap_plots
#'
#' @export
gkwfitall <- function(data,
                      method = "nlminb",
                      use_moments = FALSE,
                      profile = TRUE,
                      npoints = 20,
                      plot = TRUE) {
  # Validate data
  if (missing(data) || !is.numeric(data)) {
    stop("'data' must be a numeric vector")
  }

  if (length(data) < 5) { # 5 is the max number of parameters for the GKw family
    stop("Data must have at least 5 observations to fit all GKw family models")
  }

  # Check data bounds
  epsilon <- .Machine$double.eps^0.5
  if (any(data <= epsilon | data >= (1 - epsilon))) {
    warning("Data contains values near or at the boundaries (0 or 1). Results may be unreliable.")
  }

  # List of all family distributions
  families <- c("gkw", "bkw", "kkw", "ekw", "mc", "kw", "beta")

  # Initialize list to store fitted models
  fits <- list()

  # Fit all family distributions
  # message("Fitting all GKw family distributions...")
  for (family in families) {
    # message(paste0("  Fitting ", family, " distribution..."))

    tryCatch(
      {
        fits[[family]] <- gkwfit(
          data = data,
          family = family,
          method = method,
          use_moments = use_moments,
          profile = profile,
          npoints = npoints,
          plot = FALSE, # We'll create our own combined plots
          silent = TRUE # No messages from individual fits
        )
      },
      error = function(e) {
        warning(paste0("Failed to fit ", family, " distribution: ", e$message))
        fits[[family]] <- NULL
      }
    )
  }

  # If no models could be fitted, stop
  if (length(fits) == 0) {
    stop("Could not fit any GKw family distributions to the data.")
  }

  # Create comparison table
  comparison <- create_comparison_table(fits)

  # Create comparison plots if requested
  if (plot) {
    plots <- create_comparison_plots(fits, data)
  } else {
    plots <- NULL
  }

  # Return results
  result <- list(
    fits = fits,
    comparison = comparison,
    plots = plots
  )

  class(result) <- "gkwfitall"
  return(result)
}

#' Create comparison table of fit statistics
#'
#' @param fits List of fitted gkwfit models
#' @return Data frame with comparison statistics
#' @keywords internal
create_comparison_table <- function(fits) {
  # Create data frame with model names as rows
  model_names <- names(fits)
  n_models <- length(model_names)

  # Columns for comparison
  comparison <- data.frame(
    Family = model_names,
    Parameters = numeric(n_models),
    LogLik = numeric(n_models),
    AIC = numeric(n_models),
    BIC = numeric(n_models),
    AICc = numeric(n_models),
    KS_stat = numeric(n_models),
    KS_pvalue = numeric(n_models),
    Mean_obs = numeric(n_models),
    Mean_fit = numeric(n_models),
    Var_obs = numeric(n_models),
    Var_fit = numeric(n_models),
    Convergence = logical(n_models),
    # Additional coefficient columns
    alpha_coef = numeric(n_models),
    beta_coef = numeric(n_models),
    gama_coef = numeric(n_models), # Using 'gama' as requested instead of 'gamma'
    delta_coef = numeric(n_models),
    lambda_coef = numeric(n_models),
    # Additional standard error columns
    alpha_se = numeric(n_models),
    beta_se = numeric(n_models),
    gama_se = numeric(n_models), # Using 'gama' as requested instead of 'gamma'
    delta_se = numeric(n_models),
    lambda_se = numeric(n_models),
    row.names = model_names,
    stringsAsFactors = FALSE
  )

  # Fill in comparison statistics
  for (i in seq_along(model_names)) {
    family <- model_names[i]
    fit <- fits[[family]]

    # Basic fit statistics
    comparison$Parameters[i] <- fit$df
    comparison$LogLik[i] <- fit$loglik
    comparison$AIC[i] <- fit$AIC
    comparison$BIC[i] <- fit$BIC
    comparison$AICc[i] <- fit$AICc
    comparison$Convergence[i] <- fit$convergence

    # Goodness-of-fit statistics (if available)
    if (!is.null(fit$gof) && !is.null(fit$gof$ks)) {
      comparison$KS_stat[i] <- fit$gof$ks$statistic
      comparison$KS_pvalue[i] <- fit$gof$ks$p.value
    } else {
      comparison$KS_stat[i] <- NA
      comparison$KS_pvalue[i] <- NA
    }

    # Diagnostic statistics (if available)
    if (!is.null(fit$diagnostics)) {
      comparison$Mean_obs[i] <- fit$diagnostics$mean_obs
      comparison$Mean_fit[i] <- fit$diagnostics$mean_fitted
      comparison$Var_obs[i] <- fit$diagnostics$var_obs
      comparison$Var_fit[i] <- fit$diagnostics$var_fitted
    } else {
      comparison$Mean_obs[i] <- mean(fit$data)
      comparison$Var_obs[i] <- stats::var(fit$data)
      comparison$Mean_fit[i] <- NA
      comparison$Var_fit[i] <- NA
    }

    # Set all coefficient and SE columns to NA by default
    coef_columns <- c("alpha_coef", "beta_coef", "gama_coef", "delta_coef", "lambda_coef")
    se_columns <- c("alpha_se", "beta_se", "gama_se", "delta_se", "lambda_se")
    for (col in c(coef_columns, se_columns)) {
      comparison[[col]][i] <- NA
    }

    # Map requested column names to actual parameter names in the model
    param_mapping <- c(
      alpha_coef = "alpha",
      beta_coef = "beta",
      gama_coef = "gamma", # Map 'gama_coef' to 'gamma' parameter in the model
      delta_coef = "delta",
      lambda_coef = "lambda"
    )

    # Extract coefficients and standard errors
    for (col_name in names(param_mapping)) {
      param_name <- param_mapping[[col_name]]
      se_col_name <- gsub("_coef", "_se", col_name)

      # Check if parameter exists in coefficients
      if (!is.null(fit$coefficients) && param_name %in% names(fit$coefficients)) {
        comparison[[col_name]][i] <- fit$coefficients[param_name]

        # Check if standard error is available
        if (!is.null(fit$std.errors) && param_name %in% names(fit$std.errors)) {
          comparison[[se_col_name]][i] <- fit$std.errors[param_name]
        }
      }
    }
  }

  # Sort by AIC (ascending)
  comparison <- comparison[order(comparison$AIC), ]

  return(comparison)
}

#' Create comparison plots of all fitted distributions
#'
#' @param fits List of fitted gkwfit models
#' @param data Original data vector
#' @return A ggplot2 object with multiple panels
#' @importFrom ggplot2 ggplot aes geom_histogram geom_line geom_point geom_abline
#'   labs theme_minimal facet_wrap scale_color_brewer
#' @importFrom patchwork wrap_plots
#' @keywords internal
create_comparison_plots <- function(fits, data) {
  # Check if required packages are available
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
    !requireNamespace("patchwork", quietly = TRUE)) {
    warning("Packages 'ggplot2' and 'patchwork' are required for plotting but not installed.")
    return(NULL)
  }

  # Extract model families
  families <- names(fits)

  # Calculate density values on a grid for all models
  x_grid <- seq(0.001, 0.999, length.out = 200)
  density_data <- data.frame(x = numeric(0), density = numeric(0), family = character(0))

  for (family in families) {
    fit <- fits[[family]]

    # Use the appropriate density function based on family
    density_func <- get_density_function(fit)

    # Calculate density values
    density_values <- tryCatch(
      {
        sapply(x_grid, density_func)
      },
      error = function(e) {
        warning(paste0("Error calculating density values for ", family, ": ", e$message))
        rep(NA, length(x_grid))
      }
    )

    # Add to data frame
    density_data <- rbind(
      density_data,
      data.frame(
        x = x_grid,
        density = density_values,
        family = family
      )
    )
  }

  # Prepare data for P-P and Q-Q plots
  pp_data <- data.frame(Empirical = numeric(0), Theoretical = numeric(0), family = character(0))
  qq_data <- data.frame(Theoretical = numeric(0), Empirical = numeric(0), family = character(0))

  for (family in families) {
    fit <- fits[[family]]

    # Get empirical values
    ecdf_vals <- stats::ecdf(data)(sort(data))

    # Get theoretical CDF values
    cdf_func <- get_cdf_function(fit)
    theor_cdf <- tryCatch(
      {
        sapply(sort(data), cdf_func)
      },
      error = function(e) {
        warning(paste0("Error calculating theoretical CDF for ", family, ": ", e$message))
        rep(NA, length(data))
      }
    )

    # Add to P-P data
    pp_data <- rbind(
      pp_data,
      data.frame(
        Empirical = ecdf_vals,
        Theoretical = theor_cdf,
        family = family
      )
    )

    # Get theoretical quantiles for Q-Q plot
    quant_func <- get_quantile_function(fit)
    theor_quant <- tryCatch(
      {
        sapply(stats::ppoints(length(data)), quant_func)
      },
      error = function(e) {
        warning(paste0("Error calculating theoretical quantiles for ", family, ": ", e$message))
        rep(NA, length(data))
      }
    )

    # Add to Q-Q data
    qq_data <- rbind(
      qq_data,
      data.frame(
        Theoretical = theor_quant,
        Empirical = sort(data),
        family = family
      )
    )
  }

  # Create histogram with density curve
  p1 <- ggplot2::ggplot() +
    ggplot2::geom_histogram(
      ggplot2::aes(x = data, y = ggplot2::after_stat(density)),
      bins = min(30, ceiling(sqrt(length(data)))),
      fill = "lightblue", color = "black", alpha = 0.7
    ) +
    ggplot2::geom_line(
      data = density_data,
      ggplot2::aes(x = x, y = density, color = family),
      size = 1
    ) +
    ggplot2::labs(
      title = "(A) Histogram with Fitted Densities",
      x = "Data", y = "Density"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::scale_color_brewer(palette = "Dark2")

  # Create P-P plot
  p2 <- ggplot2::ggplot(pp_data, ggplot2::aes(x = Theoretical, y = Empirical, color = family)) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
    ggplot2::labs(
      title = "(B) P-P Plot",
      x = "Theoretical Probability", y = "Empirical Probability"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::scale_color_brewer(palette = "Dark2")

  # Create Q-Q plot
  p3 <- ggplot2::ggplot(qq_data, ggplot2::aes(x = Theoretical, y = Empirical, color = family)) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
    ggplot2::labs(
      title = "(C) Q-Q Plot",
      x = "Theoretical Quantiles", y = "Empirical Quantiles"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::scale_color_brewer(palette = "Dark2")

  # Create residual plot
  p4 <- create_residual_plot(fits, data)

  # Combine plots
  combined_plots <- patchwork::wrap_plots(p1, p2, p3, p4) +
    patchwork::plot_annotation(
      title = "Comparison of GKw Family Distributions",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 16, face = "bold")
      )
    ) &
    ggplot2::theme(legend.position = "bottom")

  return(combined_plots)
}

#' Create residual plot for all fitted models
#'
#' @param fits List of fitted gkwfit models
#' @param data Original data vector
#' @return A ggplot2 object with the residual plot
#' @keywords internal
create_residual_plot <- function(fits, data) {
  # Prepare data frame for residuals
  residual_data <- data.frame(Empirical = numeric(0), Residual = numeric(0), family = character(0))

  for (family in names(fits)) {
    fit <- fits[[family]]

    # Get theoretical CDF values
    cdf_func <- get_cdf_function(fit)
    theor_cdf <- tryCatch(
      {
        sapply(sort(data), cdf_func)
      },
      error = function(e) {
        warning(paste0("Error calculating theoretical CDF for residuals (", family, "): ", e$message))
        rep(NA, length(data))
      }
    )

    # Calculate residuals
    ecdf_vals <- stats::ecdf(data)(sort(data))
    residuals <- ecdf_vals - theor_cdf

    # Add to data frame
    residual_data <- rbind(
      residual_data,
      data.frame(
        Empirical = sort(data),
        Residual = residuals,
        family = family
      )
    )
  }

  # Create residual plot
  p <- ggplot2::ggplot(residual_data, ggplot2::aes(x = Empirical, y = Residual, color = family)) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    ggplot2::labs(
      title = "(D) Residual Plot (ECDF - CDF)",
      x = "Empirical Data", y = "Difference"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::scale_color_brewer(palette = "Dark2")

  return(p)
}

#' Get the appropriate density function for a fitted model
#'
#' @param fit A gkwfit object
#' @return A function that calculates density values
#' @keywords internal
get_density_function <- function(fit) {
  # Extract parameters
  params <- fit$coefficients
  family <- fit$family

  # Define density function based on family
  switch(family,
    "gkw" = function(x) {
      dgkw(
        x, params["alpha"], params["beta"],
        params["gamma"], params["delta"],
        params["lambda"]
      )
    },
    "bkw" = function(x) {
      dbkw(
        x, params["alpha"], params["beta"],
        params["gamma"], params["delta"]
      )
    },
    "kkw" = function(x) {
      dkkw(
        x, params["alpha"], params["beta"],
        params["delta"], params["lambda"]
      )
    },
    "ekw" = function(x) {
      dekw(
        x, params["alpha"], params["beta"],
        params["lambda"]
      )
    },
    "mc" = function(x) {
      dmc(
        x, params["gamma"], params["delta"],
        params["lambda"]
      )
    },
    "kw" = function(x) dkw(x, params["alpha"], params["beta"]),
    "beta" = function(x) dbeta_(x, params["gamma"], params["delta"])
  )
}

#' Get the appropriate CDF function for a fitted model
#'
#' @param fit A gkwfit object
#' @return A function that calculates CDF values
#' @keywords internal
get_cdf_function <- function(fit) {
  # Extract parameters
  params <- fit$coefficients
  family <- fit$family

  # Define CDF function based on family
  switch(family,
    "gkw" = function(x) {
      pgkw(
        x, params["alpha"], params["beta"],
        params["gamma"], params["delta"],
        params["lambda"]
      )
    },
    "bkw" = function(x) {
      pbkw(
        x, params["alpha"], params["beta"],
        params["gamma"], params["delta"]
      )
    },
    "kkw" = function(x) {
      pkkw(
        x, params["alpha"], params["beta"],
        params["delta"], params["lambda"]
      )
    },
    "ekw" = function(x) {
      pekw(
        x, params["alpha"], params["beta"],
        params["lambda"]
      )
    },
    "mc" = function(x) {
      pmc(
        x, params["gamma"], params["delta"],
        params["lambda"]
      )
    },
    "kw" = function(x) pkw(x, params["alpha"], params["beta"]),
    "beta" = function(x) pbeta_(x, params["gamma"], params["delta"])
  )
}

#' Get the appropriate quantile function for a fitted model
#'
#' @param fit A gkwfit object
#' @return A function that calculates quantile values
#' @keywords internal
get_quantile_function <- function(fit) {
  # Extract parameters
  params <- fit$coefficients
  family <- fit$family

  # Define quantile function based on family
  switch(family,
    "gkw" = function(p) {
      qgkw(
        p, params["alpha"], params["beta"],
        params["gamma"], params["delta"],
        params["lambda"]
      )
    },
    "bkw" = function(p) {
      qbkw(
        p, params["alpha"], params["beta"],
        params["gamma"], params["delta"]
      )
    },
    "kkw" = function(p) {
      qkkw(
        p, params["alpha"], params["beta"],
        params["delta"], params["lambda"]
      )
    },
    "ekw" = function(p) {
      qekw(
        p, params["alpha"], params["beta"],
        params["lambda"]
      )
    },
    "mc" = function(p) {
      qmc(
        p, params["gamma"], params["delta"],
        params["lambda"]
      )
    },
    "kw" = function(p) qkw(p, params["alpha"], params["beta"]),
    "beta" = function(p) qbeta_(p, params["gamma"], params["delta"])
  )
}

#' Print method for gkwfitall objects
#'
#' @param x An object of class \code{"gkwfitall"}
#' @param ... Additional arguments (currently ignored)
#' @return Invisibly returns the input object
#' @export
print.gkwfitall <- function(x, ...) {
  cat("Fit Comparison of Generalized Kumaraswamy Family Distributions\n\n")

  cat("Number of families fitted:", length(x$fits), "\n")
  cat("Families:", paste(names(x$fits), collapse = ", "), "\n\n")

  cat("Comparison table (ordered by AIC):\n")
  print(x$comparison[, c("Family", "Parameters", "LogLik", "AIC", "BIC", "KS_stat", "KS_pvalue", "Convergence")])

  cat("\nBest fitting model:", x$comparison$Family[1], "\n")

  invisible(x)
}

#' Summary method for gkwfitall objects
#'
#' @param object An object of class \code{"gkwfitall"}
#' @param ... Additional arguments (currently ignored)
#' @return A summarized version of the gkwfitall object
#' @export
summary.gkwfitall <- function(object, ...) {
  # Create a summary object
  result <- list(
    comparison = object$comparison,
    best_model = object$fits[[object$comparison$Family[1]]],
    n_models = length(object$fits),
    model_names = names(object$fits)
  )

  class(result) <- "summary.gkwfitall"
  return(result)
}

#' Print method for summary.gkwfitall objects
#'
#' @param x An object of class \code{"summary.gkwfitall"}
#' @param ... Additional arguments (currently ignored)
#' @return Invisibly returns the input object
#' @export
print.summary.gkwfitall <- function(x, ...) {
  cat("Summary of Generalized Kumaraswamy Family Distributions Fit Comparison\n\n")

  cat("Number of families fitted:", x$n_models, "\n")
  cat("Families:", paste(x$model_names, collapse = ", "), "\n\n")

  cat("Comparison table (ordered by AIC):\n")
  print(x$comparison[, c("Family", "Parameters", "LogLik", "AIC", "BIC", "KS_stat", "KS_pvalue")])

  cat("\nBest fitting model:", x$comparison$Family[1], "\n\n")

  cat("Summary of best fitting model:\n")
  print(summary(x$best_model))

  invisible(x)
}

#' Plot method for gkwfitall objects
#'
#' @param x An object of class \code{"gkwfitall"}
#' @param ... Additional arguments (currently ignored)
#' @return Invisibly returns the input object
#' @export
plot.gkwfitall <- function(x, ...) {
  if (is.null(x$plots)) {
    stop("No plots available. The 'gkwfitall' object was created with plot = FALSE.")
  }

  print(x$plots)
  invisible(x)
}



# require(gkwreg)
# data("ReadingSkills", package = "betareg")
# data("ImpreciseTask", package = "betareg")
#
# # y <- ReadingSkills$accuracy
# y <- ImpreciseTask$location
#
#
# # Fit all GKw family distributions and compare
# results <- gkwfitall(y, method = "nlminb")
#
# # View comparison statistics
# results$comparison
#
# # View comparison plots
# print(results$plots)
#
# # Examine the best-fitting model
# best_model <- results$fits[[which.min(results$comparison$AIC)]]
# summary(best_model)
