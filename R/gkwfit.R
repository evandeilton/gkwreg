# Modularized implementation of the GKw family fitting functions
# This file contains the main gkwfit function and all its helper functions

#' Convert family string to numeric code for TMB
#'
#' @param family Character string specifying the family.
#' @return Integer code for TMB. 0="gkw", 1="bkw", 2="kkw", 3="ekw", 4="mc", 5="kw", 6="beta"
#' @keywords internal
.family_to_code <- function(family) {
  codes <- c(gkw = 0, bkw = 1, kkw = 2, ekw = 3, mc = 4, kw = 5, beta = 6)
  return(codes[family])
}

#' Get family parameter information
#'
#' @param family Character string specifying the family.
#' @return List with parameter names and count.
#' @keywords internal
.get_family_param_info <- function(family) {
  family_params <- list(
    gkw = list(names = c("alpha", "beta", "gamma", "delta", "lambda"), n = 5),
    bkw = list(names = c("alpha", "beta", "gamma", "delta"), n = 4),
    kkw = list(names = c("alpha", "beta", "delta", "lambda"), n = 4),
    ekw = list(names = c("alpha", "beta", "lambda"), n = 3),
    mc = list(names = c("gamma", "delta", "lambda"), n = 3),
    kw = list(names = c("alpha", "beta"), n = 2),
    beta = list(names = c("gamma", "delta"), n = 2)
  )

  return(family_params[[family]])
}

#' Get default fixed parameters for a family
#'
#' @param family Character string specifying the family.
#' @return Named list of fixed parameters.
#' @keywords internal
.get_default_fixed <- function(family) {
  family_fixed <- list()

  if (family != "gkw") {
    # Set fixed parameters based on family
    if (family == "bkw") {
      family_fixed$lambda <- 1 # BKw: λ = 1 fixed
    } else if (family == "kkw") {
      family_fixed$gamma <- 1 # KKw: γ = 1 fixed
    } else if (family == "ekw") {
      family_fixed$gamma <- 1 # EKw: γ = 1, δ = 0 fixed
      family_fixed$delta <- 0
    } else if (family == "mc") {
      family_fixed$alpha <- 1 # Mc: α = 1, β = 1 fixed
      family_fixed$beta <- 1
    } else if (family == "kw") {
      family_fixed$gamma <- 1 # Kw: γ = 1, δ = 0, λ = 1 fixed
      family_fixed$delta <- 0
      family_fixed$lambda <- 1
    } else if (family == "beta") {
      family_fixed$alpha <- 1 # Beta: α = 1, β = 1, λ = 1 fixed
      family_fixed$beta <- 1
      family_fixed$lambda <- 1
    }
  }

  return(family_fixed)
}

#' Get default start values for a family
#'
#' @param family Character string specifying the family.
#' @return Named list of starting values.
#' @keywords internal
.get_default_start <- function(family) {
  if (family == "gkw") {
    start <- list(alpha = 2, beta = 2, gamma = 1, delta = 0.5, lambda = 2)
  } else if (family == "bkw") {
    start <- list(alpha = 2, beta = 2, gamma = 1, delta = 0.5) # λ = 1 fixed
  } else if (family == "kkw") {
    start <- list(alpha = 2, beta = 2, delta = 0.5, lambda = 2) # γ = 1 fixed
  } else if (family == "ekw") {
    start <- list(alpha = 2, beta = 2, lambda = 2) # γ = 1, δ = 0 fixed
  } else if (family == "mc") {
    start <- list(gamma = 1, delta = 0.5, lambda = 2) # α = 1, β = 1 fixed
  } else if (family == "kw") {
    start <- list(alpha = 2, beta = 2) # γ = 1, δ = 0, λ = 1 fixed
  } else if (family == "beta") {
    start <- list(gamma = 1, delta = 0.5) # α = 1, β = 1, λ = 1 fixed
  }

  return(start)
}

#' Validate data for GKw family distributions
#'
#' @param data Numeric vector with values in the (0, 1) interval.
#' @param n_params Number of parameters to estimate.
#' @return Validated and possibly adjusted data vector.
#' @keywords internal
.validate_data <- function(data, n_params) {
  # Data validation
  if (!is.numeric(data)) {
    stop("Data must be numeric")
  }

  if (any(is.na(data))) {
    warning("Missing values removed from data")
    data <- data[!is.na(data)]
  }

  if (length(data) < n_params) {
    stop(paste0("At least ", n_params, " observations are required to fit a ", n_params, "-parameter model"))
  }

  # Check data bounds with a small tolerance for numerical precision
  epsilon <- .Machine$double.eps^0.5
  if (any(data <= epsilon | data >= (1 - epsilon))) {
    # Try to automatically adjust boundary values
    orig_length <- length(data)
    boundary_vals <- data <= epsilon | data >= (1 - epsilon)

    if (sum(boundary_vals) > 0.1 * orig_length) {
      # If more than 10% of data is at boundaries, warn and stop
      stop(
        "Too many data values (", sum(boundary_vals),
        ") are at or beyond the boundaries of (0, 1). Please preprocess your data."
      )
    } else {
      # Adjust boundary values with a warning
      data[data <= epsilon] <- epsilon + epsilon * 10
      data[data >= (1 - epsilon)] <- (1 - epsilon) - epsilon * 10
      warning("Adjusted ", sum(boundary_vals), " data values that were at or beyond boundaries to be within (0, 1)")
    }
  }

  return(data)
}

#' Validate parameters for GKw family distributions
#'
#' @param start List with initial parameter values.
#' @param fixed List of parameters to be held fixed (not estimated).
#' @param param_names Character vector of parameter names needed for the model.
#' @return List containing validated start and fixed parameter lists.
#' @keywords internal
.validate_parameters <- function(start, fixed, param_names) {
  # Check that all required parameters for the selected family are provided
  if (!is.null(start)) {
    missing_params <- setdiff(param_names, names(start))

    if (length(missing_params) > 0) {
      stop("Missing parameters in 'start': ", paste(missing_params, collapse = ", "))
    }

    # Check for valid parameter values
    invalid_params <- names(start)[start <= 0]
    if (length(invalid_params) > 0) {
      warning(
        "Initial values for parameters must be positive. Adjusting: ",
        paste(invalid_params, collapse = ", ")
      )
      for (param in invalid_params) {
        start[[param]] <- 1.0 # Set to a reasonable default
      }
    }
  }

  # Apply fixed parameters if any
  if (!is.null(fixed)) {
    # Validate fixed parameters
    fixed_invalid <- names(fixed)[fixed <= 0 & names(fixed) != "delta"]
    if (length(fixed_invalid) > 0) {
      stop("Fixed parameters must be positive: ", paste(fixed_invalid, collapse = ", "))
    }

    # For delta, ensure it's at least 0
    if ("delta" %in% names(fixed) && fixed$delta < 0) {
      stop("Fixed parameter 'delta' must be non-negative")
    }

    for (param in names(fixed)) {
      if (param %in% names(start)) {
        start[[param]] <- fixed[[param]]
      } else if (!param %in% c("alpha", "beta", "gamma", "delta", "lambda")) {
        warning("Ignoring unknown fixed parameter: ", param)
      }
    }
  }

  return(list(start = start, fixed = fixed))
}

#' Determine initial parameter values
#'
#' @param data Numeric vector with values in the (0, 1) interval.
#' @param start Optional list with initial parameter values.
#' @param fixed Optional list of parameters to be held fixed.
#' @param family Character string specifying the distribution family.
#' @param use_moments Logical; if TRUE, uses method of moments for initial values.
#' @param silent Logical; if TRUE, suppresses messages.
#' @return List with initial parameter values.
#' @keywords internal
.determine_start_values <- function(data, start, fixed, family, use_moments, silent) {
  param_info <- .get_family_param_info(family)
  param_names <- param_info$names

  # Get default fixed parameters for the family
  family_fixed <- .get_default_fixed(family)

  # Merge user-provided fixed parameters with family-specific fixed parameters
  if (!is.null(fixed)) {
    fixed <- c(fixed, family_fixed[!names(family_fixed) %in% names(fixed)])
  } else {
    fixed <- family_fixed
  }

  # Generate initial values using method of moments if requested
  if (use_moments) {
    if (!silent) {
      message("Computing starting values using method of moments...")
    }

    # Try to use gkwgetstartvalues with a fallback
    moment_start <- tryCatch(
      {
        gkwgetstartvalues(data, n_starts = 10)
      },
      error = function(e) {
        warning(
          "Error in method of moments estimation: ", e$message,
          ". Using default starting values."
        )
        c(alpha = 2, beta = 2, gamma = 1, delta = 0.5, lambda = 2)
      },
      warning = function(w) {
        warning("Warning in method of moments estimation: ", w$message)
        # Continue execution but still return the result
        NULL
      }
    )

    if (is.null(start) && !is.null(moment_start)) {
      start <- as.list(moment_start)
    } else if (!silent) {
      message("Using provided start values instead of method of moments estimates.")
    }
  }

  # If start is still NULL, use default starting values
  if (is.null(start)) {
    # Get default starting values for the family
    start <- .get_default_start(family)

    if (!silent) {
      message("Using default starting values for ", family, " family parameters")
    }
  }

  # Validate parameters
  valid <- .validate_parameters(start, fixed, param_names)
  start <- valid$start
  fixed <- valid$fixed

  return(list(start = start, fixed = fixed))
}

#' Fit GKw family distributions using Newton-Raphson
#'
#' @param data Numeric vector with values in the (0, 1) interval.
#' @param family Character string specifying the distribution family.
#' @param start List with initial parameter values.
#' @param fixed List of parameters to be held fixed.
#' @param hessian Logical; if TRUE, computes standard errors and covariance matrix.
#' @param conf.level Confidence level for intervals.
#' @param optimizer.control List of control parameters for the optimizer.
#' @param silent Logical; if TRUE, suppresses messages.
#' @return List containing fit results.
#' @keywords internal
.fit_nr <- function(data, family, start, fixed, hessian, conf.level, optimizer.control, silent) {
  # Create full parameter vector with all parameters (fixed and free)
  full_params <- c(1, 1, 1, 0, 1) # Default values for GKw (alpha, beta, gamma, delta, lambda)
  names(full_params) <- c("alpha", "beta", "gamma", "delta", "lambda")

  # Update with start values
  for (param in names(start)) {
    full_params[param] <- start[[param]]
  }

  # Apply fixed parameters
  if (!is.null(fixed)) {
    for (param in names(fixed)) {
      full_params[param] <- fixed[[param]]
    }
  }

  # Extract only the parameters needed for the specific family
  param_names <- .get_family_param_info(family)$names
  start_vec <- full_params[param_names]

  # Set up Newton-Raphson default control parameters
  nr_defaults <- list(
    tol = 1e-6,
    max_iter = 100,
    verbose = !silent,
    use_hessian = hessian,
    step_size = 1.0,
    enforce_bounds = TRUE,
    min_param_val = 1e-5,
    max_param_val = 1e5,
    get_num_hess = !hessian
  )

  # Merge user controls with defaults, with user controls taking precedence
  nr_control <- modifyList(nr_defaults, optimizer.control)

  # Extract control parameters
  tol <- nr_control$tol
  max_iter <- nr_control$max_iter
  verbose <- nr_control$verbose
  use_hessian <- nr_control$use_hessian
  step_size <- nr_control$step_size
  enforce_bounds <- nr_control$enforce_bounds
  min_param_val <- nr_control$min_param_val
  max_param_val <- nr_control$max_param_val
  get_num_hess <- nr_control$get_num_hess

  # Run Newton-Raphson optimization
  nr_result <- tryCatch(
    {
      nrgkw(
        start = start_vec, # Corrected from start_params to start
        data = data,
        family = family,
        tol = tol,
        max_iter = max_iter,
        verbose = verbose,
        use_hessian = use_hessian,
        step_size = step_size,
        enforce_bounds = enforce_bounds,
        min_param_val = min_param_val,
        max_param_val = max_param_val,
        get_num_hess = get_num_hess
      )
    },
    error = function(e) {
      stop("Newton-Raphson optimization failed: ", e$message)
    }
  )

  # Extract and process results
  params <- nr_result$parameters
  names(params) <- param_names

  # We only need the parameters specific to this family
  # No need to convert to a full parameter vector
  filtered_coefficients <- params

  # Process standard errors (will include NAs for fixed parameters)
  filtered_std_errors <- rep(NA, length(param_names))
  names(filtered_std_errors) <- param_names

  if (!is.null(nr_result$std_errors)) {
    for (i in seq_along(param_names)) {
      if (i <= length(nr_result$std_errors)) {
        filtered_std_errors[i] <- nr_result$std_errors[i]
      }
    }
  }

  # Create z-values and p-values only for the relevant parameters
  filtered_z_values <- rep(NA, length(param_names))
  filtered_p_values <- rep(NA, length(param_names))
  names(filtered_z_values) <- names(filtered_p_values) <- param_names

  # Update with z-values and p-values where available
  if (!is.null(nr_result$z_values)) {
    for (i in seq_along(param_names)) {
      if (i <= length(nr_result$z_values) && !is.na(nr_result$z_values[i])) {
        filtered_z_values[i] <- nr_result$z_values[i]
        filtered_p_values[i] <- nr_result$p_values[i]
      }
    }
  }

  # Create coefficient summary only for relevant parameters
  coef_summary <- data.frame(
    Estimate = filtered_coefficients,
    `Std. Error` = filtered_std_errors,
    `z value` = filtered_z_values,
    `Pr(>|z|)` = filtered_p_values,
    row.names = param_names,
    check.names = FALSE
  )

  # Calculate confidence intervals if valid standard errors are available
  conf_int <- NULL
  if (any(!is.na(filtered_std_errors))) {
    z_value <- qnorm(1 - (1 - conf.level) / 2)
    conf_int_params <- character()
    conf_int_estimates <- numeric()
    conf_int_se <- numeric()
    conf_int_lower <- numeric()
    conf_int_upper <- numeric()

    for (i in seq_along(param_names)) {
      param <- param_names[i]
      if (!is.na(filtered_std_errors[i])) {
        conf_int_params <- c(conf_int_params, param)
        conf_int_estimates <- c(conf_int_estimates, filtered_coefficients[i])
        conf_int_se <- c(conf_int_se, filtered_std_errors[i])
        lower <- max(filtered_coefficients[i] - z_value * filtered_std_errors[i], .Machine$double.eps)
        upper <- filtered_coefficients[i] + z_value * filtered_std_errors[i]
        conf_int_lower <- c(conf_int_lower, lower)
        conf_int_upper <- c(conf_int_upper, upper)
      }
    }

    if (length(conf_int_params) > 0) {
      conf_int <- data.frame(
        parameter = conf_int_params,
        estimate = conf_int_estimates,
        std.error = conf_int_se,
        lower = conf_int_lower,
        upper = conf_int_upper,
        row.names = NULL, check.names = FALSE
      )
    }
  }

  # Create comprehensive result object for Newton-Raphson
  result <- list(
    coefficients = filtered_coefficients, # Only the parameters for this family
    std.errors = filtered_std_errors,
    coef_summary = coef_summary,
    vcov = nr_result$hessian, # Note: this is actually the Hessian, not the inverse
    loglik = nr_result$loglik,
    AIC = nr_result$aic,
    BIC = nr_result$bic,
    AICc = nr_result$aic + (2 * length(param_names) * (length(param_names) + 1)) / (length(data) - length(param_names) - 1),
    data = data,
    nobs = length(data),
    df = length(param_names),
    convergence = nr_result$converged,
    message = nr_result$status,
    method = "nr",
    optimizer_method = "newton-raphson",
    conf.int = conf_int,
    conf.level = conf.level,
    optimizer = nr_result,
    fixed = fixed,
    iterations = nr_result$iterations,
    param_history = nr_result$param_history,
    loglik_history = nr_result$loglik_history,
    gradient = nr_result$gradient
  )

  return(result)
}

#' Fit GKw family distributions using TMB
#'
#' @param data Numeric vector with values in the (0, 1) interval.
#' @param family Character string specifying the distribution family.
#' @param start List with initial parameter values.
#' @param fixed List of parameters to be held fixed.
#' @param method Optimization method: "nlminb" or "optim".
#' @param hessian Logical; if TRUE, computes standard errors and covariance matrix.
#' @param conf.level Confidence level for intervals.
#' @param optimizer.control List of control parameters for the optimizer.
#' @param silent Logical; if TRUE, suppresses messages.
#' @return List containing fit results.
#' @keywords internal
.fit_tmb <- function(data, family, start, fixed, method, hessian, conf.level, optimizer.control, silent) {
  if (!requireNamespace("TMB", quietly = TRUE)) {
    stop("Package 'TMB' is required but not installed. Please install it with install.packages('TMB')")
  }

  # dll_name <- "gkwmletmb"

  # Check and compile TMB code
  tryCatch(
    {
      .check_and_compile_TMB_code("gkwmletmb", verbose = silent)
    },
    error = function(e) {
      stop("Failed to compile TMB code: ", e$message)
    }
  )

  if (!silent) {
    message("Checking TMB model compilation...")
  }

  # Get the numeric code for family
  family_code <- .family_to_code(family)

  # Get the parameter names for this family
  param_names <- .get_family_param_info(family)$names

  # Create full parameter list for TMB
  full_start_tmb <- list(
    log_alpha = log(1), # Default values, will be overwritten if needed
    log_beta = log(1),
    log_gamma = log(1),
    log_delta = log(0.0001), # Small value for delta to avoid log(0)
    log_lambda = log(1)
  )

  # Update with user-provided start values
  for (param in names(start)) {
    full_start_tmb[[paste0("log_", param)]] <- log(start[[param]])
  }

  # Apply fixed parameters for TMB
  map <- list()

  # Add family-specific fixed parameters to the map
  all_fixed_params <- c(names(fixed), setdiff(
    c("alpha", "beta", "gamma", "delta", "lambda"),
    c(param_names, names(fixed))
  ))

  for (param in all_fixed_params) {
    param_name <- paste0("log_", param)
    if (param %in% names(fixed)) {
      # User or family provided fixed value
      full_start_tmb[[param_name]] <- log(fixed[[param]])
    } else if (param == "alpha" || param == "beta") {
      # Default α, β = 1 when not specified
      full_start_tmb[[param_name]] <- log(1)
    } else if (param == "gamma") {
      # Default γ = 1 when not specified
      full_start_tmb[[param_name]] <- log(1)
    } else if (param == "delta") {
      # Default δ = 0 when not specified (use small value to avoid log(0))
      full_start_tmb[[param_name]] <- log(0.0001)
    } else if (param == "lambda") {
      # Default λ = 1 when not specified
      full_start_tmb[[param_name]] <- log(1)
    }

    # Add to map to keep parameter fixed
    map[[param_name]] <- factor(NA)
  }

  # Prepare data for TMB
  tmb_data <- list(
    x = data,
    family = family_code
  )

  # Create TMB object
  obj <- tryCatch(
    {
      TMB::MakeADFun(
        data = tmb_data,
        parameters = full_start_tmb,
        map = if (length(map) > 0) map else NULL,
        DLL = "gkwmletmb",
        silent = silent
      )
    },
    error = function(e) {
      stop("Error creating TMB model: ", e$message)
    }
  )

  # Set up optimizer controls based on the chosen method
  if (method == "nlminb") {
    control_defaults <- list(eval.max = 500, iter.max = 300, trace = ifelse(silent, 0, 1))
  } else if (method == "optim") {
    control_defaults <- list(maxit = 500, trace = ifelse(silent, 0, 1))
  } else {
    stop("Unknown optimizer method: ", method)
  }

  # Merge user controls with defaults, with user controls taking precedence
  control <- modifyList(control_defaults, optimizer.control)

  # Run optimization
  if (method == "nlminb") {
    opt <- tryCatch(
      {
        nlminb(
          start = obj$par,
          objective = obj$fn,
          gradient = obj$gr,
          control = control
        )
      },
      error = function(e) {
        stop("Optimization with nlminb failed: ", e$message)
      }
    )

    opt$convergence <- opt$convergence == 0
    opt$message <- opt$message
  } else if (method == "optim") {
    opt <- tryCatch(
      {
        optim(
          par = obj$par,
          fn = obj$fn,
          gr = obj$gr,
          method = "BFGS",
          control = control
        )
      },
      error = function(e) {
        stop("Optimization with optim failed: ", e$message)
      }
    )

    opt$convergence <- opt$convergence == 0
    opt$message <- if (opt$convergence) "Successful convergence" else "Optimization failed to converge"
  }

  if (!opt$convergence) {
    warning("Model did not converge: ", opt$message)
  }

  # Get parameter estimates for all parameters
  all_params <- exp(opt$par)
  names(all_params) <- sub("log_", "", names(all_params))

  # Filter parameters to include only those relevant for the family
  filtered_coefficients <- all_params[paste0("log_", param_names) %in% names(opt$par)]
  names(filtered_coefficients) <- param_names[param_names %in% sub("log_", "", names(opt$par))]

  # If some parameters are missing (fixed), add them from the fixed list
  for (param in param_names) {
    if (!param %in% names(filtered_coefficients)) {
      if (param %in% names(fixed)) {
        filtered_coefficients[param] <- fixed[[param]]
      } else {
        # Use default values for missing parameters
        if (param == "alpha" || param == "beta" || param == "gamma" || param == "lambda") {
          filtered_coefficients[param] <- 1.0
        } else if (param == "delta") {
          filtered_coefficients[param] <- 0.0
        }
      }
    }
  }

  # Ensure the order matches param_names
  filtered_coefficients <- filtered_coefficients[param_names]

  # Calculate standard errors and Hessian if requested
  filtered_std_errors <- rep(NA, length(param_names))
  names(filtered_std_errors) <- param_names

  cov_matrix <- NULL
  coef_summary <- NULL
  conf_int <- NULL

  if (hessian) {
    sd_report <- tryCatch(
      {
        TMB::sdreport(obj)
      },
      error = function(e) {
        warning("Error calculating standard errors: ", e$message)
        NULL
      }
    )

    if (!is.null(sd_report) && !is.character(sd_report)) {
      # Extract parameter estimates, SEs, and covariance matrix
      cov_matrix <- sd_report$cov.fixed
      std_errors_log <- sqrt(diag(cov_matrix))

      # Get the SEs only for non-fixed parameters
      for (i in seq_along(param_names)) {
        param <- param_names[i]
        log_param <- paste0("log_", param)

        if (log_param %in% names(std_errors_log)) {
          # Use the Delta method to transform SEs to original scale
          filtered_std_errors[i] <- std_errors_log[log_param] * filtered_coefficients[i]
        }
      }

      # Create coefficient summary
      coef_summary <- data.frame(
        Estimate = filtered_coefficients,
        `Std. Error` = filtered_std_errors,
        `z value` = filtered_coefficients / filtered_std_errors,
        `Pr(>|z|)` = 2 * pnorm(abs(filtered_coefficients / filtered_std_errors), lower.tail = FALSE),
        row.names = param_names,
        check.names = FALSE
      )

      # Calculate confidence intervals
      z_value <- qnorm(1 - (1 - conf.level) / 2)
      conf_int_params <- character()
      conf_int_estimates <- numeric()
      conf_int_se <- numeric()
      conf_int_lower <- numeric()
      conf_int_upper <- numeric()

      for (i in seq_along(param_names)) {
        if (!is.na(filtered_std_errors[i])) {
          param <- param_names[i]
          conf_int_params <- c(conf_int_params, param)
          conf_int_estimates <- c(conf_int_estimates, filtered_coefficients[i])
          conf_int_se <- c(conf_int_se, filtered_std_errors[i])
          lower <- max(filtered_coefficients[i] - z_value * filtered_std_errors[i], .Machine$double.eps)
          upper <- filtered_coefficients[i] + z_value * filtered_std_errors[i]
          conf_int_lower <- c(conf_int_lower, lower)
          conf_int_upper <- c(conf_int_upper, upper)
        }
      }

      if (length(conf_int_params) > 0) {
        conf_int <- data.frame(
          parameter = conf_int_params,
          estimate = conf_int_estimates,
          std.error = conf_int_se,
          lower = conf_int_lower,
          upper = conf_int_upper,
          row.names = NULL, check.names = FALSE
        )
      }
    } else {
      warning("Hessian calculation failed, standard errors not available")

      # Create coefficient summary without SEs
      coef_summary <- data.frame(
        Estimate = filtered_coefficients,
        `Std. Error` = filtered_std_errors,
        `z value` = rep(NA, length(param_names)),
        `Pr(>|z|)` = rep(NA, length(param_names)),
        row.names = param_names,
        check.names = FALSE
      )
    }
  } else {
    # If hessian = FALSE, only report parameter estimates
    coef_summary <- data.frame(
      Estimate = filtered_coefficients,
      `Std. Error` = filtered_std_errors,
      `z value` = rep(NA, length(param_names)),
      `Pr(>|z|)` = rep(NA, length(param_names)),
      row.names = param_names,
      check.names = FALSE
    )
  }

  # Calculate log likelihood, AIC, BIC
  loglik <- -opt$objective
  n <- length(data)
  k <- sum(!is.na(opt$par)) # Count non-fixed parameters
  aic <- -2 * loglik + 2 * k
  bic <- -2 * loglik + log(n) * k
  aicc <- aic + (2 * k * (k + 1)) / (n - k - 1)

  # Create comprehensive result object for TMB
  result <- list(
    coefficients = filtered_coefficients, # Only the parameters for this family
    std.errors = filtered_std_errors,
    coef_summary = coef_summary,
    vcov = cov_matrix,
    loglik = loglik,
    AIC = aic,
    BIC = bic,
    AICc = aicc,
    data = data,
    nobs = n,
    df = k,
    convergence = opt$convergence,
    message = opt$message,
    method = "tmb",
    optimizer_method = method,
    conf.int = conf_int,
    conf.level = conf.level,
    optimizer = opt,
    obj = obj,
    fixed = fixed
  )

  return(result)
}

#' Calculate profile likelihoods
#'
#' @param result Fit result from TMB or Newton-Raphson.
#' @param data Numeric vector with values in the (0, 1) interval.
#' @param family Character string specifying the distribution family.
#' @param fixed List of parameters to be held fixed.
#' @param fit Estimation method: "tmb" or "nr".
#' @param method Optimization method (for TMB): "nlminb" or "optim".
#' @param npoints Number of points in profile.
#' @param silent Logical; if TRUE, suppresses messages.
#' @return List of profile likelihoods.
#' @keywords internal
.calculate_profiles <- function(result, data, family, fixed, fit, method, npoints, silent) {
  prof_list <- list()

  # Filter out fixed parameters
  param_info <- .get_family_param_info(family)
  params_to_profile <- setdiff(param_info$names, names(fixed))

  # Configure progress reporting
  total_profiles <- length(params_to_profile)
  if (!silent && total_profiles > 0) {
    message("Computing profile likelihoods for ", total_profiles, " parameters...")
  }

  # Process each parameter
  for (param_idx in seq_along(params_to_profile)) {
    param <- params_to_profile[param_idx]

    if (!silent) {
      message("  Computing profile for ", param, " (", param_idx, "/", total_profiles, ")...")
    }

    # Get current parameter estimate
    est_value <- result$coefficients[param]

    # Determine appropriate range for profiling
    if (!is.null(result$std.errors) && !is.na(result$std.errors[param])) {
      # Use standard error to determine range if available
      se <- result$std.errors[param]
      # Calculate range ensuring positive values and reasonable breadth
      min_value <- max(est_value - 3 * se, .Machine$double.eps * 10)
      max_value <- est_value + 3 * se

      profile_range <- seq(min_value, max_value, length.out = npoints)
    } else {
      # If no standard error, use a proportional range
      min_value <- max(est_value * 0.2, .Machine$double.eps * 10)
      max_value <- est_value * 2.0

      profile_range <- seq(min_value, max_value, length.out = npoints)
    }

    # Calculate log-likelihood at each point in the profile
    prof_ll <- numeric(length(profile_range))

    for (i in seq_along(profile_range)) {
      if (fit == "tmb") {
        # For TMB, modify log parameters
        tmb_par <- result$optimizer$par
        log_param <- paste0("log_", param)
        tmb_par[log_param] <- log(profile_range[i])

        # Evaluate negative log-likelihood
        prof_ll[i] <- -result$obj$fn(tmb_par)
      } else {
        # For Newton-Raphson, directly modify the parameter vector
        mod_params <- result$coefficients
        mod_params[param] <- profile_range[i]

        # Transform to the family-specific parameter vector
        param_names <- .get_family_param_info(family)$names
        family_mod_params <- mod_params[param_names]

        # Use the appropriate log-likelihood function based on family
        ll_func <- switch(family,
          "gkw" = llgkw,
          "bkw" = llbkw,
          "kkw" = llkkw,
          "ekw" = llekw,
          "mc" = llbp,
          "kw" = llkw,
          "beta" = llbeta
        )

        # Use ll_func directly (already negated for consistency)
        prof_ll[i] <- ll_func(family_mod_params, data)
      }
    }

    # Create profile data frame
    prof_list[[param]] <- data.frame(
      parameter = param,
      value = profile_range,
      loglik = prof_ll
    )
  }

  return(prof_list)
}

#' Fit submodels for comparison
#'
#' @param data Numeric vector with values in the (0, 1) interval.
#' @param result Main fit result.
#' @param fit Estimation method: "tmb" or "nr".
#' @param method Optimization method (for TMB): "nlminb" or "optim".
#' @param hessian Logical; if TRUE, computes standard errors and covariance matrix.
#' @param optimizer.control List of control parameters for the optimizer.
#' @param silent Logical; if TRUE, suppresses messages.
#' @return List containing submodel fits and LRT results.
#' @keywords internal
.fit_submodels <- function(data, result, fit, method, hessian, optimizer.control, silent) {
  # Only fit submodels if current model is GKw (full model)
  if (result$family != "gkw") {
    if (!silent) {
      message("Submodel fitting is only available when family = 'gkw'. Skipping.")
    }
    return(NULL)
  }

  if (!silent) {
    message("Fitting submodels for comparison...")
  }

  submodel_results <- list()
  submodel_families <- c("bkw", "kkw", "ekw", "mc", "kw", "beta")

  for (submodel in submodel_families) {
    if (!silent) {
      message("  Fitting ", submodel, " submodel...")
    }

    tryCatch(
      {
        # Use current parameter estimates as starting values for submodel
        submodel_start <- as.list(result$coefficients)

        # Subset to include only parameters needed for this family
        submodel_param_names <- .get_family_param_info(submodel)$names
        submodel_start <- submodel_start[submodel_param_names]

        submodel_result <- gkwfit(
          data = data,
          family = submodel,
          start = submodel_start,
          fit = fit,
          method = method,
          use_moments = FALSE,
          hessian = hessian,
          profile = FALSE,
          plot = FALSE,
          silent = TRUE,
          optimizer.control = optimizer.control
        )

        submodel_results[[submodel]] <- submodel_result
      },
      error = function(e) {
        warning("Failed to fit ", submodel, " submodel: ", e$message)
      }
    )
  }

  # Perform likelihood ratio tests
  lrt_results <- list()

  for (submodel in names(submodel_results)) {
    lrt_result <- data.frame(
      Model = paste(submodel, "vs GKw"),
      df = result$df - submodel_results[[submodel]]$df,
      loglik_full = result$loglik,
      loglik_reduced = submodel_results[[submodel]]$loglik,
      LR_statistic = 2 * (result$loglik - submodel_results[[submodel]]$loglik),
      p_value = 1 - pchisq(
        2 * (result$loglik - submodel_results[[submodel]]$loglik),
        result$df - submodel_results[[submodel]]$df
      )
    )
    lrt_results[[submodel]] <- lrt_result
  }

  return(list(submodels = submodel_results, lrt = lrt_results))
}

#' Generate diagnostic plots
#'
#' @param result Fit result from TMB or Newton-Raphson.
#' @param data Numeric vector with values in the (0, 1) interval.
#' @param family Character string specifying the distribution family.
#' @param silent Logical; if TRUE, suppresses messages.
#' @return List of ggplot objects.
#' @keywords internal
.generate_plots <- function(result, data, family, silent) {
  plots <- list()

  # Check if ggplot2 is available
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("Package 'ggplot2' is required for plotting but not installed. Plotting will be disabled.")
    return(NULL)
  } else {
    # Load ggplot2
    require(ggplot2)

    if (!silent) {
      message("Generating diagnostic plots...")
    }

    # Calculate density values on a grid
    x_grid <- seq(0.001, 0.999, length.out = 200)

    # Use the appropriate density function based on family
    density_func <- switch(family,
      "gkw" = function(x) {
        dgkw(
          x, result$coefficients["alpha"], result$coefficients["beta"],
          result$coefficients["gamma"], result$coefficients["delta"],
          result$coefficients["lambda"]
        )
      },
      "bkw" = function(x) {
        dbkw(
          x, result$coefficients["alpha"], result$coefficients["beta"],
          result$coefficients["gamma"], result$coefficients["delta"]
        )
      },
      "kkw" = function(x) {
        dkkw(
          x, result$coefficients["alpha"], result$coefficients["beta"],
          result$coefficients["delta"], result$coefficients["lambda"]
        )
      },
      "ekw" = function(x) {
        dekw(
          x, result$coefficients["alpha"], result$coefficients["beta"],
          result$coefficients["lambda"]
        )
      },
      "mc" = function(x) {
        dbp(
          x, result$coefficients["gamma"], result$coefficients["delta"],
          result$coefficients["lambda"]
        )
      },
      "kw" = function(x) dkw(x, result$coefficients["alpha"], result$coefficients["beta"]),
      "beta" = function(x) dbeta_(x, result$coefficients["gamma"], result$coefficients["delta"])
    )

    density_values <- tryCatch(
      {
        sapply(x_grid, density_func)
      },
      error = function(e) {
        warning("Error calculating density values for plot: ", e$message)
        rep(NA, length(x_grid))
      }
    )

    # Histogram with fitted density overlay
    plot_data <- data.frame(x = x_grid, density = density_values)

    # Create histogram with density curve
    p1 <- ggplot() +
      geom_histogram(aes(x = data, y = after_stat(density)),
        bins = min(30, ceiling(sqrt(length(data)))),
        fill = "lightblue", color = "black", alpha = 0.7
      ) +
      geom_line(
        data = plot_data, aes(x = x, y = density),
        color = "red", size = 1
      ) +
      labs(
        title = paste("Fitted", toupper(family), "Distribution"),
        x = "Data", y = "Density"
      ) +
      theme_minimal()

    plots$histogram <- p1

    # P-P plot (Probability-Probability plot)
    ecdf_vals <- ecdf(data)(sort(data))

    # Use the appropriate CDF function based on family
    cdf_func <- switch(family,
      "gkw" = function(x) {
        pgkw(
          x, result$coefficients["alpha"], result$coefficients["beta"],
          result$coefficients["gamma"], result$coefficients["delta"],
          result$coefficients["lambda"]
        )
      },
      "bkw" = function(x) {
        pbkw(
          x, result$coefficients["alpha"], result$coefficients["beta"],
          result$coefficients["gamma"], result$coefficients["delta"]
        )
      },
      "kkw" = function(x) {
        pkkw(
          x, result$coefficients["alpha"], result$coefficients["beta"],
          result$coefficients["delta"], result$coefficients["lambda"]
        )
      },
      "ekw" = function(x) {
        pekw(
          x, result$coefficients["alpha"], result$coefficients["beta"],
          result$coefficients["lambda"]
        )
      },
      "mc" = function(x) {
        pbp(
          x, result$coefficients["gamma"], result$coefficients["delta"],
          result$coefficients["lambda"]
        )
      },
      "kw" = function(x) pkw(x, result$coefficients["alpha"], result$coefficients["beta"]),
      "beta" = function(x) pbeta_(x, result$coefficients["gamma"], result$coefficients["delta"])
    )

    # Calculate theoretical CDF values
    theor_cdf <- tryCatch(
      {
        sapply(sort(data), cdf_func)
      },
      error = function(e) {
        warning("Error calculating theoretical CDF for P-P plot: ", e$message)
        rep(NA, length(data))
      }
    )

    # Create P-P plot data frame
    pp_data <- data.frame(Empirical = ecdf_vals, Theoretical = theor_cdf)

    # Create P-P plot
    p2 <- ggplot(pp_data, aes(x = Theoretical, y = Empirical)) +
      geom_point() +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
      labs(
        title = "P-P Plot",
        x = "Theoretical Probability", y = "Empirical Probability"
      ) +
      theme_minimal()

    plots$pp_plot <- p2

    # Q-Q plot (Quantile-Quantile plot)
    # Use the appropriate quantile function based on family
    quant_func <- switch(family,
      "gkw" = function(p) {
        qgkw(
          p, result$coefficients["alpha"], result$coefficients["beta"],
          result$coefficients["gamma"], result$coefficients["delta"],
          result$coefficients["lambda"]
        )
      },
      "bkw" = function(p) {
        qbkw(
          p, result$coefficients["alpha"], result$coefficients["beta"],
          result$coefficients["gamma"], result$coefficients["delta"]
        )
      },
      "kkw" = function(p) {
        qkkw(
          p, result$coefficients["alpha"], result$coefficients["beta"],
          result$coefficients["delta"], result$coefficients["lambda"]
        )
      },
      "ekw" = function(p) {
        qekw(
          p, result$coefficients["alpha"], result$coefficients["beta"],
          result$coefficients["lambda"]
        )
      },
      "mc" = function(p) {
        qbp(
          p, result$coefficients["gamma"], result$coefficients["delta"],
          result$coefficients["lambda"]
        )
      },
      "kw" = function(p) qkw(p, result$coefficients["alpha"], result$coefficients["beta"]),
      "beta" = function(p) qbeta_(p, result$coefficients["gamma"], result$coefficients["delta"])
    )

    # Calculate theoretical quantiles
    theor_quant <- tryCatch(
      {
        sapply(ppoints(length(data)), quant_func)
      },
      error = function(e) {
        warning("Error calculating theoretical quantiles for Q-Q plot: ", e$message)
        rep(NA, length(data))
      }
    )

    # Create Q-Q plot data frame
    qq_data <- data.frame(Theoretical = theor_quant, Empirical = sort(data))

    # Create Q-Q plot
    p3 <- ggplot(qq_data, aes(x = Theoretical, y = Empirical)) +
      geom_point() +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
      labs(
        title = "Q-Q Plot",
        x = "Theoretical Quantiles", y = "Empirical Quantiles"
      ) +
      theme_minimal()

    plots$qq_plot <- p3

    # Add profile likelihood plots if available
    if (!is.null(result$profile) && length(result$profile) > 0) {
      for (param in names(result$profile)) {
        prof_data <- result$profile[[param]]

        # Calculate reference line at max - qchisq(0.95, 1)/2 for 95% confidence
        ref_level <- max(prof_data$loglik, na.rm = TRUE) - qchisq(0.95, 1) / 2

        # Create profile likelihood plot
        p <- ggplot(prof_data, aes(x = value, y = loglik)) +
          geom_line(size = 1) +
          geom_vline(
            xintercept = result$coefficients[param],
            linetype = "dashed", color = "red"
          ) +
          geom_hline(
            yintercept = ref_level,
            linetype = "dotted", color = "blue"
          ) +
          labs(
            title = paste("Profile Likelihood for", param),
            x = param, y = "Log-likelihood"
          ) +
          theme_minimal()

        plots[[paste0("profile_", param)]] <- p
      }
    }
  }

  return(plots)
}

#' Calculate goodness-of-fit statistics
#'
#' @param result Fit result from TMB or Newton-Raphson.
#' @param data Numeric vector with values in the (0, 1) interval.
#' @param family Character string specifying the distribution family.
#' @param silent Logical; if TRUE, suppresses messages.
#' @return List of goodness-of-fit statistics.
#' @keywords internal
.calculate_gof <- function(result, data, family, silent) {
  if (!silent) {
    message("Calculating goodness-of-fit statistics...")
  }

  gof <- list()

  # Use the appropriate CDF function based on family for Kolmogorov-Smirnov test
  cdf_func <- switch(family,
    "gkw" = function(x) {
      pgkw(
        x, result$coefficients["alpha"], result$coefficients["beta"],
        result$coefficients["gamma"], result$coefficients["delta"],
        result$coefficients["lambda"]
      )
    },
    "bkw" = function(x) {
      pbkw(
        x, result$coefficients["alpha"], result$coefficients["beta"],
        result$coefficients["gamma"], result$coefficients["delta"]
      )
    },
    "kkw" = function(x) {
      pkkw(
        x, result$coefficients["alpha"], result$coefficients["beta"],
        result$coefficients["delta"], result$coefficients["lambda"]
      )
    },
    "ekw" = function(x) {
      pekw(
        x, result$coefficients["alpha"], result$coefficients["beta"],
        result$coefficients["lambda"]
      )
    },
    "mc" = function(x) {
      pbp(
        x, result$coefficients["gamma"], result$coefficients["delta"],
        result$coefficients["lambda"]
      )
    },
    "kw" = function(x) pkw(x, result$coefficients["alpha"], result$coefficients["beta"]),
    "beta" = function(x) pbeta_(x, result$coefficients["gamma"], result$coefficients["delta"])
  )

  # Kolmogorov-Smirnov test
  ks_test <- suppressWarnings(
    ks.test(data, cdf_func)
  )
  gof$ks <- ks_test

  # Use the appropriate density function based on family
  density_func <- switch(family,
    "gkw" = function(x) {
      dgkw(
        x, result$coefficients["alpha"], result$coefficients["beta"],
        result$coefficients["gamma"], result$coefficients["delta"],
        result$coefficients["lambda"]
      )
    },
    "bkw" = function(x) {
      dbkw(
        x, result$coefficients["alpha"], result$coefficients["beta"],
        result$coefficients["gamma"], result$coefficients["delta"]
      )
    },
    "kkw" = function(x) {
      dkkw(
        x, result$coefficients["alpha"], result$coefficients["beta"],
        result$coefficients["delta"], result$coefficients["lambda"]
      )
    },
    "ekw" = function(x) {
      dekw(
        x, result$coefficients["alpha"], result$coefficients["beta"],
        result$coefficients["lambda"]
      )
    },
    "mc" = function(x) {
      dbp(
        x, result$coefficients["gamma"], result$coefficients["delta"],
        result$coefficients["lambda"]
      )
    },
    "kw" = function(x) dkw(x, result$coefficients["alpha"], result$coefficients["beta"]),
    "beta" = function(x) dbeta_(x, result$coefficients["gamma"], result$coefficients["delta"])
  )

  # Calculate diagnostic statistics
  diagnostics <- list(
    mean_obs = mean(data),
    var_obs = var(data),
    mean_fitted = tryCatch(
      {
        integrate(function(x) x * density_func(x), 0, 1)$value
      },
      error = function(e) NA
    ),
    var_fitted = tryCatch(
      {
        m1 <- integrate(function(x) x * density_func(x), 0, 1)$value
        m2 <- integrate(function(x) x^2 * density_func(x), 0, 1)$value
        m2 - m1^2
      },
      error = function(e) NA
    )
  )

  return(list(gof = gof, diagnostics = diagnostics))
}


#' Fit Generalized Kumaraswamy Distribution via Maximum Likelihood Estimation
#'
#' @description
#' Fits any distribution from the Generalized Kumaraswamy (GKw) family to data using maximum
#' likelihood estimation via TMB (Template Model Builder) or Newton-Raphson optimization.
#' The function supports all seven submodels of the GKw family.
#'
#' @param data A numeric vector with values in the (0, 1) interval.
#' @param family A character string specifying the distribution family. One of: "gkw" (default),
#'        "bkw", "kkw", "ekw", "mc", "kw", or "beta". See Details for parameter specifications.
#' @param start Optional list with initial parameter values. If \code{NULL}, reasonable starting
#'        values will be determined from the data.
#' @param fixed Optional list of parameters to be held fixed (not estimated).
#' @param fit Estimation method to be used: \code{"tmb"} (default) for TMB or \code{"nr"} for Newton-Raphson.
#' @param method (Only for \code{fit = "tmb"}) Optimization method to be used: \code{"nlminb"} (default) or \code{"optim"}.
#' @param use_moments Logical; if \code{TRUE}, uses the method of moments for initial values.
#' @param hessian Logical; if \code{TRUE}, computes standard errors and the covariance matrix.
#' @param profile Logical; if \code{TRUE}, computes likelihood profiles for parameters.
#' @param npoints Number of points to use in profile likelihood calculations.
#' @param plot Logical; if \code{TRUE}, generates diagnostic plots.
#' @param conf.level Confidence level for intervals (default: 0.95).
#' @param optimizer.control List of control parameters passed to the optimizer. Default values
#'        are set internally based on the chosen method.
#' @param submodels Logical; if \code{TRUE}, fits nested submodels for comparison.
#' @param silent Logical; if \code{TRUE}, suppresses messages.
#' @param ... Additional arguments passed to internal functions.
#'
#' @return An object of class \code{"gkwfit"} containing the fitted model. The names of the objects
#' in the output list are preserved for compatibility.
#'
#' @details
#' The \code{gkwfit} function implements fitting for all seven distributions in the Generalized Kumaraswamy family:
#' \itemize{
#'   \item \strong{GKw} (Generalized Kumaraswamy): 5 parameters (α, β, γ, δ, λ)
#'   \item \strong{BKw} (Beta-Kumaraswamy): 4 parameters (α, β, γ, δ), with λ = 1 fixed
#'   \item \strong{KKw} (Kumaraswamy-Kumaraswamy): 4 parameters (α, β, δ, λ), with γ = 1 fixed
#'   \item \strong{EKw} (Exponentiated Kumaraswamy): 3 parameters (α, β, λ), with γ = 1, δ = 0 fixed
#'   \item \strong{Mc} (McDonald/Beta Power): 3 parameters (γ, δ, λ), with α = 1, β = 1 fixed
#'   \item \strong{Kw} (Kumaraswamy): 2 parameters (α, β), with γ = 1, δ = 0, λ = 1 fixed
#'   \item \strong{Beta}: 2 parameters (γ, δ), with α = 1, β = 1, λ = 1 fixed
#' }
#'
#' The function offers two estimation methods:
#' \itemize{
#'   \item \code{fit = "tmb"}: Uses Template Model Builder (TMB) for fitting, with the option to choose
#'         the optimization method via the \code{method} argument, which can be \code{"nlminb"} or \code{"optim"}.
#'   \item \code{fit = "nr"}: Uses the Newton-Raphson method for fitting. In this case, the \code{method}
#'         argument is ignored.
#' }
#'
#' Default values for optimizer control parameters (\code{optimizer.control}) are:
#' \itemize{
#'   \item For \code{"nlminb"}: \code{list(eval.max = 500, iter.max = 300)}.
#'   \item For \code{"optim"}: \code{list(maxit = 500)}.
#'   \item For Newton-Raphson (\code{fit = "nr"}): \code{list(tol = 1e-6, max_iter = 100, step_size = 1.0,
#'         enforce_bounds = TRUE, min_param_val = 1e-5, max_param_val = 1e5)}.
#' }
#'
#' @references
#' Kumaraswamy, P. (1980). A generalized probability density function for double-bounded
#' random processes. Journal of Hydrology, 46(1-2), 79-88.
#'
#' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized distributions.
#' Journal of Statistical Computation and Simulation, 81(7), 883-898.
#'
#' @export
gkwfit <- function(data,
                   family = c("gkw", "bkw", "kkw", "ekw", "mc", "kw", "beta"),
                   start = NULL,
                   fixed = NULL,
                   fit = c("nr", "tmb"),
                   method = c("nlminb", "optim"),
                   use_moments = FALSE,
                   hessian = TRUE,
                   profile = FALSE,
                   npoints = 20,
                   plot = TRUE,
                   conf.level = 0.95,
                   optimizer.control = list(),
                   submodels = FALSE,
                   silent = TRUE,
                   ...) {
  # Match arguments
  family <- match.arg(family, choices = c("gkw", "bkw", "kkw", "ekw", "mc", "kw", "beta"))
  fit <- match.arg(fit, choices = c("nr", "tmb"))
  method <- match.arg(method, choices = c("nlminb", "optim"))
  call <- match.call()

  # Parameter validation
  if (!is.numeric(conf.level) || conf.level <= 0 || conf.level >= 1) {
    stop("conf.level must be a number between 0 and 1")
  }

  if (!is.numeric(npoints) || npoints < 5) {
    stop("npoints must be a positive integer >= 5")
  }
  npoints <- as.integer(npoints)

  # Get parameter information for the family
  param_info <- .get_family_param_info(family)
  n_params <- param_info$n
  param_names <- param_info$names

  # Validate data
  data <- .validate_data(data, n_params)

  # Initialize result list
  result <- list(call = call, family = family)

  # Determine initial values and fixed parameters
  start_info <- .determine_start_values(data, start, fixed, family, use_moments, silent)
  start <- start_info$start
  fixed <- start_info$fixed

  # Fit model based on selected method
  if (fit == "tmb") {
    fit_result <- .fit_tmb(data, family, start, fixed, method, hessian, conf.level, optimizer.control, silent)
  } else {
    fit_result <- .fit_nr(data, family, start, fixed, hessian, conf.level, optimizer.control, silent)
  }

  # Merge fit results into the main result object
  result <- c(result, fit_result)

  # Calculate profile likelihoods if requested
  if (profile) {
    prof_list <- .calculate_profiles(result, data, family, fixed, fit, method, npoints, silent)
    result$profile <- prof_list
  }

  # Fit submodels if requested
  if (submodels) {
    submodel_results <- .fit_submodels(data, result, fit, method, hessian, optimizer.control, silent)
    if (!is.null(submodel_results)) {
      result$submodels <- submodel_results$submodels
      result$lrt <- submodel_results$lrt
    }
  }

  # Generate diagnostic plots if requested
  if (plot) {
    plots <- .generate_plots(result, data, family, silent)
    if (!is.null(plots)) {
      result$plots <- plots
    }
  }

  # Calculate goodness of fit tests and diagnostics
  gof_results <- .calculate_gof(result, data, family, silent)
  result$gof <- gof_results$gof
  result$diagnostics <- gof_results$diagnostics

  # Set the class for S3 methods
  class(result) <- c("gkwfit", "list")

  # Return the final result
  return(result)
}





#' Print method for gkwfit objects
#'
#' @param x A gkwfit object
#' @param digits Number of significant digits to display
#' @param ... Additional arguments passed to print methods
#'
#' @export
print.gkwfit <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("Generalized Kumaraswamy Distribution Fit\n\n")

  cat("Call:\n")
  print(x$call)

  cat("\nCoefficients:\n")
  print.default(format(x$coefficients, digits = digits),
    print.gap = 2,
    quote = FALSE
  )

  cat("\nLog-likelihood:", formatC(x$loglik, digits = digits, format = "f"))
  cat("\nAIC:", formatC(x$AIC, digits = digits, format = "f"))
  cat("\nBIC:", formatC(x$BIC, digits = digits, format = "f"))
  cat("\nNumber of observations:", x$nobs)

  if (!x$convergence) {
    cat("\n\nWarning: Model did not converge\n")
  }

  invisible(x)
}

#' Summary method for gkwfit objects
#'
#' @param object A gkwfit object
#' @param ... Additional arguments
#'
#' @export
summary.gkwfit <- function(object, ...) {
  # Return the object with class "summary.gkwfit"
  class(object) <- c("summary.gkwfit", "gkwfit", "list")
  return(object)
}

#' Print method for summary.gkwfit objects
#'
#' @param x A summary.gkwfit object
#' @param digits Number of significant digits to display
#' @param ... Additional arguments
#'
#' @export
print.summary.gkwfit <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("Generalized Kumaraswamy Distribution Fit\n\n")

  cat("Call:\n")
  print(x$call)

  cat("\nCoefficients:\n")
  printCoefmat(x$coef_summary,
    digits = digits, P.values = TRUE,
    has.Pvalue = TRUE
  )

  cat("\nLog-likelihood:", formatC(x$loglik, digits = digits, format = "f"))
  cat("\nAIC:", formatC(x$AIC, digits = digits, format = "f"))
  cat("\nBIC:", formatC(x$BIC, digits = digits, format = "f"))
  cat("\nAICc:", formatC(x$AICc, digits = digits, format = "f"))
  cat("\nNumber of observations:", x$nobs)

  if (!is.null(x$gof)) {
    cat("\n\nGoodness-of-fit test:\n")
    cat(
      "  Kolmogorov-Smirnov test: D =",
      formatC(x$gof$ks$statistic, digits = digits),
      ", p-value =", formatC(x$gof$ks$p.value, digits = digits)
    )
  }

  if (!is.null(x$conf.int)) {
    cat("\n\nConfidence intervals (", x$conf.level * 100, "%):\n", sep = "")
    ci_table <- x$conf.int[, c("parameter", "estimate", "lower", "upper")]
    rownames(ci_table) <- ci_table$parameter
    ci_table$parameter <- NULL
    colnames(ci_table) <- c("Estimate", "Lower", "Upper")
    print(ci_table, digits = digits)
  }

  if (!is.null(x$lrt) && length(x$lrt) > 0) {
    cat("\nLikelihood ratio tests:\n")
    for (test_name in names(x$lrt)) {
      test <- x$lrt[[test_name]]
      cat("  ", test$Model, ": chi^2 = ",
        formatC(test$LR_statistic, digits = digits),
        ", df = ", test$df,
        ", p-value = ", formatC(test$p_value, digits = digits),
        "\n",
        sep = ""
      )
    }
  }

  if (!x$convergence) {
    cat("\nWarning: Model did not converge\n")
  }

  invisible(x)
}

#' Plot method for gkwfit objects
#'
#' @param x A gkwfit object
#' @param which Which plots to display
#' @param ask Whether to ask before displaying each plot
#' @param ... Additional arguments
#'
#' @export
plot.gkwfit <- function(x,
                        which = c("histogram", "pp", "qq", "profile"),
                        ask = (which != "all" && length(which) > 1 && dev.interactive()),
                        ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting")
  }

  # Determine which plots to show
  if (identical(which, "all")) {
    which <- names(x$plots)
  } else {
    # Map shorthand names to full names
    name_map <- c(
      histogram = "histogram",
      pp = "pp_plot",
      qq = "qq_plot"
    )

    which_full <- character(0)
    for (w in which) {
      if (w == "profile") {
        which_full <- c(which_full, grep("^profile_", names(x$plots), value = TRUE))
      } else if (w %in% names(name_map)) {
        which_full <- c(which_full, name_map[w])
      } else {
        which_full <- c(which_full, w)
      }
    }

    which <- intersect(which_full, names(x$plots))

    if (length(which) == 0) {
      warning("No plots to display")
      return(invisible())
    }
  }

  # Set up plot paging
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }

  # Display the plots
  for (plot_name in which) {
    print(x$plots[[plot_name]])
  }

  invisible()
}

#' Extract log-likelihood from a gkwfit object
#'
#' @param object A gkwfit object
#' @param ... Additional arguments
#'
#' @export
logLik.gkwfit <- function(object, ...) {
  val <- object$loglik
  attr(val, "df") <- object$df
  attr(val, "nobs") <- object$nobs
  class(val) <- "logLik"
  return(val)
}

#' Extract coefficients from a gkwfit object
#'
#' @param object A gkwfit object
#' @param ... Additional arguments
#'
#' @export
coef.gkwfit <- function(object, ...) {
  object$coefficients
}

#' Extract variance-covariance matrix from a gkwfit object
#'
#' @param object A gkwfit object
#' @param ... Additional arguments
#'
#' @export
vcov.gkwfit <- function(object, ...) {
  object$vcov
}

#' Extract confidence intervals from a gkwfit object
#'
#' @param object A gkwfit object
#' @param parm Parameters to extract confidence intervals for
#' @param level Confidence level
#' @param ... Additional arguments
#'
#' @export
confint.gkwfit <- function(object, parm, level = 0.95, ...) {
  if (is.null(object$conf.int)) {
    stop("Confidence intervals not available")
  }

  # Use original conf.int if level matches
  if (abs(object$conf.level - level) < .Machine$double.eps^0.5) {
    ci <- object$conf.int
  } else {
    # Recalculate for different confidence level
    z_value <- qnorm(1 - (1 - level) / 2)

    ci <- object$conf.int
    ci$lower <- pmax(object$coefficients - z_value * object$std.errors, .Machine$double.eps)
    ci$upper <- object$coefficients + z_value * object$std.errors
  }

  # Filter for requested parameters
  if (!missing(parm)) {
    if (is.character(parm)) {
      ci <- ci[ci$parameter %in% parm, ]
    } else if (is.numeric(parm)) {
      ci <- ci[parm, ]
    }
  }

  # Format for output
  result <- ci[, c("parameter", "lower", "upper")]
  rownames(result) <- result$parameter
  result$parameter <- NULL
  colnames(result) <- c("lower", "upper")

  return(result)
}

#' Calculate predicted density values for a gkwfit object
#'
#' @param object A gkwfit object
#' @param newdata Optional new data points
#' @param type Type of prediction: "density" or "cdf"
#' @param ... Additional arguments
#'
#' @export
predict.gkwfit <- function(object, newdata = NULL, type = c("density", "cdf"), ...) {
  type <- match.arg(type)

  # Use original data if newdata not provided
  if (is.null(newdata)) {
    newdata <- object$data
  }

  # Check bounds
  if (any(newdata <= 0 | newdata >= 1)) {
    stop("All values must be in the open interval (0, 1)")
  }

  # Extract parameters
  params <- object$coefficients

  # Calculate density or cdf
  if (type == "density") {
    result <- dgkw(
      newdata,
      params["alpha"],
      params["beta"],
      params["gamma"],
      params["delta"],
      params["lambda"]
    )
  } else {
    result <- pgkw(
      newdata,
      params["alpha"],
      params["beta"],
      params["gamma"],
      params["delta"],
      params["lambda"]
    )
  }

  return(result)
}

#' Calculate AIC for a gkwfit object
#'
#' @param object A gkwfit object
#' @param ... Additional models
#' @param k Penalty per parameter
#'
#' @export
AIC.gkwfit <- function(object, ..., k = 2) {
  objects <- list(object, ...)
  values <- sapply(objects, function(obj) -2 * obj$loglik + k * obj$df)

  if (length(objects) == 1) {
    return(values)
  } else {
    models <- sapply(objects, function(obj) {
      if (is.call(obj$call)) {
        deparse(obj$call)
      } else {
        "Model"
      }
    })

    result <- data.frame(
      df = sapply(objects, function(obj) obj$df),
      AIC = values,
      row.names = models
    )

    # Sort by AIC
    result <- result[order(result$AIC), ]
    return(result)
  }
}

#' Calculate BIC for a gkwfit object
#'
#' @param object A gkwfit object
#' @param ... Additional models
#'
#' @export
BIC.gkwfit <- function(object, ...) {
  objects <- list(object, ...)
  values <- sapply(objects, function(obj) -2 * obj$loglik + log(obj$nobs) * obj$df)

  if (length(objects) == 1) {
    return(values)
  } else {
    models <- sapply(objects, function(obj) {
      if (is.call(obj$call)) {
        deparse(obj$call)
      } else {
        "Model"
      }
    })

    result <- data.frame(
      df = sapply(objects, function(obj) obj$df),
      BIC = values,
      row.names = models
    )

    # Sort by BIC
    result <- result[order(result$BIC), ]
    return(result)
  }
}
