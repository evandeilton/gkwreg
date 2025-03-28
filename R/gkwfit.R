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
#'
#' @importFrom utils modifyList
#'
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
    z_value <- stats::qnorm(1 - (1 - conf.level) / 2)
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
#'
#' @importFrom stats nlminb optim pnorm
#' @importFrom utils modifyList
#'
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
      z_value <- stats::qnorm(1 - (1 - conf.level) / 2)
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
          "mc" = llmc,
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
#'
#' @importFrom stats pchisq
#'
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



#' Generate diagnostic plots for distribution models
#'
#' This internal function creates a set of diagnostic plots for evaluating
#' the fit of various distribution families to bounded data in the (0, 1) interval.
#' It generates histograms with fitted density overlays, P-P plots, Q-Q plots,
#' and profile likelihood plots when available.
#'
#' @param result A list containing model fit results from TMB or Newton-Raphson,
#'        must include a 'coefficients' element with named parameters.
#' @param data Numeric vector with values in the (0, 1) interval.
#' @param family Character string specifying the distribution family.
#'        Supported values: "gkw", "bkw", "kkw", "ekw", "mc", "kw", "beta".
#' @param silent Logical; if TRUE, suppresses messages. Default is FALSE.
#'
#' @return A list of ggplot2 objects:
#'   \item{histogram}{Histogram with fitted density overlay}
#'   \item{pp_plot}{Probability-Probability plot}
#'   \item{qq_plot}{Quantile-Quantile plot}
#'   \item{profile_*}{Profile likelihood plots for each parameter (if available)}
#'
#' @importFrom stats ecdf ppoints qchisq density
#'
#' @keywords internal
.generate_plots <- function(result, data, family, silent = FALSE) {
  plots <- list()

  # Check if ggplot2 is available
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("Package 'ggplot2' is required for plotting but not installed. Plotting will be disabled.")
    return(NULL)
  }

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
      dmc(
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

  # Create data frame for plotting
  plot_data <- data.frame(x = x_grid, density = density_values)

  # Create histogram with density curve
  p1 <- ggplot2::ggplot() +
    ggplot2::geom_histogram(ggplot2::aes(x = data, y = ggplot2::after_stat(density)),
      bins = min(30, ceiling(sqrt(length(data)))),
      fill = "lightblue", color = "black", alpha = 0.7
    ) +
    ggplot2::geom_line(
      data = plot_data, ggplot2::aes(x = x, y = density),
      color = "red", size = 1
    ) +
    ggplot2::labs(
      title = paste("Fitted", toupper(family), "Distribution"),
      x = "Data", y = "Density"
    ) +
    ggplot2::theme_minimal()

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
      pmc(
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
  p2 <- ggplot2::ggplot(pp_data, ggplot2::aes(x = Theoretical, y = Empirical)) +
    ggplot2::geom_point() +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    ggplot2::labs(
      title = "P-P Plot",
      x = "Theoretical Probability", y = "Empirical Probability"
    ) +
    ggplot2::theme_minimal()

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
      qmc(
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
  p3 <- ggplot2::ggplot(qq_data, ggplot2::aes(x = Theoretical, y = Empirical)) +
    ggplot2::geom_point() +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
    ggplot2::labs(
      title = "Q-Q Plot",
      x = "Theoretical Quantiles", y = "Empirical Quantiles"
    ) +
    ggplot2::theme_minimal()

  plots$qq_plot <- p3

  # Add profile likelihood plots if available
  if (!is.null(result$profile) && length(result$profile) > 0) {
    for (param in names(result$profile)) {
      prof_data <- result$profile[[param]]

      # Calculate reference line at max - qchisq(0.95, 1)/2 for 95% confidence
      ref_level <- max(prof_data$loglik, na.rm = TRUE) - qchisq(0.95, 1) / 2

      # Create profile likelihood plot
      p <- ggplot2::ggplot(prof_data, ggplot2::aes(x = value, y = loglik)) +
        ggplot2::geom_line(size = 1) +
        ggplot2::geom_vline(
          xintercept = result$coefficients[param],
          linetype = "dashed", color = "red"
        ) +
        ggplot2::geom_hline(
          yintercept = ref_level,
          linetype = "dotted", color = "blue"
        ) +
        ggplot2::labs(
          title = paste("Profile Likelihood for", param),
          x = param, y = "Log-likelihood"
        ) +
        ggplot2::theme_minimal()

      plots[[paste0("profile_", param)]] <- p
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
#' @importFrom stats ks.test var integrate
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
      pmc(
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
      dmc(
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


#' @title Fit Generalized Kumaraswamy Distribution via Maximum Likelihood Estimation
#'
#' @description
#' Fits any distribution from the Generalized Kumaraswamy (GKw) family to data using maximum
#' likelihood estimation via TMB (Template Model Builder) or Newton-Raphson optimization.
#' The function supports all seven submodels of the GKw family.
#'
#' @param data A numeric vector with values in the (0, 1) interval.
#' @param family A character string specifying the distribution family. One of: \code{"gkw"} (default),
#'   \code{"bkw"}, \code{"kkw"}, \code{"ekw"}, \code{"mc"}, \code{"kw"}, or \code{"beta"}. See Details for parameter specifications.
#' @param start Optional list with initial parameter values. If \code{NULL}, reasonable starting
#'   values will be determined from the data.
#' @param fixed Optional list of parameters to be held fixed (not estimated), e.g., \code{list(lambda = 1)}.
#' @param fit Estimation method to be used: \code{"tmb"} (recommended) for TMB or \code{"nr"} for Newton-Raphson.
#' @param method (Only for \code{fit = "tmb"}) Optimization method to be used by \code{\link[stats]{optim}} or \code{\link[stats]{nlminb}}: \code{"nlminb"} (default) or \code{"optim"}.
#' @param use_moments Logical; if \code{TRUE}, uses the method of moments for initial values.
#' @param hessian Logical; if \code{TRUE}, computes standard errors and the covariance matrix using the observed Hessian matrix.
#' @param profile Logical; if \code{TRUE}, computes likelihood profiles for parameters.
#' @param npoints Integer; number of points to use in profile likelihood calculations (minimum 5).
#' @param plot Logical; if \code{TRUE}, generates diagnostic plots (histogram with fitted density, QQ-plot).
#' @param conf.level Confidence level for intervals (default: 0.95).
#' @param optimizer.control List of control parameters passed to the optimizer (\code{\link[stats]{nlminb}} or \code{\link[stats]{optim}} for \code{fit = "tmb"}, or the internal Newton-Raphson). Default values
#'   are set internally based on the chosen method (see Details).
#' @param submodels Logical; if \code{TRUE}, fits nested submodels for comparison via likelihood ratio tests.
#' @param silent Logical; if \code{TRUE}, suppresses messages during fitting.
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class \code{"gkwfit"} (inheriting from \code{"list"}) containing the fitted model results. Key components include:
#' \item{coefficients}{Named vector of estimated parameters.}
#' \item{vcov}{Variance-covariance matrix of the estimates.}
#' \item{logLik}{Log-likelihood value at the maximum.}
#' \item{aic}{Akaike Information Criterion.}
#' \item{bic}{Bayesian Information Criterion.}
#' \item{convergence}{Convergence status code from the optimizer.}
#' \item{hessian}{Observed Hessian matrix (if \code{hessian = TRUE}).}
#' \item{data}{The input data vector.}
#' \item{family}{The specified distribution family.}
#' \item{call}{The matched function call.}
#' \item{plots}{A list containing ggplot objects for diagnostics (if \code{plot = TRUE}).}
#' \item{profile}{A list containing likelihood profile results (if \code{profile = TRUE}).}
#' \item{submodels}{A list of fitted submodels (if \code{submodels = TRUE}).}
#' \item{lrt}{A list of likelihood ratio test results comparing nested models (if \code{submodels = TRUE}).}
#' \item{gof}{Goodness-of-fit statistics (e.g., AD, CvM).}
#' \item{diagnostics}{Diagnostic information related to GOF tests.}
#' The names of the objects in the output list are preserved for compatibility with potential previous versions.
#'
#' @details
#' The \code{gkwfit} function implements fitting for all seven distributions in the Generalized Kumaraswamy family:
#' \itemize{
#'   \item \strong{GKw} (Generalized Kumaraswamy): 5 parameters (\eqn{\alpha, \beta, \gamma, \delta, \lambda})
#'   \item \strong{BKw} (Beta-Kumaraswamy): 4 parameters (\eqn{\alpha, \beta, \gamma, \delta}), with \eqn{\lambda = 1} fixed
#'   \item \strong{KKw} (Kumaraswamy-Kumaraswamy): 4 parameters (\eqn{\alpha, \beta, \delta, \lambda}), with \eqn{\gamma = 1} fixed
#'   \item \strong{EKw} (Exponentiated Kumaraswamy): 3 parameters (\eqn{\alpha, \beta, \lambda}), with \eqn{\gamma = 1, \delta = 0} fixed
#'   \item \strong{Mc} (McDonald / Beta Power): 3 parameters (\eqn{\gamma, \delta, \lambda}), with \eqn{\alpha = 1, \beta = 1} fixed
#'   \item \strong{Kw} (Kumaraswamy): 2 parameters (\eqn{\alpha, \beta}), with \eqn{\gamma = 1, \delta = 0, \lambda = 1} fixed
#'   \item \strong{Beta}: 2 parameters (\eqn{\gamma, \delta}), with \eqn{\alpha = 1, \beta = 1, \lambda = 1} fixed (equivalent to standard Beta distribution parameters)
#' }
#' All parameters are restricted to be positive.
#'
#' The function offers two estimation methods via the \code{fit} argument:
#' \itemize{
#'   \item \code{fit = "tmb"}: Uses Template Model Builder (TMB) for robust and efficient fitting, leveraging automatic differentiation. The optimization backend (\code{\link[stats]{nlminb}} or \code{\link[stats]{optim}}) can be selected via the \code{method} argument. This is generally the recommended method.
#'   \item \code{fit = "nr"}: Uses a custom Newton-Raphson implementation. In this case, the \code{method} argument is ignored. This might be faster for simple cases but potentially less stable than TMB.
#' }
#'
#' Default values for \code{optimizer.control} depend on the fitting method:
#' \itemize{
#'   \item For \code{fit = "tmb", method = "nlminb"}: \code{list(eval.max = 500, iter.max = 300)}
#'   \item For \code{fit = "tmb", method = "optim"}: \code{list(maxit = 500)} (uses \code{method = "BFGS"} by default within optim)
#'   \item For \code{fit = "nr"}: \code{list(tol = 1e-6, max_iter = 100, step_size = 1.0, enforce_bounds = TRUE, min_param_val = 1e-5, max_param_val = 1e5)}
#' }
#' Users can override these defaults by providing their own list to \code{optimizer.control}.
#'
#' @examples
#' \dontrun{
#' # Ensure the package is loaded (if not already)
#' # library(gkwreg) # Or your package name
#'
#' ## Example 1: Basic Kumaraswamy distribution fitting
#' # Generate sample data from a Kumaraswamy distribution
#' set.seed(123)
#' n <- 200 # Reduced size for faster example
#' # Assuming rkw is part of your package or loaded
#' kw_data <- rkw(n, alpha = 2.5, beta = 1.5)
#'
#' # Fit the Kumaraswamy distribution to the data using TMB
#' kw_fit <- gkwfit(data = kw_data, family = "kw", fit = "tmb")
#'
#' # Display summary of the fitted model (requires summary.gkwfit method)
#' summary(kw_fit)
#'
#' # Show diagnostic plots (requires plot.gkwfit method or access internal list)
#' # print(kw_fit$plots[[1]]) # Example: print histogram plot
#' # print(kw_fit$plots[[2]]) # Example: print QQ plot
#'
#' ## Example 2: Fitting a Generalized Kumaraswamy distribution
#' # Generate sample data from a GKw distribution
#' set.seed(456)
#' # Assuming rgkw is part of your package or loaded
#' gkw_data <- rgkw(n,
#'   alpha = 2.0, beta = 3.0,
#'   gamma = 1.5, delta = 2.5, lambda = 0.8
#' )
#'
#' # Fit the GKw distribution using TMB
#' gkw_fit <- gkwfit(data = gkw_data, family = "gkw", fit = "tmb")
#'
#' # Display parameter estimates with confidence intervals (requires confint.gkwfit)
#' confint(gkw_fit)
#'
#' ## Example 3: Comparing different estimation methods
#' # Generate sample data from a Beta-Kumaraswamy distribution
#' set.seed(789) # Corrected seed from 7809 to 789
#' # Assuming rbkw is part of your package or loaded
#' bkw_data <- rbkw(n, alpha = 1.8, beta = 2.2, gamma = 0.9, delta = 1.2)
#'
#' # Fit using TMB with nlminb optimizer
#' bkw_fit_tmb <- gkwfit(
#'   data = bkw_data, family = "bkw",
#'   fit = "tmb", method = "nlminb"
#' )
#'
#' # Fit using Newton-Raphson method
#' bkw_fit_nr <- gkwfit(data = bkw_data, family = "bkw", fit = "nr")
#'
#' # Compare parameter estimates (requires coef.gkwfit method)
#' print(cbind(TMB = coef(bkw_fit_tmb), NR = coef(bkw_fit_nr)))
#'
#' ## Example 4: Fixing parameters during estimation
#' # Generate data from a Kumaraswamy-Kumaraswamy distribution
#' set.seed(101)
#' # Assuming rkkw is part of your package or loaded
#' kkw_data <- rkkw(n, alpha = 2.5, beta = 1.8, delta = 1.5, lambda = 0.7)
#'
#' # Fit the model with lambda fixed at 0.7
#' kkw_fit <- gkwfit(
#'   data = kkw_data, family = "kkw",
#'   fixed = list(lambda = 0.7)
#' )
#'
#' # Display results (requires coef.gkwfit method)
#' coef(kkw_fit)
#'
#' ## Example 5: Profile likelihoods
#' # Generate data from McDonald distribution (Beta Power)
#' set.seed(202)
#' # Assuming rbp is part of your package or loaded (or maybe rmc?)
#' mc_data <- rmc(n, gamma = 2.0, delta = 1.5, lambda = 0.9) # Assuming rmc exists
#'
#' # Fit the model with profile likelihoods (may take time)
#' mc_fit <- gkwfit(
#'   data = mc_data, family = "mc",
#'   profile = TRUE, npoints = 10 # Reduced points for speed
#' )
#'
#' # Plot profile likelihoods (requires plot.profile.gkwfit method or similar)
#' # plot(mc_fit$profile) # Or specific plotting function for profiles
#'
#' ## Example 6: Fitting and comparing submodels
#' # Generate data from Exponentiated Kumaraswamy distribution
#' set.seed(303)
#' # Assuming rekw is part of your package or loaded
#' ekw_data <- rekw(n, alpha = 1.7, beta = 2.3, lambda = 1.2)
#'
#' # Fit the GKw model and all its submodels (may take significant time)
#' # Consider reducing 'n' for this example
#' ekw_fit_all <- gkwfit(data = ekw_data, family = "gkw", fit = "tmb", submodels = TRUE)
#'
#' # Display likelihood ratio tests comparing GKw to nested models
#' print(do.call(rbind, ekw_fit_all$lrt))
#'
#' ## Example 7: Using method of moments for initialization
#' # Generate data from Beta distribution
#' set.seed(404)
#' beta_data <- stats::rbeta(n, shape1 = 2.0, shape2 = 3.0)
#'
#' # Fit using method of moments for starting values
#' beta_fit <- gkwfit(data = beta_data, family = "beta", use_moments = TRUE)
#'
#' # Check final estimates (requires coef.gkwfit method)
#' coef(beta_fit)
#'
#' ## Example 8: Custom optimizer control
#' # Generate data from Kumaraswamy distribution
#' set.seed(505)
#' kw_data2 <- rkw(n, alpha = 0.8, beta = 1.2)
#'
#' # Fit with custom optimizer settings for more iterations
#' kw_fit2 <- gkwfit(
#'   data = kw_data2, family = "kw", fit = "tmb",
#'   optimizer.control = list(eval.max = 1000, iter.max = 600)
#' )
#'
#' # Check convergence status
#' print(kw_fit2$convergence)
#'
#' ## Example 9: Goodness-of-fit assessment
#' # Generate data from GKw distribution
#' set.seed(606)
#' gkw_data2 <- rgkw(n,
#'   alpha = 1.5, beta = 2.0,
#'   gamma = 1.2, delta = 0.8, lambda = 1.1
#' )
#'
#' # Fit the model
#' gkw_fit2 <- gkwfit(data = gkw_data2, family = "gkw")
#'
#' # Check goodness-of-fit statistics
#' print(gkw_fit2$gof)
#'
#' # Perform a Kolmogorov-Smirnov test against the fitted distribution
#' # Requires pgkw function and coef method
#' ks_result <- stats::ks.test(gkw_data2, function(x) {
#'   pgkw(x,
#'     alpha = coef(gkw_fit2)["alpha"],
#'     beta = coef(gkw_fit2)["beta"],
#'     gamma = coef(gkw_fit2)["gamma"],
#'     delta = coef(gkw_fit2)["delta"],
#'     lambda = coef(gkw_fit2)["lambda"]
#'   )
#' })
#' print(ks_result)
#'
#' ## Example 10: Dealing with data outside (0,1) via transformation
#' set.seed(707)
#' original_data <- stats::rnorm(n, mean = 10, sd = 2)
#'
#' # Transform data to (0,1) interval (common practice)
#' # Add small epsilon to avoid exact 0 or 1 after transformation
#' epsilon <- 1e-6
#' orig_min <- min(original_data)
#' orig_range <- max(original_data) - orig_min
#' transformed_data <- (original_data - orig_min) / orig_range * (1 - 2 * epsilon) + epsilon
#' transformed_data <- pmax(epsilon, pmin(1 - epsilon, transformed_data)) # Ensure bounds
#'
#' # Fit a suitable model (e.g., GKw) to transformed data
#' trans_fit <- gkwfit(data = transformed_data, family = "gkw")
#'
#' # Example: Predict quantiles and back-transform to original scale
#' # Requires qgkw function and coef method
#' fitted_coefs <- coef(trans_fit)
#' quantiles_01 <- qgkw(
#'   p = c(0.1, 0.25, 0.5, 0.75, 0.9),
#'   alpha = fitted_coefs["alpha"],
#'   beta = fitted_coefs["beta"],
#'   gamma = fitted_coefs["gamma"],
#'   delta = fitted_coefs["delta"],
#'   lambda = fitted_coefs["lambda"]
#' )
#' # Back-transform quantiles
#' original_quantiles <- (quantiles_01 - epsilon) / (1 - 2 * epsilon) * orig_range + orig_min
#' print(original_quantiles)
#' }
#' @references
#' Kumaraswamy, P. (1980). A generalized probability density function for double-bounded
#' random processes. *Journal of Hydrology*, 46(1-2), 79-88. \doi{10.1016/0022-1694(80)90036-0}
#'
#' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized distributions.
#' *Journal of Statistical Computation and Simulation*, 81(7), 883-898. \doi{10.1080/00949650903530745}
#'
#' Bourguignon, M., Silva, R. B., & Cordeiro, G. M. (2014). The Weibull-G family of probability distributions. *Journal of Data Science*, 12(1), 53-68. (Mentioned as GKw is related to this framework)
#'
#' @seealso \code{\link{dgkw}}, \code{\link{pgkw}}, \code{\link{qgkw}}, \code{\link{rgkw}}, \code{\link{summary.gkwfit}}, \code{\link{print.gkwfit}}, \code{\link{plot.gkwfit}}, \code{\link{coef.gkwfit}}, \code{\link{vcov.gkwfit}}, \code{\link{logLik.gkwfit}}, \code{\link{confint.gkwfit}}
#' @keywords distribution models mle hplot
#' @author Lopes, J. E.
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



#' @title Print Method for gkwfit Objects
#'
#' @description
#' Prints a concise summary of a model fitted by the \code{\link{gkwfit}} function,
#' displaying the call, estimated coefficients, log-likelihood, AIC, BIC,
#' number of observations, and a convergence warning if applicable.
#'
#' @param x An object of class \code{"gkwfit"}, typically the result of a call to \code{\link{gkwfit}}.
#' @param digits Integer; the minimum number of significant digits to be printed in values.
#'   Defaults to \code{max(3, getOption("digits") - 3)}.
#' @param ... Additional arguments passed to underlying print methods (currently unused).
#'
#' @return Invisibly returns the original input object \code{x}. Called for its side effect of printing to the console.
#'
#' @seealso \code{\link{gkwfit}}, \code{\link{summary.gkwfit}}
#'
#' @examples
#' \dontrun{
#' # Assume 'rkw' function exists to generate Kumaraswamy data
#' # Assume 'gkwfit' function exists to fit the model
#'
#' set.seed(123)
#' kw_data_sample <- rkw(50, alpha = 2.5, beta = 1.5) # Small n for quick example
#'
#' # Fit the model (suppressing plots and messages for cleaner example)
#' fit_object <- gkwfit(
#'   data = kw_data_sample, family = "kw",
#'   plot = FALSE, silent = TRUE, hessian = FALSE
#' ) # Turn off hessian for speed
#'
#' # Print the fitted object summary
#' print(fit_object)
#'
#' # Alternatively, auto-printing works when the object is returned:
#' fit_object
#' }
#'
#' @keywords methods internal
#' @author Lopes, J. E.
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



#' @title Summary Method for gkwfit Objects
#'
#' @description
#' Calculates and prepares a detailed summary of a model fitted by \code{\link{gkwfit}}.
#' This includes coefficients, standard errors, test statistics (z-values), p-values,
#' log-likelihood, information criteria (AIC, BIC), and convergence status.
#'
#' @param object An object of class \code{"gkwfit"}, typically the result of a call to \code{\link{gkwfit}}.
#' @param correlation Logical; if \code{TRUE}, the correlation matrix of the estimated
#'   parameters is computed and included in the summary. Defaults to \code{FALSE}.
#' @param ... Additional arguments (currently unused).
#'
#' @details
#' This function computes standard errors, z-values (\eqn{Estimate / SE}), and
#' p-values (two-tailed test based on the standard normal distribution) for the
#' estimated model parameters, provided that the variance-covariance matrix was
#' successfully computed and is available in the \code{object} (typically requires
#' \code{hessian = TRUE} in the original \code{\link{gkwfit}} call).
#'
#' The returned object is of class \code{"summary.gkwfit"}, and its printing is
#' handled by \code{\link{print.summary.gkwfit}}.
#'
#' @return An object of class \code{"summary.gkwfit"}, which is a list containing:
#' \item{call}{The original function call.}
#' \item{family}{The specified distribution family.}
#' \item{coefficients}{A matrix of estimates, standard errors, z-values, and p-values.}
#' \item{logLik}{The maximized log-likelihood value.}
#' \item{aic}{Akaike Information Criterion.}
#' \item{bic}{Bayesian Information Criterion.}
#' \item{nobs}{Number of observations used in fitting.}
#' \item{convergence_code}{The convergence code returned by the optimizer.}
#' \item{se_available}{Logical indicating if standard errors could be computed.}
#' \item{correlation}{The correlation matrix of coefficients (if \code{correlation = TRUE} and calculable).}
#' \item{fixed}{A list of parameters that were held fixed during estimation, or \code{NULL}.}
#'
#' @seealso \code{\link{gkwfit}}, \code{\link{print.summary.gkwfit}}, \code{\link{coef.gkwfit}}, \code{\link{vcov.gkwfit}}, \code{\link{logLik.gkwfit}}
#'
#' @examples
#' \dontrun{
#' # Assume 'rkw' and 'gkwfit' functions exist
#'
#' set.seed(123)
#' kw_data_sample <- rkw(100, alpha = 2.5, beta = 1.5) # n=100
#'
#' # Fit the model, ensuring hessian=TRUE (default) for SEs
#' fit_object <- gkwfit(
#'   data = kw_data_sample, family = "kw",
#'   plot = FALSE, silent = TRUE
#' )
#'
#' # Calculate the summary object
#' summary_info <- summary(fit_object)
#'
#' # Print the summary (uses print.summary.gkwfit)
#' print(summary_info)
#'
#' # Calculate summary including correlation matrix
#' summary_info_corr <- summary(fit_object, correlation = TRUE)
#' print(summary_info_corr)
#'
#' # Access specific parts of the summary object
#' print(summary_info$coefficients)
#' print(summary_info$logLik)
#' }
#'
#' @keywords methods
#' @author Lopes, J. E.
#' @export
summary.gkwfit <- function(object, correlation = FALSE, ...) {
  if (!inherits(object, "gkwfit")) {
    stop("Input 'object' must be of class 'gkwfit'")
  }

  # Use accessor methods if they exist, otherwise direct access
  coefs <- if (exists("coef.gkwfit", mode = "function")) stats::coef(object) else object$coefficients
  vcov_matrix <- if (exists("vcov.gkwfit", mode = "function")) stats::vcov(object) else object$vcov
  logLik_val <- if (exists("logLik.gkwfit", mode = "function")) stats::logLik(object) else object$logLik
  nobs_val <- if (exists("nobs.gkwfit", mode = "function")) stats::nobs(object) else object$nobs
  aic_val <- if (exists("AIC.gkwfit", mode = "function")) stats::AIC(object) else object$aic
  bic_val <- if (exists("BIC.gkwfit", mode = "function")) stats::BIC(object) else object$bic
  conv_status <- object$convergence # Assuming direct access is okay here
  fixed_params <- object$fixed

  se_available <- FALSE
  coef_table <- NULL

  valid_vcov <- !is.null(vcov_matrix) &&
    is.matrix(vcov_matrix) &&
    !is.null(coefs) && length(coefs) > 0 && # Ensure coefficients exist
    all(dim(vcov_matrix) == length(coefs)) &&
    !anyNA(diag(vcov_matrix)) && # Check for NAs on diagonal
    all(diag(vcov_matrix) >= 0) # Check for non-negative variances

  if (valid_vcov) {
    # Ensure coefs is a named vector if it's not already
    if (is.null(names(coefs)) && length(coefs) == nrow(vcov_matrix) && !is.null(rownames(vcov_matrix))) {
      names(coefs) <- rownames(vcov_matrix)
    } else if (is.null(names(coefs))) {
      # Fallback if names cannot be inferred
      warning("Coefficient names are missing; using generic names.")
      names(coefs) <- paste0("param", seq_along(coefs))
    }

    std_err <- sqrt(diag(vcov_matrix))
    z_vals <- coefs / std_err
    # Handle cases where std_err might be zero -> infinite z -> NA
    z_vals[!is.finite(z_vals)] <- NA
    p_vals <- 2 * stats::pnorm(abs(z_vals), lower.tail = FALSE)

    coef_table <- cbind(
      Estimate = coefs,
      `Std. Error` = std_err,
      `z value` = z_vals,
      `Pr(>|z|)` = p_vals
    )
    se_available <- TRUE
    # Ensure row names match coefficient names
    rownames(coef_table) <- names(coefs)
  } else {
    # Ensure coefs is named even if SEs fail
    if (!is.null(coefs) && is.null(names(coefs))) {
      names(coefs) <- paste0("param", seq_along(coefs))
    } else if (is.null(coefs)) {
      coefs <- numeric(0) # Handle case where coefs might be NULL
    }
    # Create table with NAs if SEs cannot be calculated
    coef_table <- cbind(
      Estimate = coefs,
      `Std. Error` = NA_real_,
      `z value` = NA_real_,
      `Pr(>|z|)` = NA_real_
    )
    rownames(coef_table) <- names(coefs)
  }

  cor_matrix <- NULL
  if (correlation && se_available) {
    # Ensure vcov has non-zero diagonal > machine epsilon for cov2cor stability
    if (all(diag(vcov_matrix) > .Machine$double.eps^0.5)) { # Use sqrt(eps) as tolerance
      cor_matrix <- stats::cov2cor(vcov_matrix)
      # Ensure dimnames match coef names
      dimnames(cor_matrix) <- list(names(coefs), names(coefs))
    } else {
      warning("Cannot compute correlations: near-zero or zero variance estimates found.")
    }
  }

  summary_list <- list(
    call = object$call,
    family = object$family,
    coefficients = coef_table,
    logLik = as.numeric(logLik_val), # Ensure raw number
    aic = aic_val,
    bic = bic_val,
    nobs = as.integer(nobs_val), # Ensure integer
    convergence_code = conv_status,
    se_available = se_available,
    correlation = cor_matrix,
    fixed = fixed_params
  )

  class(summary_list) <- c("summary.gkwfit", "list")
  return(summary_list)
}


#' @title S3 method for class 'summary.gkwfit'
#'
#' @description
#' Prints the formatted summary of a \code{gkwfit} model object, as generated by
#' \code{summary.gkwfit}. It displays the call, family, coefficient table
#' (with standard errors, z-values, p-values, and significance stars), fit statistics
#' (LogLik, AIC, BIC, N), fixed parameters, convergence status, and optionally the
#' correlation matrix of coefficients.
#'
#' @param x An object of class \code{"summary.gkwfit"}, usually the result of \code{summary(gkwfit_object)}.
#' @param digits Integer; the minimum number of significant digits to display in the output.
#' @param signif.stars Logical; if \code{TRUE}, p-values are additionally encoded visually using
#'   "significance stars". Defaults to \code{getOption("show.signif.stars")}.
#' @param ... Additional arguments passed to \code{\link[stats]{printCoefmat}}.
#'
#' @return Invisibly returns the original input object \code{x}. Called primarily for its side effect of printing to the console.
#'
#' @seealso \code{\link{summary.gkwfit}}, \code{\link{gkwfit}}, \code{\link[stats]{printCoefmat}}
#'
#' @keywords methods internal
#' @author Lopes, J. E.
#' @export
print.summary.gkwfit <- function(x, digits = max(3L, getOption("digits") - 3L),
                                 signif.stars = getOption("show.signif.stars", TRUE),
                                 ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Family:", x$family, "\n\n")
  cat("Coefficients:\n")

  stats::printCoefmat(x$coefficients,
    digits = digits, signif.stars = signif.stars,
    na.print = "NA", P.values = x$se_available, has.Pvalue = x$se_available,
    cs.ind = 1:2, tst.ind = 3, zap.ind = integer(), ...
  )

  if (!x$se_available) {
    cat("\n(Standard Errors & Tests cannot be computed or were not requested)\n")
  }

  cat("\n--- Fit Statistics ---\n")
  # Use format() for consistent digit handling respecting user options
  cat("Log-likelihood:", format(x$logLik, digits = digits), "\n")
  cat("AIC:", format(x$aic, digits = digits), ", BIC:", format(x$bic, digits = digits), "\n")
  cat("Number of observations:", x$nobs, "\n")

  if (!is.null(x$fixed) && length(x$fixed) > 0) {
    fixed_text <- paste(names(x$fixed), "=", format(unlist(x$fixed), digits = digits), collapse = ", ")
    cat("Fixed parameters:", fixed_text, "\n")
  }

  conv_code <- x$convergence_code
  if (is.null(conv_code) || is.na(conv_code)) {
    conv_status_text <- "Convergence status unknown"
  } else if (conv_code == 0) {
    # Add check for potential false convergence if applicable? Maybe too complex here.
    conv_status_text <- "Successful convergence"
  } else {
    conv_status_text <- paste("Optimizer potentially did NOT converge (code:", conv_code, ")")
    # Suggest checking optimizer output/message if available in the original fit object
  }
  cat("Convergence:", conv_status_text, "\n")

  if (!is.null(x$correlation)) {
    cat("\nCorrelation of Coefficients:\n")
    cor_mat <- x$correlation
    # Ensure lower triangle only for printing clarity if matrix is large
    # But printCoefmat style often prints full symmetric matrix nicely
    print(cor_mat, digits = digits, ...) # Let print handle it directly
  }

  cat("\n") # Final newline
  invisible(x)
}


#' @title Plot Diagnostics for a gkwfit Object
#'
#' @description
#' Creates a panel of diagnostic plots for assessing the fit of a model estimated
#' by \code{\link{gkwfit}}. It displays a histogram of the data overlaid with the
#' fitted density, a Probability-Probability (P-P) plot, a Quantile-Quantile (Q-Q)
#' plot, and profile likelihood plots for each parameter if they were computed
#' during the fit (i.e., if \code{profile = TRUE} was used in \code{\link{gkwfit}}).
#'
#' @details
#' This function utilizes \code{ggplot2} for creating the plots and \code{patchwork}
#' for arranging them into a single figure. All plots use \code{ggplot2::theme_minimal()}.
#'
#' If the plots were already generated during the original \code{\link{gkwfit}} call
#' (because \code{plot = TRUE}), they are retrieved from the fitted object.
#' Otherwise, this function will attempt to generate the plots on the fly,
#' which requires the \code{ggplot2} package and the necessary distribution
#' functions (like \code{dgkw}, \code{pgkw}, \code{qgkw}, etc.) for the specific
#' \code{family} to be available.
#'
#' The arrangement of plots is handled automatically by \code{patchwork::wrap_plots}.
#' No user interaction (like menu selection) is required.
#'
#' @param x An object of class \code{"gkwfit"}, typically the result of a call to \code{\link{gkwfit}}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns the original input object \code{x}. This function is called for its side effect of producing a plot.
#'
#' @seealso \code{\link{gkwfit}}, \code{\link{summary.gkwfit}}
#' @examples
#' \dontrun{
#' # Assume 'rkw' and 'gkwfit' functions exist
#'
#' set.seed(123)
#' kw_data_sample <- rkw(150, alpha = 2.5, beta = 1.5)
#'
#' # Fit the model (with plot=FALSE initially)
#' fit_object_no_plots <- gkwfit(
#'   data = kw_data_sample, family = "kw",
#'   plot = FALSE, silent = TRUE
#' )
#'
#' # Generate plots using the plot method
#' plot(fit_object_no_plots)
#'
#' # Fit the model including profile likelihoods (with plot=TRUE)
#' fit_object_profiles <- gkwfit(
#'   data = kw_data_sample, family = "kw",
#'   profile = TRUE, npoints = 10, # Fewer points for example speed
#'   plot = TRUE, silent = TRUE
#' )
#'
#' # Display the pre-generated plots using the plot method
#' plot(fit_object_profiles)
#' }
#'
#' @keywords hplot methods
#' @author Lopes, J. E.
#' @export
plot.gkwfit <- function(x, ...) {
  if (!inherits(x, "gkwfit")) {
    stop("Input 'x' must be of class 'gkwfit'")
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. Please install it.")
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required to arrange plots. Please install it.")
  }

  plot_list <- x$plots

  if (is.null(plot_list)) {
    # Attempt to generate plots if not found in the object
    message("Plots not found in the fitted object, generating them now...")

    # Check if necessary components are present in the gkwfit object
    required_comps <- c("coefficients", "data", "family")
    if (!all(required_comps %in% names(x))) {
      stop(
        "Cannot generate plots: the 'gkwfit' object is missing required components ",
        "(e.g., coefficients, data, family). Re-run gkwfit perhaps?"
      )
    }

    # Call the internal plotting function (defined previously)
    # Need to ensure all d/p/q functions are available in the environment
    plot_list <- tryCatch(
      {
        .generate_plots(
          result = x, # Pass the whole object as 'result' argument
          data   = x$data,
          family = x$family,
          silent = FALSE # Show messages from .generate_plots
        )
      },
      error = function(e) {
        warning("Failed to generate plots: ", e$message)
        NULL
      }
    )

    # Handle potential failure in .generate_plots
    if (is.null(plot_list)) {
      warning("Plot generation failed or produced no plots.")
      return(invisible(x)) # Return early
    }
  }

  if (!is.list(plot_list) || length(plot_list) == 0) {
    warning("No plots available in the 'gkwfit' object.")
    return(invisible(x))
  }

  # Filter out any non-ggplot objects or NULLs just in case
  valid_plots <- Filter(function(p) inherits(p, "ggplot"), plot_list)

  if (length(valid_plots) == 0) {
    warning("No valid ggplot objects found to plot.")
    return(invisible(x))
  }

  # Use patchwork::wrap_plots for automatic arrangement in a grid
  # It handles determining the number of columns/rows reasonably well
  combined_plot <- tryCatch(
    {
      patchwork::wrap_plots(valid_plots) +
        patchwork::plot_annotation(
          title = paste("Diagnostic Plots for Fitted", toupper(x$family), "Model")
        )
    },
    error = function(e) {
      warning("Failed to arrange plots using patchwork: ", e$message)
      NULL
    }
  )

  # Print the combined plot if successful
  if (!is.null(combined_plot)) {
    print(combined_plot)
  } else {
    # Fallback if arrangement failed: print plots individually
    warning("Plot arrangement failed. Printing available plots individually.")
    for (i in seq_along(valid_plots)) {
      try(print(valid_plots[[i]]), silent = TRUE)
    }
  }

  invisible(x)
}



#' @title Extract Log-Likelihood from a gkwfit Object
#'
#' @description
#' Extracts the maximized log-likelihood value from a model fitted by \code{\link{gkwfit}}.
#' It returns an object of class \code{"logLik"}, which includes attributes for the
#' degrees of freedom (\code{"df"}) and the number of observations (\code{"nobs"}) used in the fit.
#'
#' @param object An object of class \code{"gkwfit"}, typically the result of a call to \code{\link{gkwfit}}.
#' @param ... Additional arguments (currently ignored).
#'
#' @details
#' This method provides compatibility with standard R functions that operate on
#' log-likelihood values, such as \code{\link[stats]{AIC}}, \code{\link[stats]{BIC}},
#' and likelihood ratio tests. It retrieves the log-likelihood stored during the
#' model fitting process (in \code{object$loglik}) and attaches the required
#' attributes (\code{object$df} for the number of estimated parameters and
#' \code{object$nobs} for the number of observations).
#'
#' @return An object of class \code{"logLik"}. This is the numeric log-likelihood value
#'   with the following attributes:
#'   \item{df}{The number of estimated parameters in the model.}
#'   \item{nobs}{The number of observations used for fitting the model.}
#'
#' @seealso \code{\link{gkwfit}}, \code{\link[stats]{AIC}}, \code{\link[stats]{BIC}}
#'
#' @examples
#' \dontrun{
#' # Assume 'rkw' and 'gkwfit' functions exist
#'
#' set.seed(123)
#' kw_data_sample <- rkw(100, alpha = 2.5, beta = 1.5)
#'
#' # Fit the model
#' fit_obj <- gkwfit(data = kw_data_sample, family = "kw", silent = TRUE)
#'
#' # Extract the log-likelihood object
#' ll <- logLik(fit_obj)
#'
#' # Print the value
#' print(ll)
#'
#' # Show attributes
#' attributes(ll)
#'
#' # Use with standard functions
#' AIC(fit_obj) # Equivalent to AIC(ll)
#' BIC(fit_obj) # Equivalent to BIC(ll)
#'
#' # Extract components directly
#' logLik_value <- as.numeric(ll)
#' degrees_freedom <- attr(ll, "df")
#' num_observations <- attr(ll, "nobs")
#' cat("LogLik:", logLik_value, "\n")
#' cat("df:", degrees_freedom, "\n")
#' cat("nobs:", num_observations, "\n")
#' }
#'
#' @keywords methods models
#' @author Lopes, J. E.
#' @export
logLik.gkwfit <- function(object, ...) {
  # Ensure components exist before accessing
  if (!all(c("loglik", "df", "nobs") %in% names(object))) {
    stop("The 'gkwfit' object is missing required components: 'loglik', 'df', or 'nobs'.")
  }

  val <- object$loglik
  # Basic validation of attributes
  df_val <- object$df
  nobs_val <- object$nobs
  if (!is.numeric(val) || length(val) != 1) stop("'loglik' component must be a single numeric value.")
  if (!is.numeric(df_val) || length(df_val) != 1 || df_val < 0) stop("'df' component must be a single non-negative numeric value.")
  if (!is.numeric(nobs_val) || length(nobs_val) != 1 || nobs_val < 0) stop("'nobs' component must be a single non-negative numeric value.")

  attr(val, "df") <- as.integer(df_val) # Store df as integer
  attr(val, "nobs") <- as.integer(nobs_val) # Store nobs as integer
  class(val) <- "logLik"
  return(val)
}


#' @title Extract Model Coefficients from a gkwfit Object
#'
#' @description
#' Extracts the estimated coefficients for the parameters of a model fitted by
#' \code{\link{gkwfit}}. This is an S3 method for the generic \code{\link[stats]{coef}} function.
#'
#' @param object An object of class \code{"gkwfit"}, typically the result of a call to \code{\link{gkwfit}}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A named numeric vector containing the estimated coefficients for the
#'   parameters of the specified GKw family distribution. The names correspond
#'   to the parameter names (e.g., \code{"alpha"}, \code{"beta"}, etc.).
#'
#' @seealso \code{\link{gkwfit}}, \code{\link[stats]{coef}}, \code{\link{vcov.gkwfit}}, \code{\link{logLik.gkwfit}}
#'
#' @examples
#' \dontrun{
#' # Assume 'rkw' and 'gkwfit' functions exist
#'
#' set.seed(123)
#' kw_data_sample <- rkw(100, alpha = 2.5, beta = 1.5)
#'
#' # Fit the model
#' fit_obj <- gkwfit(data = kw_data_sample, family = "kw", silent = TRUE)
#'
#' # Extract coefficients using the method
#' coefficients_vector <- coef(fit_obj)
#' print(coefficients_vector)
#'
#' # Can also use the generic function
#' stats::coef(fit_obj)
#' }
#'
#' @keywords methods models
#' @author Lopes, J. E.
#' @export
coef.gkwfit <- function(object, ...) {
  # Basic check for existence
  if (is.null(object$coefficients)) {
    warning("Component 'coefficients' not found in the 'gkwfit' object.")
    return(NULL)
  }
  object$coefficients
}


#' @title Extract Variance-Covariance Matrix from a gkwfit Object
#'
#' @description
#' Extracts the variance-covariance matrix of the estimated parameters from a model
#' fitted by \code{\link{gkwfit}}. This matrix is typically derived from the
#' inverse of the observed Hessian matrix calculated during fitting (requires
#' \code{hessian = TRUE} in the \code{\link{gkwfit}} call). This is an S3 method
#' for the generic \code{\link[stats]{vcov}} function.
#'
#' @param object An object of class \code{"gkwfit"}, typically the result of a call to \code{\link{gkwfit}}.
#' @param ... Additional arguments (currently ignored).
#'
#' @return A numeric matrix representing the variance-covariance matrix of the
#'   estimated model parameters. Row and column names correspond to the parameter
#'   names (e.g., \code{"alpha"}, \code{"beta"}, etc.). Returns \code{NULL} or
#'   raises a warning/error if the matrix is not available (e.g., if \code{hessian=FALSE}
#'   was used or if the Hessian computation failed).
#'
#' @seealso \code{\link{gkwfit}}, \code{\link[stats]{vcov}}, \code{\link{coef.gkwfit}}, \code{\link{logLik.gkwfit}}
#'
#' @examples
#' \dontrun{
#' # Assume 'rkw' and 'gkwfit' functions exist
#'
#' set.seed(123)
#' kw_data_sample <- rkw(100, alpha = 2.5, beta = 1.5)
#'
#' # Fit the model, ensuring hessian is computed (default TRUE)
#' fit_obj <- gkwfit(data = kw_data_sample, family = "kw", silent = TRUE)
#'
#' # Extract variance-covariance matrix using the method
#' vcov_matrix <- vcov(fit_obj)
#' print(vcov_matrix)
#'
#' # Can also use the generic function
#' stats::vcov(fit_obj)
#'
#' # Extract standard errors from the diagonal
#' std_errors <- sqrt(diag(vcov_matrix))
#' print(std_errors)
#' }
#'
#' @keywords methods models
#' @author Lopes, J. E.
#' @export
vcov.gkwfit <- function(object, ...) {
  # Basic check for existence
  vcov_mat <- object$vcov
  if (is.null(vcov_mat)) {
    warning(
      "Component 'vcov' (variance-covariance matrix) not found in the 'gkwfit' object.",
      "\nWas the model fitted with hessian=TRUE?"
    )
    return(NULL)
  }
  # Optional: Further checks if vcov_mat is a valid matrix?
  if (!is.matrix(vcov_mat) || !is.numeric(vcov_mat)) {
    warning("'vcov' component is not a numeric matrix.")
    return(NULL)
  }
  # Ensure dimnames match coefficients if possible
  if (!is.null(object$coefficients)) {
    coef_names <- names(object$coefficients)
    if (length(coef_names) == nrow(vcov_mat) && length(coef_names) == ncol(vcov_mat)) {
      dimnames(vcov_mat) <- list(coef_names, coef_names)
    }
  }
  vcov_mat
}


#' @title Compute Confidence Intervals for gkwfit Parameters
#'
#' @description
#' Computes confidence intervals for one or more parameters in a model fitted by
#' \code{\link{gkwfit}}. It uses the Wald method based on the estimated coefficients
#' and their standard errors derived from the variance-covariance matrix.
#'
#' @param object An object of class \code{"gkwfit"}, typically the result of a call to \code{\link{gkwfit}}.
#'   The object must contain valid coefficient estimates and a corresponding
#'   variance-covariance matrix (usually requires fitting with \code{hessian = TRUE}).
#' @param parm A specification of which parameters are to be given confidence intervals,
#'   either a vector of numbers (indices) or a vector of names. If missing,
#'   confidence intervals are computed for all parameters that have a valid
#'   standard error available. Parameter indices refer to the order of parameters
#'   for which standard errors could be calculated.
#' @param level The confidence level required (default: 0.95).
#' @param ... Additional arguments (currently ignored).
#'
#' @details
#' This function calculates confidence intervals using the Wald method:
#' \eqn{Estimate \pm z \times SE}, where \eqn{z} is the appropriate quantile
#' from the standard normal distribution for the given confidence `level`.
#'
#' It relies on the results from \code{\link{coef.gkwfit}} and \code{\link{vcov.gkwfit}}
#' (or directly accesses \code{object$coefficients} and \code{object$vcov} if those
#' methods aren't defined). It checks for the validity of the variance-covariance
#' matrix before proceeding.
#'
#' Since all parameters of the GKw family distributions are constrained to be positive,
#' the lower bound of the confidence interval is truncated at a small positive value
#' (\code{.Machine$double.eps^0.5}) if the calculated lower bound is non-positive.
#'
#' If `parm` is specified, it selects the parameters for which to compute intervals.
#' Numeric indices in `parm` refer to the parameters that have calculable standard
#' errors, not necessarily all parameters in the model (if some were fixed or had
#' estimation issues).
#'
#' @return A matrix with columns giving lower and upper confidence limits for each
#'   parameter specified in `parm`. The columns are labeled with quantile percentages
#'   (e.g., \code{"2.5 %"} and \code{"97.5 %"} for \code{level = 0.95}). Row names
#'   are taken from the parameter names. Returns \code{NULL} or stops with an error
#'   if coefficients or a valid variance-covariance matrix cannot be extracted.
#'
#' @seealso \code{\link{gkwfit}}, \code{\link[stats]{confint}}, \code{\link{coef.gkwfit}}, \code{\link{vcov.gkwfit}}
#'
#' @examples
#' \dontrun{
#' # Assume 'rkw' and 'gkwfit' functions exist
#'
#' set.seed(123)
#' kw_data_sample <- rkw(100, alpha = 2.5, beta = 1.5)
#'
#' # Fit the model, ensuring hessian is computed (default TRUE)
#' fit_obj <- gkwfit(data = kw_data_sample, family = "kw", silent = TRUE)
#'
#' # Calculate default 95% confidence intervals for all parameters
#' ci_default <- confint(fit_obj)
#' print(ci_default)
#'
#' # Calculate 90% confidence intervals
#' ci_90 <- confint(fit_obj, level = 0.90)
#' print(ci_90)
#'
#' # Calculate interval for a specific parameter by name
#' ci_alpha <- confint(fit_obj, parm = "alpha")
#' print(ci_alpha)
#'
#' # Calculate interval for specific parameters by index
#' # (Indices refer to parameters with SEs, here alpha=1, beta=2)
#' ci_indexed <- confint(fit_obj, parm = 1:2) # Should be same as default for Kw
#' print(ci_indexed)
#'
#' # Fit a model with a fixed parameter (example for KKw)
#' kkw_data_sample <- rkkw(100, alpha = 2, beta = 3, delta = 1.5, lambda = 0.8) # Assume rkkw
#' fit_kkw_fixed <- gkwfit(kkw_data_sample, family = "kkw", fixed = list(lambda = 0.8))
#' ci_kkw_fixed <- confint(fit_kkw_fixed) # Should only show CIs for alpha, beta, delta
#' print(ci_kkw_fixed)
#' }
#'
#' @keywords methods models
#' @author Lopes, J. E.
#' @export
confint.gkwfit <- function(object, parm, level = 0.95, ...) {
  if (!inherits(object, "gkwfit")) {
    stop("Input 'object' must be of class 'gkwfit'.")
  }
  if (!is.numeric(level) || length(level) != 1 || level <= 0 || level >= 1) {
    stop("'level' must be a single numeric value between 0 and 1.")
  }

  # Prefer methods if they exist and work
  cf <- tryCatch(stats::coef(object), error = function(e) object$coefficients)
  vc <- tryCatch(stats::vcov(object), error = function(e) object$vcov)

  if (is.null(cf)) {
    stop("Could not extract coefficients ('coefficients') from the object.")
  }
  if (is.null(vc)) {
    stop(
      "Could not extract variance-covariance matrix ('vcov') from the object.",
      "\nWas the model fitted with hessian=TRUE?"
    )
  }

  # --- Validate VCOV and Calculate SE ---
  coef_names <- names(cf)
  if (is.null(coef_names)) {
    stop("Coefficients must be named.")
  }
  if (!is.matrix(vc) || nrow(vc) != ncol(vc) || nrow(vc) != length(cf) ||
    is.null(rownames(vc)) || is.null(colnames(vc)) || # Check names exist
    !all(rownames(vc) == coef_names) || !all(colnames(vc) == coef_names)) {
    stop(
      "Variance-covariance matrix ('vcov') has unexpected dimensions or names",
      " relative to coefficients."
    )
  }
  variances <- diag(vc)
  names(variances) <- coef_names # Ensure names consistency
  if (anyNA(variances)) {
    warning(
      "NA values found in the diagonal of the vcov matrix.",
      " Parameters with NA variance will have NA intervals."
    )
  }
  if (any(variances[!is.na(variances)] < 0)) {
    warning(
      "Negative variances found in vcov matrix.",
      " Affected parameters will have NA intervals."
    )
    variances[variances < 0] <- NA # Set negative variances to NA
  }
  ses <- sqrt(variances) # Will produce NAs where variances were NA or negative

  # --- Determine Parameters for Intervals ---
  # Identify parameters with valid (non-NA) standard errors
  param_names_with_se <- coef_names[!is.na(ses)]

  if (length(param_names_with_se) == 0) {
    warning(
      "No parameters with valid standard errors available ",
      "for confidence interval calculation."
    )
    # Return an empty matrix with correct column names
    a <- (1 - level) / 2
    a <- c(a, 1 - a)
    pct <- paste(format(100 * a, trim = TRUE, scientific = FALSE, digits = 3), "%")
    return(matrix(NA_real_, nrow = 0, ncol = 2, dimnames = list(NULL, pct)))
  }

  # Select parameters based on 'parm' argument
  if (missing(parm)) {
    params_to_compute <- param_names_with_se # Default: all with valid SE
  } else {
    if (is.numeric(parm)) {
      # Check numeric indices against the list of parameters *with SEs*
      if (any(parm <= 0) || any(parm > length(param_names_with_se))) {
        stop(
          "Numeric 'parm' indices are out of bounds (1 to ",
          length(param_names_with_se), ") for parameters with available SEs."
        )
      }
      params_to_compute <- param_names_with_se[parm]
    } else if (is.character(parm)) {
      unknown_parm <- setdiff(parm, coef_names) # Check against all coef names
      if (length(unknown_parm) > 0) {
        stop(
          "Unknown parameter(s) requested in 'parm': ",
          paste(dQuote(unknown_parm), collapse = ", ")
        )
      }
      # Filter requested parameters to only those with valid SEs
      params_to_compute <- intersect(parm, param_names_with_se)
      if (length(params_to_compute) == 0) {
        # Check if the requested parameters existed but just lacked SEs
        requested_but_no_se <- intersect(parm, setdiff(coef_names, param_names_with_se))
        if (length(requested_but_no_se) > 0) {
          stop(
            "None of the requested parameters have available standard errors. ",
            "Problematic parameters might include: ",
            paste(dQuote(requested_but_no_se), collapse = ", ")
          )
        } else {
          # Should not happen given the unknown_parm check, but as fallback
          stop("No valid parameters selected for confidence intervals.")
        }
      }
    } else {
      stop("'parm' must be missing, a numeric vector, or a character vector.")
    }
  }

  # --- Calculate Confidence Intervals ---
  a <- (1 - level) / 2
  a <- c(a, 1 - a)
  pct <- paste(format(100 * a, trim = TRUE, scientific = FALSE, digits = 3), "%") # Format column names

  # Z-critical value for normal distribution
  z_crit <- stats::qnorm(a[2])

  # Subset coefficients and SEs for the selected parameters
  cf_sub <- cf[params_to_compute]
  ses_sub <- ses[params_to_compute] # NAs handled by calculation below

  # Calculate intervals matrix: Estimate +/- z * SE
  ci <- cf_sub + ses_sub %o% c(-z_crit, z_crit)

  # Ensure parameters are positive (all GKw parameters should be)
  # Using small positive value instead of exact zero
  # Apply only where the original estimate was positive, otherwise result is NA anyway
  ci <- ifelse(!is.na(ci) & ci < .Machine$double.eps^0.5, .Machine$double.eps^0.5, ci)

  # --- Format Output Matrix ---
  dimnames(ci) <- list(params_to_compute, pct) # Assign row and column names

  return(ci)
}



#' @title Predict Method for gkwfit Objects
#'
#' @description
#' Computes various types of predictions based on a model fitted by the
#' \code{\link{gkwfit}} function. This includes predicted density (pdf),
#' cumulative probability (CDF), quantiles, or the theoretical mean response
#' of the fitted distribution.
#'
#' @param object An object of class \code{"gkwfit"}, typically the result of a
#'   call to \code{\link{gkwfit}}.
#' @param newdata An optional numeric vector of values strictly between 0 and 1
#'   at which to compute predictions for \code{type = "density"} or
#'   \code{type = "cdf"}. If \code{NULL} (the default) for these types, an
#'   error is raised as specific points are required. For \code{type = "response"},
#'   if \code{newdata} is provided, the calculated mean is replicated to match
#'   the length of \code{newdata}; otherwise (if \code{newdata = NULL}), a
#'   single mean value is returned. For \code{type = "quantile"}, this argument
#'   is ignored.
#' @param type A character string specifying the type of prediction required. Valid options are:
#'   \itemize{
#'     \item \code{"density"}: Computes the probability density function (pdf) values.
#'       Requires the \code{newdata} argument.
#'     \item \code{"cdf"}: Computes the cumulative distribution function (CDF) values.
#'       Requires the \code{newdata} argument.
#'     \item \code{"quantile"}: Computes the quantiles for given probabilities.
#'       Requires the \code{p} argument. Ignores \code{newdata}.
#'     \item \code{"response"}: Computes the theoretical mean of the fitted distribution.
#'       Does not strictly require \code{newdata}, but if provided, the result is
#'       replicated to match its length.
#'   }
#' @param p A numeric vector of probabilities (values between 0 and 1 inclusive)
#'   for which to compute quantiles. Required and only used when \code{type = "quantile"}.
#' @param ... Currently unused arguments. Included for consistency with the
#'   generic \code{\link[stats]{predict}} method.
#'
#' @details
#' This function acts as an interface to the underlying distribution functions
#' corresponding to the fitted model's family (specified in \code{object$family}).
#' It dynamically selects the appropriate density function (e.g., \code{dgkw}, \code{dkw}),
#' CDF function (e.g., \code{pgkw}, \code{pkw}), quantile function (e.g., \code{qgkw}, \code{qkw}),
#' or potentially a mean function (e.g., \code{meangkw}, \code{meankw}) based on the
#' requested \code{type} and the fitted \code{family}.
#'
#' The function retrieves the estimated parameters using \code{coef(object)} and
#' passes them along with \code{newdata} or \code{p} to the relevant underlying
#' distribution function.
#'
#' **Important:** For this function to work correctly, the corresponding `d`, `p`, `q`
#' (and `mean` if `type = "response"` is used, unless it's the standard Beta)
#' functions for the specific \code{family} fitted (e.g., `dkw`, `pkw`, `qkw` for the
#' `"kw"` family) must be available in the environment where \code{predict.gkwfit}
#' is called (typically meaning they should be defined and exported by the same package).
#'
#' Input values in \code{newdata} must be strictly within the open interval (0, 1).
#' Probabilities in \code{p} must be within the closed interval (0, 1).
#'
#' @return The type of return value depends on the \code{type} argument:
#' \itemize{
#'   \item If \code{type = "density"} or \code{type = "cdf"}: A numeric vector of
#'     density or CDF values, respectively, corresponding to each value in \code{newdata}.
#'     The length matches the length of \code{newdata}.
#'   \item If \code{type = "quantile"}: A numeric vector of quantile values
#'     corresponding to each probability in \code{p}. The length matches the
#'     length of \code{p}.
#'   \item If \code{type = "response"}: A single numeric value representing the
#'     theoretical mean of the fitted distribution if \code{newdata = NULL}. If
#'     \code{newdata} was provided, a numeric vector repeating the theoretical mean,
#'     with length matching the length of \code{newdata}.
#' }
#'
#' @seealso \code{\link{gkwfit}}, \code{\link[stats]{predict}}, and the specific
#'   distribution functions used (e.g., \code{dgkw}, \code{pgkw}, \code{qgkw},
#'   \code{dkw}, \code{pkw}, \code{qkw}, etc.).
#'
#' @examples
#' \dontrun{
#' # Assume necessary functions like 'rkw', 'gkwfit',
#' # 'dkw', 'pkw', 'qkw', and 'meankw' (or equivalent logic) exist.
#'
#' set.seed(123)
#' kw_data_sample <- rkw(100, alpha = 2.5, beta = 1.5)
#' fit_obj <- gkwfit(data = kw_data_sample, family = "kw", silent = TRUE)
#'
#' # Example points for prediction
#' new_pts <- c(0.1, 0.25, 0.5, 0.75, 0.9)
#'
#' # --- Predict Density ---
#' pred_d <- predict(fit_obj, newdata = new_pts, type = "density")
#' print(pred_d)
#' # plot(new_pts, pred_d, type = "l", main = "Predicted Density")
#'
#' # --- Predict CDF ---
#' pred_p <- predict(fit_obj, newdata = new_pts, type = "cdf")
#' print(pred_p)
#' # plot(new_pts, pred_p, type = "l", main = "Predicted CDF")
#'
#' # --- Predict Quantiles ---
#' probabilities <- c(0.01, 0.05, 0.5, 0.95, 0.99)
#' pred_q <- predict(fit_obj, p = probabilities, type = "quantile")
#' print(pred_q)
#'
#' # --- Predict Mean Response ---
#' # Get the single theoretical mean value
#' pred_mean_single <- predict(fit_obj, type = "response")
#' print(pred_mean_single)
#'
#' # Get the mean replicated for the length of newdata
#' pred_mean_vector <- predict(fit_obj, newdata = new_pts, type = "response")
#' print(pred_mean_vector)
#'
#' # --- Example with a different family (assuming setup exists) ---
#' # set.seed(456)
#' # gkw_data <- rgkw(100, alpha=2, beta=3, gamma=1.5, delta=2.5, lambda=0.8)
#' # fit_gkw <- gkwfit(gkw_data, family="gkw", silent=TRUE)
#' # pred_gkw_cdf <- predict(fit_gkw, newdata = new_pts, type = "cdf")
#' # print(pred_gkw_cdf)
#' }
#'
#' @keywords methods models prediction
#' @author Lopes, J. E.
#' @export
predict.gkwfit <- function(object, newdata = NULL,
                           type = c("density", "cdf", "quantile", "response"),
                           p = NULL, # Probability for quantile type
                           ...) {
  if (!inherits(object, "gkwfit")) {
    stop("Input 'object' must be of class 'gkwfit'.")
  }
  type <- match.arg(type)

  family <- object$family
  params <- stats::coef(object) # Use coef extractor
  if (is.null(params) || length(params) == 0) {
    stop("Could not extract coefficients from the object.")
  }
  param_names_fit <- names(params) # Names of parameters actually fitted

  result <- NULL

  switch(type,
    "density" = {
      if (is.null(newdata)) stop("'newdata' must be provided for type = 'density'.")
      if (!is.numeric(newdata)) stop("'newdata' must be numeric.")
      # Validate bounds slightly relaxed to avoid floating point issues at edges
      eps <- .Machine$double.eps^0.5
      if (any(newdata <= eps | newdata >= (1 - eps))) {
        warning("Some 'newdata' values are very close to or outside (0, 1). Results may be NA/NaN/Inf or unreliable.")
        # Allow calculation but warn
      }

      density_func_name <- paste0("d", family)
      # Special case for beta? Assuming dbeta_ exists
      if (family == "beta") density_func_name <- "dbeta_"

      if (!exists(density_func_name, mode = "function")) {
        stop("Density function '", density_func_name, "' not found for family '", family, "'.")
      }
      density_func <- get(density_func_name, mode = "function")

      # Prepare arguments, matching fitted params to function args
      func_args <- list(newdata) # First argument is usually the data/x
      # Add only the parameters relevant for this family
      func_args <- c(func_args, params)

      result <- tryCatch(do.call(density_func, func_args),
        error = function(e) {
          stop("Error calling ", density_func_name, ": ", e$message)
        }
      )
    },
    "cdf" = {
      if (is.null(newdata)) stop("'newdata' must be provided for type = 'cdf'.")
      if (!is.numeric(newdata)) stop("'newdata' must be numeric.")
      eps <- .Machine$double.eps^0.5
      if (any(newdata <= eps | newdata >= (1 - eps))) {
        warning("Some 'newdata' values are very close to or outside (0, 1). Results may be less reliable at boundaries.")
        # Clamp newdata for CDF calculation robustness
        newdata <- pmax(eps, pmin(1 - eps, newdata))
      }


      cdf_func_name <- paste0("p", family)
      if (family == "beta") cdf_func_name <- "pbeta_" # Assuming pbeta_

      if (!exists(cdf_func_name, mode = "function")) {
        stop("CDF function '", cdf_func_name, "' not found for family '", family, "'.")
      }
      cdf_func <- get(cdf_func_name, mode = "function")

      func_args <- list(newdata)
      func_args <- c(func_args, params)

      result <- tryCatch(do.call(cdf_func, func_args),
        error = function(e) {
          stop("Error calling ", cdf_func_name, ": ", e$message)
        }
      )
    },
    "quantile" = {
      if (is.null(p)) stop("'p' (probabilities) must be provided for type = 'quantile'.")
      if (!is.numeric(p) || any(p < 0 | p > 1)) {
        stop("'p' must be numeric values between 0 and 1.")
      }

      quantile_func_name <- paste0("q", family)
      if (family == "beta") quantile_func_name <- "qbeta_" # Assuming qbeta_

      if (!exists(quantile_func_name, mode = "function")) {
        stop("Quantile function '", quantile_func_name, "' not found for family '", family, "'.")
      }
      quantile_func <- get(quantile_func_name, mode = "function")

      func_args <- list(p)
      func_args <- c(func_args, params)

      result <- tryCatch(do.call(quantile_func, func_args),
        error = function(e) {
          stop("Error calling ", quantile_func_name, ": ", e$message)
        }
      )
    },
    "response" = {
      # Calculate the theoretical mean of the distribution
      # Assuming functions like meangkw, meanbkw, meankw, meanbeta_ exist
      mean_func_name <- paste0("mean", family)
      if (family == "beta") mean_func_name <- "meanbeta_" # Or calculate directly: g/(g+d)

      if (!exists(mean_func_name, mode = "function")) {
        # Fallback for beta if meanbeta_ doesn't exist
        if (family == "beta" && all(c("gamma", "delta") %in% names(params))) {
          result <- params["gamma"] / (params["gamma"] + params["delta"])
        } else {
          stop(
            "Mean function '", mean_func_name, "' not found for family '", family,
            "' and type = 'response'. Cannot calculate predicted mean."
          )
        }
      } else {
        mean_func <- get(mean_func_name, mode = "function")
        # Prepare arguments - mean functions usually just take parameters
        func_args <- params

        result <- tryCatch(do.call(mean_func, func_args),
          error = function(e) {
            stop("Error calling ", mean_func_name, ": ", e$message)
          }
        )
      }
      # Mean is typically a single value for the distribution
      # If newdata was provided, maybe replicate the mean? Or just return the single value?
      # Standard predict methods often return predictions corresponding to newdata.
      # Let's return a vector matching newdata length if newdata is provided.
      if (!is.null(newdata)) {
        if (!is.numeric(newdata)) stop("'newdata' must be numeric if provided for type='response'.")
        result <- rep(result, length.out = length(newdata))
      } else {
        # If newdata is NULL, return the single mean value
      }
    },
    # Default case for switch
    stop("Invalid 'type' specified.")
  )

  return(result)
}


#' @title Calculate AIC or BIC for gkwfit Objects
#'
#' @description
#' Computes the Akaike Information Criterion (AIC) or variants like the Bayesian
#' Information Criterion (BIC) for one or more fitted model objects of class \code{"gkwfit"}.
#'
#' @param object An object of class \code{"gkwfit"}, typically the result of a call to \code{\link{gkwfit}}.
#' @param ... Optionally, more fitted model objects of class \code{"gkwfit"}.
#' @param k Numeric scalar specifying the penalty per parameter. The default \code{k = 2}
#'   corresponds to the traditional AIC. Use \code{k = log(n)} (where n is the number
#'   of observations) for the BIC (Bayesian Information Criterion).
#'
#' @details
#' This function calculates an information criterion based on the formula
#' \eqn{-2 \times \log Likelihood + k \times df}, where \eqn{df} represents the
#' number of estimated parameters in the model (degrees of freedom).
#'
#' It relies on the \code{\link{logLik.gkwfit}} method to extract the log-likelihood
#' and the degrees of freedom for each model.
#'
#' When comparing multiple models fitted to the **same data**, the model with the
#' lower AIC or BIC value is generally preferred. The function returns a sorted
#' data frame to facilitate this comparison when multiple objects are provided.
#'
#' @return
#' \itemize{
#'   \item If only one \code{object} is provided: A single numeric value representing the calculated criterion (AIC or BIC).
#'   \item If multiple objects are provided: A \code{data.frame} with rows corresponding
#'     to the models and columns for the degrees of freedom (\code{df}) and the
#'     calculated criterion value (named \code{AIC}, regardless of the value of \code{k}).
#'     The data frame is sorted in ascending order based on the criterion values.
#'     Row names are derived from the deparsed calls of the fitted models.
#' }
#'
#' @seealso \code{\link{gkwfit}}, \code{\link[stats]{AIC}}, \code{\link{logLik.gkwfit}}, \code{\link{BIC.gkwfit}}
#'
#' @examples
#' \dontrun{
#' # Assume necessary functions (rkw, rgkw, gkwfit, logLik.gkwfit) exist
#'
#' set.seed(123)
#' sample_data <- rkw(100, alpha = 2.5, beta = 1.5)
#'
#' # Fit different models to the same data
#' fit1_kw <- gkwfit(sample_data, family = "kw", silent = TRUE)
#' fit2_bkw <- gkwfit(sample_data, family = "bkw", silent = TRUE) # Assuming bkw fits
#' fit3_gkw <- gkwfit(sample_data, family = "gkw", silent = TRUE) # Assuming gkw fits
#'
#' # Calculate AIC for a single model
#' aic1 <- AIC(fit1_kw)
#' print(aic1)
#'
#' # Compare AIC values for multiple models
#' aic_comparison <- AIC(fit1_kw, fit2_bkw, fit3_gkw)
#' print(aic_comparison) # Sorted by AIC
#'
#' # Compare BIC values for multiple models
#' # Get nobs from one of the models (must be fitted to same data)
#' n_obs <- attr(logLik(fit1_kw), "nobs")
#' bic_comparison <- AIC(fit1_kw, fit2_bkw, fit3_gkw, k = log(n_obs))
#' print(bic_comparison)
#' # Note: Column is still named 'AIC', but values represent BIC
#' # Optional: rename column for clarity
#' # colnames(bic_comparison)[colnames(bic_comparison) == "AIC"] <- "BIC"
#' # print(bic_comparison)
#' }
#'
#' @keywords models methods
#' @author Lopes, J. E.
#' @export
AIC.gkwfit <- function(object, ..., k = 2) {
  # --- Input Validation ---
  if (!inherits(object, "gkwfit")) {
    stop("Input 'object' must be of class 'gkwfit'.")
  }
  objects <- list(object, ...)
  obj_classes <- sapply(objects, class)
  if (any(obj_classes != "gkwfit")) {
    warning("All objects passed to AIC should ideally be of class 'gkwfit'.")
  }
  if (!is.numeric(k) || length(k) != 1 || k < 0) {
    stop("'k' must be a single non-negative numeric value.")
  }


  # --- Use logLik method to get value and df attribute consistently ---
  lls <- lapply(objects, function(obj) {
    ll <- tryCatch(logLik(obj), error = function(e) NULL)
    if (is.null(ll) || !inherits(ll, "logLik")) {
      stop("Could not extract valid 'logLik' object for one of the models.")
    }
    if (is.null(attr(ll, "df")) || is.null(attr(ll, "nobs"))) {
      stop("The 'logLik' object is missing 'df' or 'nobs' attribute.")
    }
    ll
  })

  # --- Check if nobs are consistent for model comparison ---
  if (length(lls) > 1) {
    nobs_vals <- sapply(lls, attr, "nobs")
    if (length(unique(nobs_vals)) > 1) {
      warning(
        "Models were not all fitted to the same number of observations.",
        "\nAIC/BIC comparison might be problematic."
      )
    }
  }

  # --- Calculate Criterion ---
  vals <- sapply(lls, function(ll) -2 * as.numeric(ll) + k * attr(ll, "df"))
  dfs <- sapply(lls, function(ll) attr(ll, "df"))

  # --- Format Output ---
  if (length(objects) == 1) {
    return(vals)
  } else {
    # Try to get model names from call
    calls <- lapply(objects, function(obj) {
      # Attempt to retrieve the call safely
      obj_call <- tryCatch(obj$call, error = function(e) NULL)
      if (!is.call(obj_call)) obj_call <- NULL # Ensure it's a call or NULL
      obj_call
    })
    # Use deparse1 for concise name, provide default if call missing
    mnames <- sapply(calls, function(cc) {
      if (!is.null(cc)) deparse1(cc, collapse = " ") else paste0("Model", seq_along(calls))
    })
    # Handle potential duplicate names
    if (anyDuplicated(mnames)) mnames <- make.unique(mnames, sep = ".")

    result <- data.frame(
      df = dfs,
      AIC = vals # Column name remains "AIC" regardless of k, per stats::AIC convention
    )
    rownames(result) <- mnames

    # Sort by the criterion value
    result <- result[order(result$AIC), ]
    return(result)
  }
}



#' @title Calculate Bayesian Information Criterion (BIC) for gkwfit Objects
#'
#' @description
#' Computes the Bayesian Information Criterion (BIC), sometimes called the
#' Schwarz criterion (SIC), for one or more fitted model objects of class \code{"gkwfit"}.
#'
#' @param object An object of class \code{"gkwfit"}, typically the result of a call to \code{\link{gkwfit}}.
#' @param ... Optionally, more fitted model objects of class \code{"gkwfit"}.
#'
#' @details
#' This function calculates the BIC based on the formula
#' \eqn{-2 \times \log Likelihood + \log(n) \times df}, where \eqn{n} is the number
#' of observations and \eqn{df} represents the number of estimated parameters in the
#' model (degrees of freedom).
#'
#' It relies on the \code{\link{logLik.gkwfit}} method to extract the log-likelihood,
#' the degrees of freedom (\code{df}), and the number of observations (\code{nobs})
#' for each model.
#'
#' When comparing multiple models fitted to the **same data**, the model with the
#' lower BIC value is generally preferred, as BIC tends to penalize model complexity
#' more heavily than AIC for larger sample sizes. The function returns a sorted
#' data frame to facilitate this comparison when multiple objects are provided.
#'
#' @return
#' \itemize{
#'   \item If only one \code{object} is provided: A single numeric value representing the calculated BIC.
#'   \item If multiple objects are provided: A \code{data.frame} with rows corresponding
#'     to the models and columns for the degrees of freedom (\code{df}) and the
#'     calculated BIC value (named \code{BIC}). The data frame is sorted in
#'     ascending order based on the BIC values. Row names are derived from the
#'     deparsed calls of the fitted models.
#' }
#'
#' @seealso \code{\link{gkwfit}}, \code{\link[stats]{BIC}}, \code{\link{logLik.gkwfit}}, \code{\link{AIC.gkwfit}}
#'
#' @examples
#' \dontrun{
#' # Assume necessary functions (rkw, rgkw, gkwfit, logLik.gkwfit) exist
#'
#' set.seed(123)
#' sample_data <- rkw(100, alpha = 2.5, beta = 1.5)
#'
#' # Fit different models to the same data
#' fit1_kw <- gkwfit(sample_data, family = "kw", silent = TRUE)
#' fit2_bkw <- gkwfit(sample_data, family = "bkw", silent = TRUE) # Assuming bkw fits
#' fit3_gkw <- gkwfit(sample_data, family = "gkw", silent = TRUE) # Assuming gkw fits
#'
#' # Calculate BIC for a single model
#' bic1 <- BIC(fit1_kw)
#' print(bic1)
#'
#' # Compare BIC values for multiple models
#' bic_comparison <- BIC(fit1_kw, fit2_bkw, fit3_gkw)
#' print(bic_comparison) # Sorted by BIC
#' }
#'
#' @keywords models methods
#' @author Lopes, J. E.
#' @export
BIC.gkwfit <- function(object, ...) {
  # --- Input Validation ---
  if (!inherits(object, "gkwfit")) {
    stop("Input 'object' must be of class 'gkwfit'.")
  }
  objects <- list(object, ...)
  obj_classes <- sapply(objects, class)
  if (any(!sapply(obj_classes, function(cls) inherits(cls, "gkwfit")))) {
    warning("All objects passed to BIC should ideally inherit from 'gkwfit'.")
  }

  # --- Use logLik method to get value, df, and nobs attribute consistently ---
  lls <- lapply(objects, function(obj) {
    ll <- tryCatch(logLik(obj), error = function(e) NULL)
    if (is.null(ll) || !inherits(ll, "logLik")) {
      stop("Could not extract valid 'logLik' object for model: ", deparse1(substitute(obj)))
    }
    req_attrs <- c("df", "nobs")
    if (!all(req_attrs %in% names(attributes(ll)))) {
      stop("The 'logLik' object is missing 'df' or 'nobs' attribute for model: ", deparse1(substitute(obj)))
    }
    # Check nobs > 0 for log()
    n_obs_val <- attr(ll, "nobs")
    if (is.null(n_obs_val) || !is.numeric(n_obs_val) || length(n_obs_val) != 1 || n_obs_val <= 0) {
      stop("Number of observations ('nobs') must be a single positive value for BIC calculation for model: ", deparse1(substitute(obj)))
    }
    ll
  })

  # --- Check if nobs are consistent for model comparison ---
  nobs_vals <- sapply(lls, attr, "nobs")
  if (length(lls) > 1) {
    if (length(unique(nobs_vals)) > 1) {
      warning(
        "Models were not all fitted to the same number of observations.",
        "\nBIC comparison might be problematic."
      )
    }
  }

  # --- Calculate BIC ---
  # BIC = -2*logLik + log(nobs) * df
  vals <- sapply(lls, function(ll) -2 * as.numeric(ll) + log(attr(ll, "nobs")) * attr(ll, "df"))
  dfs <- sapply(lls, function(ll) attr(ll, "df"))

  # --- Format Output ---
  if (length(objects) == 1) {
    return(vals)
  } else {
    # Try to get model names from call
    calls <- lapply(objects, function(obj) {
      obj_call <- tryCatch(obj$call, error = function(e) NULL)
      if (!is.call(obj_call)) obj_call <- NULL
      obj_call
    })
    # Use deparse1 for concise name, provide default if call missing
    mnames <- character(length(calls))
    for (i in seq_along(calls)) {
      if (!is.null(calls[[i]])) {
        mnames[i] <- deparse1(calls[[i]], collapse = " ")
      } else {
        # Try getting name from the list(...) call if possible
        deparsed_arg <- deparse1(substitute(list(object, ...))[[i + 1]])
        mnames[i] <- if (deparsed_arg == ".") paste0("Model", i) else deparsed_arg
      }
    }

    if (anyDuplicated(mnames)) mnames <- make.unique(mnames, sep = ".")

    result <- data.frame(
      df = dfs,
      BIC = vals # Column name is BIC
    )
    rownames(result) <- mnames

    # Sort by the criterion value
    result <- result[order(result$BIC), ]
    return(result)
  }
}


#' @title Compare Fitted gkwfit Models using Likelihood Ratio Tests
#'
#' @description
#' Computes Likelihood Ratio Tests (LRT) to compare two or more nested models
#' fitted using \code{\link{gkwfit}}. It produces a table summarizing the models
#' and the test statistics.
#'
#' @param object An object of class \code{"gkwfit"}, representing the first fitted model.
#' @param ... One or more additional objects of class \code{"gkwfit"}, representing
#'   subsequent fitted models, assumed to be nested within each other or the first model.
#'
#' @details
#' This function performs pairwise likelihood ratio tests between consecutively ordered
#' models (ordered by their degrees of freedom). It assumes the models are nested
#' and are fitted to the same dataset. A warning is issued if the number of
#' observations differs between models.
#'
#' The Likelihood Ratio statistic is calculated as \eqn{LR = 2 \times (\log L_{complex} - \log L_{simple})}.
#' This statistic is compared to a Chi-squared distribution with degrees of freedom
#' equal to the difference in the number of parameters between the two models
#' (\eqn{\Delta df = df_{complex} - df_{simple}}).
#'
#' The output table includes the number of parameters (`N.Par`), AIC, BIC, log-likelihood (`LogLik`),
#' the test description (`Test`), the LR statistic (`LR stat.`), and the p-value (`Pr(>Chi)`).
#' Models are ordered by increasing complexity (number of parameters).
#'
#' Warnings are issued if models do not appear correctly nested based on degrees of
#' freedom or if the log-likelihood decreases for a more complex model, as the LRT
#' results may not be meaningful in such cases.
#'
#' The function relies on a working \code{\link{logLik.gkwfit}} method to extract
#' necessary information (log-likelihood, df, nobs).
#'
#' @return An object of class \code{c("anova.gkwfit", "anova", "data.frame")}.
#'   This data frame contains rows for each model and columns summarizing the fit
#'   and the pairwise likelihood ratio tests. It includes:
#'   \item{N.Par}{Number of estimated parameters (degrees of freedom).}
#'   \item{AIC}{Akaike Information Criterion.}
#'   \item{BIC}{Bayesian Information Criterion.}
#'   \item{LogLik}{Log-likelihood value.}
#'   \item{Test}{Description of the pairwise comparison (e.g., "1 vs 2").}
#'   \item{LR stat.}{Likelihood Ratio test statistic.}
#'   \item{Pr(>Chi)}{P-value from the Chi-squared test.}
#'   The table is printed using a method that mimics \code{print.anova}.
#'
#' @seealso \code{\link{gkwfit}}, \code{\link{logLik.gkwfit}}, \code{\link{AIC.gkwfit}}, \code{\link{BIC.gkwfit}}, \code{\link[stats]{anova}}
#'
#' @examples
#' \dontrun{
#' # Assume necessary functions (rkw, rbkw, gkwfit, logLik.gkwfit) exist
#'
#' set.seed(123)
#' sample_data <- rkw(150, alpha = 2.5, beta = 1.5) # Data from Kw
#'
#' # Fit nested models (e.g., Kw is nested within BKw, which is nested within GKw)
#' fit_kw <- gkwfit(sample_data, family = "kw", silent = TRUE) # 2 params
#' fit_bkw <- gkwfit(sample_data, family = "bkw", silent = TRUE) # 4 params
#' fit_gkw <- gkwfit(sample_data, family = "gkw", silent = TRUE) # 5 params
#'
#' # Compare the nested models using LRT
#' lrt_results <- anova(fit_kw, fit_bkw, fit_gkw)
#' print(lrt_results)
#'
#' # Example comparing only two models
#' anova(fit_kw, fit_gkw)
#' }
#' @keywords models methods regression htest
#' @author Lopes, J. E.
#' @export anova.gkwfit
anova.gkwfit <- function(object, ...) {
  # --- Gather objects and perform initial validation ---
  objects <- list(object, ...)
  is_gkwfit <- sapply(objects, inherits, "gkwfit")
  if (!all(is_gkwfit)) {
    stop("All objects provided must inherit from class 'gkwfit'.")
  }
  nmodels <- length(objects)
  if (nmodels < 2) {
    stop("Need at least two models to compare using anova().")
  }

  # Use substitute to capture original variable names if possible
  object_names_call <- match.call()
  object_names <- vapply(as.list(object_names_call[-1L])[seq_len(nmodels)], deparse1, "")

  lls <- vector("list", nmodels)
  for (i in seq_len(nmodels)) {
    lls[[i]] <- tryCatch(logLik(objects[[i]]), error = function(e) {
      stop("Could not extract valid 'logLik' object for model ", object_names[i], ": ", e$message)
    })
    if (!inherits(lls[[i]], "logLik") || is.null(attr(lls[[i]], "df")) || is.null(attr(lls[[i]], "nobs"))) {
      stop(
        "Invalid 'logLik' object returned for model ", object_names[i],
        " (missing class or 'df'/'nobs' attribute)."
      )
    }
    if (attr(lls[[i]], "nobs") <= 0) {
      stop("Number of observations ('nobs') must be positive for model ", object_names[i])
    }
  }

  dfs <- sapply(lls, attr, "df")
  nobs_vals <- sapply(lls, attr, "nobs")
  loglik_vals <- sapply(lls, as.numeric)

  if (length(unique(nobs_vals)) > 1) {
    warning(
      "Models were not all fitted to the same number of observations.\n",
      "Likelihood ratio tests assume comparison on the same dataset."
    )
  }
  n_obs <- nobs_vals[1] # Use first for AIC/BIC calculation, assumes consistency or accepts warning

  mnames <- object_names
  if (anyDuplicated(mnames)) mnames <- make.unique(mnames, sep = ".")

  ord <- order(dfs)
  objects <- objects[ord]
  lls <- lls[ord]
  dfs <- dfs[ord]
  loglik_vals <- loglik_vals[ord]
  mnames <- mnames[ord]

  aics <- -2 * loglik_vals + 2 * dfs
  bics <- -2 * loglik_vals + log(n_obs) * dfs

  # --- Perform pairwise LRTs ---
  lr_stat <- rep(NA_real_, nmodels)
  pr_chi <- rep(NA_real_, nmodels)
  delta_df <- rep(NA_integer_, nmodels)
  test_desc <- rep("", nmodels)

  for (i in 2:nmodels) {
    df_diff <- dfs[i] - dfs[i - 1]

    # Check if degrees of freedom increased (necessary for meaningful LRT)
    if (df_diff <= 0) {
      warning("Model ", mnames[i], " (", dfs[i], " df)",
        " does not have more parameters than model ", mnames[i - 1], " (", dfs[i - 1], " df).",
        "\nLRT requires models to be nested with increasing complexity.",
        call. = FALSE
      )
      next # Skip LRT calculation for this pair
    }

    # Check if log-likelihood increased (as expected for nested models)
    loglik_diff <- loglik_vals[i] - loglik_vals[i - 1]
    if (loglik_diff < -1e-6) { # Allow for small numerical errors
      warning("Log-likelihood decreased unexpectedly for the more complex model (",
        mnames[i], " vs ", mnames[i - 1], ").\nCheck model convergence or nesting.",
        call. = FALSE
      )
      lr <- NA_real_ # LRT statistic is not meaningful
    } else {
      lr <- max(0, 2 * loglik_diff) # Ensure non-negative LR statistic
    }

    delta_df[i] <- df_diff
    lr_stat[i] <- lr
    # Calculate p-value only if LR is not NA
    pr_chi[i] <- if (!is.na(lr)) stats::pchisq(lr, df_diff, lower.tail = FALSE) else NA_real_
    test_desc[i] <- paste(i - 1, "vs", i)
  }

  result_table <- data.frame(
    `N.Par` = dfs,
    AIC = aics,
    BIC = bics,
    LogLik = loglik_vals,
    Test = test_desc,
    `LR stat.` = lr_stat,
    `Pr(>Chi)` = pr_chi,
    row.names = mnames,
    check.names = FALSE # Prevent conversion of Pr(>Chi) etc.
  )

  heading <- c("Likelihood Ratio Test Comparison\n")
  attr(result_table, "heading") <- heading
  class(result_table) <- c("anova.gkwfit", "anova", "data.frame")

  return(result_table)
}

#' @title S3 method for class 'anova.gkwfit'
#' @param x An object of class \code{"anova.gkwfit"}.
#' @param digits Minimum number of significant digits to print.
#' @param signif.stars Logical; if TRUE, add significance stars.
#' @param ... Other args passed to
#' @export print.anova.gkwfit
print.anova.gkwfit <- function(x, digits = max(getOption("digits") - 2L, 3L),
                               signif.stars = getOption("show.signif.stars", TRUE), ...) {
  if (!inherits(x, "anova")) stop("x not of class anova") # Basic check

  # Print header if it exists
  if (!is.null(heading <- attr(x, "heading"))) {
    cat(heading, "\n")
  }

  # Use printCoefmat for standard ANOVA table formatting
  # Identify the P-value column index
  pval_col_idx <- which(colnames(x) == "Pr(>Chi)")
  # Identify the test statistic column index
  test_col_idx <- which(colnames(x) == "LR stat.")
  if (length(pval_col_idx) == 0) pval_col_idx <- NULL # Handle if column name changes
  if (length(test_col_idx) == 0) test_col_idx <- NULL

  # stats:::printCoefmat is not exported, use stats::printCoefmat directly
  # We need to tell it which column contains p-values
  stats::printCoefmat(x,
    digits = digits, signif.stars = signif.stars,
    has.Pvalue = !is.null(pval_col_idx), # Does it have P-values?
    P.values = !is.null(pval_col_idx), # Are they P-values? (for formatting)
    cs.ind = integer(), # No specific estimate/SE cols
    tst.ind = test_col_idx, # Which column has test stats?
    zap.ind = integer(),
    na.print = "", ...
  ) # Don't print NA for test on first line
  cat("\n")
  invisible(x)
}
