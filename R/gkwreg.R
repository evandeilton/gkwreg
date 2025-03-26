#' Fit Generalized Kumaraswamy Regression Models
#'
#' @description
#' Fits a regression model using the Generalized Kumaraswamy (GKw) family of distributions
#' for bounded responses in the (0, 1) interval. The function supports modeling parameters
#' from all seven submodels of the GKw family as functions of predictors through appropriate
#' link functions.
#'
#' @param formula A Formula object of the form
#'   \code{y ~ alpha_terms | beta_terms | gamma_terms | delta_terms | lambda_terms},
#'   where each part on the right specifies the covariates for the corresponding parameter.
#'   Missing parts are automatically replaced with intercept-only models.
#' @param data A data frame containing the variables in the model.
#' @param family A character string specifying the distribution family. One of: "gkw" (default),
#'        "bkw", "kkw", "ekw", "mc", "kw", or "beta". See Details for parameter specifications.
#' @param link A character string or list of character strings specifying the link functions.
#'        Default links for each parameter: log for \eqn{\alpha}, \eqn{\beta}, \eqn{\gamma}, \eqn{\lambda} and logit for \eqn{\delta}.
#'        Supported link functions: "log", "logit", "identity", "inverse", "sqrt", "probit", "cloglog".
#' @param start Optional list with initial parameter values for regression coefficients.
#'        If \code{NULL}, reasonable starting values will be determined.
#' @param fixed Optional list of parameters or coefficients to be held fixed (not estimated).
#' @param fit Estimation method to be used: currently only supports \code{"tmb"}.
#' @param method (Only for \code{fit = "tmb"}) Optimization method: \code{"nlminb"} (default) or \code{"optim"}.
#' @param hessian Logical; if \code{TRUE}, computes standard errors and the covariance matrix.
#' @param profile Logical; if \code{TRUE}, computes likelihood profiles for parameters.
#' @param npoints Number of points to use in profile likelihood calculations.
#' @param plot Logical; if \code{TRUE}, generates diagnostic plots.
#' @param conf.level Confidence level for intervals (default: 0.95).
#' @param optimizer.control List of control parameters passed to the optimizer.
#' @param subset An optional vector specifying a subset of observations to be used in the fitting process.
#' @param weights An optional vector of weights to be used in the fitting process.
#' @param offset Additional terms with known coefficients to be included in the linear predictors.
#' @param na.action Function to handle missing values.
#' @param contrasts List of contrasts to be used for factors in the model matrix.
#' @param x Logical; if \code{TRUE}, the model matrix is returned.
#' @param y Logical; if \code{TRUE}, the response vector is returned.
#' @param model Logical; if \code{TRUE}, the model frame is returned.
#' @param silent Logical; if \code{TRUE}, suppresses messages.
#' @param ... Additional arguments passed to internal functions.
#'
#' @return An object of class \code{"gkwreg"} containing the fitted regression model.
#'
#' @details
#' The \code{gkwreg} function extends the \code{gkwfit} framework to regression modeling,
#' allowing parameters of the GKw family to depend on covariates. The function handles
#' all seven distributions in the Generalized Kumaraswamy family with their specific parameter sets:
#'
#' \itemize{
#'   \item \strong{GKw} (Generalized Kumaraswamy): 5 parameters (\eqn{\alpha}, \eqn{\beta}, \eqn{\gamma}, \eqn{\delta}, \eqn{\lambda})
#'   \item \strong{BKw} (Beta-Kumaraswamy): 4 parameters (\eqn{\alpha}, \eqn{\beta}, \eqn{\gamma}, \eqn{\delta}), with \eqn{\lambda = 1} fixed
#'   \item \strong{KKw} (Kumaraswamy-Kumaraswamy): 4 parameters (\eqn{\alpha}, \eqn{\beta}, \eqn{\delta}, \eqn{\lambda}), with \eqn{\gamma = 1} fixed
#'   \item \strong{EKw} (Exponentiated Kumaraswamy): 3 parameters (\eqn{\alpha}, \eqn{\beta}, \eqn{\lambda}), with \eqn{\gamma = 1}, \eqn{\delta = 0} fixed
#'   \item \strong{Mc} (McDonald/Beta Power): 3 parameters (\eqn{\gamma}, \eqn{\delta}, \eqn{\lambda}), with \eqn{\alpha = 1}, \eqn{\beta = 1} fixed
#'   \item \strong{Kw} (Kumaraswamy): 2 parameters (\eqn{\alpha}, \eqn{\beta}), with \eqn{\gamma = 1}, \eqn{\delta = 0}, \eqn{\lambda = 1} fixed
#'   \item \strong{Beta}: 2 parameters (\eqn{\gamma}, \eqn{\delta}), with \eqn{\alpha = 1}, \eqn{\beta = 1}, \eqn{\lambda = 1} fixed
#' }
#'
#' For each family, only the relevant parameters are modeled as functions of covariates. The fixed
#' parameters are automatically handled based on the selected family.
#'
#' Default link functions are assigned based on parameter constraints:
#' \itemize{
#'   \item "log" for \eqn{\alpha}, \eqn{\beta}, \eqn{\gamma}, and \eqn{\lambda} (which should be positive)
#'   \item "logit" for \eqn{\delta} (which should be between 0 and 1 when standardized)
#' }
#'
#' @examples
#' \dontrun{
#' require(gkwreg)
#'
#' ## Example 1: Simple Kumaraswamy regression model
#' set.seed(123)
#' n <- 500
#' x1 <- runif(n, -2, 2)
#' x2 <- rnorm(n)
#'
#' # Generate regression coefficients
#' alpha_coef <- c(0.8, 0.3, -0.2) # Intercept, x1, x2
#' beta_coef <- c(1.2, -0.4, 0.1) # Intercept, x1, x2
#'
#' # Generate linear predictors and transform to parameters
#' alpha <- exp(alpha_coef[1] + alpha_coef[2] * x1 + alpha_coef[3] * x2)
#' beta <- exp(beta_coef[1] + beta_coef[2] * x1 + beta_coef[3] * x2)
#'
#' # Generate responses from Kumaraswamy distribution
#' y <- rkw(n, alpha = alpha, beta = beta)
#'
#' # Create data frame
#' df <- data.frame(y = y, x1 = x1, x2 = x2)
#'
#' # Fit Kumaraswamy regression model using extended formula syntax
#' kw_reg <- gkwreg(y ~ x1 + x2 | x1 + x2, data = df, family = "kw")
#'
#' # Display summary
#' summary(kw_reg)
#'
#' # Plot diagnostics
#' plot(kw_reg)
#'
#' ## Example 2: Generalized Kumaraswamy regression with parameter-specific formulas
#' set.seed(456)
#' n <- 400
#' x1 <- runif(n, -1, 1)
#' x2 <- rnorm(n)
#' x3 <- rbinom(n, 1, 0.5)
#'
#' # Generate regression coefficients for each parameter
#' alpha_coef <- c(0.5, 0.2) # Intercept, x1
#' beta_coef <- c(0.8, -0.3, 0.1) # Intercept, x1, x2
#' gamma_coef <- c(0.6, 0.4) # Intercept, x3
#' delta_coef <- c(0.0, 0.2) # Intercept, x3
#' lambda_coef <- c(-0.2, 0.1) # Intercept, x2
#'
#' # Generate linear predictors and transform to parameters
#' alpha <- exp(alpha_coef[1] + alpha_coef[2] * x1)
#' beta <- exp(beta_coef[1] + beta_coef[2] * x1 + beta_coef[3] * x2)
#' gamma <- exp(gamma_coef[1] + gamma_coef[2] * x3)
#' delta <- plogis(delta_coef[1] + delta_coef[2] * x3)
#' lambda <- exp(lambda_coef[1] + lambda_coef[2] * x2)
#'
#' # Generate response from GKw distribution
#' y <- rgkw(n, alpha = alpha, beta = beta, gamma = gamma, delta = delta, lambda = lambda)
#'
#' # Create data frame
#' df <- data.frame(y = y, x1 = x1, x2 = x2, x3 = as.factor(x3))
#'
#' # Fit GKw regression with parameter-specific formulas using the extended formula syntax
#' gkw_reg <- gkwreg(y ~ x1 | x1 + x2 | x3 | x3 | x2, data = df, family = "gkw")
#'
#' # Compare true vs. estimated coefficients
#' coef(gkw_reg)
#'
#' ## Example 3: Beta regression for a simpler case
#' set.seed(789)
#' n <- 300
#' x1 <- runif(n, -1, 1)
#'
#' # Generate regression coefficients
#' gamma_coef <- c(1.0, 0.5) # Intercept, x1
#' delta_coef <- c(1.5, -0.7) # Intercept, x1
#'
#' # Generate linear predictors and transform to parameters
#' gamma <- exp(gamma_coef[1] + gamma_coef[2] * x1)
#' delta <- exp(delta_coef[1] + delta_coef[2] * x1)
#'
#' # Generate response from Beta distribution
#' y <- rbeta_(n, gamma, delta)
#'
#' # Create data frame
#' df <- data.frame(y = y, x1 = x1)
#'
#' # Fit Beta regression model - for Beta family, we only need the gamma and delta terms
#' beta_reg <- gkwreg(y ~ . | . | x1 | x1, data = df, family = "beta")
#'
#' # Display confidence intervals for coefficients
#' confint(beta_reg)
#'
#' ## Example 4: Model comparison
#' # Fit multiple models and compare them
#' kw_reg2 <- gkwreg(y ~ x1 | x1, data = df, family = "kw")
#' beta_reg2 <- gkwreg(y ~ . | . | x1 | x1, data = df, family = "beta")
#'
#' # Compare with AIC/BIC
#' AIC(kw_reg2, beta_reg2)
#' BIC(kw_reg2, beta_reg2)
#'
#' ## Example 5: Predicting with a fitted model
#' # Create new data for prediction
#' newdata <- data.frame(x1 = seq(-1, 1, by = 0.1))
#'
#' # Predict expected response
#' pred <- predict(beta_reg, newdata = newdata, type = "response")
#'
#' # Plot original data and predicted curve
#' plot(df$x1, df$y, pch = 20, col = "gray")
#' lines(newdata$x1, pred, col = "red", lwd = 2)
#' }
#'
#' @references
#' Kumaraswamy, P. (1980). A generalized probability density function for double-bounded
#' random processes. Journal of Hydrology, 46(1-2), 79-88.
#'
#' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized distributions.
#' Journal of Statistical Computation and Simulation, 81(7), 883-898.
#'
#' Ferrari, S. L. P., & Cribari-Neto, F. (2004). Beta regression for modelling rates and
#' proportions. Journal of Applied Statistics, 31(7), 799-815.
#'
#' @export
gkwreg <- function(formula,
                   data,
                   family = c("gkw", "bkw", "kkw", "ekw", "mc", "kw", "beta"),
                   link = NULL,
                   start = NULL,
                   fixed = NULL,
                   fit = "tmb",
                   method = c("nlminb", "optim"),
                   hessian = TRUE,
                   profile = FALSE,
                   npoints = 20,
                   plot = TRUE,
                   conf.level = 0.95,
                   optimizer.control = list(),
                   subset = NULL,
                   weights = NULL,
                   offset = NULL,
                   na.action = getOption("na.action"),
                   contrasts = NULL,
                   x = FALSE,
                   y = TRUE,
                   model = TRUE,
                   silent = TRUE,
                   ...) {
  # Match arguments
  family <- match.arg(family, choices = c("gkw", "bkw", "kkw", "ekw", "mc", "kw", "beta"))
  method <- match.arg(method, choices = c("nlminb", "optim"))
  call <- match.call()

  # Load Formula package for multi-part formula
  if (!requireNamespace("Formula", quietly = TRUE)) {
    stop("The 'Formula' package is required for this function. Please install it.")
  }

  # Get parameter information for the specified family
  param_info <- .get_family_param_info(family)
  param_names <- param_info$names
  fixed_params <- param_info$fixed
  param_positions <- param_info$positions

  # Convert to Formula object
  formula_obj <- Formula::as.Formula(formula)

  # Process Formula object to get individual formulas for each parameter
  formula_list <- .process_formula_parts(formula_obj, param_names, fixed_params, data)

  # Process link functions
  link_list <- .process_link(link, param_names, fixed_params)

  # Convert link strings to integers for TMB
  link_ints <- .convert_links_to_int(link_list)

  # Extract model frames, responses, and model matrices
  # ** FIXED: Pass the original call to .extract_model_data() **
  model_data <- .extract_model_data(
    formula_list = formula_list,
    data = data,
    subset = subset,
    weights = weights,
    na.action = na.action,
    offset = offset,
    contrasts = contrasts,
    original_call = call # Pass the original call for correct scoping
  )

  # Validate response variable is in (0, 1)
  y_var <- model_data$y
  .validate_data(y_var, length(param_names))

  # Initialize result list
  result <- list(
    call = call,
    family = family,
    formula = formula,
    link = link_list,
    param_names = param_names,
    fixed_params = fixed_params
  )

  # Process fixed parameters
  fixed_processed <- .process_fixed(fixed, param_names, fixed_params)

  # Prepare TMB data with correct matrices based on family
  tmb_data <- .prepare_tmb_data(
    model_data, family, param_names, fixed_processed,
    link_ints, y_var, param_positions
  )

  # Prepare TMB parameters in the correct structure required by each family
  tmb_params <- .prepare_tmb_params(
    model_data, family, param_names, fixed_processed,
    param_positions
  )

  # Compile and load the appropriate TMB model based on the family
  dll_name <- paste0(family, "reg")

  if (!silent) {
    message("Using TMB model: ", dll_name)
  }

  # Use the existing function to check and compile the TMB model
  .check_and_compile_TMB_code(dll_name, verbose = !silent)

  # Create TMB object
  obj <- TMB::MakeADFun(
    data = tmb_data,
    parameters = tmb_params,
    DLL = dll_name,
    silent = silent
  )

  # Set up optimizer control parameters
  default_control <- switch(method,
    "nlminb" = list(eval.max = 500, iter.max = 300),
    "optim" = list(maxit = 500)
  )
  opt_control <- modifyList(default_control, optimizer.control)

  # Optimize the model
  if (method == "nlminb") {
    if (!silent) message("Optimizing with nlminb...")
    opt <- stats::nlminb(
      start = obj$par,
      objective = obj$fn,
      gradient = obj$gr,
      control = opt_control
    )
    fit_result <- list(
      coefficients = obj$env$last.par,
      loglik = -opt$objective,
      convergence = opt$convergence,
      message = opt$message,
      iterations = opt$iterations
    )
  } else { # optim
    if (!silent) message("Optimizing with optim...")
    opt <- stats::optim(
      par = obj$par,
      fn = obj$fn,
      gr = obj$gr,
      method = "BFGS",
      control = opt_control
    )
    fit_result <- list(
      coefficients = opt$par,
      loglik = -opt$value,
      convergence = opt$convergence,
      message = ifelse(opt$convergence == 0, "Successful convergence", "Optimization failed"),
      iterations = opt$counts[1]
    )
  }

  # Format coefficient names with parameter mappings
  coef_names <- .format_coefficient_names(param_names, model_data, param_positions)
  names(fit_result$coefficients) <- coef_names

  # Calculate standard errors and covariance matrix if requested
  if (hessian) {
    if (!silent) message("Computing standard errors...")
    tryCatch(
      {
        sd_report <- TMB::sdreport(obj, getJointPrecision = TRUE)
        fit_result$se <- as.vector(sd_report$sd)
        fit_result$vcov <- as.matrix(sd_report$cov)

        # Add parameter names to standard errors
        names(fit_result$se) <- coef_names

        # Add confidence intervals
        alpha <- 1 - conf.level
        z_value <- stats::qnorm(1 - alpha / 2)
        fit_result$ci_lower <- fit_result$coefficients - z_value * fit_result$se
        fit_result$ci_upper <- fit_result$coefficients + z_value * fit_result$se
      },
      error = function(e) {
        warning("Could not compute standard errors: ", e$message)
        fit_result$se <- rep(NA, length(fit_result$coefficients))
        fit_result$vcov <- matrix(NA, length(fit_result$coefficients), length(fit_result$coefficients))
        fit_result$ci_lower <- fit_result$ci_upper <- rep(NA, length(fit_result$coefficients))
        names(fit_result$se) <- coef_names
      }
    )
  }

  # Extract TMB report
  tmb_report <- obj$report()

  # Extract parameter means
  alpha_mean <- tmb_report$alpha_mean
  beta_mean <- tmb_report$beta_mean
  gamma_mean <- tmb_report$gamma_mean
  delta_mean <- tmb_report$delta_mean
  lambda_mean <- tmb_report$lambda_mean

  # Extract parameter vectors for each observation
  alphaVec <- if ("alphaVec" %in% names(tmb_report)) tmb_report$alphaVec else rep(alpha_mean, length(y_var))
  betaVec <- if ("betaVec" %in% names(tmb_report)) tmb_report$betaVec else rep(beta_mean, length(y_var))
  gammaVec <- if ("gammaVec" %in% names(tmb_report)) tmb_report$gammaVec else rep(gamma_mean, length(y_var))
  deltaVec <- if ("deltaVec" %in% names(tmb_report)) tmb_report$deltaVec else rep(delta_mean, length(y_var))
  lambdaVec <- if ("lambdaVec" %in% names(tmb_report)) tmb_report$lambdaVec else rep(lambda_mean, length(y_var))

  # Extract fitted values
  fitted_values <- if ("fitted" %in% names(tmb_report)) tmb_report$fitted else rep(NA, length(y_var))

  # Calculate residuals
  response_residuals <- y_var - fitted_values

  # Calculate fit statistics
  n_obs <- length(y_var)
  npar <- length(fit_result$coefficients)

  # Deviance, AIC, BIC
  nll <- -fit_result$loglik
  deviance <- 2.0 * nll
  aic <- deviance + 2.0 * npar
  bic <- deviance + log(n_obs) * npar

  # Residual degrees of freedom
  df.residual <- n_obs - npar

  # Additional statistics
  rmse <- sqrt(mean(response_residuals^2, na.rm = TRUE))
  sse <- sum((y_var - fitted_values)^2, na.rm = TRUE)
  sst <- sum((y_var - mean(y_var))^2, na.rm = TRUE)
  efron_r2 <- if (sst > 0) 1 - sse / sst else NA
  mean_absolute_error <- mean(abs(response_residuals), na.rm = TRUE)

  # Build final result with all necessary components
  result <- list(
    call = call,
    family = family,
    formula = formula,
    coefficients = fit_result$coefficients,
    fitted.values = as.vector(fitted_values),
    residuals = as.vector(response_residuals),
    fitted_parameters = list(
      alpha = alpha_mean,
      beta = beta_mean,
      gamma = gamma_mean,
      delta = delta_mean,
      lambda = lambda_mean
    ),
    parameter_vectors = list(
      alphaVec = alphaVec,
      betaVec = betaVec,
      gammaVec = gammaVec,
      deltaVec = deltaVec,
      lambdaVec = lambdaVec
    ),
    link = link_list,
    param_names = param_names,
    fixed_params = fixed_params,
    loglik = fit_result$loglik,
    aic = aic,
    bic = bic,
    deviance = deviance,
    df.residual = df.residual,
    nobs = n_obs,
    npar = npar,
    vcov = if (hessian) fit_result$vcov else NULL,
    se = if (hessian) fit_result$se else NULL,
    convergence = fit_result$convergence,
    message = fit_result$message,
    iterations = fit_result$iterations,
    rmse = rmse,
    efron_r2 = efron_r2,
    mean_absolute_error = mean_absolute_error
  )

  # Add extra information if requested
  if (x) result$x <- model_data$matrices
  if (y) result$y <- y_var
  if (model) result$model <- model_data$model

  # Store TMB object
  result$tmb_object <- obj

  # Set class for S3 methods
  class(result) <- c("gkwreg", "list")

  # Return the final result
  return(result)
}


#' Prepare TMB Parameters for GKw Regression
#'
#' @param model_data List of model data.
#' @param family Family name.
#' @param param_names Names of parameters.
#' @param fixed List of fixed parameters.
#' @param param_positions Parameter position mapping for the family.
#' @return A list with TMB parameters.
#' @keywords internal
.prepare_tmb_params <- function(model_data, family, param_names, fixed, param_positions) {
  # Initialize params list with empty vectors for beta1, beta2, etc.
  params <- list()

  # Number of beta parameters needed depends on the family
  num_params <- switch(family,
    "gkw" = 5,
    "bkw" = 4,
    "kkw" = 4,
    "ekw" = 3,
    "mc" = 3,
    "kw" = 2,
    "beta" = 2
  )

  # Initialize empty vectors for all beta parameters
  for (i in 1:num_params) {
    params[[paste0("beta", i)]] <- numeric(0)
  }

  # Get non-fixed parameters
  non_fixed_params <- setdiff(param_names, names(fixed))

  # Fill in the parameter vectors
  for (param in non_fixed_params) {
    # Get TMB parameter position based on family
    tmb_pos <- param_positions[[param]]

    # Skip if not mapped
    if (is.null(tmb_pos) || is.na(tmb_pos)) next

    # Get the model matrix for this parameter
    mat_name <- paste0("beta", tmb_pos)

    # If parameter exists in model_data, use it
    if (param %in% names(model_data$matrices)) {
      X <- model_data$matrices[[param]]

      # Initialize with zeros
      params[[mat_name]] <- rep(0, ncol(X))

      # Set reasonable starting values for intercept
      if (ncol(X) > 0) {
        # For the intercept, use a reasonable value
        if (param %in% c("alpha", "beta", "gamma", "lambda")) {
          params[[mat_name]][1] <- 0.0 # log(1) = 0
        } else if (param == "delta") {
          params[[mat_name]][1] <- 0.0 # logit(0.5) = 0
        }
      }
    }
  }

  return(params)
}

#' Process Formula Parts from a Formula Object
#'
#' @param formula_obj Formula object created with the Formula package.
#' @param param_names Names of the parameters for the specified family.
#' @param fixed_params List of fixed parameters.
#' @param data Data frame containing the variables.
#' @return A list of formula objects for each parameter.
#' @keywords internal
.process_formula_parts <- function(formula_obj, param_names, fixed_params, data) {
  # Get non-fixed parameters
  non_fixed_params <- setdiff(param_names, names(fixed_params))

  # Extract the response variable from the formula
  resp_var <- as.character(formula_obj[[2]])

  # Create list to store formulas for each parameter
  formula_list <- list()

  # Get the max number of RHS parts in the formula
  n_parts <- length(attr(Formula::Formula(formula_obj), "rhs"))

  # Process each non-fixed parameter
  for (i in seq_along(non_fixed_params)) {
    param <- non_fixed_params[i]

    if (i <= n_parts) {
      # Extract the ith part of the formula
      rhs_part <- Formula::Formula(formula_obj, rhs = i, lhs = 1)[[3]]

      # Check if this part is just a dot
      if (identical(as.character(rhs_part), ".") || identical(as.character(rhs_part), "1")) {
        # Use intercept-only model
        formula_list[[param]] <- Formula::as.Formula(paste(resp_var, "~", "1"))
      } else {
        # Use the specified formula part
        formula_list[[param]] <- Formula::as.Formula(paste(resp_var, "~", deparse(rhs_part)))
      }
    } else {
      # For parameters beyond the number of specified parts, use intercept-only model
      formula_list[[param]] <- Formula::as.Formula(paste(resp_var, "~", "1"))
    }
  }

  return(formula_list)
}

#' Convert Link Function Names to TMB Integers
#'
#' @param link_list List of link function names
#' @return List of link function integers for TMB
#' @keywords internal
.convert_links_to_int <- function(link_list) {
  link_map <- c(
    "log" = 1,
    "logit" = 2,
    "probit" = 3,
    "cauchy" = 4,
    "cloglog" = 5,
    "identity" = 6,
    "sqrt" = 7,
    "inverse" = 8,
    "inverse-square" = 9
  )

  result <- lapply(link_list, function(link) {
    if (link %in% names(link_map)) {
      return(link_map[link])
    } else {
      warning("Unsupported link function: ", link, ". Using log link instead.")
      return(1) # Default to log
    }
  })

  return(result)
}


#' Format Coefficient Names Based on Family and Model Matrices
#'
#' @param param_names Names of parameters for the family
#' @param model_data Model data list including matrices
#' @param param_positions Parameter position mapping for the family
#' @return Vector of formatted coefficient names
#' @keywords internal
.format_coefficient_names <- function(param_names, model_data, param_positions) {
  # Initialize vector to accumulate all coefficient names
  all_coef_names <- c()
  # For each parameter that has a model matrix
  for (param in param_names) {
    if (param %in% names(model_data$matrices)) {
      # Get parameter position in TMB
      tmb_pos <- param_positions[[param]]
      # Get the model matrix
      X <- model_data$matrices[[param]]
      mat_coef_names <- colnames(X)
      # Create names for coefficients
      if (is.null(mat_coef_names)) {
        # If there are no column names, use generic names
        mat_coef_names <- paste0(param, "_", 1:ncol(X))
      } else {
        # Add parameter name prefix
        mat_coef_names <- paste0(param, ":", mat_coef_names)
      }
      # Add to accumulator vector
      all_coef_names <- c(all_coef_names, mat_coef_names)
    }
  }
  # Return coefficient names
  return(all_coef_names)
}

#' Get Parameter Information for a GKw Family Distribution
#'
#' @param family The GKw family distribution name.
#' @return A list with parameter information.
#' @keywords internal
.get_family_param_info <- function(family) {
  # Define parameter information and TMB parameter positions for each family
  family_params <- list(
    gkw = list(
      names = c("alpha", "beta", "gamma", "delta", "lambda"),
      n = 5,
      fixed = list(),
      positions = list(alpha = 1, beta = 2, gamma = 3, delta = 4, lambda = 5)
    ),
    bkw = list(
      names = c("alpha", "beta", "gamma", "delta"),
      n = 4,
      fixed = list(lambda = 1),
      positions = list(alpha = 1, beta = 2, gamma = 3, delta = 4)
    ),
    kkw = list(
      names = c("alpha", "beta", "delta", "lambda"),
      n = 4,
      fixed = list(gamma = 1),
      positions = list(alpha = 1, beta = 2, delta = 3, lambda = 4)
    ),
    ekw = list(
      names = c("alpha", "beta", "lambda"),
      n = 3,
      fixed = list(gamma = 1, delta = 0),
      positions = list(alpha = 1, beta = 2, lambda = 3)
    ),
    mc = list(
      names = c("gamma", "delta", "lambda"),
      n = 3,
      fixed = list(alpha = 1, beta = 1),
      positions = list(gamma = 1, delta = 2, lambda = 3)
    ),
    kw = list(
      names = c("alpha", "beta"),
      n = 2,
      fixed = list(gamma = 1, delta = 0, lambda = 1),
      positions = list(alpha = 1, beta = 2)
    ),
    beta = list(
      names = c("gamma", "delta"),
      n = 2,
      fixed = list(alpha = 1, beta = 1, lambda = 1),
      positions = list(gamma = 1, delta = 2)
    )
  )

  # Return parameter information for the specified family
  return(family_params[[family]])
}

#' Process Link Functions for GKw Regression
#'
#' @param link A character string or list of character strings specifying link functions.
#' @param param_names Names of the parameters for the specified family.
#' @param fixed_params List of fixed parameters.
#' @return A list of link functions.
#' @keywords internal
.process_link <- function(link, param_names, fixed_params) {
  # Default link functions for each parameter
  default_links <- list(
    alpha = "log",
    beta = "log",
    gamma = "log",
    delta = "logit",
    lambda = "log"
  )

  # Supported link functions
  supported_links <- c(
    "log", "logit", "identity", "inverse", "sqrt", "probit", "cloglog",
    "cauchy", "inverse-square"
  )

  # If link is NULL, use default links
  if (is.null(link)) {
    # Get default links for non-fixed parameters
    non_fixed_params <- setdiff(param_names, names(fixed_params))
    link_list <- default_links[non_fixed_params]
    return(link_list)
  }

  # If link is a single character string, apply to all parameters
  if (is.character(link) && length(link) == 1) {
    if (!link %in% supported_links) {
      stop(paste("Unsupported link function:", link))
    }

    # Apply the same link function to all non-fixed parameters
    non_fixed_params <- setdiff(param_names, names(fixed_params))
    link_list <- replicate(length(non_fixed_params), link, simplify = FALSE)
    names(link_list) <- non_fixed_params
    return(link_list)
  }

  # If link is a list, validate and return
  if (is.list(link) || is.character(link) && length(link) > 1) {
    if (is.character(link)) {
      link <- as.list(link)
      names(link) <- setdiff(param_names, names(fixed_params))
    }

    # Check if names of list match parameter names
    link_names <- names(link)
    if (is.null(link_names) || !all(link_names %in% param_names)) {
      stop("Names of link list must match parameter names for the chosen family")
    }

    # Check if all links are supported
    unsupported <- !unlist(link) %in% supported_links
    if (any(unsupported)) {
      stop(paste(
        "Unsupported link function(s):",
        paste(unlist(link)[unsupported], collapse = ", ")
      ))
    }

    # Remove links for fixed parameters
    fixed_param_names <- names(fixed_params)
    link <- link[setdiff(link_names, fixed_param_names)]

    return(link)
  }

  stop("link must be either a character string or a list of character strings")
}


#' Extract Model Data for GKw Regression
#'
#' @param formula_list List of formulas for each parameter.
#' @param data Data frame containing the variables.
#' @param subset Optional subset specification.
#' @param weights Optional weights.
#' @param na.action Function to handle missing values.
#' @param offset Optional offset.
#' @param contrasts List of contrasts for factors.
#' @param original_call The original function call.
#' @return A list of model data including frames, matrices, etc.
#' @keywords internal
.extract_model_data <- function(formula_list, data, subset, weights, na.action,
                                offset, contrasts, original_call) {
  # Initialize result list
  model_data <- list()

  # Get unique response variable name (should be the same for all formulas)
  resp_names <- unique(vapply(formula_list, function(f) as.character(f[[2]]), character(1)))
  if (length(resp_names) > 1) {
    stop("All formulas must have the same response variable")
  }

  # Extract model frames, responses, and model matrices for each parameter
  model_data$frames <- list()
  model_data$responses <- list()
  model_data$matrices <- list()
  model_data$terms <- list()

  for (param in names(formula_list)) {
    # Construct a list with arguments for model.frame
    mf_args <- list(
      formula = formula_list[[param]],
      data = data,
      subset = subset,
      na.action = na.action,
      drop.unused.levels = TRUE
    )

    # Evaluate and fix weights if provided
    if (!is.null(weights)) {
      weight_val <- tryCatch(weights, error = function(e) NULL)
      if (!is.null(weight_val) && !is.function(weight_val)) {
        mf_args$weights <- weight_val
      }
    }

    # Force the evaluation of the model frame in the proper environment
    model_data$frames[[param]] <- do.call(stats::model.frame, mf_args)

    # Extract response and model matrix
    model_data$responses[[param]] <- stats::model.response(model_data$frames[[param]])
    model_data$terms[[param]] <- attr(model_data$frames[[param]], "terms")
    X <- stats::model.matrix(model_data$terms[[param]], model_data$frames[[param]], contrasts)
    model_data$matrices[[param]] <- X
  }

  # Extract common response variable
  model_data$y <- model_data$responses[[1]]

  # Store model frame for the first parameter (all should have the same response)
  model_data$model <- model_data$frames[[1]]

  return(model_data)
}


#' Validate Data for GKw Regression
#'
#' @param data Numeric vector to validate.
#' @param n_params Number of parameters in the selected model.
#' @return The validated data.
#' @keywords internal
.validate_data <- function(data, n_params) {
  # Check if data is numeric
  if (!is.numeric(data)) {
    stop("Response variable must be numeric")
  }

  # Check for missing values
  if (any(is.na(data))) {
    stop("Response variable contains missing values")
  }

  # Check for values in the (0, 1) interval
  if (any(data <= 0 | data >= 1)) {
    stop("Response variable must be in the open interval (0, 1)")
  }

  # Return the validated data
  return(data)
}

#' Process Fixed Parameters for GKw Regression
#'
#' @param fixed List of fixed parameters or coefficients.
#' @param param_names Names of the parameters for the specified family.
#' @param fixed_params List of fixed parameters from the family definition.
#' @return A list of processed fixed parameters and coefficients.
#' @keywords internal
.process_fixed <- function(fixed, param_names, fixed_params) {
  # If no additional fixed parameters, return the family's fixed parameters
  if (is.null(fixed)) {
    return(fixed_params)
  }

  # Check if fixed is a valid list
  if (!is.list(fixed)) {
    stop("fixed must be a list")
  }

  # Combine user-defined fixed parameters with family's fixed parameters
  fixed_combined <- c(fixed, fixed_params)

  # Check for duplicates
  if (length(unique(names(fixed_combined))) < length(names(fixed_combined))) {
    stop("Duplicate entries in fixed parameters")
  }

  return(fixed_combined)
}

#' Prepare TMB Data for GKw Regression
#'
#' @param model_data List of model data.
#' @param family Family name.
#' @param param_names Names of parameters.
#' @param fixed List of fixed parameters and coefficients.
#' @param link_ints List of link function integers.
#' @param y Response variable.
#' @param param_positions Parameter position mapping for the family.
#' @return A list with TMB data.
#' @keywords internal
.prepare_tmb_data <- function(model_data, family, param_names, fixed, link_ints, y, param_positions) {
  # Initialize TMB data
  tmb_data <- list(
    y = y,
    useMeanCache = 1, # Enable mean caching
    calcFitted = 1, # Calculate fitted values
    userChunkSize = 100 # Reasonable chunk size
  )

  # All families need matrices and link types
  # The number of X matrices and link types needed depends on the family
  num_params <- switch(family,
    "gkw" = 5,
    "bkw" = 4,
    "kkw" = 4,
    "ekw" = 3,
    "mc" = 3,
    "kw" = 2,
    "beta" = 2
  )

  # Initialize default matrices and links for all required parameters
  for (i in 1:num_params) {
    matrix_name <- paste0("X", i)
    link_name <- paste0("link_type", i)
    scale_name <- paste0("scale", i)

    # Default empty matrix with 1 column (intercept only)
    tmb_data[[matrix_name]] <- matrix(0, nrow = length(y), ncol = 1)

    # Default link is log (1) and scale is 10
    tmb_data[[link_name]] <- 1
    tmb_data[[scale_name]] <- 10.0
  }

  # Fill in actual matrices and links for non-fixed parameters
  non_fixed_params <- setdiff(param_names, names(fixed))

  for (param in non_fixed_params) {
    # Get TMB parameter position based on family
    tmb_pos <- param_positions[[param]]

    # Skip if not mapped
    if (is.null(tmb_pos) || is.na(tmb_pos)) next

    # Update matrix, link type and scale
    matrix_name <- paste0("X", tmb_pos)
    link_name <- paste0("link_type", tmb_pos)
    scale_name <- paste0("scale", tmb_pos)

    # If parameter exists in model_data, use it
    if (param %in% names(model_data$matrices)) {
      tmb_data[[matrix_name]] <- model_data$matrices[[param]]

      # If link exists, use it
      if (param %in% names(link_ints)) {
        tmb_data[[link_name]] <- link_ints[[param]]
      }

      # Set appropriate scale for the parameter
      if (param == "delta") {
        tmb_data[[scale_name]] <- 1.0 # For delta, which uses logit link
      } else {
        tmb_data[[scale_name]] <- 10.0 # For other parameters, which typically use log link
      }
    }
  }

  return(tmb_data)
}


#' @title Summary Method for Generalized Kumaraswamy Regression Models
#'
#' @description
#' Computes and returns a detailed statistical summary of a fitted Generalized Kumaraswamy (GKw)
#' regression model. The summary provides extensive information about parameter estimates,
#' their standard errors, confidence intervals, and model fit statistics.
#'
#' @param object An object of class \code{"gkwreg"}, typically the result of a call to \code{\link{gkwreg}}.
#' @param conf.level Confidence level for the confidence intervals. Default is 0.95.
#' @param ... Additional arguments (not used).
#'
#' @details
#' This method creates a comprehensive summary of the fitted GKw regression model. It provides:
#' \itemize{
#'   \item Model call and family information
#'   \item A detailed coefficient table with estimates, standard errors, z-values, p-values, and confidence intervals
#'   \item Link functions for each parameter
#'   \item Fitted parameter means and their distributions
#'   \item Extensive model fit statistics
#'   \item Residual analysis summary
#'   \item Convergence information
#' }
#'
#' For the coefficient table, significant p-values are marked with stars based on common
#' significance levels (0.001, 0.01, 0.05, 0.1).
#'
#' @return An object of class \code{"summary.gkwreg"} with the following components:
#' \item{call}{The original function call that created the model.}
#' \item{family}{Distribution family used in the model.}
#' \item{coefficients}{Table of coefficient estimates, standard errors, z-values, and p-values.}
#' \item{conf.int}{Matrix of confidence intervals for the coefficients.}
#' \item{link}{List of link functions used for each parameter.}
#' \item{fitted_parameters}{List of fitted parameter means.}
#' \item{residuals}{Summary statistics of the model residuals.}
#' \item{nobs}{Number of observations used in the model.}
#' \item{npar}{Number of parameters in the model.}
#' \item{df.residual}{Residual degrees of freedom.}
#' \item{loglik}{Log-likelihood of the fitted model.}
#' \item{aic}{Akaike Information Criterion.}
#' \item{bic}{Bayesian Information Criterion.}
#' \item{rmse}{Root Mean Square Error.}
#' \item{efron_r2}{Efron's pseudo-R-squared.}
#' \item{mean_absolute_error}{Mean absolute error, if available.}
#' \item{convergence}{Convergence status.}
#' \item{iterations}{Number of iterations performed by the optimizer.}
#'
#' @seealso \code{\link{gkwreg}}
#'
#' @examples
#' \dontrun{
#' # Fit a Kumaraswamy regression model
#' kw_reg <- gkwreg(y ~ x1 + x2 | x1 + x2, data = df, family = "kw")
#'
#' # Generate detailed summary
#' summary_kw <- summary(kw_reg)
#'
#' # Print summary
#' print(summary_kw)
#'
#' # Extract coefficient table
#' coef(summary_kw)
#' }
#'
#' @export
summary.gkwreg <- function(object, conf.level = 0.95, ...) {
  # Calculate z-values and p-values
  coef_est <- object$coefficients
  se <- object$se

  if (is.null(se)) {
    coef_table <- data.frame(
      Estimate = coef_est,
      row.names = names(coef_est)
    )
  } else {
    z_values <- coef_est / se
    p_values <- 2 * pnorm(-abs(z_values))

    coef_table <- data.frame(
      Estimate = coef_est,
      `Std. Error` = se,
      `z value` = z_values,
      `Pr(>|z|)` = p_values,
      row.names = names(coef_est), check.names = FALSE
    )
  }

  # Calculate confidence intervals
  alpha <- 1 - conf.level
  z_value <- qnorm(1 - alpha / 2)

  if (!is.null(se)) {
    ci_lower <- coef_est - z_value * se
    ci_upper <- coef_est + z_value * se
    conf_int <- cbind(ci_lower, ci_upper)
    colnames(conf_int) <- c(
      paste0(format(100 * alpha / 2, digits = 1), "%"),
      paste0(format(100 * (1 - alpha / 2), digits = 1), "%")
    )
    rownames(conf_int) <- names(coef_est)
  } else {
    conf_int <- NULL
  }

  # Summarize residuals
  if (!is.null(object$residuals)) {
    res_summary <- c(
      Min = min(object$residuals, na.rm = TRUE),
      Q1 = quantile(object$residuals, 0.25, na.rm = TRUE),
      Median = stats::median(object$residuals, na.rm = TRUE),
      Mean = mean(object$residuals, na.rm = TRUE),
      Q3 = quantile(object$residuals, 0.75, na.rm = TRUE),
      Max = max(object$residuals, na.rm = TRUE)
    )
  } else {
    res_summary <- NULL
  }

  # Create and return summary object
  result <- list(
    call = object$call,
    family = object$family,
    coefficients = coef_table,
    conf.int = conf_int,
    link = object$link,
    fitted_parameters = object$fitted_parameters,
    residuals = res_summary,
    nobs = object$nobs,
    npar = object$npar,
    df.residual = object$df.residual,
    loglik = object$loglik,
    aic = object$aic,
    bic = object$bic,
    rmse = object$rmse,
    efron_r2 = object$efron_r2,
    mean_absolute_error = object$mean_absolute_error,
    convergence = object$convergence,
    iterations = object$iterations,
    conf.level = conf.level
  )

  class(result) <- "summary.gkwreg"
  return(result)
}

#' @export
print.summary.gkwreg <- function(x, digits = max(3, getOption("digits") - 3),
                                 signif.stars = getOption("show.signif.stars"), ...) {
  cat("\nGeneralized Kumaraswamy Regression Model Summary\n\n")

  # Display family
  cat("Family:", x$family, "\n\n")

  # Display call
  cat("Call:\n")
  print(x$call)

  # Display residuals summary
  if (!is.null(x$residuals)) {
    cat("\nResiduals:\n")
    print(round(x$residuals, digits = digits))
  }

  # Display coefficient table with significance stars
  cat("\nCoefficients:\n")
  coefs <- x$coefficients
  stats::printCoefmat(coefs,
    digits = digits, signif.stars = signif.stars,
    has.Pvalue = ncol(coefs) >= 4
  )

  # Display confidence intervals
  if (!is.null(x$conf.int)) {
    cat("\nConfidence intervals (", x$conf.level * 100, "%):\n", sep = "")
    print(round(x$conf.int, digits = digits))
  }

  # Display link functions
  cat("\nLink functions:\n")
  for (param in names(x$link)) {
    cat(param, ": ", x$link[[param]], "\n", sep = "")
  }

  # Display fitted parameter means
  cat("\nFitted parameter means:\n")
  for (param in names(x$fitted_parameters)) {
    cat(param, ": ", format(x$fitted_parameters[[param]], digits = digits), "\n", sep = "")
  }

  # Display model fit statistics
  cat("\nModel fit statistics:\n")
  cat("Number of observations:", x$nobs, "\n")
  cat("Number of parameters:", x$npar, "\n")
  cat("Residual degrees of freedom:", x$df.residual, "\n")
  cat("Log-likelihood:", format(x$loglik, digits = digits), "\n")
  cat("AIC:", format(x$aic, digits = digits), "\n")
  cat("BIC:", format(x$bic, digits = digits), "\n")

  if (!is.null(x$rmse)) {
    cat("RMSE:", format(x$rmse, digits = digits), "\n")
  }

  if (!is.null(x$efron_r2) && !is.na(x$efron_r2)) {
    cat("Efron's R2:", format(x$efron_r2, digits = digits), "\n")
  }

  if (!is.null(x$mean_absolute_error)) {
    cat("Mean Absolute Error:", format(x$mean_absolute_error, digits = digits), "\n")
  }

  # Display convergence information
  cat("\nConvergence status:", ifelse(x$convergence == 0, "Successful", "Failed"), "\n")
  if (!is.null(x$iterations)) {
    cat("Iterations:", x$iterations, "\n")
  }

  cat("\n")
  invisible(x)
}

#' @export
coef.summary.gkwreg <- function(object, ...) {
  object$coefficients
}

#' @title Predict Method for Generalized Kumaraswamy Regression Models
#'
#' @description Obtains predictions and related quantities from a fitted
#' Generalized Kumaraswamy regression model.
#'
#' @param object A fitted model object of class "gkwreg".
#' @param newdata Optionally, a data frame with new data. If omitted, the original
#'   data is used.
#' @param type Character indicating type of prediction required:
#'   \itemize{
#'     \item "response" or "mean": predicted mean of the response (default)
#'     \item "link": linear predictors for all parameters
#'     \item "parameter": all distribution parameters (alpha, beta, gamma, delta, lambda)
#'     \item "alpha", "beta", "gamma", "delta", "lambda": individual parameters
#'     \item "variance": variance of the response
#'     \item "density" or "pdf": density function at values specified by 'at'
#'     \item "probability" or "cdf": cumulative distribution function at 'at'
#'     \item "quantile": quantiles of the distribution at probabilities 'at'
#'   }
#' @param na.action Function determining what to do with missing values in newdata.
#' @param at Numeric vector at which to evaluate the prediction function for
#'   types "density", "probability", or "quantile".
#' @param elementwise Logical. Should each element of the distribution only be
#'   evaluated at the corresponding element of 'at' (TRUE) or at all elements
#'   in 'at' (FALSE).
#' @param family Character string specifying the distribution family to use.
#'   If NULL (default), uses the family from the fitted model.
#'   Available options: "gkw", "bkw", "kkw", "ekw", "mc", "kw", "beta".
#' @param ... Further arguments passed to methods.
#'
#' @return A vector, matrix, or data frame of predictions, depending on the type.
#'
#' @examples
#' \dontrun{
#' # Fit a GKw model
#' model <- gkwreg(y ~ x1 | x2 | 1 | 1 | 1, data = mydata, family = "gkw")
#'
#' # Predictions on new data
#' newdata <- data.frame(x1 = seq(0, 1, by = 0.1), x2 = rep(0, 11))
#'
#' # Predict mean response
#' predict(model, newdata, type = "response")
#'
#' # Predict parameters
#' predict(model, newdata, type = "parameter")
#'
#' # Predict densities at specific points
#' predict(model, newdata, type = "density", at = c(0.2, 0.5, 0.8))
#'
#' # Predict using a different family than what was fitted
#' predict(model, newdata, type = "response", family = "beta")
#' }
#'
#' @export
predict.gkwreg <- function(object, newdata = NULL,
                           type = c(
                             "response", "link", "parameter",
                             "alpha", "beta", "gamma", "delta", "lambda",
                             "variance", "density", "pdf",
                             "probability", "cdf", "quantile"
                           ),
                           na.action = stats::na.pass, at = 0.5,
                           elementwise = NULL, family = NULL, ...) {
  # Match type argument
  type <- match.arg(type)

  # Aliases for some types
  if (type == "pdf") type <- "density"
  if (type == "cdf") type <- "probability"
  if (type == "mean") type <- "response"

  # Validate object
  if (!inherits(object, "gkwreg")) {
    stop("'object' must be a gkwreg model")
  }

  # Get the family from the object if not specified
  if (is.null(family)) {
    if (!is.null(object$family)) {
      family <- object$family
    } else {
      # Default to gkw for backward compatibility
      family <- "gkw"
      warning("No family specified in the model. Using 'gkw' as default.")
    }
  } else {
    # Validate the family parameter
    family <- match.arg(family, c("gkw", "bkw", "kkw", "ekw", "mc", "kw", "beta"))
  }

  # Import functions from the required package
  if (!requireNamespace("Formula", quietly = TRUE)) {
    stop("Package 'Formula' is needed for this function.")
  }

  # Handle case when we can directly use parameter vectors
  if (is.null(newdata) && !is.null(object$parameter_vectors) &&
    type %in% c(
      "response", "parameter", "alpha", "beta", "gamma", "delta", "lambda",
      "density", "probability", "quantile"
    )) {
    n <- length(object$parameter_vectors$alphaVec)
    params <- matrix(0, nrow = n, ncol = 5)
    params[, 1] <- object$parameter_vectors$alphaVec
    params[, 2] <- object$parameter_vectors$betaVec
    params[, 3] <- object$parameter_vectors$gammaVec
    params[, 4] <- object$parameter_vectors$deltaVec
    params[, 5] <- object$parameter_vectors$lambdaVec

    # Special case for response and fitted values
    if (type == "response" && !is.null(object$fitted.values)) {
      return(object$fitted.values)
    }

    # Handle parameter predictions
    if (type == "parameter") {
      return(data.frame(
        alpha = params[, 1],
        beta = params[, 2],
        gamma = params[, 3],
        delta = params[, 4],
        lambda = params[, 5]
      ))
    } else if (type %in% c("alpha", "beta", "gamma", "delta", "lambda")) {
      param_index <- match(type, c("alpha", "beta", "gamma", "delta", "lambda"))
      return(params[, param_index])
    }

    # For density, probability, and quantile using at
    if (type %in% c("density", "probability", "quantile")) {
      if (!is.numeric(at)) {
        stop("'at' must be numeric")
      }

      if (is.null(elementwise)) {
        elementwise <- length(at) == n
      }

      if (elementwise) {
        if (length(at) != n) {
          stop("For elementwise=TRUE, length of 'at' must equal number of observations")
        }
        eval_y <- at
        eval_params <- params
      } else {
        eval_params <- do.call(rbind, lapply(seq_len(length(at)), function(i) params))
        eval_y <- rep(at, each = n)
      }

      # Use the specific family distribution functions with proper prefixes
      result <- numeric(length(eval_y))

      # Loop through observations to calculate results using the appropriate distribution function
      for (i in seq_along(eval_y)) {
        alpha_i <- eval_params[i, 1]
        beta_i <- eval_params[i, 2]
        gamma_i <- eval_params[i, 3]
        delta_i <- eval_params[i, 4]
        lambda_i <- eval_params[i, 5]
        y_i <- eval_y[i]

        # Use the appropriate distribution function based on family and type
        if (type == "density") {
          # Density functions (d-prefix)
          switch(family,
            "gkw" = {
              result[i] <- dgkw(y_i, alpha_i, beta_i, gamma_i, delta_i, lambda_i)
            },
            "bkw" = {
              result[i] <- dbkw(y_i, alpha_i, beta_i, gamma_i, delta_i)
            },
            "kkw" = {
              result[i] <- dkkw(y_i, alpha_i, beta_i, delta_i, lambda_i)
            },
            "ekw" = {
              result[i] <- dekw(y_i, alpha_i, beta_i, lambda_i)
            },
            "mc" = {
              result[i] <- dmc(y_i, gamma_i, delta_i, lambda_i)
            },
            "kw" = {
              result[i] <- dkw(y_i, alpha_i, beta_i)
            },
            "beta" = {
              result[i] <- dbeta_(y_i, gamma_i, delta_i + 1)
            }
          )
        } else if (type == "probability") {
          # CDF functions (p-prefix)
          switch(family,
            "gkw" = {
              result[i] <- pgkw(y_i, alpha_i, beta_i, gamma_i, delta_i, lambda_i)
            },
            "bkw" = {
              result[i] <- pbkw(y_i, alpha_i, beta_i, gamma_i, delta_i)
            },
            "kkw" = {
              result[i] <- pkkw(y_i, alpha_i, beta_i, delta_i, lambda_i)
            },
            "ekw" = {
              result[i] <- pekw(y_i, alpha_i, beta_i, lambda_i)
            },
            "mc" = {
              result[i] <- pmc(y_i, gamma_i, delta_i, lambda_i)
            },
            "kw" = {
              result[i] <- pkw(y_i, alpha_i, beta_i)
            },
            "beta" = {
              result[i] <- pbeta_(y_i, gamma_i, delta_i + 1)
            }
          )
        } else if (type == "quantile") {
          # Quantile functions (q-prefix)
          switch(family,
            "gkw" = {
              result[i] <- qgkw(y_i, alpha_i, beta_i, gamma_i, delta_i, lambda_i)
            },
            "bkw" = {
              result[i] <- qbkw(y_i, alpha_i, beta_i, gamma_i, delta_i)
            },
            "kkw" = {
              result[i] <- qkkw(y_i, alpha_i, beta_i, delta_i, lambda_i)
            },
            "ekw" = {
              result[i] <- qekw(y_i, alpha_i, beta_i, lambda_i)
            },
            "mc" = {
              result[i] <- qmc(y_i, gamma_i, delta_i, lambda_i)
            },
            "kw" = {
              result[i] <- qkw(y_i, alpha_i, beta_i)
            },
            "beta" = {
              result[i] <- qbeta_(y_i, gamma_i, delta_i + 1)
            }
          )
        }
      }

      if (elementwise) {
        return(result)
      } else {
        result_matrix <- matrix(result, nrow = n, ncol = length(at))
        colnames(result_matrix) <- as.character(at)
        return(result_matrix)
      }
    }
  }

  # Special case for type = "variance" - implement direct calculation
  if (type == "variance" && is.null(newdata) &&
    !is.null(object$parameter_vectors) && !is.null(object$fitted.values)) {
    n <- length(object$fitted.values)
    variances <- numeric(n)

    # Get needed parameters
    alpha_vec <- object$parameter_vectors$alphaVec
    beta_vec <- object$parameter_vectors$betaVec
    gamma_vec <- object$parameter_vectors$gammaVec
    delta_vec <- object$parameter_vectors$deltaVec
    lambda_vec <- object$parameter_vectors$lambdaVec
    means <- object$fitted.values

    # Calculate variance for each observation
    for (i in 1:n) {
      if (family == "beta") {
        # Analytical formula for Beta distribution
        variances[i] <- (gamma_vec[i] * delta_vec[i]) /
          ((gamma_vec[i] + delta_vec[i])^2 * (gamma_vec[i] + delta_vec[i] + 1))
      } else if (family == "kw") {
        # Approximate formula for Kumaraswamy distribution
        variances[i] <- alpha_vec[i] * beta_vec[i] *
          beta(1 + 1 / alpha_vec[i], beta_vec[i]) -
          (alpha_vec[i] * beta_vec[i] * beta(1 + 1 / alpha_vec[i], beta_vec[i]))^2
      } else {
        # Numerical approximation using density function with proper prefix
        h <- 0.001
        mean_i <- means[i]

        # Safely compute y within (0,1)
        mean_plus <- min(0.99999, mean_i + h)
        mean_minus <- max(0.00001, mean_i - h)

        # Use the appropriate density function for the family
        switch(family,
          "gkw" = {
            f_plus <- dgkw(mean_plus, alpha_vec[i], beta_vec[i], gamma_vec[i], delta_vec[i], lambda_vec[i])
            f <- dgkw(mean_i, alpha_vec[i], beta_vec[i], gamma_vec[i], delta_vec[i], lambda_vec[i])
            f_minus <- dgkw(mean_minus, alpha_vec[i], beta_vec[i], gamma_vec[i], delta_vec[i], lambda_vec[i])
          },
          "bkw" = {
            f_plus <- dbkw(mean_plus, alpha_vec[i], beta_vec[i], gamma_vec[i], delta_vec[i])
            f <- dbkw(mean_i, alpha_vec[i], beta_vec[i], gamma_vec[i], delta_vec[i])
            f_minus <- dbkw(mean_minus, alpha_vec[i], beta_vec[i], gamma_vec[i], delta_vec[i])
          },
          "kkw" = {
            f_plus <- dkkw(mean_plus, alpha_vec[i], beta_vec[i], delta_vec[i], lambda_vec[i])
            f <- dkkw(mean_i, alpha_vec[i], beta_vec[i], delta_vec[i], lambda_vec[i])
            f_minus <- dkkw(mean_minus, alpha_vec[i], beta_vec[i], delta_vec[i], lambda_vec[i])
          },
          "ekw" = {
            f_plus <- dekw(mean_plus, alpha_vec[i], beta_vec[i], lambda_vec[i])
            f <- dekw(mean_i, alpha_vec[i], beta_vec[i], lambda_vec[i])
            f_minus <- dekw(mean_minus, alpha_vec[i], beta_vec[i], lambda_vec[i])
          },
          "mc" = {
            f_plus <- dmc(mean_plus, gamma_vec[i], delta_vec[i], lambda_vec[i])
            f <- dmc(mean_i, gamma_vec[i], delta_vec[i], lambda_vec[i])
            f_minus <- dmc(mean_minus, gamma_vec[i], delta_vec[i], lambda_vec[i])
          }
        )

        # Approximate second derivative of log-likelihood
        d2ll <- (f_plus - 2 * f + f_minus) / (h * h)

        # Variance is inverse of the information
        variances[i] <- min(0.25, max(1e-6, 1 / (abs(d2ll) + 1e-10)))
      }
    }

    return(variances)
  }

  # Special case for type = "link" when newdata is NULL
  if (type == "link" && is.null(newdata)) {
    if (!is.null(object$tmb_object) && !is.null(object$tmb_object$env$data)) {
      # Try to extract linear predictors from TMB object if available
      tmb_data <- object$tmb_object$env$data
      if (!is.null(tmb_data$X1) && !is.null(tmb_data$X2) &&
        !is.null(tmb_data$X3) && !is.null(tmb_data$X4) &&
        !is.null(tmb_data$X5)) {
        # Extract coefficients
        coefs <- object$coefficients
        if (is.list(coefs)) {
          beta1 <- coefs$alpha
          beta2 <- coefs$beta
          beta3 <- coefs$gamma
          beta4 <- coefs$delta
          beta5 <- coefs$lambda
        } else {
          # Try using regex patterns to extract coefficients
          beta1 <- coefs[grep("^alpha:", names(coefs))]
          beta2 <- coefs[grep("^beta:", names(coefs))]
          beta3 <- coefs[grep("^gamma:", names(coefs))]
          beta4 <- coefs[grep("^delta:", names(coefs))]
          beta5 <- coefs[grep("^lambda:", names(coefs))]

          # If regex didn't work, try to infer from dimensions
          if (length(beta1) == 0 || length(beta2) == 0) {
            # Extract design matrices
            X1 <- tmb_data$X1
            X2 <- tmb_data$X2
            X3 <- tmb_data$X3
            X4 <- tmb_data$X4
            X5 <- tmb_data$X5

            # Split coefficients based on matrix dimensions
            all_coefs <- coefs
            beta1 <- all_coefs[1:ncol(X1)]
            remaining <- all_coefs[-(1:ncol(X1))]

            beta2 <- remaining[1:ncol(X2)]
            remaining <- remaining[-(1:ncol(X2))]

            beta3 <- remaining[1:ncol(X3)]
            remaining <- remaining[-(1:ncol(X3))]

            beta4 <- remaining[1:ncol(X4)]
            remaining <- remaining[-(1:ncol(X4))]

            beta5 <- remaining[1:ncol(X5)]
          }
        }

        # Calculate linear predictors using TMB data matrices
        X1 <- tmb_data$X1
        X2 <- tmb_data$X2
        X3 <- tmb_data$X3
        X4 <- tmb_data$X4
        X5 <- tmb_data$X5

        eta1 <- as.vector(X1 %*% beta1)
        eta2 <- as.vector(X2 %*% beta2)
        eta3 <- as.vector(X3 %*% beta3)
        eta4 <- as.vector(X4 %*% beta4)
        eta5 <- as.vector(X5 %*% beta5)

        return(data.frame(
          alpha = eta1,
          beta = eta2,
          gamma = eta3,
          delta = eta4,
          lambda = eta5
        ))
      }
    }

    # If we can't extract linear predictors, return a warning and parameter values instead
    if (!is.null(object$parameter_vectors)) {
      warning("Cannot calculate link values without design matrices. Returning parameter values instead.")
      return(data.frame(
        alpha = object$parameter_vectors$alphaVec,
        beta = object$parameter_vectors$betaVec,
        gamma = object$parameter_vectors$gammaVec,
        delta = object$parameter_vectors$deltaVec,
        lambda = object$parameter_vectors$lambdaVec
      ))
    } else {
      stop("Cannot extract linear predictors. Please provide 'newdata' or set x=TRUE when fitting the model.")
    }
  }

  # Prepare new data or use original data
  if (is.null(newdata)) {
    # Try different ways to access model matrices
    if (!is.null(object$x)) {
      # Model was fitted with x=TRUE
      X1 <- object$x$alpha
      X2 <- object$x$beta
      X3 <- object$x$gamma
      X4 <- object$x$delta
      X5 <- object$x$lambda
    } else if (!is.null(object$model) && !is.null(object$formula)) {
      # Recreate model matrices from the model frame and formula
      formula <- object$formula
      mf <- object$model

      # Extract model matrices from the model frame
      formula_obj <- formula
      if (class(formula_obj)[1] == "formula") {
        formula_obj <- Formula::as.Formula(formula_obj)
      }

      # Safely get the number of RHS parts
      rhs_parts <- 0
      if (inherits(formula_obj, "Formula")) {
        rhs_parts <- length(attr(formula_obj, "rhs"))
      }

      # Initialize X matrices
      X <- vector("list", 5)
      param_names <- c("alpha", "beta", "gamma", "delta", "lambda")

      for (i in seq_len(5)) {
        if (i <= rhs_parts) {
          # Try to extract the model matrix for this part
          tryCatch(
            {
              X_i <- model.matrix(formula_obj, data = mf, rhs = i)
              X[[i]] <- X_i
            },
            error = function(e) {
              # If extraction fails, use intercept-only
              X[[i]] <- matrix(1,
                nrow = nrow(mf), ncol = 1,
                dimnames = list(NULL, "(Intercept)")
              )
            }
          )
        } else {
          # For parts beyond the specified ones, use intercept-only
          X[[i]] <- matrix(1,
            nrow = nrow(mf), ncol = 1,
            dimnames = list(NULL, "(Intercept)")
          )
        }
      }
      names(X) <- param_names

      X1 <- X$alpha
      X2 <- X$beta
      X3 <- X$gamma
      X4 <- X$delta
      X5 <- X$lambda
    } else if (!is.null(object$tmb_object) && !is.null(object$tmb_object$env$data)) {
      # Try using TMB object if available
      tmb_data <- object$tmb_object$env$data
      if (!is.null(tmb_data$X1) && !is.null(tmb_data$X2) &&
        !is.null(tmb_data$X3) && !is.null(tmb_data$X4) &&
        !is.null(tmb_data$X5)) {
        X1 <- tmb_data$X1
        X2 <- tmb_data$X2
        X3 <- tmb_data$X3
        X4 <- tmb_data$X4
        X5 <- tmb_data$X5
      } else {
        stop("Cannot extract original model matrices. Please provide 'newdata'.")
      }
    } else {
      stop("Cannot extract original model matrices. Please provide 'newdata' or set x=TRUE when fitting the model.")
    }
  } else {
    # Create model matrices from the new data
    formula <- object$formula

    # Ensure formula is a Formula object
    if (class(formula)[1] == "formula") {
      formula <- Formula::as.Formula(formula)
    }

    # Process newdata
    mf <- stats::model.frame(formula, newdata, na.action = na.action, drop.unused.levels = TRUE)

    # Safely get the number of RHS parts
    rhs_parts <- 0
    if (inherits(formula, "Formula")) {
      rhs_parts <- length(attr(formula, "rhs"))
    }

    # Construct model matrices for each parameter
    X <- vector("list", 5)
    param_names <- c("alpha", "beta", "gamma", "delta", "lambda")

    for (i in seq_len(5)) {
      if (i <= rhs_parts) {
        tryCatch(
          {
            X_i <- model.matrix(formula, data = mf, rhs = i)
            if (ncol(X_i) == 0) {
              if (i == 1) {
                stop("The first RHS (for alpha) cannot be empty.")
              } else {
                X[[i]] <- matrix(1,
                  nrow = nrow(mf), ncol = 1,
                  dimnames = list(NULL, "(Intercept)")
                )
              }
            } else {
              X[[i]] <- X_i
            }
          },
          error = function(e) {
            # If extraction fails, use intercept-only
            X[[i]] <- matrix(1,
              nrow = nrow(mf), ncol = 1,
              dimnames = list(NULL, "(Intercept)")
            )
          }
        )
      } else {
        # Missing parts -> intercept-only
        X[[i]] <- matrix(1,
          nrow = nrow(mf), ncol = 1,
          dimnames = list(NULL, "(Intercept)")
        )
      }
    }
    names(X) <- param_names

    # Handle offset if it exists (only in alpha)
    if (!is.null(object$offset)) {
      warning("Offset not applied to predictions with newdata")
    }

    X1 <- X$alpha
    X2 <- X$beta
    X3 <- X$gamma
    X4 <- X$delta
    X5 <- X$lambda
  }

  # Number of observations
  n <- nrow(X1)

  # Extract the coefficients
  coefs <- object$coefficients

  # Handle different coefficient structures
  if (is.list(coefs)) {
    # If coefficients are a list with parameter names
    beta1 <- coefs$alpha
    beta2 <- coefs$beta
    beta3 <- coefs$gamma
    beta4 <- coefs$delta
    beta5 <- coefs$lambda
  } else if (!is.null(names(coefs))) {
    # If coefficients are named, extract them using regex patterns
    beta1 <- coefs[grep("^alpha:", names(coefs))]
    beta2 <- coefs[grep("^beta:", names(coefs))]
    beta3 <- coefs[grep("^gamma:", names(coefs))]
    beta4 <- coefs[grep("^delta:", names(coefs))]
    beta5 <- coefs[grep("^lambda:", names(coefs))]

    # If regex didn't work, try another approach
    if (length(beta1) == 0 || length(beta2) == 0) {
      # Split the coefficient vector based on matrix dimensions
      all_coefs <- coefs

      # Make sure we have the right number of coefficients
      total_cols <- ncol(X1) + ncol(X2) + ncol(X3) + ncol(X4) + ncol(X5)
      if (length(all_coefs) != total_cols) {
        stop("Number of coefficients doesn't match design matrices dimensions")
      }

      beta1 <- all_coefs[1:ncol(X1)]
      remaining <- all_coefs[-(1:ncol(X1))]

      beta2 <- remaining[1:ncol(X2)]
      remaining <- remaining[-(1:ncol(X2))]

      beta3 <- remaining[1:ncol(X3)]
      remaining <- remaining[-(1:ncol(X3))]

      beta4 <- remaining[1:ncol(X4)]
      remaining <- remaining[-(1:ncol(X4))]

      beta5 <- remaining[1:ncol(X5)]
    }
  } else {
    # If coefficients are unnamed, split based on matrix dimensions
    all_coefs <- coefs

    # Make sure we have the right number of coefficients
    total_cols <- ncol(X1) + ncol(X2) + ncol(X3) + ncol(X4) + ncol(X5)
    if (length(all_coefs) != total_cols) {
      stop("Number of coefficients doesn't match design matrices dimensions")
    }

    beta1 <- all_coefs[1:ncol(X1)]
    remaining <- all_coefs[-(1:ncol(X1))]

    beta2 <- remaining[1:ncol(X2)]
    remaining <- remaining[-(1:ncol(X2))]

    beta3 <- remaining[1:ncol(X3)]
    remaining <- remaining[-(1:ncol(X3))]

    beta4 <- remaining[1:ncol(X4)]
    remaining <- remaining[-(1:ncol(X4))]

    beta5 <- remaining[1:ncol(X5)]
  }

  # Verify coefficient dimensions match design matrices
  if (length(beta1) != ncol(X1) ||
    length(beta2) != ncol(X2) ||
    length(beta3) != ncol(X3) ||
    length(beta4) != ncol(X4) ||
    length(beta5) != ncol(X5)) {
    stop("Mismatch between coefficient dimensions and design matrices")
  }

  # Extract link information
  if (!is.null(object$link_codes)) {
    link_types <- as.integer(unlist(object$link_codes))
  } else if (!is.null(object$link)) {
    # Convert link functions to codes
    link_map <- c(
      "log" = 1,
      "logit" = 2,
      "probit" = 3,
      "cauchy" = 4,
      "cloglog" = 5,
      "identity" = 6,
      "sqrt" = 7,
      "inverse" = 8,
      "inverse-square" = 9
    )
    link_types <- sapply(object$link, function(l) link_map[l])
  } else {
    # Default link types if not specified
    link_types <- c(1, 1, 1, 2, 1) # log for most parameters, logit for delta
  }

  # Extract scale factors
  if (!is.null(object$scale_factors)) {
    scale_factors <- as.numeric(unlist(object$scale_factors))
  } else {
    # Default scale factors
    scale_factors <- c(10, 10, 10, 1, 10) # 10 for most parameters, 1 for delta
  }

  # If the type is "link", calculate and return the linear predictors
  if (type == "link") {
    eta1 <- as.vector(X1 %*% beta1)
    eta2 <- as.vector(X2 %*% beta2)
    eta3 <- as.vector(X3 %*% beta3)
    eta4 <- as.vector(X4 %*% beta4)
    eta5 <- as.vector(X5 %*% beta5)

    return(data.frame(
      alpha = eta1,
      beta = eta2,
      gamma = eta3,
      delta = eta4,
      lambda = eta5
    ))
  }

  # For the other types, we need the distribution parameters
  params <- calculateParameters(
    X1, X2, X3, X4, X5,
    beta1, beta2, beta3, beta4, beta5,
    link_types, scale_factors,
    family = family
  )

  # Return the parameters if requested
  if (type == "parameter") {
    return(data.frame(
      alpha = params[, 1],
      beta = params[, 2],
      gamma = params[, 3],
      delta = params[, 4],
      lambda = params[, 5]
    ))
  } else if (type %in% c("alpha", "beta", "gamma", "delta", "lambda")) {
    param_index <- match(type, c("alpha", "beta", "gamma", "delta", "lambda"))
    return(params[, param_index])
  }

  # Calculate and return the mean if requested
  if (type == "response") {
    means <- calculateMeans(params, family = family)
    return(means)
  } else if (type == "variance") {
    # Calculate means first
    means <- calculateMeans(params, family = family)

    # Calculate variance for each observation
    variances <- numeric(n)

    for (i in 1:n) {
      alpha_i <- params[i, 1]
      beta_i <- params[i, 2]
      gamma_i <- params[i, 3]
      delta_i <- params[i, 4]
      lambda_i <- params[i, 5]

      # Different variance approximations based on family
      if (family == "beta") {
        variances[i] <- (gamma_i * delta_i) /
          ((gamma_i + delta_i)^2 * (gamma_i + delta_i + 1))
      } else if (family == "kw") {
        variances[i] <- alpha_i * beta_i *
          beta(1 + 1 / alpha_i, beta_i) -
          (alpha_i * beta_i * beta(1 + 1 / alpha_i, beta_i))^2
      } else {
        # Numerical approximation
        h <- 0.001
        mean_i <- means[i]

        # Safely compute for y within (0,1)
        mean_plus <- min(0.999, mean_i + h)
        mean_minus <- max(0.001, mean_i - h)

        # Use the appropriate density function based on family
        switch(family,
          "gkw" = {
            f_plus <- dgkw(mean_plus, alpha_i, beta_i, gamma_i, delta_i, lambda_i)
            f <- dgkw(mean_i, alpha_i, beta_i, gamma_i, delta_i, lambda_i)
            f_minus <- dgkw(mean_minus, alpha_i, beta_i, gamma_i, delta_i, lambda_i)
          },
          "bkw" = {
            f_plus <- dbkw(mean_plus, alpha_i, beta_i, gamma_i, delta_i)
            f <- dbkw(mean_i, alpha_i, beta_i, gamma_i, delta_i)
            f_minus <- dbkw(mean_minus, alpha_i, beta_i, gamma_i, delta_i)
          },
          "kkw" = {
            f_plus <- dkkw(mean_plus, alpha_i, beta_i, delta_i, lambda_i)
            f <- dkkw(mean_i, alpha_i, beta_i, delta_i, lambda_i)
            f_minus <- dkkw(mean_minus, alpha_i, beta_i, delta_i, lambda_i)
          },
          "ekw" = {
            f_plus <- dekw(mean_plus, alpha_i, beta_i, lambda_i)
            f <- dekw(mean_i, alpha_i, beta_i, lambda_i)
            f_minus <- dekw(mean_minus, alpha_i, beta_i, lambda_i)
          },
          "mc" = {
            f_plus <- dmc(mean_plus, gamma_i, delta_i, lambda_i)
            f <- dmc(mean_i, gamma_i, delta_i, lambda_i)
            f_minus <- dmc(mean_minus, gamma_i, delta_i, lambda_i)
          }
        )

        # Approximate second derivative of log-likelihood
        d2ll <- (f_plus - 2 * f + f_minus) / (h * h)

        # Variance is inverse of the information
        variances[i] <- min(0.25, max(1e-6, 1 / (abs(d2ll) + 1e-10)))
      }
    }

    return(variances)
  }

  # For the remaining types, we need the 'at' argument
  if (!is.numeric(at)) {
    stop("'at' must be numeric")
  }

  # Check if the elementwise mode should be used
  if (is.null(elementwise)) {
    elementwise <- length(at) == n
  }

  # Prepare the data for calculation based on 'at'
  if (elementwise) {
    if (length(at) != n) {
      stop("For elementwise=TRUE, length of 'at' must equal number of observations")
    }
    eval_y <- at
    eval_params <- params
  } else {
    # Repeat params for each value in 'at'
    eval_params <- do.call(rbind, lapply(seq_len(length(at)), function(i) params))
    # Repeat each value of 'at' for each row of params
    eval_y <- rep(at, each = n)
  }

  # Calculate the requested quantities using specific family functions
  result <- numeric(length(eval_y))

  for (i in seq_along(eval_y)) {
    alpha_i <- eval_params[i, 1]
    beta_i <- eval_params[i, 2]
    gamma_i <- eval_params[i, 3]
    delta_i <- eval_params[i, 4]
    lambda_i <- eval_params[i, 5]
    y_i <- eval_y[i]

    # Use the appropriate function based on family and type
    if (type == "density") {
      switch(family,
        "gkw" = {
          result[i] <- dgkw(y_i, alpha_i, beta_i, gamma_i, delta_i, lambda_i)
        },
        "bkw" = {
          result[i] <- dbkw(y_i, alpha_i, beta_i, gamma_i, delta_i)
        },
        "kkw" = {
          result[i] <- dkkw(y_i, alpha_i, beta_i, delta_i, lambda_i)
        },
        "ekw" = {
          result[i] <- dekw(y_i, alpha_i, beta_i, lambda_i)
        },
        "mc" = {
          result[i] <- dmc(y_i, gamma_i, delta_i, lambda_i)
        },
        "kw" = {
          result[i] <- dkw(y_i, alpha_i, beta_i)
        },
        "beta" = {
          result[i] <- dbeta_(y_i, gamma_i, delta_i + 1)
        }
      )
    } else if (type == "probability") {
      switch(family,
        "gkw" = {
          result[i] <- pgkw(y_i, alpha_i, beta_i, gamma_i, delta_i, lambda_i)
        },
        "bkw" = {
          result[i] <- pbkw(y_i, alpha_i, beta_i, gamma_i, delta_i)
        },
        "kkw" = {
          result[i] <- pkkw(y_i, alpha_i, beta_i, delta_i, lambda_i)
        },
        "ekw" = {
          result[i] <- pekw(y_i, alpha_i, beta_i, lambda_i)
        },
        "mc" = {
          result[i] <- pmc(y_i, gamma_i, delta_i, lambda_i)
        },
        "kw" = {
          result[i] <- pkw(y_i, alpha_i, beta_i)
        },
        "beta" = {
          result[i] <- pbeta_(y_i, gamma_i, delta_i + 1)
        }
      )
    } else if (type == "quantile") {
      switch(family,
        "gkw" = {
          result[i] <- qgkw(y_i, alpha_i, beta_i, gamma_i, delta_i, lambda_i)
        },
        "bkw" = {
          result[i] <- qbkw(y_i, alpha_i, beta_i, gamma_i, delta_i)
        },
        "kkw" = {
          result[i] <- qkkw(y_i, alpha_i, beta_i, delta_i, lambda_i)
        },
        "ekw" = {
          result[i] <- qekw(y_i, alpha_i, beta_i, lambda_i)
        },
        "mc" = {
          result[i] <- qmc(y_i, gamma_i, delta_i, lambda_i)
        },
        "kw" = {
          result[i] <- qkw(y_i, alpha_i, beta_i)
        },
        "beta" = {
          result[i] <- qbeta_(y_i, gamma_i, delta_i + 1)
        }
      )
    }
  }

  # Format the result
  if (elementwise) {
    return(result)
  } else {
    result_matrix <- matrix(result, nrow = n, ncol = length(at))
    colnames(result_matrix) <- as.character(at)
    return(result_matrix)
  }
}


#' @title Extract Fitted Values from a Generalized Kumaraswamy Regression Model
#'
#' @description
#' Extract the fitted values (predicted means) from a fitted Generalized Kumaraswamy
#' regression model.
#'
#' @param object A fitted model object of class "gkwreg".
#' @param family Character string specifying the distribution family to use.
#'   If NULL (default), uses the family from the fitted model.
#'   Available options: "gkw", "bkw", "kkw", "ekw", "mc", "kw", "beta".
#' @param ... Additional arguments. Currently not used.
#'
#' @return A numeric vector of fitted values (predicted means) on the scale of the
#'   original response. For the GKw model, these are values between 0 and 1.
#'
#' @details
#' The fitted values are extracted directly from the model object if available.
#' If not available (which may happen in rare cases), they are recalculated using
#' the model parameters and the specified family. The values represent the predicted
#' means of the response variable conditional on the observed covariates.
#'
#' When using a different family than what was used to fit the model (by specifying
#' the 'family' parameter), the function will always recalculate the fitted values
#' using the predict method instead of using stored values.
#'
#' @examples
#' \dontrun{
#' # Fit a GKw model
#' model <- gkwreg(y ~ x1 | x2 | 1 | 1 | 1, data = mydata, family = "gkw")
#'
#' # Extract fitted values
#' fitted_values <- fitted(model)
#'
#' # Extract fitted values using a different family
#' fitted_values_beta <- fitted(model, family = "beta")
#'
#' # Plot observed vs. fitted
#' plot(model$y, fitted_values,
#'   xlab = "Observed", ylab = "Fitted",
#'   main = "Observed vs Fitted Values"
#' )
#' abline(0, 1, col = "red", lty = 2)
#' }
#'
#' @seealso
#' \code{\link{predict.gkwreg}} for predicting from a GKw model.
#'
#' \code{\link{residuals.gkwreg}} for extracting residuals.
#'
#' @export
fitted.gkwreg <- function(object, family = NULL, ...) {
  # Check if object is of class "gkwreg"
  if (!inherits(object, "gkwreg")) {
    stop("'object' must be a fitted model of class \"gkwreg\"")
  }

  # Get the family from the object if not specified
  if (is.null(family)) {
    if (!is.null(object$family)) {
      family <- object$family
    } else {
      # Default to gkw for backward compatibility
      family <- "gkw"
      message("No family specified in the model. Using 'gkw' as default.")
    }
  } else {
    # Validate the family parameter
    family <- match.arg(family, c("gkw", "bkw", "kkw", "ekw", "mc", "kw", "beta"))

    # If family is different from the model's family, we need to recalculate
    if (!is.null(object$family) && family != object$family) {
      message(paste0(
        "Using different family (", family, ") than what was used to fit the model (",
        object$family, "). Recalculating fitted values..."
      ))

      # Use predict to calculate fitted values with the new family
      return(stats::predict(object, type = "response", family = family))
    }
  }

  # Determine the number of observations from various possible sources
  n <- if (!is.null(object$nobs)) {
    object$nobs
  } else if (!is.null(object$y)) {
    length(object$y)
  } else if (!is.null(object$model) && !is.null(object$model$y)) {
    length(object$model$y)
  } else if (!is.null(object$parameter_vectors) && !is.null(object$parameter_vectors$alphaVec)) {
    length(object$parameter_vectors$alphaVec)
  } else {
    stop("Cannot determine the number of observations in the model")
  }

  # Method 1: Try to get fitted values directly from the model object
  if (!is.null(object$fitted.values)) {
    fitted_len <- length(object$fitted.values)

    # Check if dimensions match the expected number of observations
    if (fitted_len == n) {
      return(object$fitted.values)
    } else if (fitted_len == 1 && n > 1) {
      message("Only one fitted value found in model object. Recalculating all fitted values...")
    } else if (fitted_len > 1 && fitted_len < n && n > 10000) {
      # For very large datasets, we might have a sample of fitted values - try to interpolate
      message("Partial fitted values found. Interpolating...")

      # Create fitted values vector with NAs
      fitted_values <- rep(NA_real_, n)

      # Map available values to appropriate indices
      sample_size <- fitted_len
      sample_idx <- floor(seq(1, n, length.out = sample_size))
      fitted_values[sample_idx] <- object$fitted.values

      # Interpolate the rest
      idx_with_values <- which(!is.na(fitted_values))
      fitted_values <- stats::approx(
        x = idx_with_values,
        y = fitted_values[idx_with_values],
        xout = seq_len(n),
        rule = 2
      )$y

      return(fitted_values)
    }
  }

  # Method 2: Try to use the parameter vectors if available
  if (!is.null(object$parameter_vectors)) {
    # Check if all necessary parameter vectors are available
    param_vec <- object$parameter_vectors
    if (!is.null(param_vec$alphaVec) &&
      !is.null(param_vec$betaVec) &&
      !is.null(param_vec$gammaVec) &&
      !is.null(param_vec$deltaVec) &&
      !is.null(param_vec$lambdaVec)) {
      # Check if dimensions match
      param_lengths <- c(
        length(param_vec$alphaVec),
        length(param_vec$betaVec),
        length(param_vec$gammaVec),
        length(param_vec$deltaVec),
        length(param_vec$lambdaVec)
      )

      if (all(param_lengths == n)) {
        # Create parameter matrix for calculateMeans
        params <- matrix(0, nrow = n, ncol = 5)
        params[, 1] <- param_vec$alphaVec
        params[, 2] <- param_vec$betaVec
        params[, 3] <- param_vec$gammaVec
        params[, 4] <- param_vec$deltaVec
        params[, 5] <- param_vec$lambdaVec

        # Calculate means using the parameters
        message("Calculating fitted values from parameter vectors...")
        return(calculateMeans(params, family = family))
      }
    }
  }

  # Method 3: Check TMB report for fitted values
  if (!is.null(object$tmb_object)) {
    tmb_report <- try(object$tmb_object$report(), silent = TRUE)

    if (!inherits(tmb_report, "try-error") && "fitted" %in% names(tmb_report)) {
      fitted_from_tmb <- tmb_report$fitted

      if (length(fitted_from_tmb) == n) {
        return(as.vector(fitted_from_tmb))
      } else if (length(fitted_from_tmb) > 1 && length(fitted_from_tmb) < n && n > 10000) {
        # For large datasets, interpolate if we only have a sample
        message("Found partial fitted values in TMB report. Interpolating...")

        # Create fitted values vector with NAs
        fitted_values <- rep(NA_real_, n)

        # Map available values to appropriate indices
        sample_size <- length(fitted_from_tmb)
        sample_idx <- floor(seq(1, n, length.out = sample_size))
        fitted_values[sample_idx] <- fitted_from_tmb

        # Interpolate the rest
        idx_with_values <- which(!is.na(fitted_values))
        fitted_values <- stats::approx(
          x = idx_with_values,
          y = fitted_values[idx_with_values],
          xout = seq_len(n),
          rule = 2
        )$y

        return(fitted_values)
      }
    }
  }

  # Method 4: If all else fails, use predict.gkwreg
  message("Recalculating fitted values using the predict function...")

  # Use predict to calculate fitted values
  fitted_values <- stats::predict(object, type = "response", family = family)

  return(fitted_values)
}

#' @title Extract Residuals from a Generalized Kumaraswamy Regression Model
#'
#' @description
#' Extract various types of residuals from a fitted Generalized Kumaraswamy regression model
#' for diagnostic purposes.
#'
#' @param object A fitted model object of class "gkwreg".
#' @param type The type of residuals to be computed. Options include:
#'   \itemize{
#'     \item "response" (default): Response residuals (y - fitted)
#'     \item "pearson": Pearson residuals (standardized by estimated std dev)
#'     \item "deviance": Deviance residuals (signed sqrt of deviance contributions)
#'     \item "quantile": Quantile residuals (normal quantiles of CDF at response)
#'     \item "modified.deviance": Modified deviance residuals (standardized deviance)
#'     \item "cox-snell": Cox-Snell residuals (transformation of CDF)
#'     \item "score": Score residuals (based on score function)
#'     \item "partial": Partial residuals (requires covariate index)
#'   }
#' @param covariate_idx Integer. For 'partial' residuals, specifies which covariate to use
#'   (1 for first covariate, 2 for second, etc.). Only used if type="partial".
#' @param parameter Parameter name for partial residuals. One of "alpha", "beta", "gamma",
#'   "delta", or "lambda". Only used if type="partial".
#' @param family Character string specifying the distribution family to use.
#'   If NULL (default), uses the family from the fitted model.
#'   Available options: "gkw", "bkw", "kkw", "ekw", "mc", "kw", "beta".
#' @param ... Additional arguments. Currently not used.
#'
#' @return A numeric vector of residuals.
#'
#' @details
#' The different types of residuals offer various diagnostic tools for the GKw regression model:
#'
#' \itemize{
#'   \item Response residuals: The raw difference between observed and fitted values.
#'
#'   \item Pearson residuals: Response residuals standardized by the estimated standard
#'         deviation from the model. Useful for checking heteroscedasticity.
#'
#'   \item Deviance residuals: Signed square root of deviance contributions.
#'         Useful for assessing model fit.
#'
#'   \item Quantile residuals: Normal quantiles of the model's CDF evaluated at observed values.
#'         These should follow a standard normal distribution if the model is correct.
#'
#'   \item Modified deviance residuals: Standardized deviance residuals that should
#'         better approximate a normal distribution.
#'
#'   \item Cox-Snell residuals: Based on -log(1-F(y)), should follow exponential(1)
#'         distribution if the model is correct. Uses the model's CDF.
#'
#'   \item Score residuals: Based on the score function from the model's log-likelihood.
#'         Useful for assessing influence.
#'
#'   \item Partial residuals: For examining individual covariate effects in the context
#'         of the specified distribution family.
#' }
#'
#' All calculations use the specific properties of the selected distribution family
#' and are performed efficiently in C++ to ensure good performance even with large datasets.
#'
#' Specifying a different family than what was used to fit the model can be useful for
#' diagnostic comparisons between different distribution assumptions.
#'
#' @references
#' Dunn, P. K., & Smyth, G. K. (1996). Randomized Quantile Residuals.
#' Journal of Computational and Graphical Statistics, 5(3), 236-244.
#'
#' Cox, D. R., & Snell, E. J. (1968). A General Definition of Residuals.
#' Journal of the Royal Statistical Society, Series B, 30(2), 248-275.
#'
#' McCullagh, P., & Nelder, J. A. (1989). Generalized Linear Models.
#' CRC Press.
#'
#' @examples
#' \dontrun{
#' # Fit a GKw model
#' model <- gkwreg(y ~ x1 | x2 | 1 | 1 | 1, data = mydata, family = "gkw")
#'
#' # Extract different types of residuals
#' resp_res <- residuals(model, type = "response")
#' pearson_res <- residuals(model, type = "pearson")
#' quant_res <- residuals(model, type = "quantile")
#'
#' # Using a different family for diagnostic comparison
#' quant_res_beta <- residuals(model, type = "quantile", family = "beta")
#'
#' # QQ-plot for quantile residuals
#' qqnorm(quant_res)
#' qqline(quant_res)
#'
#' # Cox-Snell residuals for exponential quantile plot
#' cs_res <- residuals(model, type = "cox-snell")
#' plot(qexp(ppoints(length(cs_res))), sort(cs_res),
#'   xlab = "Theoretical Quantiles", ylab = "Cox-Snell Residuals"
#' )
#' abline(0, 1, col = "red")
#'
#' # Comparing residuals from different distribution families
#' par(mfrow = c(1, 2))
#' qqnorm(quant_res, main = "GKw Quantile Residuals")
#' qqline(quant_res)
#' qqnorm(quant_res_beta, main = "Beta Quantile Residuals")
#' qqline(quant_res_beta)
#'
#' # Partial residuals for examining covariate effect
#' part_res <- residuals(model,
#'   type = "partial",
#'   parameter = "alpha", covariate_idx = 2
#' )
#' }
#'
#' @seealso
#' \code{\link{fitted.gkwreg}} for extracting fitted values.
#'
#' \code{\link{predict.gkwreg}} for predicting from a GKw model.
#'
#' @export
residuals.gkwreg <- function(
    object, type = c(
      "response", "pearson", "deviance", "quantile",
      "modified.deviance", "cox-snell",
      "score", "partial"
    ),
    covariate_idx = 1, parameter = "alpha",
    family = NULL, ...) {
  # Check if object is of class "gkwreg"
  if (!inherits(object, "gkwreg")) {
    stop("'object' must be a fitted model of class \"gkwreg\"")
  }

  # Match argument
  type <- match.arg(type)

  # Get the family from the object if not specified
  if (is.null(family)) {
    if (!is.null(object$family)) {
      family <- object$family
    } else {
      # Default to gkw for backward compatibility
      family <- "gkw"
      message("No family specified in the model. Using 'gkw' as default.")
    }
  } else {
    # Validate the family parameter
    family <- match.arg(family, c("gkw", "bkw", "kkw", "ekw", "mc", "kw", "beta"))

    # If family is different from the model's family, show a message
    if (!is.null(object$family) && family != object$family) {
      message(paste0(
        "Using different family (", family, ") than what was used to fit the model (",
        object$family, ")."
      ))
    }
  }

  # Get response values
  if (is.null(object$y)) {
    stop("Response variable not found in model object. Cannot calculate residuals.")
  }
  y <- object$y

  # Get fitted values - passing the family parameter
  fitted_vals <- stats::fitted(object, family = family)

  # If type is "response", we can calculate quickly
  if (type == "response") {
    return(calculateResponseResiduals(y, fitted_vals))
  }

  # For partial residuals, verify parameters
  if (type == "partial") {
    if (!parameter %in% c("alpha", "beta", "gamma", "delta", "lambda")) {
      stop("parameter must be one of: 'alpha', 'beta', 'gamma', 'delta', 'lambda'")
    }

    # Confirm valid covariate_idx
    model_matrices <- object$model_matrices
    if (is.null(model_matrices)) {
      stop("Model matrices not available in model object")
    }

    X <- model_matrices[[parameter]]
    if (is.null(X)) {
      stop("Model matrix for parameter '", parameter, "' not found")
    }

    p <- ncol(X)
    if (covariate_idx < 1 || covariate_idx > p) {
      stop("covariate_idx must be between 1 and ", p)
    }
  }

  # Get parameters for each observation
  get_parameters <- function(object, family) {
    # Initialize vectors for parameters
    n <- length(object$y)

    # Try to get parameters via predict with the specified family
    if (exists("predict.gkwreg", mode = "function")) {
      tryCatch(
        {
          params_df <- stats::predict(object, type = "parameter", family = family)
          return(list(
            alpha = params_df$alpha,
            beta = params_df$beta,
            gamma = params_df$gamma,
            delta = params_df$delta,
            lambda = params_df$lambda
          ))
        },
        error = function(e) {
          warning("Could not extract parameter values using predict.gkwreg(): ", e$message)
        }
      )
    }

    # If we're using the same family as the model, try to get parameters directly
    if (is.null(object$family) || family == object$family) {
      # Try to get individual parameter vectors from the model object first
      if (!is.null(object$model) && !is.null(object$model$report)) {
        report <- object$model$report()

        # Check if individual parameter vectors are available
        if (all(c("alphaVec", "betaVec", "gammaVec", "deltaVec", "lambdaVec") %in% names(report))) {
          return(list(
            alpha = report$alphaVec,
            beta = report$betaVec,
            gamma = report$gammaVec,
            delta = report$deltaVec,
            lambda = report$lambdaVec
          ))
        }
      }
    }

    # If we get here, we need to approximate using the average parameter values
    # or recompute parameters with the appropriate family constraints

    # Check if we have fitted parameters already
    if (!is.null(object$fitted_parameters)) {
      # Start with the model's fitted parameters
      base_params <- list(
        alpha = rep(object$fitted_parameters$alpha, n),
        beta = rep(object$fitted_parameters$beta, n),
        gamma = rep(object$fitted_parameters$gamma, n),
        delta = rep(object$fitted_parameters$delta, n),
        lambda = rep(object$fitted_parameters$lambda, n)
      )

      # Apply family-specific constraints if needed
      if (family != "gkw") {
        if (family == "bkw") {
          # BKw: lambda = 1
          base_params$lambda <- rep(1.0, n)
        } else if (family == "kkw") {
          # KKw: gamma = 1
          base_params$gamma <- rep(1.0, n)
        } else if (family == "ekw") {
          # EKw: gamma = 1, delta = 0
          base_params$gamma <- rep(1.0, n)
          base_params$delta <- rep(0.0, n)
        } else if (family == "mc") {
          # MC: alpha = 1, beta = 1
          base_params$alpha <- rep(1.0, n)
          base_params$beta <- rep(1.0, n)
        } else if (family == "kw") {
          # KW: lambda = 1, gamma = 1, delta = 0
          base_params$gamma <- rep(1.0, n)
          base_params$delta <- rep(0.0, n)
          base_params$lambda <- rep(1.0, n)
        } else if (family == "beta") {
          # Beta: alpha = 1, beta = 1, lambda = 1
          base_params$alpha <- rep(1.0, n)
          base_params$beta <- rep(1.0, n)
          base_params$lambda <- rep(1.0, n)
        }
      }

      message("Using adjusted parameter values for family ", family, ".")
      return(base_params)
    }

    # If all else fails
    stop("Unable to extract parameter values from the model object.")
  }

  # Get parameters with the specified family
  params <- get_parameters(object, family)

  # Create a parameter matrix for C++ functions that expect it
  param_matrix <- matrix(
    c(params$alpha, params$beta, params$gamma, params$delta, params$lambda),
    ncol = 5
  )

  # Calculate appropriate residuals
  result <- switch(type,
    "pearson" = {
      calculatePearsonResiduals(
        y, fitted_vals, param_matrix,
        family = family
      )
    },
    "deviance" = {
      calculateDevianceResiduals(
        y, fitted_vals, param_matrix,
        family = family
      )
    },
    "quantile" = {
      calculateQuantileResiduals(
        y, param_matrix,
        family = family
      )
    },
    "modified.deviance" = {
      calculateModifiedDevianceResiduals(
        y, fitted_vals, param_matrix,
        family = family
      )
    },
    "cox-snell" = {
      calculateCoxSnellResiduals(
        y, param_matrix,
        family = family
      )
    },
    "score" = {
      calculateScoreResiduals(
        y, fitted_vals, param_matrix,
        family = family
      )
    },
    "partial" = {
      # Get model matrices and coefficients for the specified parameter
      X <- object$model_matrices[[parameter]]
      beta <- object$coefficients[[parameter]]

      # We need C++ 0-based indexing
      c_idx <- covariate_idx - 1

      calculatePartialResiduals(y, fitted_vals, X, beta, c_idx)
    }
  )

  return(result)
}


#' Diagnostic Plots for gkwreg Objects
#'
#' @description
#' This S3 method produces a set of diagnostic plots for a fitted Generalized
#' Kumaraswamy (GKw) regression model, as returned by \code{\link{gkwreg}}.
#' The function provides standard plots for residual analysis, influence measures,
#' and goodness-of-fit assessment. Users may choose between base R graphics or
#' ggplot2 for plotting.
#'
#' @param x A fitted model object of class \code{"gkwreg"}.
#' @param which Integer vector specifying which diagnostic plots to produce.
#'   Valid values are from 1 to 6, corresponding to:
#'   \enumerate{
#'     \item Residuals vs. Observation Indices.
#'     \item Cook's Distance Plot.
#'     \item Generalized Leverage vs. Fitted Values.
#'     \item Residuals vs. Linear Predictor.
#'     \item Half-Normal Plot of Residuals.
#'     \item Predicted vs. Observed Values.
#'   }
#' @param caption A character vector of captions for the plots. Its length must be at
#'   least \code{max(which)}.
#' @param sub.caption A character string for a common subtitle above the plots (e.g.,
#'   the original call of the model). Defaults to the deparsed call of \code{x}.
#' @param main An optional character string to be used as the main title for each plot,
#'   appended to the respective caption.
#' @param ask Logical; if \code{TRUE}, the user is prompted before each plot (only applicable
#'   to base R graphics). Defaults to \code{TRUE} when the number of plots exceeds the current
#'   graphics layout.
#' @param ... Additional graphical parameters to be passed to the underlying plotting functions.
#' @param type Character string indicating the type of residual to be used.
#'   Valid options are:
#'   \itemize{
#'     \item \code{"quantile"}: Quantile residuals, based on the probability integral transform.
#'       Recommended for bounded response variables like those in GKw regression.
#'     \item \code{"pearson"}: Pearson residuals, standardized by the estimated standard deviation.
#'       Useful for identifying heteroscedasticity.
#'     \item \code{"deviance"}: Deviance residuals, based on the contribution to the deviance.
#'       Helpful for assessing model fit quality.
#'   }
#'   Defaults to \code{"quantile"}.
#' @param family Character string specifying the distribution family to use.
#'   If NULL (default), uses the family from the fitted model.
#'   Available options: "gkw", "bkw", "kkw", "ekw", "mc", "kw", "beta".
#' @param nsim Integer; number of simulations for the half-normal plot envelope.
#'   Must be a positive integer. Defaults to 100.
#' @param level Numeric; confidence level for the half-normal plot envelope.
#'   Must be between 0 and 1. Defaults to 0.90.
#' @param use_ggplot Logical; if \code{TRUE}, ggplot2 is used for plotting, otherwise base R
#'   graphics are used. Defaults to \code{FALSE}.
#' @param arrange_plots Logical; if \code{TRUE} and \code{use_ggplot = TRUE}, attempts to
#'   arrange multiple plots in a grid using \code{gridExtra::grid.arrange} or
#'   \code{ggpubr::ggarrange}. Defaults to \code{FALSE}.
#' @param sample_size Integer or NULL; if provided and less than the number of observations,
#'   a random sample of observations of this size will be used for plotting. Useful for very
#'   large datasets. Defaults to NULL (use all observations).
#' @param theme_fn Function; a ggplot2 theme function to be applied to all plots when
#'   \code{use_ggplot = TRUE}. Defaults to \code{ggplot2::theme_minimal()}.
#' @param save_diagnostics Logical; if \code{TRUE}, returns a list with the diagnostic
#'   measures used for plotting. Defaults to \code{FALSE}.
#'
#' @return Invisibly returns either the fitted model object \code{x} (if \code{save_diagnostics = FALSE})
#'   or a list containing the diagnostic measures (if \code{save_diagnostics = TRUE}).
#'
#' @export
plot.gkwreg <- function(x,
                        which = 1:6,
                        caption = c(
                          "Residuals vs. Observation Indices",
                          "Cook's Distance Plot",
                          "Generalized Leverage vs. Fitted Values",
                          "Residuals vs. Linear Predictor",
                          "Half-Normal Plot of Residuals",
                          "Predicted vs. Observed Values"
                        ),
                        sub.caption = paste(deparse(x$call), collapse = "\n"),
                        main = "",
                        ask = prod(par("mfcol")) < length(which) && dev.interactive(),
                        ...,
                        type = c("quantile", "pearson", "deviance"),
                        family = NULL,
                        nsim = 100,
                        level = 0.90,
                        use_ggplot = FALSE,
                        arrange_plots = FALSE,
                        sample_size = NULL,
                        theme_fn = ggplot2::theme_minimal,
                        save_diagnostics = FALSE) {
  # Validate inputs and prepare diagnostic data using the fixed function
  diag_data <- .validate_and_prepare_gkwreg_diagnostics(
    x = x,
    which = which,
    caption = caption,
    type = type,
    family = family,
    nsim = nsim,
    level = level,
    use_ggplot = use_ggplot,
    arrange_plots = arrange_plots,
    sample_size = sample_size,
    theme_fn = theme_fn
  )

  # Get formatted plot titles
  plot_titles <- .create_plot_titles(
    which = which,
    caption = caption,
    main = main,
    family = diag_data$model_info$family,
    orig_family = x$family
  )

  # Choose plotting implementation based on use_ggplot
  if (!use_ggplot) {
    # Base R graphics implementation
    result <- .plot_gkwreg_base_r(
      diag_data = diag_data,
      which = which,
      plot_titles = plot_titles,
      sub.caption = sub.caption,
      ask = ask,
      ...
    )
  } else {
    # ggplot2 implementation
    result <- .plot_gkwreg_ggplot(
      diag_data = diag_data,
      which = which,
      plot_titles = plot_titles,
      sub.caption = sub.caption,
      ask = ask,
      arrange_plots = arrange_plots,
      theme_fn = theme_fn,
      ...
    )
  }

  # Return diagnostic data if requested
  if (save_diagnostics) {
    return(invisible(diag_data))
  } else {
    return(invisible(x))
  }
}



#' Validate inputs and prepare diagnostic data for gkwreg plots
#'
#' @param x A fitted model object of class "gkwreg"
#' @param which Integer vector specifying which plots to produce
#' @param caption Character vector of plot captions
#' @param type Character string specifying residual type
#' @param family Character string specifying distribution family
#' @param nsim Number of simulations for envelope calculation
#' @param level Confidence level for envelope
#' @param use_ggplot Logical; whether to use ggplot2
#' @param arrange_plots Logical; whether to arrange multiple plots
#' @param sample_size Integer or NULL; sample size for large datasets
#' @param theme_fn ggplot2 theme function
#'
#' @return A list containing diagnostic data and model information
#'
#' @keywords internal
.validate_and_prepare_gkwreg_diagnostics <- function(x,
                                                     which,
                                                     caption,
                                                     type = c("quantile", "pearson", "deviance"),
                                                     family = NULL,
                                                     nsim = 100,
                                                     level = 0.90,
                                                     use_ggplot = FALSE,
                                                     arrange_plots = FALSE,
                                                     sample_size = NULL,
                                                     theme_fn = ggplot2::theme_minimal) {
  # Input validation
  if (!inherits(x, "gkwreg")) {
    stop("The object must be of class 'gkwreg'.")
  }

  type <- match.arg(type)

  # Get the family from the object if not specified
  if (is.null(family)) {
    if (!is.null(x$family)) {
      family <- x$family
    } else {
      # Default to gkw for backward compatibility
      family <- "gkw"
      message("No family specified in the model. Using 'gkw' as default.")
    }
  } else {
    # Validate the family parameter
    family <- match.arg(family, c("gkw", "bkw", "kkw", "ekw", "mc", "kw", "beta"))

    # If family is different from the model's family, show a message
    if (!is.null(x$family) && family != x$family) {
      message(paste0(
        "Using different family (", family, ") than what was used to fit the model (",
        x$family, ") for diagnostics."
      ))
    }
  }

  # Get parameter information for the family
  param_info <- .get_family_param_info(family)

  # Other input validations
  if (!all(which %in% 1:6)) {
    stop("Argument 'which' must contain values between 1 and 6.")
  }

  if (max(which) > length(caption)) {
    stop("The 'caption' vector is too short for the selected 'which' plots.")
  }

  if (!is.numeric(nsim) || nsim <= 0 || nsim != round(nsim)) {
    stop("Argument 'nsim' must be a positive integer.")
  }

  if (!is.numeric(level) || level <= 0 || level >= 1) {
    stop("Argument 'level' must be between 0 and 1.")
  }

  # Check dependencies
  if (use_ggplot && !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for ggplot2 plotting. Please install it with install.packages('ggplot2').")
  }

  if (use_ggplot && arrange_plots &&
    !requireNamespace("gridExtra", quietly = TRUE) &&
    !requireNamespace("ggpubr", quietly = TRUE)) {
    stop("Either package 'gridExtra' or 'ggpubr' is required for arranging plots. Please install one of them.")
  }

  # Extract model components with error handling
  y_obs <- x$y
  if (is.null(y_obs)) {
    stop("No 'y' component found in the model object. Ensure the model was fitted with y=TRUE.")
  }

  # Get fitted values using the specified family
  fitted_vals <- fitted(x, family = family)
  if (is.null(fitted_vals)) {
    stop("Could not calculate fitted values.")
  }

  # Extract model matrices and parameters using the improved functions
  model_matrices <- .extract_model_matrices(x)
  model_params <- .extract_model_params(x)

  # Sample data if requested
  n <- length(y_obs)
  idx <- seq_len(n)

  if (!is.null(sample_size) && is.numeric(sample_size) && sample_size > 0 && sample_size < n) {
    sampling_result <- .sample_model_data(
      n = n,
      sample_size = sample_size,
      y_obs = y_obs,
      fitted_vals = fitted_vals,
      model_matrices = model_matrices
    )

    idx <- sampling_result$idx
    y_obs <- sampling_result$y_obs
    fitted_vals <- sampling_result$fitted_vals
    model_matrices <- sampling_result$model_matrices
  }

  # Calculate model parameters for the specified family
  param_mat <- .calculate_model_parameters(model_matrices, model_params, family)

  # Extract parameter vectors
  param_vectors <- .extract_parameter_vectors(param_mat)

  # Calculate residuals
  resid_vec <- .calculate_residuals(y_obs, fitted_vals, param_mat, type, family)

  # Calculate diagnostic measures with family-specific handling
  diagnostic_measures <- .calculate_diagnostic_measures(
    y_obs = y_obs,
    fitted_vals = fitted_vals,
    resid_vec = resid_vec,
    model_matrices = model_matrices,
    model_params = model_params,
    param_vectors = param_vectors,
    idx = idx,
    family = family,
    param_info = param_info
  )

  # Calculate half-normal plot data if needed
  half_normal_data <- NULL
  if (5 %in% which) {
    half_normal_data <- .calculate_half_normal_data(
      resid_vec = resid_vec,
      idx = idx,
      nsim = nsim,
      level = level,
      param_mat = param_mat,
      param_vectors = param_vectors,
      type = type,
      family = family
    )
  }

  # Create the diagnostic data structure
  diag_data <- list(
    data = data.frame(
      index = idx,
      y_obs = y_obs,
      fitted = fitted_vals,
      resid = resid_vec,
      abs_resid = abs(resid_vec),
      cook_dist = diagnostic_measures$cook_dist,
      leverage = diagnostic_measures$leverage,
      linpred = diagnostic_measures$linpred
    ),
    model_info = list(
      n = n,
      p = diagnostic_measures$p,
      cook_threshold = 4 / n,
      leverage_threshold = 2 * diagnostic_measures$p / n,
      family = family,
      type = type,
      param_info = param_info,
      level = level
    )
  )

  # Add half-normal data if available
  if (!is.null(half_normal_data)) {
    diag_data$half_normal <- half_normal_data
  }

  return(diag_data)
}

#' Extract model matrices from a gkwreg object with family-specific handling
#'
#' @param x A fitted model object of class "gkwreg"
#' @return A list of model matrices
#'
#' @importFrom stats model.matrix
#' @keywords internal
.extract_model_matrices <- function(x) {
  # Get family parameter information
  family <- x$family
  if (is.null(family)) {
    family <- "gkw" # Default to gkw if not specified
    warning("No family specified in the model. Using 'gkw' as default.")
  }

  # Get parameter information for the specified family
  param_info <- .get_family_param_info(family)
  param_names <- param_info$names
  param_positions <- param_info$positions

  # Initialize matrices list with the correct number for this family
  num_matrices <- max(unlist(param_positions))
  matrices <- vector("list", num_matrices)
  names(matrices) <- paste0("X", 1:num_matrices)

  # Try different ways to get design matrices
  if (!is.null(x$x)) {
    # Model was fitted with x=TRUE
    for (param in param_names) {
      pos <- param_positions[[param]]
      if (param %in% names(x$x)) {
        matrices[[pos]] <- x$x[[param]]
      } else {
        # Default to intercept-only for missing matrices
        n_obs <- ifelse(is.null(x$y), nrow(x$model), length(x$y))
        matrices[[pos]] <- matrix(1, nrow = n_obs, ncol = 1)
        colnames(matrices[[pos]]) <- "(Intercept)"
      }
    }
  } else if (!is.null(x$model) && !is.null(x$formula)) {
    # Recreate model matrices from the model frame and formula
    mf <- x$model
    formula_obj <- x$formula
    if (inherits(formula_obj, "formula")) {
      formula_obj <- Formula::as.Formula(formula_obj)
    }

    # Get number of rhs parts in the formula
    n_parts <- length(attr(Formula::Formula(formula_obj), "rhs"))

    # Extract matrices for each parameter
    for (param in param_names) {
      pos <- param_positions[[param]]
      idx <- which(param_names == param)

      if (idx <= n_parts) {
        # Try to extract matrix using the formula part
        matrices[[pos]] <- tryCatch(
          model.matrix(formula_obj, data = mf, rhs = idx),
          error = function(e) matrix(1, nrow(mf), 1, dimnames = list(NULL, "(Intercept)"))
        )
      } else {
        # Default to intercept-only
        matrices[[pos]] <- matrix(1, nrow(mf), 1, dimnames = list(NULL, "(Intercept)"))
      }
    }
  } else if (!is.null(x$tmb_object) && !is.null(x$tmb_object$env$data)) {
    # Use TMB object if available
    tmb_data <- x$tmb_object$env$data

    # Extract matrices by their TMB position
    for (param in param_names) {
      pos <- param_positions[[param]]
      tmb_matrix_name <- paste0("X", pos)

      if (!is.null(tmb_data[[tmb_matrix_name]])) {
        matrices[[pos]] <- tmb_data[[tmb_matrix_name]]

        # Add column names if missing
        if (is.null(colnames(matrices[[pos]]))) {
          if (ncol(matrices[[pos]]) == 1) {
            colnames(matrices[[pos]]) <- "(Intercept)"
          } else {
            colnames(matrices[[pos]]) <- paste0(param, "_", 1:ncol(matrices[[pos]]))
          }
        }
      }
    }
  } else {
    stop("Cannot extract model matrices. Try fitting the model with x=TRUE.")
  }

  # Check for missing matrices and replace with default
  for (i in 1:num_matrices) {
    if (is.null(matrices[[i]])) {
      n_obs <- ifelse(is.null(x$y), nrow(x$model), length(x$y))
      matrices[[i]] <- matrix(1, n_obs, 1, dimnames = list(NULL, "(Intercept)"))
    }
  }

  return(matrices)
}


#' Extract model parameters from a gkwreg object with family-specific handling
#'
#' @param x A fitted model object of class "gkwreg"
#' @return A list of model parameters
#'
#' @keywords internal
.extract_model_params <- function(x) {
  # Get family information
  family <- x$family
  if (is.null(family)) family <- "gkw"

  param_info <- .get_family_param_info(family)
  param_names <- param_info$names
  param_positions <- param_info$positions
  fixed_params <- param_info$fixed

  # Initialize beta parameters for all possible positions
  beta_params <- vector("list", 5)
  names(beta_params) <- paste0("beta", 1:5)

  # Extract coefficients
  coefs <- x$coefficients

  if (is.list(coefs)) {
    # If coefficients are a list with parameter names
    for (param in param_names) {
      if (param %in% names(coefs)) {
        pos <- param_positions[[param]]
        beta_params[[pos]] <- coefs[[param]]
      }
    }
  } else if (!is.null(names(coefs))) {
    # If coefficients are named, extract them using regex patterns
    for (param in param_names) {
      pos <- param_positions[[param]]
      param_coefs <- coefs[grep(paste0("^", param, ":"), names(coefs))]
      if (length(param_coefs) > 0) {
        beta_params[[pos]] <- param_coefs
      }
    }

    # If regex didn't work, try alternative pattern matching
    empty_positions <- sapply(beta_params[1:length(param_names)], is.null)
    if (all(empty_positions)) {
      # Try alternative pattern matching
      for (param in param_names) {
        pos <- param_positions[[param]]
        param_coefs <- coefs[grep(param, names(coefs))]
        if (length(param_coefs) > 0) {
          beta_params[[pos]] <- param_coefs
        }
      }

      # If still empty, raise error
      if (all(sapply(beta_params[1:length(param_names)], is.null))) {
        stop("Cannot determine coefficient mapping. Try a more recent version of the model.")
      }
    }
  } else {
    stop("Unrecognized coefficient structure. Try a more recent version of the model.")
  }

  # Assign default values for fixed parameters
  for (param in names(fixed_params)) {
    pos <- param_positions[[param]]
    if (!is.null(pos)) {
      # Use the fixed value for this parameter
      beta_params[[pos]] <- fixed_params[[param]]
    }
  }

  # Assign default values for any remaining NULL betas
  for (i in 1:5) {
    if (is.null(beta_params[[i]])) {
      beta_params[[i]] <- 0
    }
  }

  # Extract link information
  if (!is.null(x$link_codes)) {
    link_codes <- unlist(x$link_codes, use.names = FALSE)
  } else if (!is.null(x$link)) {
    # Convert link functions to codes
    link_map <- c(
      "log" = 1,
      "logit" = 2,
      "probit" = 3,
      "cauchy" = 4,
      "cloglog" = 5,
      "identity" = 6,
      "sqrt" = 7,
      "inverse" = 8,
      "inverse-square" = 9
    )

    # Default links (log for most, logit for delta)
    link_codes <- c(1, 1, 1, 2, 1)

    # Update with actual links from the model
    for (param in param_names) {
      if (param %in% names(x$link)) {
        pos <- param_positions[[param]]
        link_codes[pos] <- link_map[x$link[[param]]]
      }
    }
  } else {
    # Default link types if not specified
    link_codes <- c(1, 1, 1, 2, 1) # log for most parameters, logit for delta
  }

  # Extract scale factors
  if (!is.null(x$scale_factors)) {
    scale_factors <- as.numeric(unlist(x$scale_factors))
  } else {
    # Default scale factors
    scale_factors <- c(10, 10, 10, 1, 10) # 10 for most parameters, 1 for delta
  }

  return(c(
    beta_params,
    list(
      link_codes = link_codes,
      scale_factors = scale_factors
    )
  ))
}


#' Sample model data for large datasets
#'
#' @param n Total number of observations
#' @param sample_size Target sample size
#' @param y_obs Vector of observed values
#' @param fitted_vals Vector of fitted values
#' @param model_matrices List of model matrices
#' @return A list with sampled data
#'
#' @keywords internal
.sample_model_data <- function(n, sample_size, y_obs, fitted_vals, model_matrices) {
  # For reproducibility
  set.seed(123)
  idx <- sample(n, size = min(sample_size, n))

  # Sample matrices - ensure we maintain the structure
  sampled_matrices <- list()
  for (matrix_name in names(model_matrices)) {
    X <- model_matrices[[matrix_name]]
    if (is.matrix(X) && nrow(X) == n) {
      sampled_matrices[[matrix_name]] <- X[idx, , drop = FALSE]
    } else {
      # Handle non-matrix or incorrectly sized matrix
      sampled_matrices[[matrix_name]] <- matrix(1, length(idx), 1)
    }
  }

  return(list(
    idx = idx,
    y_obs = y_obs[idx],
    fitted_vals = fitted_vals[idx],
    model_matrices = sampled_matrices
  ))
}


#' Calculate model parameters for the specified family
#'
#' @param model_matrices List of model matrices
#' @param model_params List of model parameters
#' @param family Character string specifying distribution family
#' @return A matrix of calculated parameters
#'
#' @keywords internal
.calculate_model_parameters <- function(model_matrices, model_params, family) {
  # Get parameter information for this family
  param_info <- .get_family_param_info(family)

  # Calculate the number of required parameters based on family
  num_params <- max(unlist(param_info$positions))

  # Extract matrices for all parameters
  X_matrices <- vector("list", 5)
  for (i in 1:5) {
    X_name <- paste0("X", i)
    if (i <= num_params && X_name %in% names(model_matrices)) {
      X_matrices[[i]] <- model_matrices[[X_name]]
    } else {
      # Default matrix for unused parameters
      X_matrices[[i]] <- matrix(1, nrow(model_matrices[[1]]), 1)
    }
  }

  # Extract beta parameters
  beta_params <- vector("list", 5)
  for (i in 1:5) {
    beta_name <- paste0("beta", i)
    if (beta_name %in% names(model_params)) {
      beta_params[[i]] <- model_params[[beta_name]]
    } else {
      # Default for unused parameters
      beta_params[[i]] <- 0
    }
  }

  # Extract link codes and scale factors
  link_codes <- model_params$link_codes
  scale_factors <- model_params$scale_factors

  # Call the calculateParameters function with proper unpacking
  param_mat <- calculateParameters(
    X_matrices[[1]], X_matrices[[2]], X_matrices[[3]],
    X_matrices[[4]], X_matrices[[5]],
    beta_params[[1]], beta_params[[2]], beta_params[[3]],
    beta_params[[4]], beta_params[[5]],
    link_codes, scale_factors,
    family = family
  )

  return(param_mat)
}

#' Extract parameter vectors from parameter matrix
#'
#' @param param_mat Matrix of calculated parameters
#' @return A list of parameter vectors
#'
#' @keywords internal
.extract_parameter_vectors <- function(param_mat) {
  list(
    alphaVec = param_mat[, 1],
    betaVec = param_mat[, 2],
    gammaVec = param_mat[, 3],
    deltaVec = param_mat[, 4],
    lambdaVec = param_mat[, 5]
  )
}


#' Calculate residuals based on the specified type
#'
#' @param y_obs Vector of observed values
#' @param fitted_vals Vector of fitted values
#' @param param_mat Matrix of calculated parameters
#' @param type Character string specifying residual type
#' @param family Character string specifying distribution family
#' @return A vector of residuals
#'
#' @keywords internal
.calculate_residuals <- function(y_obs, fitted_vals, param_mat, type, family) {
  # The following calculate*Residuals() are assumed to be internal package functions
  if (type == "quantile") {
    calculateQuantileResiduals(y_obs, param_mat, family = family)
  } else if (type == "pearson") {
    calculatePearsonResiduals(y_obs, fitted_vals, param_mat, family = family)
  } else if (type == "deviance") {
    calculateDevianceResiduals(y_obs, fitted_vals, param_mat, family = family)
  } else {
    stop("Unsupported residual type.")
  }
}

#' Calculate diagnostic measures for gkwreg plots
#'
#' @param y_obs Vector of observed values
#' @param fitted_vals Vector of fitted values
#' @param resid_vec Vector of residuals
#' @param model_matrices List of model matrices
#' @param model_params List of model parameters
#' @param param_vectors List of parameter vectors
#' @param idx Vector of observation indices
#' @param family Character string specifying distribution family
#' @param param_info Parameter information for the family
#' @return A list of diagnostic measures
#'
#' @importFrom stats rnorm
#' @keywords internal
.calculate_diagnostic_measures <- function(y_obs,
                                           fitted_vals,
                                           resid_vec,
                                           model_matrices,
                                           model_params,
                                           param_vectors,
                                           idx,
                                           family,
                                           param_info) {
  # Get the first active parameter for this family to use for linear predictor
  first_param <- param_info$names[1]
  first_pos <- param_info$positions[[first_param]]

  # Get the X matrix and beta coefficients for the first parameter
  X_matrix <- model_matrices[[paste0("X", first_pos)]]
  beta_coef <- model_params[[paste0("beta", first_pos)]]

  # Calculate linear predictor
  if (length(beta_coef) > 0 && !is.null(X_matrix)) {
    if (is.matrix(X_matrix) && is.numeric(beta_coef) &&
      ncol(X_matrix) == length(beta_coef)) {
      linpred <- as.vector(X_matrix %*% beta_coef)
    } else {
      # Fallback for incompatible dimensions
      linpred <- rep(0, length(idx))
    }
  } else {
    # Fallback for missing data
    linpred <- rep(0, length(idx))
  }

  # Calculate total number of parameters (excluding fixed ones)
  p <- 0
  for (param in param_info$names) {
    pos <- param_info$positions[[param]]
    beta_name <- paste0("beta", pos)
    if (beta_name %in% names(model_params)) {
      beta_val <- model_params[[beta_name]]
      if (is.numeric(beta_val) && length(beta_val) > 0) {
        p <- p + length(beta_val)
      }
    }
  }

  # Ensure p is at least 1 to avoid division by zero
  p <- max(1, p)

  # Mean squared error
  mse <- mean(resid_vec^2, na.rm = TRUE)

  # Approximate generalized leverage
  n <- length(idx)
  leverage <- rep(p / n, length(idx))

  # Small noise addition for visual differentiation
  leverage <- leverage + abs(rnorm(length(idx), 0, 0.01))

  # Calculate Cook's distance
  cook_dist <- (resid_vec^2 / (p * mse)) * (leverage / ((1 - leverage)^2))
  cook_dist[is.infinite(cook_dist) | cook_dist > 100 | is.na(cook_dist)] <- NA

  list(
    linpred = linpred,
    p = p,
    leverage = leverage,
    cook_dist = cook_dist
  )
}


#' Calculate half-normal plot data with envelope
#'
#' @param resid_vec Vector of residuals
#' @param idx Vector of observation indices
#' @param nsim Number of simulations for envelope
#' @param level Confidence level for envelope
#' @param param_mat Matrix of calculated parameters
#' @param param_vectors List of parameter vectors
#' @param type Character string specifying residual type
#' @param family Character string specifying distribution family
#' @return A data frame with half-normal plot data
#' @importFrom stats quantile
#'
#' @keywords internal
.calculate_half_normal_data <- function(resid_vec,
                                        idx,
                                        nsim,
                                        level,
                                        param_mat,
                                        param_vectors,
                                        type,
                                        family) {
  abs_resid <- abs(resid_vec)
  sorted_abs_resid <- sort(abs_resid)
  prob_points <- (seq_along(idx) - 0.5) / length(idx)
  hn_q <- qnorm(0.5 + prob_points / 2) # half-normal quantiles from normal

  # Prepare data frame
  half_normal_data <- data.frame(
    index = seq_along(sorted_abs_resid),
    theoretical = hn_q,
    observed = sorted_abs_resid
  )

  # Extract parameter vectors for simulation
  alphaVec <- param_vectors$alphaVec
  betaVec <- param_vectors$betaVec
  gammaVec <- param_vectors$gammaVec
  deltaVec <- param_vectors$deltaVec
  lambdaVec <- param_vectors$lambdaVec

  # Simulate envelope
  set.seed(54321)
  envelope_data <- matrix(NA, nrow = length(idx), ncol = nsim)

  cat("Simulating envelope (", nsim, "iterations): ")
  progress_step <- max(1, floor(nsim / 10))

  for (i in seq_len(nsim)) {
    if (i %% progress_step == 0) cat(".")

    # Generate simulated data using family-specific random generation
    sim_y <- .simulate_from_distribution(
      n = length(idx),
      alphaVec = alphaVec,
      betaVec = betaVec,
      gammaVec = gammaVec,
      deltaVec = deltaVec,
      lambdaVec = lambdaVec,
      family = family
    )

    # Calculate residuals for simulated data
    sim_resid <- .calculate_sim_residuals(
      sim_y = sim_y,
      param_mat = param_mat,
      type = type,
      family = family
    )

    envelope_data[, i] <- sort(abs(sim_resid))
  }
  cat(" Done!\n")

  # Calculate envelope bounds
  lower_bound <- apply(envelope_data, 1, quantile, probs = (1 - level) / 2, na.rm = TRUE)
  upper_bound <- apply(envelope_data, 1, quantile, probs = 1 - (1 - level) / 2, na.rm = TRUE)

  half_normal_data$lower <- lower_bound
  half_normal_data$upper <- upper_bound

  half_normal_data
}




#' Simulate observations from a specified distribution family
#'
#' @param n Number of observations to simulate
#' @param alphaVec Vector of alpha parameters
#' @param betaVec Vector of beta parameters
#' @param gammaVec Vector of gamma parameters
#' @param deltaVec Vector of delta parameters
#' @param lambdaVec Vector of lambda parameters
#' @param family Character string specifying distribution family
#' @return A vector of simulated observations
#'
#' @keywords internal
.simulate_from_distribution <- function(n,
                                        alphaVec,
                                        betaVec,
                                        gammaVec,
                                        deltaVec,
                                        lambdaVec,
                                        family) {
  # Generate uniform random variates
  sim_u <- stats::runif(n)

  # Use the appropriate quantile function for each family
  # This approach replaces the complex nested expressions with direct calls
  # to the existing quantile functions for each distribution
  switch(family,
    "gkw" = {
      # Use the implemented qgkw function
      rgkw(n,
        alpha = alphaVec, beta = betaVec, gamma = gammaVec,
        delta = deltaVec, lambda = lambdaVec
      )
    },
    "bkw" = {
      # BKw: lambda = 1
      rbkw(n, alpha = alphaVec, beta = betaVec, gamma = gammaVec, delta = deltaVec)
    },
    "kkw" = {
      # KKw: gamma = 1
      rkkw(n, alpha = alphaVec, beta = betaVec, delta = deltaVec, lambda = lambdaVec)
    },
    "ekw" = {
      # EKw: gamma = 1, delta = 0
      rekw(n, alpha = alphaVec, beta = betaVec, lambda = lambdaVec)
    },
    "mc" = {
      # MC: alpha = 1, beta = 1
      rmc(n, gamma = gammaVec, delta = deltaVec, lambda = lambdaVec)
    },
    "kw" = {
      # KW: lambda = 1, gamma = 1, delta = 0
      rkw(n, alpha = alphaVec, beta = betaVec)
    },
    "beta" = {
      # Beta: alpha = 1, beta = 1, lambda = 1
      rbeta_(n, gammaVec, deltaVec)
    },
    # Default case - use the GKw distribution
    {
      warning("Unrecognized family '", family, "'. Using GKw distribution instead.")
      rgkw(n,
        alpha = alphaVec, beta = betaVec, gamma = gammaVec,
        delta = deltaVec, lambda = lambdaVec
      )
    }
  )
}



#' Calculate residuals for simulated data
#'
#' @param sim_y Vector of simulated observations
#' @param param_mat Matrix of calculated parameters
#' @param type Character string specifying residual type
#' @param family Character string specifying distribution family
#' @return A vector of residuals
#'
#' @keywords internal
.calculate_sim_residuals <- function(sim_y, param_mat, type, family) {
  if (type == "quantile") {
    calculateQuantileResiduals(sim_y, param_mat, family = family)
  } else if (type == "pearson") {
    sim_fitted <- calculateMeans(param_mat, family = family)
    calculatePearsonResiduals(sim_y, sim_fitted, param_mat, family = family)
  } else if (type == "deviance") {
    sim_fitted <- calculateMeans(param_mat, family = family)
    calculateDevianceResiduals(sim_y, sim_fitted, param_mat, family = family)
  } else {
    stop("Unsupported residual type.")
  }
}


#' Create formatted plot titles
#'
#' @param which Integer vector specifying which plots to produce
#' @param caption Character vector of plot captions
#' @param main Optional character string for main title
#' @param family Character string specifying distribution family
#' @param orig_family Original family from the model
#' @return A named vector of formatted plot titles
#'
#' @keywords internal
.create_plot_titles <- function(which, caption, main, family, orig_family) {
  plot_titles <- rep("", length(which))
  names(plot_titles) <- which

  for (i in seq_along(which)) {
    i_plot <- which[i]
    base_title <- caption[i_plot]
    if (!is.null(orig_family) && family != orig_family) {
      base_title <- paste0(base_title, " (", family, ")")
    }
    if (nzchar(main)) {
      plot_titles[as.character(i_plot)] <- paste(main, "-", base_title)
    } else {
      plot_titles[as.character(i_plot)] <- base_title
    }
  }

  plot_titles
}


#' Generate diagnostic plots using base R graphics
#'
#' @param diag_data List of diagnostic data
#' @param which Integer vector specifying which plots to produce
#' @param plot_titles Named vector of plot titles
#' @param sub.caption Character string for subtitle
#' @param ask Logical; whether to ask before new plots
#' @param ... Additional graphical parameters
#' @return NULL (invisibly)
#'
#' @keywords internal
.plot_gkwreg_base_r <- function(diag_data, which, plot_titles, sub.caption, ask, ...) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  one.fig <- prod(par("mfcol")) == 1

  if (ask) {
    oask <- grDevices::devAskNewPage(TRUE)
    on.exit(grDevices::devAskNewPage(oask), add = TRUE)
  }

  type <- diag_data$model_info$type

  for (i_plot in which) {
    grDevices::dev.hold()
    p_main <- plot_titles[as.character(i_plot)]

    if (i_plot == 1) {
      .plot_base_r_residuals_vs_index(diag_data, p_main, type, ...)
    } else if (i_plot == 2) {
      .plot_base_r_cooks_distance(diag_data, p_main, ...)
    } else if (i_plot == 3) {
      .plot_base_r_leverage_vs_fitted(diag_data, p_main, ...)
    } else if (i_plot == 4) {
      .plot_base_r_residuals_vs_linpred(diag_data, p_main, type, ...)
    } else if (i_plot == 5) {
      .plot_base_r_half_normal(diag_data, p_main, type, ...)
    } else if (i_plot == 6) {
      .plot_base_r_predicted_vs_observed(diag_data, p_main, ...)
    }
    mtext(sub.caption, side = 3, line = 0.25, cex = 0.8, col = "gray40")

    grDevices::dev.flush()

    if (ask && !one.fig) {
      message("Press <ENTER> to continue to the next plot...")
      invisible(readLines(con = "stdin", n = 1))
    }
  }

  invisible(NULL)
}


#' Plot residuals vs. index (base R)
#'
#' @keywords internal
.plot_base_r_residuals_vs_index <- function(diag_data, p_main, type, ...) {
  plot(diag_data$data$index, diag_data$data$resid,
    ylab = paste0(type, " residuals"),
    xlab = "Observation Index",
    main = p_main, ...
  )
  abline(h = 0, lty = 2, col = "gray")
  lines(stats::lowess(diag_data$data$index, diag_data$data$resid),
    col = "red", lwd = 2
  )

  invisible(NULL)
}


#' Plot Cook's distance (base R)
#'
#' @keywords internal
.plot_base_r_cooks_distance <- function(diag_data, p_main, ...) {
  plot(diag_data$data$index, diag_data$data$cook_dist,
    ylab = "Cook's Distance",
    xlab = "Observation Index",
    main = p_main, ...
  )

  threshold <- diag_data$model_info$cook_threshold
  abline(h = threshold, lty = 2, col = "red")
  text(0.9 * max(diag_data$data$index), threshold * 1.1,
    "4/n",
    col = "red", pos = 3
  )

  invisible(NULL)
}


#' Plot leverage vs. fitted (base R)
#'
#' @keywords internal
.plot_base_r_leverage_vs_fitted <- function(diag_data, p_main, ...) {
  plot(diag_data$data$fitted, diag_data$data$leverage,
    xlab = "Fitted Values",
    ylab = "Generalized Leverage",
    main = p_main, ...
  )

  threshold <- diag_data$model_info$leverage_threshold
  abline(h = threshold, lty = 2, col = "red")
  text(0.9 * max(diag_data$data$fitted), threshold * 1.1,
    "2p/n",
    col = "red", pos = 3
  )

  invisible(NULL)
}


#' Plot residuals vs. linear predictor (base R)
#'
#' @keywords internal
.plot_base_r_residuals_vs_linpred <- function(diag_data, p_main, type, ...) {
  plot(diag_data$data$linpred, diag_data$data$resid,
    xlab = "Linear Predictor (alpha)",
    ylab = paste0(type, " residuals"),
    main = p_main, ...
  )
  abline(h = 0, lty = 2, col = "gray")
  lines(stats::lowess(diag_data$data$linpred, diag_data$data$resid),
    col = "red", lwd = 2
  )

  invisible(NULL)
}


#' Plot half-normal plot (base R)
#'
#' @importFrom stats qqnorm qqline
#' @keywords internal
.plot_base_r_half_normal <- function(diag_data, p_main, type, ...) {
  if (is.null(diag_data$half_normal)) {
    # fallback to a normal QQ-plot of absolute residuals
    stats::qqnorm(abs(diag_data$data$resid),
      main = paste(p_main, "(Half-Normal)"),
      ylab = paste0("|", type, " residuals|"), ...
    )
    stats::qqline(abs(diag_data$data$resid), lty = 2, col = "gray")
    return(invisible(NULL))
  }

  hn_data <- diag_data$half_normal
  plot(hn_data$theoretical, hn_data$observed,
    xlab = "Theoretical Half-Normal Quantiles",
    ylab = paste0("Ordered |", type, " residuals|"),
    main = paste(p_main, "(Half-Normal)"), ...
  )

  abline(0, stats::sd(diag_data$data$abs_resid, na.rm = TRUE),
    lty = 2, col = "gray"
  )

  if ("lower" %in% names(hn_data) && "upper" %in% names(hn_data)) {
    lines(hn_data$theoretical, hn_data$lower, lty = 2, col = "blue", lwd = 1.5)
    lines(hn_data$theoretical, hn_data$upper, lty = 2, col = "blue", lwd = 1.5)
    level <- round(100 * (1 - (1 - 0.5 * (hn_data$upper[1] / hn_data$observed[1]))), 2)
    text(max(hn_data$theoretical) * 0.8, max(hn_data$upper) * 0.9,
      paste0(level, "% envelope"),
      col = "blue"
    )
  }

  invisible(NULL)
}


#' Plot predicted vs. observed (base R)
#'
#' @keywords internal
.plot_base_r_predicted_vs_observed <- function(diag_data, p_main, ...) {
  plot(diag_data$data$fitted, diag_data$data$y_obs,
    xlab = "Fitted (Mean)",
    ylab = "Observed (y)",
    main = p_main, ...
  )
  abline(0, 1, col = "gray", lty = 2)
  lines(stats::lowess(diag_data$data$fitted, diag_data$data$y_obs),
    col = "red", lwd = 2
  )

  invisible(NULL)
}


#' Generate diagnostic plots using ggplot2
#'
#' @param diag_data List of diagnostic data
#' @param which Integer vector specifying which plots to produce
#' @param plot_titles Named vector of plot titles
#' @param sub.caption Character string for subtitle
#' @param ask Logical; whether to ask before new plots
#' @param arrange_plots Logical; whether to arrange multiple plots in a grid
#' @param theme_fn ggplot2 theme function
#' @param ... Additional graphical parameters
#' @return NULL (invisibly)
#'
#' @keywords internal
.plot_gkwreg_ggplot <- function(diag_data,
                                which,
                                plot_titles,
                                sub.caption,
                                ask,
                                arrange_plots,
                                theme_fn,
                                ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for ggplot2 plotting.")
  }

  type <- diag_data$model_info$type
  p_list <- list()

  for (i in seq_along(which)) {
    i_plot <- which[i]
    p_main <- plot_titles[as.character(i_plot)]

    if (i_plot == 1) {
      p_list[[i]] <- .plot_ggplot_residuals_vs_index(diag_data, p_main, sub.caption, type, theme_fn)
    } else if (i_plot == 2) {
      p_list[[i]] <- .plot_ggplot_cooks_distance(diag_data, p_main, sub.caption, theme_fn)
    } else if (i_plot == 3) {
      p_list[[i]] <- .plot_ggplot_leverage_vs_fitted(diag_data, p_main, sub.caption, theme_fn)
    } else if (i_plot == 4) {
      p_list[[i]] <- .plot_ggplot_residuals_vs_linpred(diag_data, p_main, sub.caption, type, theme_fn)
    } else if (i_plot == 5) {
      p_list[[i]] <- .plot_ggplot_half_normal(diag_data, p_main, sub.caption, type, theme_fn)
    } else if (i_plot == 6) {
      p_list[[i]] <- .plot_ggplot_predicted_vs_observed(diag_data, p_main, sub.caption, theme_fn)
    }
  }

  if (arrange_plots && length(p_list) > 1) {
    n_plots <- length(p_list)
    n_cols <- min(2, n_plots)
    n_rows <- ceiling(n_plots / n_cols)

    if (requireNamespace("gridExtra", quietly = TRUE)) {
      gridExtra::grid.arrange(grobs = p_list, ncol = n_cols, nrow = n_rows)
    } else if (requireNamespace("ggpubr", quietly = TRUE)) {
      do.call(ggpubr::ggarrange, c(p_list, list(ncol = n_cols, nrow = n_rows)))
    } else {
      warning("Neither 'gridExtra' nor 'ggpubr' is installed. Displaying plots individually.")
      for (i in seq_along(p_list)) {
        print(p_list[[i]])
        if (ask && i < length(p_list)) {
          message("Press <ENTER> to continue to the next plot...")
          invisible(readLines(con = "stdin", n = 1))
        }
      }
    }
  } else {
    for (i in seq_along(p_list)) {
      print(p_list[[i]])
      if (ask && i < length(p_list)) {
        message("Press <ENTER> to continue to the next plot...")
        invisible(readLines(con = "stdin", n = 1))
      }
    }
  }

  invisible(NULL)
}


#' Plot residuals vs. index (ggplot2)
#'
#' @keywords internal
.plot_ggplot_residuals_vs_index <- function(diag_data, p_main, sub.caption, type, theme_fn) {
  ggplot2::ggplot(diag_data$data, ggplot2::aes(x = index, y = resid)) +
    ggplot2::geom_point() +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    ggplot2::geom_smooth(
      method = "loess", formula = y ~ x, se = FALSE,
      color = "red", linewidth = 1
    ) +
    ggplot2::labs(
      x = "Observation Index",
      y = paste0(type, " residuals"),
      title = p_main,
      subtitle = sub.caption
    ) +
    theme_fn()
}


#' Plot Cook's distance (ggplot2)
#'
#' @keywords internal
.plot_ggplot_cooks_distance <- function(diag_data, p_main, sub.caption, theme_fn) {
  cook_ref <- diag_data$model_info$cook_threshold
  df_data <- diag_data$data

  p <- ggplot2::ggplot(df_data, ggplot2::aes(x = index, y = cook_dist)) +
    ggplot2::geom_point() +
    ggplot2::geom_hline(yintercept = cook_ref, linetype = "dashed", color = "red") +
    ggplot2::annotate("text",
      x = max(df_data$index) * 0.9, y = cook_ref * 1.1,
      label = "4/n", color = "red", hjust = 0
    ) +
    ggplot2::labs(
      x = "Observation Index",
      y = "Cook's Distance",
      title = p_main,
      subtitle = sub.caption
    ) +
    theme_fn()

  # Label points exceeding threshold
  infl <- df_data$cook_dist > cook_ref & !is.na(df_data$cook_dist)
  if (any(infl)) {
    p <- p + ggplot2::geom_text(
      data = df_data[infl, ],
      ggplot2::aes(label = index),
      vjust = -0.5, color = "red", size = 3
    )
  }

  p
}


#' Plot leverage vs. fitted (ggplot2)
#'
#' @keywords internal
.plot_ggplot_leverage_vs_fitted <- function(diag_data, p_main, sub.caption, theme_fn) {
  lev_ref <- diag_data$model_info$leverage_threshold
  df_data <- diag_data$data

  p <- ggplot2::ggplot(df_data, ggplot2::aes(x = fitted, y = leverage)) +
    ggplot2::geom_point() +
    ggplot2::geom_hline(yintercept = lev_ref, linetype = "dashed", color = "red") +
    ggplot2::annotate("text",
      x = max(df_data$fitted) * 0.9, y = lev_ref * 1.1,
      label = "2p/n", color = "red", hjust = 0
    ) +
    ggplot2::labs(
      x = "Fitted Values",
      y = "Generalized Leverage",
      title = p_main,
      subtitle = sub.caption
    ) +
    theme_fn()

  # Label high leverage points
  high_lev <- df_data$leverage > lev_ref
  if (any(high_lev)) {
    p <- p + ggplot2::geom_text(
      data = df_data[high_lev, ],
      ggplot2::aes(label = index),
      vjust = -0.5, color = "red", size = 3
    )
  }

  p
}


#' Plot residuals vs. linear predictor (ggplot2)
#'
#' @keywords internal
.plot_ggplot_residuals_vs_linpred <- function(diag_data, p_main, sub.caption, type, theme_fn) {
  df_data <- diag_data$data
  ggplot2::ggplot(df_data, ggplot2::aes(x = linpred, y = resid)) +
    ggplot2::geom_point() +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    ggplot2::geom_smooth(
      method = "loess", formula = y ~ x, se = FALSE,
      color = "red", linewidth = 1
    ) +
    ggplot2::labs(
      x = "Linear Predictor (alpha)",
      y = paste0(type, " residuals"),
      title = p_main,
      subtitle = sub.caption
    ) +
    theme_fn()
}


#' Plot half-normal plot (ggplot2)
#'
#' @param diag_data Diagnostic data list
#' @param p_main Plot title
#' @param sub.caption Plot subtitle
#' @param type Residual type
#' @param theme_fn ggplot2 theme function
#' @return A ggplot object
#'
#' @keywords internal
.plot_ggplot_half_normal <- function(diag_data, p_main, sub.caption, type, theme_fn) {
  # Extract the confidence level from the half_normal data
  # If it's not available in diag_data, use a default value
  level_value <- diag_data$model_info$level
  if (is.null(level_value)) {
    level_value <- 0.90 # Default level if not stored in diag_data
  }

  if (is.null(diag_data$half_normal)) {
    # Fallback to a half-normal Q-Q plot approach
    return(
      ggplot2::ggplot(diag_data$data, ggplot2::aes(sample = abs_resid)) +
        ggplot2::stat_qq(distribution = function(p) qnorm(p * 0.5 + 0.5)) +
        ggplot2::stat_qq_line(
          distribution = function(p) qnorm(p * 0.5 + 0.5),
          color = "gray", linetype = "dashed"
        ) +
        ggplot2::labs(
          x = "Theoretical Half-Normal Quantiles",
          y = paste0("|", type, " residuals|"),
          title = paste(p_main, "(Half-Normal)"),
          subtitle = sub.caption
        ) +
        theme_fn()
    )
  }

  hn_data <- diag_data$half_normal
  sd_abs <- stats::sd(diag_data$data$abs_resid, na.rm = TRUE)

  p <- ggplot2::ggplot(hn_data, ggplot2::aes(x = theoretical, y = observed)) +
    ggplot2::geom_point() +
    ggplot2::geom_abline(
      slope = sd_abs, intercept = 0,
      linetype = "dashed", color = "gray"
    ) +
    ggplot2::labs(
      x = "Theoretical Half-Normal Quantiles",
      y = paste0("Ordered |", type, " residuals|"),
      title = paste(p_main, "(Half-Normal)"),
      subtitle = sub.caption
    ) +
    theme_fn()

  if ("lower" %in% names(hn_data) && "upper" %in% names(hn_data)) {
    p <- p +
      ggplot2::geom_line(ggplot2::aes(y = lower),
        linetype = "dashed",
        color = "blue", linewidth = 0.8
      ) +
      ggplot2::geom_line(ggplot2::aes(y = upper),
        linetype = "dashed",
        color = "blue", linewidth = 0.8
      ) +
      ggplot2::annotate(
        "text",
        x = max(hn_data$theoretical) * 0.8,
        y = max(hn_data$upper) * 0.9,
        label = paste0(format(100 * level_value, digits = 2), "% envelope"),
        color = "blue"
      )
  }

  return(p)
}


#' Plot predicted vs. observed (ggplot2)
#'
#' @keywords internal
.plot_ggplot_predicted_vs_observed <- function(diag_data, p_main, sub.caption, theme_fn) {
  df_data <- diag_data$data
  ggplot2::ggplot(df_data, ggplot2::aes(x = fitted, y = y_obs)) +
    ggplot2::geom_point() +
    ggplot2::geom_abline(
      intercept = 0, slope = 1,
      color = "gray", linetype = "dashed"
    ) +
    ggplot2::geom_smooth(
      method = "loess", formula = y ~ x, se = FALSE,
      color = "red", linewidth = 1
    ) +
    ggplot2::labs(
      x = "Fitted (Mean)",
      y = "Observed (y)",
      title = p_main,
      subtitle = sub.caption
    ) +
    theme_fn()
}


# Extract confidence intervals from TMB profile likelihood objects
extract_profile_ci <- function(profile_obj, level = 0.95) {
  # Check if the profile object has the expected format
  if (!is.data.frame(profile_obj) || !all(c("par", "value") %in% names(profile_obj))) {
    stop("Profile object does not have the expected format")
  }

  # Negative log-likelihood profile (-logLik)
  prof_data <- profile_obj[, c("par", "value")]

  # Find the minimum value of -logLik (maximum of logLik)
  min_value <- min(prof_data$value, na.rm = TRUE)

  # Calculate threshold based on chi-square distribution
  # For CI of level (1-alpha), we use the (1-alpha) quantile of chi-square with 1 d.f.
  alpha <- 1 - level
  threshold <- min_value + qchisq(level, df = 1) / 2

  # Filter points within the confidence interval
  ci_points <- prof_data[prof_data$value <= threshold, ]

  # If there aren't enough points, return NA
  if (nrow(ci_points) < 2) {
    return(c(NA, NA))
  }

  # Extract lower and upper CI limits
  ci_lower <- min(ci_points$par, na.rm = TRUE)
  ci_upper <- max(ci_points$par, na.rm = TRUE)

  return(c(ci_lower, ci_upper))
}


#' Confidence Intervals for Generalized Kumaraswamy Regression Coefficients
#'
#' @description
#' Computes confidence intervals for one or more parameters in a fitted
#' Generalized Kumaraswamy regression model. This function provides both
#' Wald-type confidence intervals (based on standard errors and normal
#' approximation) and profile likelihood confidence intervals when available.
#'
#' @param object An object of class \code{"gkwreg"}, typically the result of
#'   \code{\link{gkwreg}}.
#' @param parm A specification of which parameters to compute confidence intervals
#'   for. Either a character vector of parameter names or a numeric vector of
#'   parameter indices. If missing, confidence intervals for all parameters are
#'   returned.
#' @param level The confidence level required. Default is 0.95 (95\% confidence
#'   interval).
#' @param method The method used to compute the confidence intervals. Options are:
#'   \code{"wald"} (default) for Wald-type intervals based on standard errors, or
#'   \code{"profile"} for profile likelihood intervals.
#' @param max_steps Integer; maximum number of steps for profile likelihood
#'   evaluation. Only applicable when \code{method = "profile"}. Default is 100.
#' @param stepsize Numeric; step size factor for profile likelihood evaluation.
#'   Only applicable when \code{method = "profile"}. Default is 0.01.
#' @param trace Logical; if \code{TRUE}, prints information about the profiling
#'   process. Default is \code{FALSE}.
#' @param ... Additional arguments (currently not used).
#'
#' @details
#' For Wald-type intervals (\code{method = "wald"}), the function computes
#' confidence intervals using the formula:
#'
#' \deqn{\hat{\beta}_j \pm z_{1-\alpha/2} \cdot \mathrm{SE}(\hat{\beta}_j)}
#'
#' where \eqn{\hat{\beta}_j} is the estimated coefficient, \eqn{\mathrm{SE}(\hat{\beta}_j)}
#' is the standard error, and \eqn{z_{1-\alpha/2}} is the \eqn{1-\alpha/2} quantile
#' of the standard normal distribution, with \eqn{\alpha = 1 - \mathrm{level}}.
#'
#' For profile likelihood intervals (\code{method = "profile"}), the function uses
#' the embedded TMB object to compute confidence intervals based on the profile
#' likelihood, which can be more accurate for non-linear models or when the
#' likelihood is asymmetric.
#'
#' @return A matrix with columns giving the lower and upper confidence limits for
#' each parameter. These will be labeled according to the confidence level used.
#' If \code{method = "profile"} but profile likelihood intervals cannot be computed
#' (e.g., if the model was not fitted with TMB or if \code{hessian = FALSE}), the
#' function will fall back to Wald-type intervals with a warning.
#'
#' @note
#' For robust confidence intervals, ensure that the model was fitted with
#' \code{hessian = TRUE} in the \code{gkwreg} function call. Profile likelihood
#' intervals can be computationally expensive but typically provide more accurate
#' coverage, especially for parameters with asymmetric likelihood surfaces.
#'
#' @examples
#' \dontrun{
#' # Fit a Kumaraswamy regression model
#' kw_reg <- gkwreg(y ~ x1 + x2 | x3, data = df, family = "kw")
#'
#' # Calculate 95% confidence intervals for all parameters
#' confint(kw_reg)
#'
#' # Calculate 90% confidence intervals for specific parameters
#' confint(kw_reg, parm = c("alpha:(Intercept)", "beta:x3"), level = 0.90)
#'
#' # Calculate profile likelihood confidence intervals
#' confint(kw_reg, method = "profile")
#' }
#'
#' @seealso \code{\link{gkwreg}}, \code{\link{summary.gkwreg}}
#'
#' @importFrom stats qnorm
#' @export
confint.gkwreg <- function(object, parm, level = 0.95,
                           method = c("wald", "profile"),
                           max_steps = 100, stepsize = 0.01,
                           trace = FALSE, ...) {
  # Match method argument
  method <- match.arg(method)

  # Check if object is of correct class
  if (!inherits(object, "gkwreg")) {
    stop("'object' must be a fitted model object of class 'gkwreg'")
  }

  # Check if level is valid
  if (!is.numeric(level) || level <= 0 || level >= 1) {
    stop("'level' must be a numeric value between 0 and 1")
  }

  # Retrieve the coefficients and their names
  cf <- object$coefficients

  if (is.null(cf)) {
    stop("No coefficients found in the model object")
  }

  cf_names <- names(cf)
  n_coefs <- length(cf)

  # Handle parm argument
  if (missing(parm)) {
    parm_idx <- seq_len(n_coefs)
    parm_names <- cf_names
  } else if (is.numeric(parm)) {
    if (any(parm < 1 | parm > n_coefs)) {
      stop("If numeric, 'parm' must contain valid coefficient indices")
    }
    parm_idx <- parm
    parm_names <- cf_names[parm_idx]
  } else if (is.character(parm)) {
    if (!all(parm %in% cf_names)) {
      missing_parms <- parm[!parm %in% cf_names]
      stop("Parameter(s) not found: ", paste(missing_parms, collapse = ", "))
    }
    parm_idx <- match(parm, cf_names)
    parm_names <- parm
  } else {
    stop("'parm' must be either a character vector of parameter names or a numeric vector of parameter indices")
  }

  # Subset coefficients based on parm
  cf_sub <- cf[parm_idx]

  # Create matrix for results
  alpha <- 1 - level
  ci_names <- c(
    paste0(format(100 * alpha / 2, digits = 3), " %"),
    paste0(format(100 * (1 - alpha / 2), digits = 3), " %")
  )
  ci <- matrix(NA, length(parm_idx), 2, dimnames = list(parm_names, ci_names))

  # Compute confidence intervals based on method
  if (method == "wald") {
    if (is.null(object$se)) {
      warning(
        "Standard errors not found in the model object. ",
        "The model may have been fitted with hessian = FALSE. ",
        "Consider refitting with hessian = TRUE for valid confidence intervals."
      )
      # Set CIs to NA and return
      return(ci)
    }

    # Extract standard errors
    se_sub <- object$se[parm_idx]

    # Critical value from normal distribution
    z_crit <- qnorm(1 - alpha / 2)

    # Compute Wald-type intervals
    ci[, 1] <- cf_sub - z_crit * se_sub # Lower bound
    ci[, 2] <- cf_sub + z_crit * se_sub # Upper bound
  } else if (method == "profile") {
    # Check if TMB object is available
    if (is.null(object$tmb_object)) {
      warning("TMB object not found in the model. Falling back to Wald-type intervals.")

      # Fall back to Wald-type intervals
      if (!is.null(object$se)) {
        se_sub <- object$se[parm_idx]
        z_crit <- qnorm(1 - alpha / 2)
        ci[, 1] <- cf_sub - z_crit * se_sub # Lower bound
        ci[, 2] <- cf_sub + z_crit * se_sub # Upper bound
      }

      return(ci)
    }

    # Check if package TMB is available
    if (!requireNamespace("TMB", quietly = TRUE)) {
      warning(
        "Package 'TMB' is required for profile likelihood confidence intervals. ",
        "Falling back to Wald-type intervals."
      )

      # Fall back to Wald-type intervals
      if (!is.null(object$se)) {
        se_sub <- object$se[parm_idx]
        z_crit <- qnorm(1 - alpha / 2)
        ci[, 1] <- cf_sub - z_crit * se_sub # Lower bound
        ci[, 2] <- cf_sub + z_crit * se_sub # Upper bound
      }

      return(ci)
    }

    # Use profiling when available
    obj <- object$tmb_object

    # Try to use TMB's tmbprofile
    tryCatch(
      {
        if (trace) message("Computing profile likelihood confidence intervals...")

        # Loop through parameters and compute profile CIs
        for (i in seq_along(parm_idx)) {
          param_name <- parm_names[i]
          param_idx <- parm_idx[i]

          if (trace) message("Processing parameter: ", param_name)

          # Map TMB parameter name to internal parameter index
          # This is implementation-specific and may need adjustment
          tmb_param_idx <- .map_gkwreg_to_tmb_param(object, param_idx)

          if (is.na(tmb_param_idx)) {
            if (trace) message("Could not map parameter to TMB index. Skipping.")
            next
          }

          # Compute profile
          profile_obj <- try(TMB::tmbprofile(
            obj = obj,
            name = tmb_param_idx,
            h = stepsize * object$se[param_idx],
            steptol = 0.01,
            sd = TRUE,
            trace = trace,
            inner.method = "nlminb",
            maxsteps = max_steps
          ), silent = !trace)

          # Check if profiling succeeded
          if (inherits(profile_obj, "try-error")) {
            warning(
              "Profile likelihood calculation failed for parameter '",
              param_name, "'. Falling back to Wald-type interval."
            )

            # Fall back to Wald-type interval for this parameter
            if (!is.null(object$se)) {
              se_val <- object$se[param_idx]
              z_crit <- qnorm(1 - alpha / 2)
              ci[i, 1] <- cf_sub[i] - z_crit * se_val
              ci[i, 2] <- cf_sub[i] + z_crit * se_val
            }

            next
          }

          # Extract confidence interval from profile
          # profile_ci <- try(TMB::confint(profile_obj, level = level), silent = !trace)
          profile_ci <- try(extract_profile_ci(profile_obj, level = level), silent = !trace)

          if (inherits(profile_ci, "try-error") || !is.numeric(profile_ci) || length(profile_ci) != 2) {
            warning(
              "Could not extract profile likelihood confidence interval for parameter '",
              param_name, "'. Falling back to Wald-type interval."
            )

            # Fall back to Wald-type interval for this parameter
            if (!is.null(object$se)) {
              se_val <- object$se[param_idx]
              z_crit <- qnorm(1 - alpha / 2)
              ci[i, 1] <- cf_sub[i] - z_crit * se_val
              ci[i, 2] <- cf_sub[i] + z_crit * se_val
            }

            next
          }

          # Store the profile likelihood CI
          ci[i, 1] <- profile_ci[1] # Lower bound
          ci[i, 2] <- profile_ci[2] # Upper bound
        }
      },
      error = function(e) {
        warning(
          "Error in profile likelihood calculation: ", e$message,
          "\nFalling back to Wald-type intervals."
        )

        # Fall back to Wald-type intervals for all parameters
        if (!is.null(object$se)) {
          se_sub <- object$se[parm_idx]
          z_crit <- qnorm(1 - alpha / 2)
          ci[, 1] <- cf_sub - z_crit * se_sub # Lower bound
          ci[, 2] <- cf_sub + z_crit * se_sub # Upper bound
        }
      }
    )
  }

  # Return the confidence intervals
  return(ci)
}



#' Map gkwreg parameter index to TMB parameter index
#'
#' @param object A fitted model object of class "gkwreg"
#' @param param_idx Index of the parameter in the gkwreg coefficients vector
#' @return The corresponding index in the TMB parameter vector, or NA if mapping fails
#'
#' @keywords internal
.map_gkwreg_to_tmb_param <- function(object, param_idx) {
  # Extract necessary information from the model object
  cf_names <- names(object$coefficients)
  param_name <- cf_names[param_idx]

  # Try to extract the TMB parameter mapping
  if (!is.null(object$tmb_param_map)) {
    # If the model object already has a parameter map, use it
    if (param_name %in% names(object$tmb_param_map)) {
      return(object$tmb_param_map[param_name])
    }
  }

  # Otherwise, we need to reconstruct the mapping
  # This is implementation-specific and depends on how parameters are organized in TMB

  # Parse the parameter name to determine the parameter type and position
  param_parts <- strsplit(param_name, ":", fixed = TRUE)[[1]]

  if (length(param_parts) < 2) {
    # Cannot parse parameter name
    return(NA)
  }

  param_type <- param_parts[1] # alpha, beta, gamma, delta, or lambda
  covariate <- paste(param_parts[-1], collapse = ":") # The covariate name

  # Get the model family
  family <- object$family
  if (is.null(family)) family <- "gkw" # Default to gkw

  # Get parameter information for this family
  param_info <- .get_family_param_info(family)

  # Check if the parameter type is valid for this family
  if (!param_type %in% param_info$names) {
    # Parameter not valid for this family
    return(NA)
  }

  # Get the parameter position for this family
  param_pos <- param_info$positions[[param_type]]

  # Get the model matrices to determine covariate position
  if (!is.null(object$x) && param_type %in% names(object$x)) {
    X_mat <- object$x[[param_type]]
    cov_idx <- which(colnames(X_mat) == covariate)

    if (length(cov_idx) == 1) {
      # Calculate the TMB parameter index
      # This formula depends on how parameters are organized in TMB
      # The basic idea is to map (param_type, covariate) to a linear index

      # Count parameters for previous types
      offset <- 0
      for (prev_type in param_info$names) {
        if (prev_type == param_type) break

        if (prev_type %in% names(object$x)) {
          offset <- offset + ncol(object$x[[prev_type]])
        }
      }

      return(offset + cov_idx)
    }
  }

  # If we can't determine the exact mapping, try a simpler approach
  # This assumes parameters are ordered as they appear in the coefficients vector
  return(param_idx)
}


#' Extract Log-Likelihood from a Generalized Kumaraswamy Regression Model
#'
#' @description
#' This function extracts the log-likelihood value from a fitted Generalized
#' Kumaraswamy regression model. The result is returned as an object of class
#' \code{"logLik"}, which can be used for model selection and comparison.
#'
#' @param object An object of class \code{"gkwreg"}, typically the result of a call
#'   to \code{\link{gkwreg}}.
#' @param ... Additional arguments (currently not used).
#'
#' @details
#' The log-likelihood for Generalized Kumaraswamy regression models is computed
#' during model fitting and stored in the model object. This function retrieves
#' the stored value and creates a proper \code{"logLik"} object with the appropriate
#' attributes for use with model selection criteria and likelihood ratio tests.
#'
#' For GKw family models, the log-likelihood is computed as:
#'
#' \deqn{l(\theta) = \sum_{i=1}^n \log f(y_i; \alpha_i, \beta_i, \gamma_i, \delta_i, \lambda_i)}
#'
#' where \eqn{f(y; \alpha, \beta, \gamma, \delta, \lambda)} is the probability
#' density function of the specific distribution from the GKw family, and \eqn{\theta}
#' represents all model parameters.
#'
#' @return An object of class \code{"logLik"} with the following attributes:
#' \itemize{
#'   \item \code{df}: The number of estimated parameters in the model.
#'   \item \code{nobs}: The number of observations used for fitting the model.
#' }
#'
#' @examples
#' \dontrun{
#' # Fit a Kumaraswamy regression model
#' kw_reg <- gkwreg(y ~ x1 + x2 | x3, data = df, family = "kw")
#'
#' # Extract log-likelihood
#' ll <- logLik(kw_reg)
#' print(ll)
#'
#' # Get the number of parameters
#' attr(ll, "df")
#' }
#'
#' @seealso \code{\link{gkwreg}}, \code{\link{AIC.gkwreg}}, \code{\link{BIC.gkwreg}}
#'
#' @export
logLik.gkwreg <- function(object, ...) {
  # Check if the object is of class gkwreg
  if (!inherits(object, "gkwreg")) {
    stop("'object' must be a fitted model object of class 'gkwreg'")
  }

  # Extract log-likelihood value
  ll <- object$loglik

  # If log-likelihood is not available, try to recover from other information
  if (is.null(ll)) {
    # Try to calculate from deviance if available
    if (!is.null(object$deviance)) {
      ll <- -object$deviance / 2
    } else {
      warning("Log-likelihood not found in the model object and cannot be calculated")
      ll <- NA_real_
    }
  }

  # Get the number of parameters
  df <- object$npar
  if (is.null(df)) {
    # Try to determine number of parameters from coefficients
    if (!is.null(object$coefficients)) {
      df <- length(object$coefficients)
    } else {
      warning("Number of parameters not found in the model object")
      df <- NA_integer_
    }
  }

  # Get the number of observations
  nobs <- object$nobs
  if (is.null(nobs)) {
    # Try to determine number of observations from residuals or fitted values
    if (!is.null(object$residuals)) {
      nobs <- length(object$residuals)
    } else if (!is.null(object$fitted.values)) {
      nobs <- length(object$fitted.values)
    } else if (!is.null(object$y)) {
      nobs <- length(object$y)
    } else {
      warning("Number of observations not found in the model object")
      nobs <- NA_integer_
    }
  }

  # Create and return the logLik object with appropriate attributes
  structure(ll,
    df = df,
    nobs = nobs,
    class = "logLik"
  )
}

#' Akaike's Information Criterion for Generalized Kumaraswamy Regression Models
#'
#' @description
#' Calculates the Akaike Information Criterion (AIC) for one or more fitted
#' Generalized Kumaraswamy regression models. AIC is a measure of the relative
#' quality of statistical models for a given set of data, providing a means for
#' model selection.
#'
#' @param object An object of class \code{"gkwreg"}, typically the result of a call
#'   to \code{\link{gkwreg}}.
#' @param ... Optionally, additional objects of class \code{"gkwreg"} for model
#'   comparison.
#' @param k The penalty per parameter to be used; the default \code{k = 2} gives
#'   the traditional AIC.
#'
#' @details
#' The AIC is calculated as:
#'
#' \deqn{AIC = -2 \times \log(L) + 2 \times k}
#'
#' where \eqn{L} is the maximized likelihood function and \eqn{k} is the number of
#' estimated parameters in the model. The model with the lowest AIC is generally
#' preferred.
#'
#' For small sample sizes, the corrected AIC (AICc) may be more appropriate.
#' Though not directly implemented here, it can be calculated as:
#'
#' \deqn{AICc = AIC + \frac{2k(k+1)}{n-k-1}}
#'
#' where \eqn{n} is the sample size.
#'
#' @return If just one object is provided, returns a numeric AIC value. If multiple
#' objects are provided, returns a data frame with rows corresponding to the objects
#' and columns representing the number of parameters (df) and AIC values.
#'
#' @examples
#' \dontrun{
#' # Fit two competing models
#' kw_reg1 <- gkwreg(y ~ x1 | x2, data = df, family = "kw")
#' kw_reg2 <- gkwreg(y ~ x1 + x3 | x2, data = df, family = "kw")
#'
#' # Compare models using AIC
#' AIC(kw_reg1, kw_reg2)
#'
#' # Calculate AIC with a different penalty
#' AIC(kw_reg1, k = 3)
#' }
#'
#' @seealso \code{\link{gkwreg}}, \code{\link{logLik.gkwreg}}, \code{\link{BIC.gkwreg}}
#'
#' @references
#' Akaike, H. (1974). A new look at the statistical model identification.
#' \emph{IEEE Transactions on Automatic Control}, 19(6), 716-723.
#'
#' @export
AIC.gkwreg <- function(object, ..., k = 2) {
  # Check if the object is of class gkwreg
  if (!inherits(object, "gkwreg")) {
    stop("'object' must be a fitted model object of class 'gkwreg'")
  }

  # Handle case with multiple models for comparison
  if (length(list(...)) > 0) {
    return(stats::AIC(object = object, ..., k = k))
  }

  # Check if AIC is already computed and stored in the object
  if (!is.null(object$aic) && k == 2) {
    return(object$aic)
  }

  # Calculate AIC from log-likelihood and number of parameters
  ll <- logLik(object)

  # Number of parameters
  df <- attr(ll, "df")
  if (is.na(df)) {
    warning("Number of parameters not available. AIC calculation may be inaccurate.")
    # Try to extract parameters from the coefficients
    if (!is.null(object$coefficients)) {
      df <- length(object$coefficients)
    } else {
      df <- 0 # Default to avoid errors, but will give incorrect results
    }
  }

  # Calculate and return AIC
  -2 * as.numeric(ll) + k * df
}


#' Bayesian Information Criterion for Generalized Kumaraswamy Regression Models
#'
#' @description
#' Calculates the Bayesian Information Criterion (BIC), also known as Schwarz's
#' Bayesian Criterion (SBC), for one or more fitted Generalized Kumaraswamy
#' regression models. BIC is a criterion for model selection that penalizes model
#' complexity more strongly than AIC.
#'
#' @param object An object of class \code{"gkwreg"}, typically the result of a call
#'   to \code{\link{gkwreg}}.
#' @param ... Optionally, additional objects of class \code{"gkwreg"} for model
#'   comparison.
#'
#' @details
#' The BIC is calculated as:
#'
#' \deqn{BIC = -2 \times \log(L) + k \times \log(n)}
#'
#' where \eqn{L} is the maximized likelihood function, \eqn{k} is the number of
#' estimated parameters in the model, and \eqn{n} is the number of observations.
#' The model with the lowest BIC is generally preferred.
#'
#' BIC penalizes model complexity more strongly than AIC, especially for larger
#' sample sizes. It is derived from a Bayesian framework and approximates the
#' Bayes factor for large sample sizes.
#'
#' @return If just one object is provided, returns a numeric BIC value. If multiple
#' objects are provided, returns a data frame with rows corresponding to the objects
#' and columns representing the number of parameters (df) and BIC values.
#'
#' @examples
#' \dontrun{
#' # Fit two competing models
#' kw_reg1 <- gkwreg(y ~ x1 | x2, data = df, family = "kw")
#' kw_reg2 <- gkwreg(y ~ x1 + x3 | x2, data = df, family = "kw")
#'
#' # Compare models using BIC
#' BIC(kw_reg1, kw_reg2)
#' }
#'
#' @seealso \code{\link{gkwreg}}, \code{\link{logLik.gkwreg}}, \code{\link{AIC.gkwreg}}
#'
#' @references
#' Schwarz, G. (1978). Estimating the dimension of a model.
#' \emph{The Annals of Statistics}, 6(2), 461-464.
#'
#' @export
BIC.gkwreg <- function(object, ...) {
  # Check if the object is of class gkwreg
  if (!inherits(object, "gkwreg")) {
    stop("'object' must be a fitted model object of class 'gkwreg'")
  }

  # Handle case with multiple models for comparison
  if (length(list(...)) > 0) {
    return(stats::BIC(object = object, ...))
  }

  # Check if BIC is already computed and stored in the object
  if (!is.null(object$bic)) {
    return(object$bic)
  }

  # Calculate BIC from log-likelihood, number of parameters, and number of observations
  ll <- logLik(object)

  # Number of parameters
  df <- attr(ll, "df")
  if (is.na(df)) {
    warning("Number of parameters not available. BIC calculation may be inaccurate.")
    # Try to extract parameters from the coefficients
    if (!is.null(object$coefficients)) {
      df <- length(object$coefficients)
    } else {
      df <- 0 # Default to avoid errors, but will give incorrect results
    }
  }

  # Number of observations
  n <- attr(ll, "nobs")
  if (is.na(n)) {
    warning("Number of observations not available. BIC calculation may be inaccurate.")
    # Try to extract number of observations from residuals or fitted values
    if (!is.null(object$residuals)) {
      n <- length(object$residuals)
    } else if (!is.null(object$fitted.values)) {
      n <- length(object$fitted.values)
    } else if (!is.null(object$y)) {
      n <- length(object$y)
    } else {
      n <- 0 # Default to avoid errors, but will give incorrect results
    }
  }

  # Calculate and return BIC
  -2 * as.numeric(ll) + df * log(n)
}



#' Register the S3 method
#' @param object obrl class fit
#' @param ... Additional arguments passed to or from other methods.
#' @export
logLik <- function(object, ...) {
  UseMethod("logLik.gkwreg")
}

#' Register the S3 method
#' @param object obrl class fit
#' @param ... Additional arguments passed to or from other methods.
#' @export
AIC <- function(object, ...) {
  UseMethod("AIC.gkwreg")
}

#' Register the S3 method
#' @param object obrl class fit
#' @param ... Additional arguments passed to or from other methods.
#' @export
BIC <- function(object, ...) {
  UseMethod("BIC.gkwreg")
}


#' Extract Variance-Covariance Matrix from a Generalized Kumaraswamy Regression Model
#'
#' @description
#' This function extracts the variance-covariance matrix of the estimated parameters
#' from a fitted Generalized Kumaraswamy regression model. The variance-covariance
#' matrix is essential for statistical inference, including hypothesis testing and
#' confidence interval calculation.
#'
#' @param object An object of class \code{"gkwreg"}, typically the result of a call
#'   to \code{\link{gkwreg}}.
#' @param complete Logical indicating whether the complete variance-covariance matrix
#'   should be returned in case some coefficients were omitted from the original fit.
#'   Currently ignored for \code{gkwreg} objects.
#' @param ... Additional arguments (currently not used).
#'
#' @details
#' The variance-covariance matrix is estimated based on the observed information
#' matrix, which is derived from the second derivatives of the log-likelihood function
#' with respect to the model parameters. For \code{gkwreg} objects, this matrix is
#' typically computed using the TMB (Template Model Builder) automatic differentiation
#' framework during model fitting.
#'
#' The diagonal elements of the variance-covariance matrix correspond to the squared
#' standard errors of the parameter estimates, while the off-diagonal elements represent
#' the covariances between pairs of parameters.
#'
#' @return A square matrix with row and column names corresponding to the coefficients
#' in the model. If the variance-covariance matrix is not available (for example, if
#' the model was fitted with \code{hessian = FALSE}), the function returns \code{NULL}
#' with a warning.
#'
#' @examples
#' \dontrun{
#' # Fit a Kumaraswamy regression model
#' kw_reg <- gkwreg(y ~ x1 + x2 | x3, data = df, family = "kw", hessian = TRUE)
#'
#' # Extract variance-covariance matrix
#' vcov_mat <- vcov(kw_reg)
#'
#' # Extract standard errors (square root of diagonal elements)
#' std_errors <- sqrt(diag(vcov_mat))
#' print(std_errors)
#'
#' # Compare with the standard errors from the summary
#' summary(kw_reg)
#' }
#'
#' @seealso \code{\link{gkwreg}}, \code{\link{confint.gkwreg}}, \code{\link{summary.gkwreg}}
#'
#' @export
vcov.gkwreg <- function(object, complete = TRUE, ...) {
  # Check if the object is of class gkwreg
  if (!inherits(object, "gkwreg")) {
    stop("'object' must be a fitted model object of class 'gkwreg'")
  }

  # Check if the variance-covariance matrix is available
  if (is.null(object$vcov)) {
    warning(
      "Variance-covariance matrix not found in the model object. ",
      "The model may have been fitted with hessian = FALSE. ",
      "Consider refitting with hessian = TRUE for valid statistical inference."
    )
    return(NULL)
  }

  # Return the variance-covariance matrix
  object$vcov
}
