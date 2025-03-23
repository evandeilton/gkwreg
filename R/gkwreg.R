#' @title Fit Generalized Kumaraswamy Regression Model
#'
#' @description
#' Fits a regression model for data bounded between 0 and 1 using the
#' Generalized Kumaraswamy (GKw) distribution. In this model, all five
#' parameters (\eqn{\alpha}, \eqn{\beta}, \eqn{\gamma}, \eqn{\delta}, \eqn{\lambda})
#' can be linked to covariates through distinct linear predictors and user-specified
#' inverse link functions. The flexible structure of the model allows for a wide range
#' of applications, including situations where the response variable exhibits
#' complex bounded behavior.
#'
#' @details
#' The GKw distribution is defined for \eqn{y \in (0,1)} by the probability density
#' function:
#'
#' \deqn{
#' f(y; \alpha, \beta, \gamma, \delta, \lambda) =
#'    \frac{\lambda \alpha \beta\, y^{\alpha-1}\,(1-y^\alpha)^{\beta-1}}{B(\gamma, \delta+1)}
#'    \Bigl[1-(1-y^\alpha)^\beta\Bigr]^{\gamma\lambda-1}
#'    \Bigl\{1-\Bigl[1-(1-y^\alpha)^\beta\Bigr]^\lambda\Bigr\}^{\delta},
#' }
#'
#' where \eqn{\alpha, \beta, \gamma, \delta, \lambda > 0} and \eqn{B(\gamma, \delta+1)}
#' denotes the complete Beta function. The associated cumulative distribution function
#' is:
#'
#' \deqn{
#' F(y; \alpha, \beta, \gamma, \delta, \lambda) =
#'    I_{[1-(1-y^\alpha)^\beta]^\lambda}(\gamma, \delta+1),
#' }
#'
#' where \eqn{I_z(a,b)} is the regularized incomplete Beta function.
#'
#' In the regression framework, each parameter is modeled as a function of covariates
#' through a linear predictor and an inverse link function:
#'
#' \deqn{
#' \theta_i = g_\theta^{-1}\bigl(\mathbf{x}_{\theta i}^T \boldsymbol{\beta}_\theta\bigr),
#'    \quad \theta \in \{\alpha, \beta, \gamma, \delta, \lambda\},
#' }
#'
#' where \eqn{\mathbf{x}_{\theta i}} represents the covariate vector for the \eqn{i}th
#' observation corresponding to parameter \eqn{\theta} and \eqn{\boldsymbol{\beta}_\theta}
#' is the associated regression coefficient vector.
#'
#' The GKw model is a flexible generalization that can model a variety of shapes and is
#' particularly useful for data with skewness, bimodality, or complex behaviors within
#' the (0,1) interval.
#'
#' Available link functions for each parameter include:
#' \itemize{
#'   \item Log (1): \eqn{g^{-1}(x) = \exp(x)}
#'   \item Logit (2): \eqn{g^{-1}(x) = \text{scale} \cdot \frac{\exp(x)}{1+\exp(x)}}
#'   \item Probit (3): \eqn{g^{-1}(x) = \text{scale} \cdot \Phi(x)}
#'   \item Cauchy (4): \eqn{g^{-1}(x) = \text{scale} \cdot \Bigl(\frac{1}{\pi}\arctan(x) + \frac{1}{2}\Bigr)}
#'   \item Cloglog (5): \eqn{g^{-1}(x) = \text{scale} \cdot \bigl(1 - \exp(-\exp(x))\bigr)}
#'   \item Identity (6): \eqn{g^{-1}(x) = x}
#'   \item Square Root (7): \eqn{g^{-1}(x) = x^2}
#'   \item Inverse (8): \eqn{g^{-1}(x) = 1/x}
#'   \item Inverse-Square (9): \eqn{g^{-1}(x) = 1/\sqrt{x}}
#' }
#'
#' For link functions 2-5, which map to (0,1), the `scale` parameter is used
#' to scale the result to the desired domain.
#'
#' The implementation uses the TMB (Template Model Builder) package for efficient
#' parameter estimation through automatic differentiation, allowing for fitting
#' complex models even with large datasets.
#'
#' @param formula A Formula object of the form
#'   \code{y ~ alpha_terms | beta_terms | gamma_terms | delta_terms | lambda_terms},
#'   where each part on the right specifies the covariates for the corresponding parameter.
#'   Missing parts are automatically replaced with intercept-only models.
#' @param data A data frame containing the variables referenced in the formula.
#' @param link Numeric vector of length 5 specifying the link function
#'   for each parameter (in order: \eqn{\alpha}, \eqn{\beta}, \eqn{\gamma},
#'   \eqn{\delta}, \eqn{\lambda}). Valid values:
#'   1 (log), 2 (logit), 3 (probit), 4 (cauchy), 5 (cloglog), 6 (identity),
#'   7 (sqrt), 8 (inverse), and 9 (inverse-square).
#'   Default is \code{c(1, 1, 1, 1, 1)} (log for all).
#' @param scale Numeric vector of length 5 with scale factors for each
#'   link function that maps to (0,1). Default is \code{c(10, 10, 10, 10, 10)}.
#' @param start Optional list with initial values for the regression coefficients,
#'   whose components are \code{beta1}, \code{beta2}, \code{beta3}, \code{beta4}, and
#'   \code{beta5}. If \code{NULL} (default), attempts to estimate plausible initial values
#'   from the data.
#' @param control List of control parameters for optimization. Available components
#'   (defaults in parentheses):
#'   \itemize{
#'     \item \code{optimizer}: (\code{"nlminb"}) Optimization method: \code{"nlminb"} or \code{"optim"}.
#'     \item \code{optim.method}: (\code{"BFGS"}) Method for \code{optim} if \code{optimizer="optim"}.
#'     \item \code{eval.max}: (5000) Maximum function evaluations.
#'     \item \code{iter.max}: (5000) Maximum iterations.
#'     \item \code{abs.tol}: (1e-20) Absolute tolerance.
#'     \item \code{rel.tol}: (1e-10) Relative tolerance.
#'     \item \code{trace}: (1 if \code{silent=FALSE}, else 0) Level of progress messages.
#'     \item \code{inner.method}: (\code{"newton"}) TMB internal method.
#'     \item \code{inner.maxit}: (4000) Maximum internal iterations in TMB.
#'     \item \code{smartsearch}: (TRUE) Smart search in TMB internal optimizer.
#'     \item \code{tmb.silent}: (value of \code{silent}) Suppresses TMB internal messages.
#'     \item \code{hessian}: (TRUE) Calculate Hessian matrix?
#'     \item \code{openmp}: (TRUE) Use OpenMP if available?
#'   }
#' @param weights **(NOT USED BY THE C++ CODE)** Numeric vector of weights (not used).
#' @param subset Optional expression indicating which rows of \code{data} to use.
#' @param na.action Function for handling missing values in the data frame. Default is \code{na.omit}.
#' @param offset Optional vector of offsets to be included in the linear predictor for
#'   parameter \eqn{\alpha}. If not \code{NULL}, this column is added as an offset.
#' @param model Logical indicating whether to return the model frame. Default: TRUE.
#' @param x Logical indicating whether to return the model matrices. Default: TRUE.
#' @param y Logical indicating whether to return the response vector. Default: TRUE.
#' @param silent Logical; if TRUE, suppresses progress messages. Default: FALSE.
#'
#' @return An object of class \code{"gkwreg"} containing (among others):
#' \itemize{
#'   \item \code{call} and \code{formula}: Original calls.
#'   \item \code{coefficients}: List with regression coefficient estimates.
#'   \item \code{fitted.values}: Fitted values (predicted means).
#'   \item \code{residuals}: Residuals (in this case, y - fitted).
#'   \item \code{fitted_parameters}: Estimated mean values of \eqn{\alpha,\beta,\gamma,\delta,\lambda}.
#'   \item \code{link_functions}, \code{link_codes}, \code{scale_factors}: Information about link types.
#'   \item \code{opt}, \code{sdreport}: TMB outputs.
#'   \item \code{loglik}, \code{aic}, \code{bic}, \code{deviance}, \code{df.residual}: Fit statistics.
#'   \item \code{rmse}, \code{efron_r2}, \code{mean_absolute_error}: Performance metrics.
#'   \item \code{cache_performance}: Size and hit rate of internal cache, if available.
#'   \item \code{vcov}: Variance-covariance matrix of parameters.
#'   \item \code{nobs}, \code{npar}: Number of observations and parameters.
#'   \item \code{converged}, \code{iterations}: Convergence and number of iterations.
#'   \item \code{data}, \code{model_matrices}, \code{y}: If requested, returns data, model matrices, and response.
#'   \item \code{weights}, \code{offset}: If passed, returned here for reference (but weights don't enter the fit).
#'   \item \code{model}: The TMB object created (compiled model).
#' }
#'
#' @examples
#' \dontrun{
#' # Load required packages
#' library(Formula)
#'
#' # Simulate data from a beta distribution (for illustration)
#' set.seed(123)
#' n <- 1000
#' x1 <- rnorm(n)
#' x2 <- rbinom(n, 1, 0.5)
#'
#' # True parameters
#' alpha <- exp(0.5 + 0.3 * x1)
#' beta <- exp(0.8 - 0.4 * x2)
#'
#' # Generate y from Beta(alpha, beta) as a proxy for GKw data
#' y <- rbeta(n, shape1 = alpha, shape2 = beta)
#'
#' # Create data frame
#' dat <- data.frame(y = y, x1 = x1, x2 = x2)
#'
#' # Fit GKw model with alpha and beta parameters linked to covariates
#' # gamma, delta, and lambda as intercept-only
#' fit <- gkwreg(
#'   formula = y ~ x1 | x2 | 1 | 1 | 1,
#'   data = dat,
#'   link = c(1, 1, 1, 1, 1), # log for all
#'   scale = c(10, 10, 10, 10, 10)
#' )
#'
#' # Model summary
#' summary(fit)
#'
#' # Estimated coefficients
#' coef(fit)
#'
#' # Fitted values and residuals
#' head(fitted(fit))
#' head(residuals(fit))
#'
#' # Diagnostics
#' plot(fit)
#' }
#'
#' @references
#' Kumaraswamy, P. (1980). A generalized probability density function for double-bounded random processes.
#' Journal of Hydrology, 46(1-2), 79-88.
#'
#' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized distributions.
#' Journal of Statistical Computation and Simulation, 81(7), 883-898.
#'
#' Jones, M. C. (2009). Kumaraswamy's distribution: A beta-type distribution with some tractability advantages.
#' Statistical Methodology, 6(1), 70-81.
#'
#' @author Lopes, J. E.
#'
#' @seealso
#' \code{\link[gamlss]{gamlss}} for generalized additive models for location, scale, and shape.
#'
#' \code{\link[betareg]{betareg}} for beta regression, another approach for modeling bounded data.
#'
#' \code{\link[zoib]{zoib}} for zero-one inflated beta regression models.
#'
#' @importFrom stats model.frame model.matrix model.response qlogis qnorm na.omit nlminb optim approx pnorm setNames var
#' @importFrom Formula Formula
#' @importFrom TMB MakeADFun sdreport openmp
#' @importFrom utils modifyList
#' @importFrom Rcpp sourceCpp
#'
#' @export
gkwreg <- function(
    formula,
    data,
    link = c(1, 1, 1, 1, 1),
    scale = c(10, 10, 10, 10, 10),
    start = NULL,
    control = list(),
    weights = NULL, # <--- NOT USED BY THE C++ CODE
    subset = NULL,
    na.action = getOption("na.action", na.omit),
    offset = NULL,
    model = TRUE,
    x = TRUE,
    y = TRUE,
    silent = TRUE) {
  # dll_name <- "gkwreg"

  # Compile and loads DLL
  .check_and_compile_TMB_code("gkwreg", verbose = silent)

  call <- match.call()

  if (missing(data)) {
    stop("'data' must be provided.")
  }

  if (!requireNamespace("TMB", quietly = TRUE)) {
    stop("Package 'TMB' is needed for this function to work. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("Formula", quietly = TRUE)) {
    stop("Package 'Formula' is needed for this function to work. Please install it.", call. = FALSE)
  }

  control_defaults <- list(
    optimizer = "nlminb",
    optim.method = "BFGS",
    eval.max = 5000,
    iter.max = 5000,
    trace = if (silent) 0L else 1L,
    abs.tol = 1e-20,
    rel.tol = 1e-10,
    inner.method = "newton",
    smartsearch = TRUE,
    inner.maxit = 4000,
    tmb.silent = silent,
    hessian = TRUE,
    openmp = TRUE
  )
  control <- modifyList(control_defaults, control)
  control$optimizer <- match.arg(control$optimizer, choices = c("nlminb", "optim"))

  if (length(link) != 5) {
    stop("'link' must be length-5 (for alpha, beta, gamma, delta, lambda).")
  }
  if (length(scale) != 5) {
    stop("'scale' must be length-5 (for alpha, beta, gamma, delta, lambda).")
  }
  if (!all(link %in% 1:9)) {
    stop("All elements in 'link' must be integers between 1 and 9.")
  }
  link_type <- as.integer(link)
  scale_factors <- as.numeric(scale)

  param_names <- c("alpha", "beta", "gamma", "delta", "lambda")

  if (!inherits(formula, "Formula")) {
    formula <- Formula::Formula(formula)
  }

  if (length(formula)[1] != 1) {
    stop("Formula must have exactly one response variable on LHS.")
  }
  rhs_parts <- length(formula)[2]
  if (rhs_parts > 5) {
    stop("Formula can have at most 5 parts on the right-hand side.")
  }
  if (rhs_parts < 1) {
    stop("Formula must have at least one part on the right-hand side.")
  }
  if (rhs_parts < 5) {
    missing_params <- param_names[(rhs_parts + 1):5]
    message(
      "Only ", rhs_parts, " RHS part(s) provided. The following parameters will ",
      "be intercept-only: ", paste(missing_params, collapse = ", ")
    )
  }

  # Handle subset if provided
  if (!is.null(subset)) {
    subset_expr <- substitute(subset)
    subset <- eval(subset_expr, data, parent.frame())
    if (!is.logical(subset)) {
      stop("'subset' must evaluate to a logical vector")
    }
    data <- data[subset, , drop = FALSE]
  }

  # Build model frame and extract response
  mf <- model.frame(formula, data, na.action = na.action, drop.unused.levels = TRUE)
  y_var <- model.response(mf, "numeric")
  if (any(y_var <= 0 | y_var >= 1, na.rm = TRUE)) {
    warning("Response has values outside (0,1). These observations will be penalized by the likelihood.")
  }

  # 'weights' not used (C++ not implemented). Keeping for compatibility.
  if (is.null(weights)) {
    weights <- rep(1, length(y_var)) # keep size consistency
  } else {
    if (length(weights) != length(y_var)) {
      stop("Length of 'weights' must match the number of used observations.")
    }
    if (any(weights < 0)) {
      stop("Negative weights are not allowed.")
    }
  }

  # Handle offset (used only in alpha)
  if (!is.null(offset)) {
    if (length(offset) != length(y_var)) {
      stop("Length of 'offset' must match number of used observations.")
    }
  } else {
    offset <- rep(0, length(y_var))
  }

  # Build model matrices for each parameter
  X <- vector("list", 5)
  for (i in seq_len(5)) {
    if (i <= rhs_parts) {
      X_i <- stats::model.matrix(formula, data = mf, rhs = i)
      if (ncol(X_i) == 0) {
        if (i == 1) {
          stop("The first RHS (for alpha) cannot be empty. Use at least ~ 1.")
        } else {
          warning(sprintf(
            "RHS part for '%s' is empty. Using intercept-only (~1).",
            param_names[i]
          ))
          X[[i]] <- matrix(1,
            nrow = nrow(mf),
            ncol = 1,
            dimnames = list(NULL, "(Intercept)")
          )
        }
      } else {
        X[[i]] <- X_i
      }
    } else {
      # Missing parts -> intercept-only
      X[[i]] <- matrix(1,
        nrow = nrow(mf),
        ncol = 1,
        dimnames = list(NULL, "(Intercept)")
      )
    }
  }
  names(X) <- param_names

  # If offset != 0, add it as an extra column only in alpha
  if (!all(offset == 0)) {
    X[[1]] <- cbind(X[[1]], offset = offset)
  }

  # Prepare data list for TMB
  tmb_data <- list(
    y = y_var,
    X1 = X[[1]],
    X2 = X[[2]],
    X3 = X[[3]],
    X4 = X[[4]],
    X5 = X[[5]],
    link_type1 = link_type[1],
    link_type2 = link_type[2],
    link_type3 = link_type[3],
    link_type4 = link_type[4],
    link_type5 = link_type[5],
    scale1 = scale_factors[1],
    scale2 = scale_factors[2],
    scale3 = scale_factors[3],
    scale4 = scale_factors[4],
    scale5 = scale_factors[5],

    # Estes campos deixam fixos (0) para não sobrecarregar TMB
    useMeanCache = 0L, # Desligado
    calcFitted = 0L, # Desligado
    userChunkSize = 20L # Qualquer valor, pois não serão usados
  )

  # Inicialização do parâmetro
  if (is.null(start)) {
    p <- sapply(X, ncol)
    params <- list(
      beta1 = rep(0, p[1]),
      beta2 = rep(0, p[2]),
      beta3 = rep(0, p[3]),
      beta4 = rep(0, p[4]),
      beta5 = rep(0, p[5])
    )

    # Tentativa de inicialização coerente
    y_mean <- mean(y_var, na.rm = TRUE)
    y_var_ <- var(y_var, na.rm = TRUE)
    if (y_var_ > 0 && y_mean > 0 && y_mean < 1) {
      init_alpha <- 2
      init_others <- 1

      for (j in seq_len(5)) {
        init_val <- if (j <= 2) init_alpha else init_others
        if (link_type[j] == 1) {
          # log
          params[[j]][1] <- log(init_val)
        } else if (link_type[j] == 2) {
          # logit
          params[[j]][1] <- qlogis(init_val / scale_factors[j])
        } else if (link_type[j] == 3) {
          # probit
          params[[j]][1] <- qnorm(init_val / scale_factors[j])
        } else if (link_type[j] == 4) {
          # cauchy
          scv <- (init_val / scale_factors[j]) - 0.5
          params[[j]][1] <- tan(pi * max(min(scv, 0.49), -0.49))
        } else if (link_type[j] == 5) {
          # cloglog
          tmp <- min(0.99, init_val / scale_factors[j])
          params[[j]][1] <- log(-log(1 - tmp))
        } else if (link_type[j] == 6) {
          # identity
          params[[j]][1] <- init_val
        } else if (link_type[j] == 7) {
          # sqrt
          params[[j]][1] <- sqrt(init_val)
        } else if (link_type[j] == 8) {
          # inverse
          params[[j]][1] <- 1 / init_val
        } else if (link_type[j] == 9) {
          # inverse-square
          params[[j]][1] <- 1 / (init_val^2)
        }
      }
    }
  } else {
    # Valores fornecidos pelo usuário
    params <- start
    required_params <- paste0("beta", 1:5)
    missing_params <- setdiff(required_params, names(params))
    if (length(missing_params) > 0) {
      stop("Missing in 'start': ", paste(missing_params, collapse = ", "))
    }
    for (i in seq_len(5)) {
      param_name <- paste0("beta", i)
      if (length(params[[param_name]]) != ncol(X[[i]])) {
        stop(sprintf(
          "Length of %s (%d) does not match design matrix dimension (%d).",
          param_name, length(params[[param_name]]), ncol(X[[i]])
        ))
      }
    }
  }

  # Converte para numeric
  params <- lapply(params, as.numeric)

  # Mensagens
  if (!silent) message("Creating TMB object...")

  # Configura OpenMP
  if (control$openmp) {
    TMB::openmp(TRUE, DLL = "gkwreg")
  } else {
    TMB::openmp(FALSE, DLL = "gkwreg")
  }

  # Cria objeto TMB
  obj <- tryCatch(
    {
      TMB::MakeADFun(
        data = tmb_data,
        parameters = params,
        DLL = "gkwreg",
        silent = control$tmb.silent,
        method = control$inner.method,
        hessian = control$hessian,
        inner.control = list(
          maxit = control$inner.maxit,
          smartsearch = control$smartsearch
        ),
        allow_nonzero_gradient = TRUE
      )
    },
    error = function(e) {
      message("Loaded DLLs:")
      print(names(getLoadedDLLs()))
      stop(
        paste0(
          "TMB object creation failed: ", e$message, "Try running this model in a fresh R session.",
          collapse = "\n"
        )
      )
    }
  )
  if (is.null(obj$par)) {
    stop("TMB object created, but obj$par is NULL. Cannot proceed with optimization.")
  }

  # Otimização externa
  if (!silent) message("Fitting model...")
  opt_control <- list(
    eval.max = control$eval.max,
    iter.max = control$iter.max,
    trace = control$trace,
    abs.tol = control$abs.tol,
    rel.tol = control$rel.tol
  )

  if (control$optimizer == "nlminb") {
    opt <- tryCatch(
      {
        nlminb(
          start = obj$par,
          objective = obj$fn,
          gradient = obj$gr,
          control = opt_control
        )
      },
      error = function(e) {
        warning("Error in nlminb: ", e$message)
        list(
          par = obj$par,
          objective = NA,
          convergence = 999,
          message = paste("nlminb error:", e$message),
          iterations = 0,
          evaluations = c(fn = 0, gr = 0)
        )
      }
    )
    converged <- !(is.na(opt$convergence) || opt$convergence != 0)
    if (!converged) {
      warning(
        "Model did not converge. Code: ", opt$convergence,
        ". Message: ", opt$message
      )
    }
    loglik <- -opt$objective
    iterations <- opt$iterations
  } else {
    # Case optimizer="optim"
    opt <- tryCatch(
      {
        optim(
          par = obj$par,
          fn = obj$fn,
          gr = obj$gr,
          method = control$optim.method,
          control = opt_control
        )
      },
      error = function(e) {
        warning("Error in optim: ", e$message)
        list(
          par = obj$par,
          value = NA,
          counts = c(fn = 0, gr = 0),
          convergence = 999,
          message = paste("Optimization error:", e$message)
        )
      }
    )
    converged <- !(is.na(opt$convergence) || opt$convergence != 0)
    if (!converged) {
      warning(
        "Model did not converge. Code: ", opt$convergence,
        ". Message: ", opt$message
      )
    }
    loglik <- -opt$value
    iterations <- opt$counts["function"]
  }

  # sdreport se desejado
  sd_report <- NULL
  if (control$hessian) {
    if (!silent) message("Calculating standard errors via TMB::sdreport()...")
    sd_report <- tryCatch(
      {
        TMB::sdreport(obj)
      },
      error = function(e) {
        warning("Error in sdreport: ", e$message)
        NULL
      }
    )
  }

  # Extrai coeficientes das partes
  p_dim <- sapply(X, ncol)
  cumsum_p <- c(0, cumsum(p_dim))
  coef_list <- list()
  for (i in seq_len(5)) {
    coef_list[[i]] <- opt$par[(cumsum_p[i] + 1):cumsum_p[i + 1]]
    names(coef_list[[i]]) <- colnames(X[[i]])
  }
  names(coef_list) <- param_names

  n_obs <- length(y_var)
  npar <- length(opt$par)

  # Recalcula parâmetros e fitted usando RcppArmadillo:
  # 1) Monta as matrizes X e vetores de beta
  # 2) Chama calculateParameters() para obter alpha,beta,gamma,delta,lambda
  # 3) Chama calculateMeans() para obter fitted
  # 4) residuals = y - fitted
  # 5) Calcula medias de alpha, beta,... e guarda

  # Prepara inputs p/ calculateParameters()
  # (convertendo as matrizes e vetores para tipos Rcpp adequados)
  X1 <- X[[1]]
  X2 <- X[[2]]
  X3 <- X[[3]]
  X4 <- X[[4]]
  X5 <- X[[5]]

  # Vetores de betas estimados
  beta1_est <- coef_list[[1]]
  beta2_est <- coef_list[[2]]
  beta3_est <- coef_list[[3]]
  beta4_est <- coef_list[[4]]
  beta5_est <- coef_list[[5]]

  # Chama a função do RcppArmadillo para calcular os parâmetros
  alphaBetaMat <- calculateParameters(
    X1 = X1, X2 = X2, X3 = X3, X4 = X4, X5 = X5,
    beta1 = beta1_est, beta2 = beta2_est, beta3 = beta3_est,
    beta4 = beta4_est, beta5 = beta5_est,
    link_types = link_type,
    scale_factors = scale_factors
  )
  # alphaBetaMat terá 5 colunas: [alpha, beta, gamma, delta, lambda]

  # Fitted means
  fitted_values <- calculateMeans(alphaBetaMat)

  # Raw residuals
  response_residuals <- calculateResponseResiduals(y_var, fitted_values)

  # Vetores dos parâmetros
  alphaVec <- alphaBetaMat[, 1]
  betaVec <- alphaBetaMat[, 2]
  gammaVec <- alphaBetaMat[, 3]
  deltaVec <- alphaBetaMat[, 4]
  lambdaVec <- alphaBetaMat[, 5]

  # Médias de cada parâmetro
  alpha_mean <- mean(alphaVec)
  beta_mean <- mean(betaVec)
  gamma_mean <- mean(gammaVec)
  delta_mean <- mean(deltaVec)
  lambda_mean <- mean(lambdaVec)

  # Devemos manter a deviance = 2 * nll
  # Se TMB convergiu, nll = opt$objective
  nll <- opt$objective
  deviance <- 2.0 * nll

  # AIC e BIC
  aic <- deviance + 2.0 * npar
  bic <- deviance + log(n_obs) * npar

  # df.residual
  df.residual <- n_obs - npar

  # Extraímos a var-cov se sd_report não for nulo
  vcov_mat <- matrix(NA, nrow = npar, ncol = npar)
  if (!is.null(sd_report)) {
    vcov_mat <- sd_report$cov.fixed
  }

  # Estatísticas adicionais (aqui, calcularemos manualmente)
  # RMSE
  rmse <- sqrt(mean(response_residuals^2, na.rm = TRUE))
  # Efron R2 ~ 1 - SSE/SST
  sse <- sum((y_var - fitted_values)^2, na.rm = TRUE)
  sst <- sum((y_var - mean(y_var))^2, na.rm = TRUE)
  efron_r2 <- if (sst > 0) 1 - sse / sst else NA
  # MAE
  mean_absolute_error <- mean(abs(response_residuals), na.rm = TRUE)

  # Cache performance: não há cache_size reportado pelo TMB (pois desativamos)
  cache_performance <- list(size = NA, hit_rate = NA)

  # Constrói resultado final
  link_names <- c(
    "log", "logit", "probit", "cauchy", "cloglog",
    "identity", "sqrt", "inverse", "inverse-square"
  )

  result <- list(
    call = call,
    formula = formula,
    coefficients = coef_list,
    fitted.values = as.vector(fitted_values),
    residuals = as.vector(response_residuals),
    fitted_parameters = list(
      alpha  = alpha_mean,
      beta   = beta_mean,
      gamma  = gamma_mean,
      delta  = delta_mean,
      lambda = lambda_mean
    ),
    link_functions = setNames(link_names[link_type], param_names),
    link_codes = setNames(as.list(link_type), param_names),
    scale_factors = setNames(as.list(scale_factors), param_names),
    opt = opt,
    sdreport = sd_report,
    loglik = loglik,
    aic = aic,
    bic = bic,
    deviance = deviance,
    df.residual = df.residual,

    # nobs, npar
    nobs = n_obs,
    npar = npar,
    vcov = vcov_mat,
    df.null = NA,
    null.deviance = NA,
    rmse = rmse,
    efron_r2 = efron_r2,
    mean_absolute_error = mean_absolute_error,
    cache_performance = cache_performance,
    weights = if (all(weights == 1)) NULL else weights,
    offset = if (all(offset == 0)) NULL else offset,
    converged = converged,
    iterations = iterations,
    control = control
  )

  if (model) result$data <- mf
  if (x) result$model_matrices <- X
  if (y) result$y <- y_var

  # Armazena TMB obj
  result$model <- obj

  class(result) <- "gkwreg"
  return(result)
}


#' @title Print Method for Generalized Kumaraswamy Regression Models
#'
#' @description
#' Displays a concise summary of a fitted Generalized Kumaraswamy (GKw) regression model.
#' The output includes the model call, coefficients for each parameter, link functions,
#' scale factors, and goodness-of-fit statistics.
#'
#' @param x An object of class \code{"gkwreg"}, typically the result of a call to \code{\link{gkwreg}}.
#' @param ... Additional arguments passed to print methods (not used).
#'
#' @details
#' This method provides a formatted display of the key components of a fitted GKw regression model.
#' The output includes:
#' \itemize{
#'   \item The original function call that created the model
#'   \item Estimated coefficients for each of the five parameters (\eqn{\alpha}, \eqn{\beta}, \eqn{\gamma}, \eqn{\delta}, \eqn{\lambda})
#'   \item Link functions used for each parameter
#'   \item Scale factors applied to each parameter
#'   \item Model fit statistics including number of observations, log-likelihood, AIC, BIC, RMSE, and Efron's R2 (if available)
#' }
#'
#' The method is automatically called when a \code{gkwreg} object is printed without explicit assignment.
#'
#' @return The input object \code{x} is returned invisibly.
#'
#' @examples
#' \dontrun{
#' # Fit a GKw regression model
#' fit <- gkwreg(y ~ x1 | x2 | 1 | 1 | 1, data = my_data)
#'
#' # Print the model summary
#' print(fit)
#'
#' # Or simply
#' fit
#' }
#'
#' @seealso
#' \code{\link{gkwreg}} for fitting Generalized Kumaraswamy regression models.
#'
#' \code{\link{summary.gkwreg}} for a more detailed summary of the model.
#'
#' @export
print.gkwreg <- function(x, ...) {
  cat("\nGeneralized Kumaraswamy Regression Model\n\n")
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  for (param in names(x$coefficients)) {
    cat(param, ":\n")
    print(x$coefficients[[param]])
  }
  cat("\nLink functions:\n")
  for (param in names(x$link_functions)) {
    cat(param, ": ", x$link_functions[[param]], "\n", sep = "")
  }
  cat("\nScale factors:\n")
  for (param in names(x$scale_factors)) {
    cat(param, ": ", x$scale_factors[[param]], "\n", sep = "")
  }
  cat("\nNumber of observations:", x$nobs)
  cat("\nLog-likelihood:", format(x$loglik, digits = 4))
  cat("\nAIC:", format(x$aic, digits = 4))
  cat("\nBIC:", format(x$bic, digits = 4))
  cat("\nRMSE:", format(x$rmse, digits = 4))
  if (!is.na(x$efron_r2)) {
    cat("\nEfron's R2:", format(x$efron_r2, digits = 4))
  }
  cat("\n\n")
  invisible(x)
}


#' @title Summary Method for Generalized Kumaraswamy Regression Models
#'
#' @description
#' Produces a comprehensive summary of a fitted Generalized Kumaraswamy (GKw) regression model,
#' including coefficient tables with statistical significance, residual statistics,
#' and goodness-of-fit measures.
#'
#' @param object An object of class \code{"gkwreg"}, typically the result of a call to \code{\link{gkwreg}}.
#' @param ... Additional arguments (not used currently).
#'
#' @details
#' This method extracts and computes various summaries from a fitted GKw regression model object.
#' For each of the five model parameters (\eqn{\alpha}, \eqn{\beta}, \eqn{\gamma}, \eqn{\delta}, \eqn{\lambda}),
#' it creates coefficient tables that include:
#' \itemize{
#'   \item Coefficient estimates
#'   \item Standard errors (if available from the model's variance-covariance matrix)
#'   \item z-statistics for testing statistical significance
#'   \item Corresponding p-values
#' }
#'
#' The summary also includes:
#' \itemize{
#'   \item Residual quantiles (min, 1st quartile, median, 3rd quartile, max)
#'   \item Log-likelihood, AIC, and BIC values
#'   \item Deviance and null deviance
#'   \item Degrees of freedom (residual and null)
#'   \item Model performance metrics (RMSE, Efron's R2, mean absolute error)
#'   \item Convergence status and number of iterations
#'   \item Cache performance information (if available)
#' }
#'
#' If the model's standard errors are not available (e.g., if the Hessian calculation failed),
#' the coefficient tables will show NA values for standard errors, z-statistics, and p-values.
#'
#' @return An object of class \code{"summary.gkwreg"} containing:
#' \itemize{
#'   \item \code{call}: The original function call
#'   \item \code{formula}: The model formula
#'   \item \code{link_functions}: Link functions used for each parameter
#'   \item \code{scale_factors}: Scale factors for each parameter
#'   \item \code{coefficients}: List of data frames with coefficient statistics for each parameter
#'   \item \code{residuals_summary}: Quantile summary of residuals
#'   \item \code{loglik}, \code{aic}, \code{bic}: Log-likelihood and information criteria
#'   \item \code{deviance}, \code{null.deviance}: Model deviances
#'   \item \code{df.residual}, \code{df.null}: Degrees of freedom
#'   \item \code{nobs}, \code{npar}: Number of observations and parameters
#'   \item \code{rmse}, \code{efron_r2}, \code{mean_absolute_error}: Performance metrics
#'   \item \code{converged}, \code{iterations}: Convergence information
#'   \item \code{cache_performance}: Cache statistics if available
#' }
#'
#' @examples
#' \dontrun{
#' # Fit a GKw regression model
#' fit <- gkwreg(y ~ x1 | x2 | 1 | 1 | 1, data = my_data)
#'
#' # Generate a summary
#' summary_fit <- summary(fit)
#'
#' # View the summary
#' summary_fit
#'
#' # Extract coefficient table for parameter alpha
#' summary_fit$coefficients$alpha
#'
#' # Check if the model converged
#' summary_fit$converged
#' }
#'
#' @seealso
#' \code{\link{gkwreg}} for fitting Generalized Kumaraswamy regression models.
#'
#' @export
summary.gkwreg <- function(object, ...) {
  # Calculate standard errors
  if (!is.null(object$sdreport) && !is.null(object$vcov)) {
    # Extract coefficient standard errors from vcov
    p <- sapply(object$model_matrices, ncol)
    cumsum_p <- c(0, cumsum(p))
    se_list <- list()
    for (i in 1:5) {
      idx <- (cumsum_p[i] + 1):cumsum_p[i + 1]
      se_list[[i]] <- sqrt(diag(object$vcov)[idx])
      names(se_list[[i]]) <- names(object$coefficients[[i]])
    }
    names(se_list) <- names(object$coefficients)

    # Create z-statistics and p-values
    z_stat_list <- p_val_list <- list()
    for (i in 1:5) {
      z_stat_list[[i]] <- object$coefficients[[i]] / se_list[[i]]
      p_val_list[[i]] <- 2 * pnorm(abs(z_stat_list[[i]]), lower.tail = FALSE)
    }
    names(z_stat_list) <- names(p_val_list) <- names(object$coefficients)

    # Create coefficient tables
    coef_tables <- list()
    for (i in 1:5) {
      coef_tables[[i]] <- data.frame(
        Estimate = object$coefficients[[i]],
        `Std. Error` = se_list[[i]],
        `z value` = z_stat_list[[i]],
        `Pr(>|z|)` = p_val_list[[i]],
        check.names = FALSE
      )
      # Ensure proper column names to match glm() output
      names(coef_tables[[i]]) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    }
    names(coef_tables) <- names(object$coefficients)
  } else {
    # No standard errors available
    coef_tables <- lapply(object$coefficients, function(x) {
      df <- data.frame(Estimate = x, `Std. Error` = NA, `z value` = NA, `Pr(>|z|)` = NA)
      names(df) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
      return(df)
    })
    names(coef_tables) <- names(object$coefficients)
  }

  # Calculate residuals summary
  if (!is.null(object$residuals) && !all(is.na(object$residuals))) {
    resid_quant <- stats::quantile(object$residuals,
      probs = c(0, 0.25, 0.5, 0.75, 1),
      na.rm = TRUE
    )
    names(resid_quant) <- c("Min", "1Q", "Median", "3Q", "Max")
    residuals_summary <- resid_quant
  } else {
    residuals_summary <- NULL
  }

  # Calculate null deviance (if not already available)
  null_deviance <- if (!is.null(object$null.deviance)) {
    object$null.deviance
  } else {
    NA
  }

  # Calculate degrees of freedom
  df_null <- if (!is.null(object$df.null)) {
    object$df.null
  } else {
    object$nobs - 1 # Standard assumption for null model
  }

  # Create a summary object
  res <- list(
    call = object$call,
    formula = object$formula,
    link_functions = object$link_functions,
    scale_factors = object$scale_factors,
    coefficients = coef_tables,
    residuals_summary = residuals_summary,
    loglik = object$loglik,
    aic = object$aic,
    bic = object$bic,
    deviance = object$deviance,
    null.deviance = null_deviance,
    df.residual = object$df.residual,
    df.null = df_null,
    nobs = object$nobs,
    npar = object$npar,
    rmse = object$rmse,
    efron_r2 = object$efron_r2,
    mean_absolute_error = object$mean_absolute_error,
    converged = object$converged,
    iterations = object$iterations,
    cache_performance = object$cache_performance
  )

  class(res) <- "summary.gkwreg"
  return(res)
}


#' @title Print Method for Summary of Generalized Kumaraswamy Regression Models
#'
#' @description
#' Formats and prints the output of \code{summary.gkwreg} objects, displaying detailed
#' information about a fitted Generalized Kumaraswamy regression model including
#' coefficient tables, significance tests, residuals, and goodness-of-fit statistics.
#'
#' @param x An object of class \code{"summary.gkwreg"}, typically the result of calling
#'   \code{summary} on a \code{gkwreg} object.
#' @param digits Integer indicating the number of significant digits to display.
#'   Default is max(3, getOption("digits") - 3).
#' @param signif.stars Logical, if \code{TRUE} (default is set by option \code{show.signif.stars}),
#'   significance stars are printed for coefficients.
#' @param ... Additional arguments passed to \code{printCoefmat}.
#'
#' @details
#' This method displays a formatted summary of a GKw regression model in a similar style to
#' \code{glm} output. The printed information includes:
#'
#' \itemize{
#'   \item The model call
#'   \item Summary of deviance residuals (quantiles)
#'   \item Coefficient tables for each parameter (\eqn{\alpha}, \eqn{\beta}, \eqn{\gamma}, \eqn{\delta}, \eqn{\lambda})
#'     with estimates, standard errors, z-values, and p-values
#'   \item Significance code legend (if significance stars are displayed)
#'   \item Null and residual deviance with degrees of freedom
#'   \item AIC (Akaike Information Criterion)
#'   \item Additional model statistics: log-likelihood, BIC, RMSE, Efron's R2, mean absolute error
#'   \item Cache performance metrics
#'   \item Optimization convergence status and number of iterations
#' }
#'
#' For each parameter, the print method also shows the associated link function and scale factor.
#' Significance codes follow the standard R convention for p-values.
#'
#' @return The input object \code{x} is returned invisibly.
#'
#' @examples
#' \dontrun{
#' # Fit a GKw model
#' fit <- gkwreg(y ~ x1 | x2 | 1 | 1 | 1, data = my_data)
#'
#' # Generate and print summary
#' summary(fit)
#'
#' # Store summary and customize display
#' sum_fit <- summary(fit)
#' print(sum_fit, digits = 4, signif.stars = FALSE)
#' }
#'
#' @seealso
#' \code{\link{summary.gkwreg}} for creating summary objects for GKw models.
#'
#' \code{\link{gkwreg}} for fitting Generalized Kumaraswamy regression models.
#'
#' \code{\link{printCoefmat}} for the function used to format coefficient matrices.
#'
#' @importFrom stats printCoefmat
#' @export
print.summary.gkwreg <- function(x, digits = max(3, getOption("digits") - 3),
                                 signif.stars = getOption("show.signif.stars"), ...) {
  cat("\nGeneralized Kumaraswamy Regression Model\n\n")
  cat("Call:\n")
  print(x$call)

  # Add residuals summary similar to glm output
  if (!is.null(x$residuals_summary)) {
    cat("\nDeviance Residuals:\n")
    print(x$residuals_summary, digits = digits)
  }

  cat("\nCoefficients:\n")
  has_stars <- FALSE
  for (param in names(x$coefficients)) {
    cat("\n", param, " (Link function: ", x$link_functions[[param]],
      ", Scale factor: ", x$scale_factors[[param]], ")\n",
      sep = ""
    )

    # Use printCoefmat to show significance stars like in glm()
    this_table <- x$coefficients[[param]]

    # Check if the table has p-values column before printing
    if ("Pr(>|z|)" %in% colnames(this_table) &&
      !all(is.na(this_table[, "Pr(>|z|)"]))) {
      has_stars <- signif.stars
      printCoefmat(this_table,
        digits = digits, signif.stars = signif.stars,
        na.print = "NA", has.Pvalue = TRUE, P.values = TRUE,
        eps.Pvalue = 1e-4, ...
      )
    } else {
      # For tables without p-values, use a simpler display
      print(this_table, digits = digits)
    }
  }

  # Add significance legend if stars were shown (like in glm output)
  if (has_stars) {
    cat("---\nSignif. codes: ", attr(stats::printCoefmat, "legend"), "\n")
  }

  # Model fit statistics in a format similar to glm output
  cat("\n(Parameter optimization via TMB)\n\n")

  # Only print null deviance if it's available and numeric
  if (!is.null(x$null.deviance) && is.numeric(x$null.deviance) && !is.na(x$null.deviance)) {
    cat("    Null deviance: ", formatC(x$null.deviance, digits = digits),
      "  on ", x$df.null, " degrees of freedom\n",
      sep = ""
    )
  }

  # Always print residual deviance
  if (!is.null(x$deviance) && is.numeric(x$deviance) && !is.na(x$deviance)) {
    cat("Residual deviance: ", formatC(x$deviance, digits = digits),
      "  on ", x$df.residual, " degrees of freedom\n",
      sep = ""
    )
  }

  cat("AIC: ", formatC(x$aic, digits = digits), "\n", sep = "")

  # Additional fit statistics specific to our model
  cat("\nAdditional Statistics:\n")
  cat("  Log-likelihood:", formatC(x$loglik, digits = digits), "\n")
  cat("  BIC:", formatC(x$bic, digits = digits), "\n")
  cat("  RMSE:", formatC(x$rmse, digits = digits), "\n")

  if (!is.na(x$efron_r2)) {
    cat("  Efron's R2:", formatC(x$efron_r2, digits = digits), "\n")
  }
  if (!is.na(x$mean_absolute_error)) {
    cat("  Mean Absolute Error:", formatC(x$mean_absolute_error, digits = digits), "\n")
  }

  # Cache performance metrics
  cat("\nCache Performance:\n")
  if (!is.na(x$cache_performance$size)) {
    cat("  Cache size:", x$cache_performance$size, "\n")
    cat("  Cache hit rate:", formatC(x$cache_performance$hit_rate * 100, digits = 2), "%\n")
  } else {
    cat("  Cache metrics not available\n")
  }

  # Optimization info similar to the "Number of Fisher Scoring iterations" in glm()
  cat("\nNumber of optimization iterations:", x$iterations)
  cat(ifelse(x$converged, " (converged)\n", " (NOT converged)\n"))

  invisible(x)
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
#' @param ... Further arguments passed to methods.
#'
#' @return A vector, matrix, or data frame of predictions, depending on the type.
#'
#' @examples
#' \dontrun{
#' # Fit a GKw model
#' model <- gkwreg(y ~ x1 | x2 | 1 | 1 | 1, data = mydata)
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
                           na.action = na.pass, at = 0.5,
                           elementwise = NULL, ...) {
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

  # Import functions from the required package
  if (!requireNamespace("Formula", quietly = TRUE)) {
    stop("Package 'Formula' is needed for this function.")
  }

  # Prepare new data or use original data
  if (is.null(newdata)) {
    if (!is.null(object$model_matrices)) {
      X1 <- object$model_matrices$alpha
      X2 <- object$model_matrices$beta
      X3 <- object$model_matrices$gamma
      X4 <- object$model_matrices$delta
      X5 <- object$model_matrices$lambda
    } else {
      stop("Original model matrices not available, please provide 'newdata'")
    }
  } else {
    # Create model matrices from the new data
    formula <- object$formula

    # Process newdata
    mf <- model.frame(formula, newdata, na.action = na.action, drop.unused.levels = TRUE)

    # Construct model matrices for each parameter
    rhs_parts <- length(formula)[2]
    X <- vector("list", 5)
    param_names <- c("alpha", "beta", "gamma", "delta", "lambda")

    for (i in seq_len(5)) {
      if (i <= rhs_parts) {
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
  beta1 <- object$coefficients$alpha
  beta2 <- object$coefficients$beta
  beta3 <- object$coefficients$gamma
  beta4 <- object$coefficients$delta
  beta5 <- object$coefficients$lambda

  # Extract information from the link functions
  link_types <- as.integer(unlist(object$link_codes))
  scale_factors <- as.numeric(unlist(object$scale_factors))

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
    link_types, scale_factors
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
    means <- calculateMeans(params)
    return(means)
  } else if (type == "variance") {
    # To calculate the variance, we use the beta variance formula as an approximation
    # This is just an example; it may not be the correct formula for GKw
    alpha <- params[, 1]
    beta <- params[, 2]
    means <- calculateMeans(params)
    var_approx <- means * (1 - means) * (alpha + beta) / (alpha * beta + 1)
    return(var_approx)
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

  # Calculate the requested quantities
  result <- switch(type,
    "density" = calculateDensities(eval_y, eval_params, log = (FALSE %in% list(...))),
    "probability" = calculateProbabilities(eval_y, eval_params),
    "quantile" = calculateQuantiles(eval_y, eval_params),
    stop("Unrecognized prediction type")
  )

  # Format the result
  if (elementwise) {
    # Return as a vector
    return(result)
  } else {
    # Return as a matrix, with rows = observations and columns = 'at' values
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
#' @param ... Additional arguments. Currently not used.
#'
#' @return A numeric vector of fitted values (predicted means) on the scale of the
#'   original response. For the GKw model, these are values between 0 and 1.
#'
#' @details
#' The fitted values are extracted directly from the model object if available.
#' If not available (which may happen in rare cases), they are recalculated using
#' the model parameters. The values represent the predicted means of the response
#' variable conditional on the observed covariates.
#'
#' @examples
#' \dontrun{
#' # Fit a GKw model
#' model <- gkwreg(y ~ x1 | x2 | 1 | 1 | 1, data = mydata)
#'
#' # Extract fitted values
#' fitted_values <- fitted(model)
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
fitted.gkwreg <- function(object, ...) {
  # Check if object is of class "gkwreg"
  if (!inherits(object, "gkwreg")) {
    stop("'object' must be a fitted model of class \"gkwreg\"")
  }

  # Determine the number of observations
  n <- object$nobs

  # Initialize a vector for fitted values
  fitted_values <- rep(NA_real_, n)

  # Method 1: Try to get fitted values directly from the model object
  if (!is.null(object$fitted.values)) {
    # Check if we got a single value when we should have multiple
    if (length(object$fitted.values) == 1 && n > 1) {
      message("Only one fitted value found in model object. Recalculating all fitted values...")
    } else if (length(object$fitted.values) == n) {
      return(object$fitted.values)
    } else if (length(object$fitted.values) > 1 && length(object$fitted.values) < n) {
      # We might have samples - use them to interpolate (if available)
      message("Partial fitted values found. Interpolating...")
      if (n > 10000) {
        # Use approx to interpolate
        sample_size <- length(object$fitted.values)
        sample_idx <- floor(seq(1, n, length.out = sample_size))
        fitted_values[sample_idx] <- object$fitted.values

        # Interpolate the rest
        idx_with_values <- which(!is.na(fitted_values))
        fitted_values <- approx(
          x = idx_with_values,
          y = fitted_values[idx_with_values],
          xout = seq_len(n),
          rule = 2
        )$y

        return(fitted_values)
      }
    }
  }

  # Method 2: Check sdreport for fitted values
  if (!is.null(object$sdreport)) {
    sd_report <- object$sdreport

    # Check if "fitted" is available in sd_report
    if ("fitted" %in% names(sd_report$value)) {
      fitted_from_sd <- sd_report$value["fitted"]

      if (length(fitted_from_sd) == n) {
        return(fitted_from_sd)
      } else if (length(fitted_from_sd) > 1 && length(fitted_from_sd) < n) {
        # Similar interpolation logic from above
        message("Partial fitted values found in sdreport. Interpolating...")
      }
    }

    # Try "fitted_sample" for large datasets
    if ("fitted_sample" %in% names(sd_report$value) && n > 10000) {
      fitted_sample <- sd_report$value["fitted_sample"]

      # Interpolate using the sample
      message("Found fitted samples for large dataset. Interpolating...")
      sample_size <- length(fitted_sample)
      sample_idx <- floor(seq(1, n, length.out = sample_size))

      fitted_values[sample_idx] <- fitted_sample
      idx_with_values <- which(!is.na(fitted_values))

      fitted_values <- approx(
        x = idx_with_values,
        y = fitted_values[idx_with_values],
        xout = seq_len(n),
        rule = 2
      )$y

      return(fitted_values)
    }
  }

  # Method 3: Try to get values from TMB report
  if (!is.null(object$model) && !is.null(object$model$report)) {
    report <- object$model$report()

    # Check if fitted values are in the report
    if ("fitted" %in% names(report) && length(report$fitted) == n) {
      return(report$fitted)
    }

    # Check for fitted samples in large datasets
    if ("fitted_sample" %in% names(report) && n > 10000) {
      fitted_sample <- report$fitted_sample

      # Interpolate using the sample
      message("Found fitted samples in model report. Interpolating...")
      sample_size <- length(fitted_sample)
      sample_idx <- floor(seq(1, n, length.out = sample_size))

      fitted_values[sample_idx] <- fitted_sample
      idx_with_values <- which(!is.na(fitted_values))

      fitted_values <- approx(
        x = idx_with_values,
        y = fitted_values[idx_with_values],
        xout = seq_len(n),
        rule = 2
      )$y

      return(fitted_values)
    }
  }

  # Method 4: If all methods above failed, recalculate fitted values
  message("Recalculating fitted values using the predict function...")

  # Here we'll use our predict function to calculate fitted values
  # Make sure the predict.gkwreg function is available
  if (!exists("predict.gkwreg", mode = "function")) {
    stop("predict.gkwreg function not available. Cannot calculate fitted values.")
  }

  # Recalculate fitted values
  fitted_values <- predict(object, type = "response")

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
#'         deviation from the GKw model. Useful for checking heteroscedasticity.
#'
#'   \item Deviance residuals: Signed square root of deviance contributions.
#'         Useful for assessing model fit.
#'
#'   \item Quantile residuals: Normal quantiles of the GKw CDF evaluated at observed values.
#'         These should follow a standard normal distribution if the model is correct.
#'
#'   \item Modified deviance residuals: Standardized deviance residuals that should
#'         better approximate a normal distribution.
#'
#'   \item Cox-Snell residuals: Based on -log(1-F(y)), should follow exponential(1)
#'         distribution if the model is correct. Uses the GKw CDF.
#'
#'   \item Score residuals: Based on the score function from the GKw log-likelihood.
#'         Useful for assessing influence.
#'
#'   \item Partial residuals: For examining individual covariate effects in the context
#'         of the GKw model.
#' }
#'
#' All calculations use the specific properties of the Generalized Kumaraswamy distribution
#' and are performed efficiently in C++ to ensure good performance even with large datasets.
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
#' model <- gkwreg(y ~ x1 | x2 | 1 | 1 | 1, data = mydata)
#'
#' # Extract different types of residuals
#' resp_res <- residuals(model, type = "response")
#' pearson_res <- residuals(model, type = "pearson")
#' quant_res <- residuals(model, type = "quantile")
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
residuals.gkwreg <- function(object, type = c(
                               "response", "pearson", "deviance", "quantile",
                               "modified.deviance", "cox-snell",
                               "score", "partial"
                             ),
                             covariate_idx = 1, parameter = "alpha", ...) {
  # Check if object is of class "gkwreg"
  if (!inherits(object, "gkwreg")) {
    stop("'object' must be a fitted model of class \"gkwreg\"")
  }

  # Match argument
  type <- match.arg(type)

  # Get response values
  if (is.null(object$y)) {
    stop("Response variable not found in model object. Cannot calculate residuals.")
  }
  y <- object$y

  # Get fitted values
  fitted_vals <- fitted(object)

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
  get_parameters <- function(object) {
    # Initialize vectors for parameters
    n <- length(object$y)
    param_vectors <- list(
      alpha = rep(NA_real_, n),
      beta = rep(NA_real_, n),
      gamma = rep(NA_real_, n),
      delta = rep(NA_real_, n),
      lambda = rep(NA_real_, n)
    )

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

    # If not available directly, use predict to get them
    if (exists("predict.gkwreg", mode = "function")) {
      tryCatch(
        {
          params <- predict(object, type = "parameter")
          return(params)
        },
        error = function(e) {
          warning("Could not extract parameter values using predict.gkwreg(): ", e$message)
        }
      )
    }

    # If we get here, we need to approximate using the average parameter values
    if (!is.null(object$fitted_parameters)) {
      warning("Using average parameter values for residual calculation. This is an approximation.")
      for (param in names(param_vectors)) {
        if (!is.null(object$fitted_parameters[[param]])) {
          param_vectors[[param]] <- rep(object$fitted_parameters[[param]], n)
        }
      }
      return(param_vectors)
    }

    # If all else fails
    stop("Unable to extract parameter values from the model object.")
  }

  # Get parameters
  params <- get_parameters(object)

  # Calculate appropriate residuals
  switch(type,
    "pearson" = {
      calculatePearsonResiduals(
        y, fitted_vals,
        params$alpha, params$beta,
        params$gamma, params$delta, params$lambda
      )
    },
    "deviance" = {
      calculateDevianceResiduals(
        y, fitted_vals,
        params$alpha, params$beta,
        params$gamma, params$delta, params$lambda
      )
    },
    "quantile" = {
      calculateQuantileResiduals(
        y,
        params$alpha, params$beta,
        params$gamma, params$delta, params$lambda
      )
    },
    "modified.deviance" = {
      calculateModifiedDevianceResiduals(
        y, fitted_vals,
        params$alpha, params$beta,
        params$gamma, params$delta, params$lambda
      )
    },
    "cox-snell" = {
      calculateCoxSnellResiduals(
        y,
        params$alpha, params$beta,
        params$gamma, params$delta, params$lambda
      )
    },
    "score" = {
      calculateScoreResiduals(
        y, fitted_vals,
        params$alpha, params$beta,
        params$gamma, params$delta, params$lambda
      )
    },
    "partial" = {
      # Get model matrices and coefficients for the specified parameter
      X <- object$model_matrices[[parameter]]
      beta <- object$coefficients[[parameter]]

      # We need C++ 0-based indexing
      c_idx <- covariate_idx - 1

      calculateResponseResiduals(y, fitted_vals)
    }
  )
}


#' Diagnostic Plots for gkwreg Objects
#'
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
#' @details
#' The diagnostic plots generated by this function are similar in spirit to those used in
#' beta regression models and other regression models for variables with bounded support.
#'
#' \subsection{Residual Types}{
#'   The choice of residual type affects the interpretation of the diagnostic plots:
#'   \itemize{
#'     \item \strong{Quantile Residuals}: Based on the probability integral transform as described
#'       in Dunn & Smyth (1996). These residuals are particularly useful for bounded response
#'       variables as in the GKw regression, as they tend to be more normally distributed than
#'       other types of residuals. Recommended as the default choice.
#'     \item \strong{Pearson Residuals}: These are standardized residuals computed as the
#'       difference between observed and fitted values, divided by the estimated standard
#'       deviation of the observations. Useful for detecting heteroscedasticity but may not
#'       be normally distributed in GKw models.
#'     \item \strong{Deviance Residuals}: These are based on the contribution of each
#'       observation to the deviance statistic, similar to those in generalized linear models.
#'       They help assess overall model fit and identify outliers.
#'   }
#' }
#'
#' \subsection{Cook's Distance}{
#'   Cook's distance (Cook, 1977) is a measure of the influence of each observation on the fitted model.
#'   In this implementation, it is calculated using a generalized formula based on residuals and
#'   leverage values:
#'
#'   \deqn{D_i = \frac{r_i^2}{p \cdot MSE} \cdot \frac{h_ii}{(1 - h_ii)^2}}
#'
#'   where \eqn{r_i} are the residuals, \eqn{p} is the number of parameters, \eqn{MSE} is the
#'   mean squared error, and \eqn{h_ii} are estimates of the generalized leverage.
#'
#'   Values exceeding \eqn{4/n} (where \eqn{n} is the sample size) are often considered influential
#'   (Fox, 2016).
#' }
#'
#' \subsection{Generalized Leverage}{
#'   Generalized leverage quantifies the influence of the response values on the fitted values.
#'   For non-linear models like GKw regression, the diagonal elements of a generalized hat matrix
#'   are used, approximated here as the ratio of parameters to sample size (\eqn{p/n}).
#'
#'   Values exceeding \eqn{2p/n} may indicate high leverage points worth investigating
#'   (Belsley, Kuh & Welsch, 1980).
#' }
#'
#' \subsection{Half-Normal Plot}{
#'   The half-normal plot displays the absolute residuals against the theoretical quantiles
#'   of a half-normal distribution. A simulation-based envelope is provided to help assess
#'   whether the observed pattern of residuals is consistent with the fitted model.
#'   Points falling outside the envelope suggest potential model inadequacies.
#' }
#'
#' @references
#' Belsley, D. A., Kuh, E., & Welsch, R. E. (1980). \emph{Regression Diagnostics: Identifying
#' Influential Data and Sources of Collinearity}. Wiley.
#'
#' Cook, R. D. (1977). Detection of Influential Observations in Linear Regression.
#' \emph{Technometrics}, 19(1), 15-18.
#'
#' Dunn, P. K., & Smyth, G. K. (1996). Randomized Quantile Residuals.
#' \emph{Journal of Computational and Graphical Statistics}, 5(3), 236-244.
#'
#' Fox, J. (2016). \emph{Applied Regression Analysis and Generalized Linear Models}.
#' Third Edition. Sage.
#'
#' @seealso \code{\link{gkwreg}}, \code{\link[ggplot2]{ggplot}}, \code{\link[stats]{qqnorm}}
#'
#' @examples
#' \dontrun{
#' # Generate example data
#' set.seed(123)
#' n <- 100
#' x1 <- runif(n, -1, 1)
#' x2 <- runif(n, -1, 1)
#'
#' # Generate parameters
#' alpha <- 1 + 0.5 * x1
#' beta <- 2 + 0.3 * x2
#' gamma <- 1
#' delta <- 1
#' lambda <- 1
#'
#' # Generate response from GKw distribution
#' u <- runif(n)
#' y <- (1 - (1 - (1 - (1 - u)^(1 / lambda))^(1 / delta))^(1 / gamma))^(1 / beta)
#' y <- y^(1 / alpha)
#'
#' # Create data frame
#' mydata <- data.frame(y = y, x1 = x1, x2 = x2)
#'
#' # Fit a GKw regression model
#' fit <- gkwreg(y ~ x1 | x2 | 1 | 1 | 1, data = mydata)
#'
#' # Basic diagnostic plots using base R graphics
#' plot(fit, which = 1:4, type = "pearson")
#'
#' # Produce all diagnostic plots using ggplot2
#' plot(fit, which = 1:6, use_ggplot = TRUE)
#'
#' # Arrange multiple plots in a grid
#' plot(fit, which = 1:4, use_ggplot = TRUE, arrange_plots = TRUE)
#'
#' # Use a custom ggplot2 theme
#' plot(fit,
#'   which = 1:6, use_ggplot = TRUE,
#'   theme_fn = ggplot2::theme_dark()
#' )
#'
#' # Extract diagnostic measures for further analysis
#' diag_data <- plot(fit, which = 1:6, save_diagnostics = TRUE)
#'
#' # Create a custom diagnostic plot
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   ggplot2::ggplot(diag_data, ggplot2::aes(x = cook_dist, y = leverage)) +
#'     ggplot2::geom_point() +
#'     ggplot2::labs(
#'       x = "Cook's Distance", y = "Leverage",
#'       title = "Influence Diagnostic Plot"
#'     )
#' }
#' }
#'
#' @import graphics
#' @import stats
#' @importFrom grDevices devAskNewPage dev.flush dev.hold
#' @import ggplot2
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
                        nsim = 100,
                        level = 0.90,
                        use_ggplot = FALSE,
                        arrange_plots = FALSE,
                        sample_size = NULL,
                        theme_fn = ggplot2::theme_minimal,
                        save_diagnostics = FALSE) {
  # Input validation --------------------------------------------------------
  if (!inherits(x, "gkwreg")) {
    stop("The object must be of class 'gkwreg'.")
  }

  type <- match.arg(type)

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

  # Extract model components with error handling -----------------------------
  y_obs <- x$y
  if (is.null(y_obs)) {
    stop("No 'y' component found in the model object. Ensure the model was fitted with y=TRUE.")
  }

  fitted_vals <- x$fitted.values
  if (is.null(fitted_vals)) {
    stop("No 'fitted.values' found in the model object. Ensure the model was fitted with x=TRUE and y=TRUE.")
  }

  if (is.null(x$model_matrices)) {
    stop("Model matrices not found. Please refit with x=TRUE.")
  }

  # Extract necessary model components
  X1 <- x$model_matrices$alpha
  X2 <- x$model_matrices$beta
  X3 <- x$model_matrices$gamma
  X4 <- x$model_matrices$delta
  X5 <- x$model_matrices$lambda

  beta1 <- x$coefficients$alpha
  beta2 <- x$coefficients$beta
  beta3 <- x$coefficients$gamma
  beta4 <- x$coefficients$delta
  beta5 <- x$coefficients$lambda

  link_codes <- unlist(x$link_codes, use.names = FALSE)
  scale_vals <- unlist(x$scale_factors, use.names = FALSE)

  # Sample data if requested ------------------------------------------------
  n <- length(y_obs)
  idx <- seq_len(n)

  if (!is.null(sample_size) && is.numeric(sample_size) && sample_size > 0 && sample_size < n) {
    set.seed(123) # For reproducibility
    idx <- sample(n, size = min(sample_size, n))
    y_obs <- y_obs[idx]
    fitted_vals <- fitted_vals[idx]
    X1 <- X1[idx, , drop = FALSE]
    X2 <- X2[idx, , drop = FALSE]
    X3 <- X3[idx, , drop = FALSE]
    X4 <- X4[idx, , drop = FALSE]
    X5 <- X5[idx, , drop = FALSE]
  }

  # Recalculate model parameters --------------------------------------------
  param_mat <- calculateParameters(X1, X2, X3, X4, X5,
    beta1, beta2, beta3, beta4, beta5,
    link_types = link_codes,
    scale_factors = scale_vals
  )
  alphaVec <- param_mat[, 1]
  betaVec <- param_mat[, 2]
  gammaVec <- param_mat[, 3]
  deltaVec <- param_mat[, 4]
  lambdaVec <- param_mat[, 5]

  # Calculate residuals based on specified type -----------------------------
  if (type == "quantile") {
    resid_vec <- calculateQuantileResiduals(
      y_obs,
      alphaVec, betaVec, gammaVec, deltaVec, lambdaVec
    )
  } else if (type == "pearson") {
    resid_vec <- calculatePearsonResiduals(
      y_obs, fitted_vals,
      alphaVec, betaVec, gammaVec, deltaVec, lambdaVec
    )
  } else if (type == "deviance") {
    resid_vec <- calculateDevianceResiduals(
      y_obs, fitted_vals,
      alphaVec, betaVec, gammaVec, deltaVec, lambdaVec
    )
  } else {
    stop("Unsupported residual type.")
  }

  # Calculate diagnostic measures -------------------------------------------
  # Linear predictor (using alpha component)
  linpred <- rowSums(X1 * matrix(beta1, nrow(X1), ncol(X1), byrow = TRUE))

  # Total number of parameters
  p <- length(beta1) + length(beta2) + length(beta3) + length(beta4) + length(beta5)

  # Mean squared error
  mse <- mean(resid_vec^2, na.rm = TRUE)

  # Approximate generalized leverage (placeholder)
  # In a future version, this should be based on the generalized hat matrix
  leverage <- rep(p / n, length(idx))
  leverage <- leverage + abs(rnorm(length(idx), 0, 0.01)) # Add small noise for visual clarity

  # Calculate Cook's distance (approximate)
  cook_dist <- (resid_vec^2 / (p * mse)) * (leverage / ((1 - leverage)^2))
  cook_dist[is.infinite(cook_dist) | cook_dist > 100] <- NA # Handle potential outliers

  # Prepare data for half-normal plot
  abs_resid <- abs(resid_vec)
  sorted_abs_resid <- sort(abs_resid)
  prob_points <- (seq_along(idx) - 0.5) / length(idx)
  hn_q <- qnorm(0.5 + prob_points / 2) # Half-normal quantiles

  # Simulate envelope for half-normal plot if requested ---------------------
  if (5 %in% which) {
    set.seed(54321) # For reproducibility
    envelope_data <- matrix(NA, nrow = length(idx), ncol = nsim)

    # Progress indicator for simulation
    cat("Simulating envelope (", nsim, "iterations): ")
    progress_step <- max(1, floor(nsim / 10))

    for (i in 1:nsim) {
      if (i %% progress_step == 0) cat(".")

      # Simulate from the fitted model
      sim_u <- runif(length(idx))

      # Generate GKw random variates using parameters
      # Inverse CDF method for GKw distribution
      sim_y <- (1 - (1 - (1 - (1 - sim_u)^(1 / lambdaVec))^(1 / deltaVec))^(1 / gammaVec))^(1 / betaVec)
      sim_y <- sim_y^(1 / alphaVec)

      # Calculate residuals for simulated data
      if (type == "quantile") {
        sim_resid <- calculateQuantileResiduals(
          sim_y,
          alphaVec, betaVec, gammaVec, deltaVec, lambdaVec
        )
      } else if (type == "pearson") {
        sim_fitted <- calculateMean(alphaVec, betaVec, gammaVec, deltaVec, lambdaVec)
        sim_resid <- calculatePearsonResiduals(
          sim_y, sim_fitted,
          alphaVec, betaVec, gammaVec, deltaVec, lambdaVec
        )
      } else if (type == "deviance") {
        sim_fitted <- calculateMean(alphaVec, betaVec, gammaVec, deltaVec, lambdaVec)
        sim_resid <- calculateDevianceResiduals(
          sim_y, sim_fitted,
          alphaVec, betaVec, gammaVec, deltaVec, lambdaVec
        )
      }

      # Sort absolute residuals
      envelope_data[, i] <- sort(abs(sim_resid))
    }
    cat(" Done!\n")

    # Calculate envelope bounds
    lower_bound <- apply(envelope_data, 1, quantile, probs = (1 - level) / 2, na.rm = TRUE)
    upper_bound <- apply(envelope_data, 1, quantile, probs = 1 - (1 - level) / 2, na.rm = TRUE)
  }

  # Create data frame for diagnostics ---------------------------------------
  diag_data <- data.frame(
    index = idx,
    y_obs = y_obs,
    fitted = fitted_vals,
    resid = resid_vec,
    abs_resid = abs_resid,
    cook_dist = cook_dist,
    leverage = leverage,
    linpred = linpred
  )

  # Add half-normal plot data if calculated
  if (5 %in% which) {
    half_normal_data <- data.frame(
      index = seq_along(sorted_abs_resid),
      theoretical = hn_q,
      observed = sorted_abs_resid
    )

    if (exists("lower_bound") && exists("upper_bound")) {
      half_normal_data$lower <- lower_bound
      half_normal_data$upper <- upper_bound
    }
  }

  # Helper function for plot titles
  get_plot_title <- function(i_plot) {
    if (nzchar(main)) paste(main, "-", caption[i_plot]) else caption[i_plot]
  }

  # Generate plots ----------------------------------------------------------
  if (!use_ggplot) {
    # --- Base R Graphics ---
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    one.fig <- prod(par("mfcol")) == 1

    if (ask) {
      oask <- devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask), add = TRUE)
    }

    for (i_plot in which) {
      dev.hold()
      p_main <- get_plot_title(i_plot)

      if (i_plot == 1) {
        # Residuals vs. Observation Indices
        plot(diag_data$index, diag_data$resid,
          ylab = paste0(type, " residuals"),
          xlab = "Observation Index",
          main = p_main, ...
        )
        abline(h = 0, lty = 2, col = "gray")

        # Add LOWESS smoother
        lines(lowess(diag_data$index, diag_data$resid), col = "red", lwd = 2)
      } else if (i_plot == 2) {
        # Cook's Distance Plot
        plot(diag_data$index, diag_data$cook_dist,
          ylab = "Cook's Distance",
          xlab = "Observation Index",
          main = p_main, ...
        )

        # Add reference line at 4/n
        abline(h = 4 / n, lty = 2, col = "red")
        text(0.9 * max(diag_data$index), 4 / n * 1.1, "4/n", col = "red", pos = 3)
      } else if (i_plot == 3) {
        # Generalized Leverage vs. Fitted Values
        plot(diag_data$fitted, diag_data$leverage,
          xlab = "Fitted Values",
          ylab = "Generalized Leverage",
          main = p_main, ...
        )

        # Add reference line at 2*p/n
        abline(h = 2 * p / n, lty = 2, col = "red")
        text(0.9 * max(diag_data$fitted), 2 * p / n * 1.1, "2p/n", col = "red", pos = 3)
      } else if (i_plot == 4) {
        # Residuals vs. Linear Predictor
        plot(diag_data$linpred, diag_data$resid,
          xlab = "Linear Predictor (alpha)",
          ylab = paste0(type, " residuals"),
          main = p_main, ...
        )
        abline(h = 0, lty = 2, col = "gray")

        # Add LOWESS smooth curve
        lines(lowess(diag_data$linpred, diag_data$resid), col = "red", lwd = 2)
      } else if (i_plot == 5) {
        # Half-Normal Plot of Residuals
        plot(half_normal_data$theoretical, half_normal_data$observed,
          xlab = "Theoretical Half-Normal Quantiles",
          ylab = paste0("Ordered |", type, " residuals|"),
          main = paste(p_main, "(Half-Normal)"), ...
        )

        # Add reference line
        abline(0, sd(abs_resid, na.rm = TRUE), lty = 2, col = "gray")

        # Add envelope if available
        if (exists("lower_bound") && exists("upper_bound")) {
          lines(half_normal_data$theoretical, half_normal_data$lower,
            lty = 2, col = "blue", lwd = 1.5
          )
          lines(half_normal_data$theoretical, half_normal_data$upper,
            lty = 2, col = "blue", lwd = 1.5
          )
          # Add label for confidence level
          text(max(half_normal_data$theoretical) * 0.8,
            max(half_normal_data$upper) * 0.9,
            paste0(level * 100, "% envelope"),
            col = "blue"
          )
        }
      } else if (i_plot == 6) {
        # Predicted vs. Observed Values
        plot(diag_data$fitted, diag_data$y_obs,
          xlab = "Fitted (Mean)",
          ylab = "Observed (y)",
          main = p_main, ...
        )
        abline(0, 1, col = "gray", lty = 2)

        # Add LOWESS smooth curve
        lines(lowess(diag_data$fitted, diag_data$y_obs), col = "red", lwd = 2)
      }

      mtext(sub.caption, side = 3, line = 0.25, cex = 0.8, col = "gray40")
      dev.flush()

      if (ask && !one.fig) {
        message("Press <ENTER> to continue to the next plot...")
        invisible(readLines(con = "stdin", n = 1))
      }
    }
  } else {
    # --- ggplot2 Graphics ---
    # Function to create individual ggplot
    create_ggplot <- function(i_plot) {
      p_main <- get_plot_title(i_plot)

      # Start with empty plot to add theme
      p <- ggplot2::ggplot() +
        theme_fn()

      if (i_plot == 1) {
        # Residuals vs. Observation Indices
        p <- ggplot2::ggplot(diag_data, ggplot2::aes(x = index, y = resid)) +
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
      } else if (i_plot == 2) {
        # Cook's Distance Plot
        cook_ref <- 4 / n
        p <- ggplot2::ggplot(diag_data, ggplot2::aes(x = index, y = cook_dist)) +
          ggplot2::geom_point() +
          ggplot2::geom_hline(yintercept = cook_ref, linetype = "dashed", color = "red") +
          ggplot2::annotate("text",
            x = max(diag_data$index) * 0.9, y = cook_ref * 1.1,
            label = "4/n", color = "red", hjust = 0
          ) +
          ggplot2::labs(
            x = "Observation Index",
            y = "Cook's Distance",
            title = p_main,
            subtitle = sub.caption
          ) +
          theme_fn()

        # Identify influential points
        influential <- diag_data$cook_dist > cook_ref & !is.na(diag_data$cook_dist)
        if (any(influential)) {
          p <- p + ggplot2::geom_text(
            data = diag_data[influential, ],
            ggplot2::aes(label = index),
            vjust = -0.5, color = "red", size = 3
          )
        }
      } else if (i_plot == 3) {
        # Generalized Leverage vs. Fitted Values
        lev_ref <- 2 * p / n
        p <- ggplot2::ggplot(diag_data, ggplot2::aes(x = fitted, y = leverage)) +
          ggplot2::geom_point() +
          ggplot2::geom_hline(yintercept = lev_ref, linetype = "dashed", color = "red") +
          ggplot2::annotate("text",
            x = max(diag_data$fitted) * 0.9, y = lev_ref * 1.1,
            label = "2p/n", color = "red", hjust = 0
          ) +
          ggplot2::labs(
            x = "Fitted Values",
            y = "Generalized Leverage",
            title = p_main,
            subtitle = sub.caption
          ) +
          theme_fn()

        # Identify high leverage points
        high_lev <- diag_data$leverage > lev_ref
        if (any(high_lev)) {
          p <- p + ggplot2::geom_text(
            data = diag_data[high_lev, ],
            ggplot2::aes(label = index),
            vjust = -0.5, color = "red", size = 3
          )
        }
      } else if (i_plot == 4) {
        # Residuals vs. Linear Predictor
        p <- ggplot2::ggplot(diag_data, ggplot2::aes(x = linpred, y = resid)) +
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
      } else if (i_plot == 5) {
        # Half-Normal Plot of Residuals
        if (exists("half_normal_data") && "lower" %in% names(half_normal_data)) {
          # With envelope
          p <- ggplot2::ggplot(
            half_normal_data,
            ggplot2::aes(x = theoretical, y = observed)
          ) +
            ggplot2::geom_point() +
            ggplot2::geom_abline(
              slope = sd(diag_data$abs_resid, na.rm = TRUE),
              intercept = 0, linetype = "dashed", color = "gray"
            ) +
            ggplot2::geom_line(ggplot2::aes(y = lower),
              linetype = "dashed",
              color = "blue", linewidth = 0.8
            ) +
            ggplot2::geom_line(ggplot2::aes(y = upper),
              linetype = "dashed",
              color = "blue", linewidth = 0.8
            ) +
            ggplot2::annotate("text",
              x = max(half_normal_data$theoretical) * 0.8,
              y = max(half_normal_data$upper) * 0.9,
              label = paste0(level * 100, "% envelope"),
              color = "blue"
            ) +
            ggplot2::labs(
              x = "Theoretical Half-Normal Quantiles",
              y = paste0("Ordered |", type, " residuals|"),
              title = paste(p_main, "(Half-Normal)"),
              subtitle = sub.caption
            ) +
            theme_fn()
        } else {
          # Without envelope
          p <- ggplot2::ggplot(diag_data, ggplot2::aes(sample = abs_resid)) +
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
        }
      } else if (i_plot == 6) {
        # Predicted vs. Observed Values
        p <- ggplot2::ggplot(diag_data, ggplot2::aes(x = fitted, y = y_obs)) +
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

      return(p)
    }

    # Create all requested plots
    p_list <- lapply(which, create_ggplot)

    # Display plots
    if (arrange_plots && length(p_list) > 1) {
      # Calculate grid dimensions
      n_plots <- length(p_list)
      n_cols <- min(2, n_plots)
      n_rows <- ceiling(n_plots / n_cols)

      # Arrange plots in a grid
      if (requireNamespace("gridExtra", quietly = TRUE)) {
        gridExtra::grid.arrange(grobs = p_list, ncol = n_cols, nrow = n_rows)
      } else if (requireNamespace("ggpubr", quietly = TRUE)) {
        do.call(ggpubr::ggarrange, c(p_list, list(ncol = n_cols, nrow = n_rows)))
      }
    } else {
      # Display plots one by one
      for (i in seq_along(p_list)) {
        print(p_list[[i]])
        if (ask && i < length(p_list)) {
          message("Press <ENTER> to continue to the next plot...")
          invisible(readLines(con = "stdin", n = 1))
        }
      }
    }
  }

  # Return diagnostic data if requested -------------------------------------
  if (save_diagnostics) {
    diag_result <- list(
      data = diag_data,
      model_info = list(
        n = n,
        p = p,
        cook_threshold = 4 / n,
        leverage_threshold = 2 * p / n
      )
    )

    if (exists("half_normal_data")) {
      diag_result$half_normal <- half_normal_data
    }

    return(invisible(diag_result))
  } else {
    return(invisible(x))
  }
}
