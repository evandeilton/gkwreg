#' @title Extract Log-Likelihood from Generalized Kumaraswamy Regression Models
#'
#' @description
#' Extracts the log-likelihood value from a fitted Generalized Kumaraswamy (GKw)
#' regression model object.
#'
#' @param object An object of class \code{"gkwreg"}, typically obtained from
#'   \code{\link{gkwreg}}.
#' @param ... Currently not used.
#'
#' @details
#' The log-likelihood is extracted from the fitted model object and returned as
#' an object of class \code{"logLik"} with appropriate attributes for the number
#' of parameters (\code{df}) and observations (\code{nobs}). These attributes
#' are required for information criteria calculations.
#'
#' For a GKw regression model with parameter vector \eqn{\theta}, the log-likelihood
#' is defined as:
#' \deqn{\ell(\theta \mid y) = \sum_{i=1}^n \log f(y_i; \alpha_i, \beta_i, \gamma_i, \delta_i, \lambda_i)}
#' where \eqn{f(\cdot)} is the probability density function of the specified GKw
#' family distribution, and the parameters may depend on covariates through link
#' functions.
#'
#' @return An object of class \code{"logLik"} containing the log-likelihood value
#'   with the following attributes:
#'   \describe{
#'     \item{\code{df}}{Number of estimated parameters}
#'     \item{\code{nobs}}{Number of observations}
#'   }
#'
#' @author Lopes, J. E.
#'
#' @seealso \code{\link{gkwreg}}, \code{\link{AIC.gkwreg}}, \code{\link{BIC.gkwreg}}
#'
#' @examples
#' \donttest{
#' # Load example data
#' data(GasolineYield)
#'
#' # Fit a Kumaraswamy regression model
#' fit <- gkwreg(yield ~ batch + temp, data = GasolineYield, family = "kw")
#'
#' # Extract log-likelihood
#' ll <- logLik(fit)
#' print(ll)
#'
#' # Access attributes
#' cat("Log-likelihood:", as.numeric(ll), "\n")
#' cat("Parameters:", attr(ll, "df"), "\n")
#' cat("Observations:", attr(ll, "nobs"), "\n")
#' }
#'
#' @importFrom stats logLik
#' @method logLik gkwreg
#' @export
logLik.gkwreg <- function(object, ...) {
  if (!inherits(object, "gkwreg")) {
    stop("'object' must be of class 'gkwreg'", call. = FALSE)
  }

  # Extract log-likelihood value
  val <- object$loglik
  if (is.null(val) || !is.finite(val)) {
    stop("log-likelihood not available in fitted model object", call. = FALSE)
  }

  # Extract number of parameters
  df <- object$npar
  if (is.null(df)) {
    df <- length(object$coefficients)
  }

  # Extract number of observations
  nobs <- object$nobs
  if (is.null(nobs)) {
    nobs <- length(object$y)
  }

  # Create logLik object with attributes
  structure(
    val,
    df = as.integer(df),
    nobs = as.integer(nobs),
    class = "logLik"
  )
}


#' @title Akaike Information Criterion for GKw Regression Models
#'
#' @description
#' Calculates the Akaike Information Criterion (AIC) for fitted Generalized
#' Kumaraswamy regression models.
#'
#' @param object An object of class \code{"gkwreg"}, typically obtained from
#'   \code{\link{gkwreg}}.
#' @param ... Optionally more fitted model objects.
#' @param k Numeric, the penalty per parameter. Default is \code{k = 2} for
#'   classical AIC. Setting \code{k = log(n)} gives BIC-equivalent penalty.
#'
#' @details
#' The AIC is computed as:
#' \deqn{AIC = -2\ell(\hat{\theta}) + k \cdot p}
#' where \eqn{\ell(\hat{\theta})} is the maximized log-likelihood and \eqn{p}
#' is the number of estimated parameters.
#'
#' When multiple objects are provided, a data frame comparing all models is
#' returned. Lower AIC values indicate better models, balancing goodness-of-fit
#' against model complexity.
#'
#' For small sample sizes, consider the corrected AIC (AICc):
#' \deqn{AICc = AIC + \frac{2p(p+1)}{n-p-1}}
#' where \eqn{n} is the sample size. This correction is not automatically applied
#' but can be calculated manually.
#'
#' @return If only one object is provided, returns a numeric value with the AIC.
#'   If multiple objects are provided, returns a data frame with columns \code{df}
#'   and \code{AIC}, with rows named according to the object names in the call.
#'
#' @author Lopes, J. E.
#'
#' @references
#' Akaike, H. (1974). A new look at the statistical model identification.
#' \emph{IEEE Transactions on Automatic Control}, \strong{19}(6), 716--723.
#' \doi{10.1109/TAC.1974.1100705}
#'
#' Burnham, K. P., & Anderson, D. R. (2004). Multimodel inference: Understanding
#' AIC and BIC in model selection. \emph{Sociological Methods & Research},
#' \strong{33}(2), 261--304. \doi{10.1177/0049124104268644}
#'
#' @seealso \code{\link{gkwreg}}, \code{\link{logLik.gkwreg}}, \code{\link{BIC.gkwreg}}
#'
#' @examples
#' \donttest{
#' # Load example data
#' data(GasolineYield)
#'
#' # Fit competing models
#' fit1 <- gkwreg(yield ~ batch, data = GasolineYield, family = "kw")
#' fit2 <- gkwreg(yield ~ batch + temp, data = GasolineYield, family = "kw")
#' fit3 <- gkwreg(yield ~ temp, data = GasolineYield, family = "kw")
#'
#' # Calculate AIC for single model
#' AIC(fit1)
#'
#' # Compare multiple models (with proper names)
#' AIC(fit1, fit2, fit3)
#'
#' # Use different penalty
#' AIC(fit1, k = 4)
#' }
#'
#' @importFrom stats AIC logLik
#' @method AIC gkwreg
#' @export
AIC.gkwreg <- function(object, ..., k = 2) {
  if (!inherits(object, "gkwreg")) {
    stop("'object' must be of class 'gkwreg'", call. = FALSE)
  }

  # Handle multiple objects
  dots <- list(...)
  if (length(dots) > 0L) {
    # Get the names of the objects from the call
    object_names <- as.character(match.call(expand.dots = FALSE)$...)

    # If names are not informative (e.g., from do.call), try deparse
    if (length(object_names) == 0L || all(object_names == "")) {
      object_names <- paste0("Model", seq_along(dots))
    }

    # Get the name of the first object
    first_name <- as.character(substitute(object))
    if (length(first_name) > 1L) {
      first_name <- deparse(substitute(object))
    }

    # Combine all objects
    all_objects <- c(list(object), dots)
    all_names <- c(first_name, object_names)

    # Calculate AIC for each object
    aic_vals <- vapply(all_objects, function(obj) {
      if (inherits(obj, "gkwreg")) {
        ll <- logLik(obj)
        df <- attr(ll, "df")
        if (is.null(df) || is.na(df)) {
          return(NA_real_)
        }
        return(-2 * as.numeric(ll) + k * df)
      } else {
        # For non-gkwreg objects, use stats::AIC
        return(stats::AIC(obj, k = k))
      }
    }, FUN.VALUE = numeric(1L))

    # Get degrees of freedom for each object
    df_vals <- vapply(all_objects, function(obj) {
      if (inherits(obj, "gkwreg")) {
        ll <- logLik(obj)
        df <- attr(ll, "df")
        return(if (is.null(df) || is.na(df)) NA_integer_ else as.integer(df))
      } else {
        ll <- logLik(obj)
        df <- attr(ll, "df")
        return(if (is.null(df) || is.na(df)) NA_integer_ else as.integer(df))
      }
    }, FUN.VALUE = integer(1L))

    # Create result data frame
    result <- data.frame(
      df = df_vals,
      AIC = aic_vals,
      row.names = all_names,
      stringsAsFactors = FALSE
    )

    return(result)
  }

  # Single object: use stored value if k = 2, otherwise compute
  if (identical(k, 2) && !is.null(object$aic)) {
    return(object$aic)
  }

  # Compute AIC from log-likelihood
  ll <- logLik(object)
  df <- attr(ll, "df")

  if (is.null(df) || is.na(df)) {
    stop("number of parameters not available", call. = FALSE)
  }

  return(-2 * as.numeric(ll) + k * df)
}


#' @title Bayesian Information Criterion for GKw Regression Models
#'
#' @description
#' Calculates the Bayesian Information Criterion (BIC), also known as the
#' Schwarz Information Criterion (SIC), for fitted Generalized Kumaraswamy
#' regression models.
#'
#' @param object An object of class \code{"gkwreg"}, typically obtained from
#'   \code{\link{gkwreg}}.
#' @param ... Optionally more fitted model objects.
#'
#' @details
#' The BIC is computed as:
#' \deqn{BIC = -2\ell(\hat{\theta}) + p \cdot \log(n)}
#' where \eqn{\ell(\hat{\theta})} is the maximized log-likelihood, \eqn{p} is
#' the number of estimated parameters, and \eqn{n} is the sample size.
#'
#' When multiple objects are provided, a data frame comparing all models is
#' returned. Lower BIC values indicate better models. BIC penalizes model
#' complexity more heavily than AIC, particularly for large samples, and tends
#' to favor more parsimonious models.
#'
#' The BIC can be derived from a Bayesian perspective as an approximation to
#' the logarithm of the Bayes factor, under certain regularity conditions and
#' assuming uniform priors.
#'
#' @return If only one object is provided, returns a numeric value with the BIC.
#'   If multiple objects are provided, returns a data frame with columns \code{df}
#'   and \code{BIC}, with rows named according to the object names in the call.
#'
#' @author Lopes, J. E.
#'
#' @references
#' Schwarz, G. (1978). Estimating the dimension of a model.
#' \emph{The Annals of Statistics}, \strong{6}(2), 461--464.
#' \doi{10.1214/aos/1176344136}
#'
#' @seealso \code{\link{gkwreg}}, \code{\link{logLik.gkwreg}}, \code{\link{AIC.gkwreg}}
#'
#' @examples
#' \donttest{
#' # Load example data
#' data(GasolineYield)
#'
#' # Fit competing models
#' fit1 <- gkwreg(yield ~ batch, data = GasolineYield, family = "kw")
#' fit2 <- gkwreg(yield ~ batch + temp, data = GasolineYield, family = "kw")
#' fit3 <- gkwreg(yield ~ temp, data = GasolineYield, family = "kw")
#'
#' # Calculate BIC for single model
#' BIC(fit1)
#'
#' # Compare multiple models (with proper names)
#' BIC(fit1, fit2, fit3)
#' }
#'
#' @importFrom stats BIC logLik
#' @method BIC gkwreg
#' @export
BIC.gkwreg <- function(object, ...) {
  if (!inherits(object, "gkwreg")) {
    stop("'object' must be of class 'gkwreg'", call. = FALSE)
  }

  # Handle multiple objects
  dots <- list(...)
  if (length(dots) > 0L) {
    # Get the names of the objects from the call
    object_names <- as.character(match.call(expand.dots = FALSE)$...)

    # If names are not informative (e.g., from do.call), try deparse
    if (length(object_names) == 0L || all(object_names == "")) {
      object_names <- paste0("Model", seq_along(dots))
    }

    # Get the name of the first object
    first_name <- as.character(substitute(object))
    if (length(first_name) > 1L) {
      first_name <- deparse(substitute(object))
    }

    # Combine all objects
    all_objects <- c(list(object), dots)
    all_names <- c(first_name, object_names)

    # Calculate BIC for each object
    bic_vals <- vapply(all_objects, function(obj) {
      if (inherits(obj, "gkwreg")) {
        ll <- logLik(obj)
        df <- attr(ll, "df")
        nobs <- attr(ll, "nobs")
        if (is.null(df) || is.na(df) || is.null(nobs) || is.na(nobs)) {
          return(NA_real_)
        }
        return(-2 * as.numeric(ll) + df * log(nobs))
      } else {
        # For non-gkwreg objects, use stats::BIC
        return(stats::BIC(obj))
      }
    }, FUN.VALUE = numeric(1L))

    # Get degrees of freedom for each object
    df_vals <- vapply(all_objects, function(obj) {
      if (inherits(obj, "gkwreg")) {
        ll <- logLik(obj)
        df <- attr(ll, "df")
        return(if (is.null(df) || is.na(df)) NA_integer_ else as.integer(df))
      } else {
        ll <- logLik(obj)
        df <- attr(ll, "df")
        return(if (is.null(df) || is.na(df)) NA_integer_ else as.integer(df))
      }
    }, FUN.VALUE = integer(1L))

    # Create result data frame
    result <- data.frame(
      df = df_vals,
      BIC = bic_vals,
      row.names = all_names,
      stringsAsFactors = FALSE
    )

    return(result)
  }

  # Single object: use stored value if available, otherwise compute
  if (!is.null(object$bic)) {
    return(object$bic)
  }

  # Compute BIC from log-likelihood
  ll <- logLik(object)
  df <- attr(ll, "df")
  nobs <- attr(ll, "nobs")

  if (is.null(df) || is.na(df)) {
    stop("number of parameters not available", call. = FALSE)
  }

  if (is.null(nobs) || is.na(nobs) || nobs <= 0) {
    stop("number of observations not available", call. = FALSE)
  }

  return(-2 * as.numeric(ll) + df * log(nobs))
}












#'
#' #' @title Extract Log-Likelihood from a Generalized Kumaraswamy Regression Model
#' #'
#' #' @description
#' #' This function extracts the maximized log-likelihood value from a fitted Generalized
#' #' Kumaraswamy (GKw) regression model object (class \code{"gkwreg"}). The result is
#' #' returned as an object of class \code{"logLik"}, which includes attributes for
#' #' degrees of freedom and number of observations, suitable for use with model
#' #' selection criteria like AIC and BIC.
#' #'
#' #' @param object An object of class \code{"gkwreg"}, typically the result of a call
#' #'   to \code{\link{gkwreg}}.
#' #' @param ... Additional arguments, currently ignored by this method.
#' #'
#' #' @details
#' #' The log-likelihood value is typically computed during the model fitting process
#' #' (e.g., by \code{\link{gkwreg}}) and stored within the resulting object. This
#' #' method retrieves this stored value. If the value is not directly available, it
#' #' attempts to calculate it from the stored deviance (\eqn{logLik = -deviance / 2}).
#' #'
#' #' The log-likelihood for a GKw family model with parameters \eqn{\theta} is
#' #' generally defined as the sum of the log-density contributions for each observation:
#' #' \deqn{l(\theta | y) = \sum_{i=1}^n \log f(y_i; \alpha_i, \beta_i, \gamma_i, \delta_i, \lambda_i)}
#' #' where \eqn{f(y; \dots)} is the probability density function (PDF) of the specific
#' #' distribution from the GKw family used in the model (determined by the \code{family}
#' #' argument in \code{gkwreg}), and parameters (\eqn{\alpha_i, \dots, \lambda_i}) may
#' #' depend on covariates.
#' #'
#' #' The function also extracts the number of estimated parameters (\code{df}) and the
#' #' number of observations (\code{nobs}) used in the fit, storing them as attributes
#' #' of the returned \code{"logLik"} object, which is essential for functions like
#' #' \code{\link[stats]{AIC}} and \code{\link[stats]{BIC}}. It attempts to find \code{df}
#' #' and \code{nobs} from various components within the \code{object} if they are not
#' #' directly stored as \code{npar} and \code{nobs}.
#' #'
#' #' @return An object of class \code{"logLik"} representing the maximized
#' #'   log-likelihood value. It has the following attributes:
#' #'   \itemize{
#' #'     \item \code{df}: (numeric) The number of estimated parameters in the model
#' #'       (coefficients).
#' #'     \item \code{nobs}: (numeric) The number of observations used for fitting the model.
#' #'   }
#' #'
#' #' @author Lopes, J. E.
#' #'
#' #' @seealso \code{\link{gkwreg}}, \code{\link{AIC.gkwreg}}, \code{\link{BIC.gkwreg}},
#' #'   \code{\link[stats]{logLik}}, \code{\link[stats]{AIC}}, \code{\link[stats]{BIC}}
#' #'
#' #' @keywords log-likelihood models likelihood
#' #'
#' #' @examples
#' #' \donttest{
#' #' # Assume 'df' exists with response 'y' and predictors 'x1', 'x2', 'x3'
#' #' # and that rkw() is available and data is appropriate (0 < y < 1).
#' #' set.seed(123)
#' #' n <- 100
#' #' x1 <- runif(n)
#' #' x2 <- rnorm(n)
#' #' x3 <- factor(rbinom(n, 1, 0.4))
#' #' alpha <- exp(0.5 + 0.2 * x1)
#' #' beta <- exp(1.0 - 0.1 * x2 + 0.3 * (x3 == "1"))
#' #' y <- rkw(n, alpha = alpha, beta = beta) # Placeholder if rkw not available
#' #' y <- pmax(pmin(y, 1 - 1e-7), 1e-7)
#' #' df <- data.frame(y = y, x1 = x1, x2 = x2, x3 = x3)
#' #'
#' #' # Fit a Kumaraswamy regression model
#' #' kw_reg <- gkwreg(y ~ x1 | x2 + x3, data = df, family = "kw")
#' #'
#' #' # Extract log-likelihood object
#' #' ll <- logLik(kw_reg)
#' #'
#' #' # Print the log-likelihood value (with attributes)
#' #' print(ll)
#' #'
#' #' # Access the value directly
#' #' ll_value <- as.numeric(ll)
#' #' print(ll_value)
#' #'
#' #' # Get the number of parameters (degrees of freedom)
#' #' df_model <- attr(ll, "df")
#' #' print(paste("Number of parameters:", df_model))
#' #'
#' #' # Get the number of observations
#' #' nobs_model <- attr(ll, "nobs")
#' #' print(paste("Number of observations:", nobs_model))
#' #'
#' #' # Use with AIC/BIC
#' #' AIC(kw_reg)
#' #' BIC(kw_reg)
#' #' }
#' #'
#' #' @importFrom stats logLik
#' #' @method logLik gkwreg
#' #' @export
#' logLik.gkwreg <- function(object, ...) {
#'   # Check if the object is of class gkwreg
#'   if (!inherits(object, "gkwreg")) {
#'     stop("'object' must be a fitted model object of class 'gkwreg'")
#'   }
#'
#'   # Extract log-likelihood value
#'   ll <- object$loglik
#'
#'   # If log-likelihood is not available, try to recover from other information
#'   if (is.null(ll)) {
#'     # Try to calculate from deviance if available
#'     if (!is.null(object$deviance)) {
#'       ll <- -object$deviance / 2
#'     } else {
#'       warning("Log-likelihood not found in the model object and cannot be calculated")
#'       ll <- NA_real_
#'     }
#'   }
#'
#'   # Get the number of parameters (degrees of freedom for the model)
#'   # Use npar if available, otherwise count coefficients
#'   df <- object$npar
#'   if (is.null(df)) {
#'     # Try to determine number of parameters from coefficients
#'     if (!is.null(object$coefficients)) {
#'       df <- length(object$coefficients)
#'     } else {
#'       warning("Number of parameters ('npar') not found in the model object")
#'       df <- NA_integer_
#'     }
#'   }
#'
#'   # Get the number of observations used in the fit
#'   # Use nobs if available, otherwise infer from other components
#'   nobs <- object$nobs
#'   if (is.null(nobs)) {
#'     # Try to determine number of observations from residuals or fitted values or y
#'     if (!is.null(object$residuals)) {
#'       nobs <- length(object$residuals)
#'     } else if (!is.null(object$fitted.values)) {
#'       nobs <- length(object$fitted.values)
#'     } else if (!is.null(object$y)) {
#'       nobs <- length(object$y)
#'     } else {
#'       warning("Number of observations ('nobs') not found in the model object")
#'       nobs <- NA_integer_
#'     }
#'   }
#'
#'   # Ensure df and nobs are numeric, even if NA
#'   df <- as.numeric(df)
#'   nobs <- as.numeric(nobs)
#'
#'   # Create and return the logLik object with appropriate attributes
#'   # Use structure() to assign attributes and class simultaneously
#'   structure(ll,
#'     df = df,
#'     nobs = nobs,
#'     class = "logLik"
#'   )
#' }
#'
#'
#'
#' #' @title Akaike's Information Criterion for GKw Regression Models
#' #'
#' #' @description
#' #' Calculates the Akaike Information Criterion (AIC) for one or more fitted
#' #' Generalized Kumaraswamy (GKw) regression model objects (class \code{"gkwreg"}).
#' #' AIC is commonly used for model selection, penalizing model complexity.
#' #'
#' #' @param object An object of class \code{"gkwreg"}, typically the result of a
#' #'   call to \code{\link{gkwreg}}.
#' #' @param ... Optionally, one or more additional fitted model objects of class
#' #'   \code{"gkwreg"}, for which AIC should also be calculated.
#' #' @param k Numeric, the penalty per parameter. The default \code{k = 2} corresponds
#' #'   to the traditional AIC. Using \code{k = log(nobs)} would yield BIC (though using
#' #'   \code{\link{BIC.gkwreg}} is preferred for that).
#' #'
#' #' @details
#' #' The AIC is calculated based on the maximized log-likelihood (\eqn{L}) and the
#' #' number of estimated parameters (\eqn{p}) in the model:
#' #' \deqn{AIC = -2 \log(L) + k \times p}
#' #' This function retrieves the log-likelihood and the number of parameters (\code{df})
#' #' using the \code{\link{logLik.gkwreg}} method for the fitted \code{gkwreg} object(s).
#' #' Models with lower AIC values are generally preferred, as they indicate a better
#' #' balance between goodness of fit and model parsimony.
#' #'
#' #' When comparing multiple models passed via \code{...}, the function relies on
#' #' \code{\link[stats]{AIC}}'s default method for creating a comparison table,
#' #' which in turn calls \code{logLik} for each provided object.
#' #'
#' #' For small sample sizes relative to the number of parameters, the second-order
#' #' AIC (AICc) might be more appropriate:
#' #' \deqn{AICc = AIC + \frac{2p(p+1)}{n-p-1}}
#' #' where \eqn{n} is the number of observations. AICc is not directly computed by
#' #' this function but can be calculated manually using the returned AIC, \eqn{p}
#' #' (from \code{attr(logLik(object), "df")}), and \eqn{n}
#' #' (from \code{attr(logLik(object), "nobs")}).
#' #'
#' #' @return If just one \code{object} is provided, returns a single numeric AIC value.
#' #'   If multiple objects are provided via \code{...}, returns a \code{data.frame}
#' #'   with rows corresponding to the models and columns for the degrees of freedom
#' #'   (\code{df}) and the AIC values, sorted by AIC.
#' #'
#' #' @author Lopes, J. E.
#' #'
#' #' @references
#' #' Akaike, H. (1974). A new look at the statistical model identification.
#' #' \emph{IEEE Transactions on Automatic Control}, \strong{19}(6), 716-723.
#' #'
#' #'
#' #' Burnham, K. P., & Anderson, D. R. (2002). \emph{Model Selection and Multimodel
#' #' Inference: A Practical Information-Theoretic Approach} (2nd ed.). Springer-Verlag.
#' #'
#' #' @seealso \code{\link{gkwreg}}, \code{\link{logLik.gkwreg}}, \code{\link{BIC.gkwreg}},
#' #'   \code{\link[stats]{AIC}}
#' #'
#' #' @keywords AIC models likelihood model selection
#' #'
#' #' @examples
#' #' \donttest{
#' #' # Assume 'df' exists with response 'y' and predictors 'x1', 'x2', 'x3'
#' #' # and that rkw() is available and data is appropriate (0 < y < 1).
#' #' set.seed(123)
#' #' n <- 100
#' #' x1 <- runif(n)
#' #' x2 <- rnorm(n)
#' #' x3 <- factor(rbinom(n, 1, 0.4))
#' #' alpha <- exp(0.5 + 0.2 * x1)
#' #' beta <- exp(1.0 - 0.1 * x2 + 0.3 * (x3 == "1"))
#' #' y <- rkw(n, alpha = alpha, beta = beta) # Placeholder if rkw not available
#' #' y <- pmax(pmin(y, 1 - 1e-7), 1e-7)
#' #' df <- data.frame(y = y, x1 = x1, x2 = x2, x3 = x3)
#' #'
#' #' # Fit two competing models
#' #' kw_reg1 <- gkwreg(y ~ x1 | x2, data = df, family = "kw")
#' #' kw_reg2 <- gkwreg(y ~ x1 | x2 + x3, data = df, family = "kw") # More complex beta model
#' #' kw_reg3 <- gkwreg(y ~ 1 | x2 + x3, data = df, family = "kw") # Simpler alpha model
#' #'
#' #' # Calculate AIC for a single model
#' #' aic1 <- AIC(kw_reg1)
#' #' print(aic1)
#' #'
#' #' # Compare models using AIC (lower is better)
#' #' model_comparison_aic <- c(AIC(kw_reg1), AIC(kw_reg2), AIC(kw_reg3))
#' #' print(model_comparison_aic)
#' #'
#' #' # Calculate AIC with a different penalty (e.g., k=4)
#' #' aic1_k4 <- AIC(kw_reg1, k = 4)
#' #' print(aic1_k4)
#' #' }
#' #'
#' #' @importFrom stats AIC
#' #' @method AIC gkwreg
#' #' @export
#' AIC.gkwreg <- function(object, ..., k = 2) {
#'   # Check if the object is of class gkwreg
#'   if (!inherits(object, "gkwreg")) {
#'     stop("'object' must be a fitted model object of class 'gkwreg'")
#'   }
#'
#'   # Capture additional model objects
#'   dot_objects <- list(...)
#'
#'   # Handle case with multiple models for comparison using stats::AIC generic logic
#'   # This relies on the generic dispatching to logLik.gkwreg for each object
#'   if (length(dot_objects) > 0) {
#'     # Ensure all additional objects are also gkwreg models (or compatible)
#'     obj_list <- c(list(object), dot_objects)
#'     classes <- vapply(obj_list, function(o) class(o)[1], character(1))
#'     if (!all(classes == "gkwreg")) {
#'       # Allow comparison if logLik methods exist for other object types
#'       # stats::AIC handles this, just pass them through
#'       warning("Comparing objects of different classes.")
#'     }
#'     # Call the default stats::AIC logic which handles multiple objects
#'     # It uses logLik() on each object.
#'     return(stats::AIC(object = object, ..., k = k))
#'   }
#'
#'   # --- Handle single object case ---
#'
#'   # Check if AIC is already computed and stored in the object AND k is the default 2
#'   # Avoid recalculating standard AIC if already present
#'   if (!is.null(object$aic) && identical(k, 2)) {
#'     return(object$aic)
#'   }
#'
#'   # Calculate AIC from log-likelihood and number of parameters
#'   ll <- stats::logLik(object) # Use stats::logLik generic, dispatches to logLik.gkwreg
#'
#'   # Extract number of parameters (df) from the logLik object
#'   df <- attr(ll, "df")
#'
#'   # Check if df is valid
#'   if (is.null(df) || is.na(df) || !is.numeric(df) || df < 0) {
#'     warning("Could not extract a valid number of parameters (df) from the logLik object. AIC calculation might be incorrect.")
#'     # Attempt fallback: count coefficients if df is invalid
#'     if (!is.null(object$coefficients)) {
#'       df <- length(object$coefficients)
#'       warning("Using the count of coefficients (", df, ") as the number of parameters.")
#'     } else {
#'       df <- NA_real_ # Cannot determine df
#'       warning("Setting number of parameters to NA.")
#'     }
#'   }
#'
#'   # Check if logLik value is valid
#'   ll_val <- as.numeric(ll)
#'   if (is.null(ll_val) || is.na(ll_val) || !is.finite(ll_val)) {
#'     warning("Invalid log-likelihood value extracted. Cannot compute AIC.")
#'     return(NA_real_)
#'   }
#'
#'   # Calculate and return AIC using the formula -2*logLik + k*df
#'   aic_val <- -2 * ll_val + k * df
#'
#'   return(aic_val)
#' }
#'
#'
#'
#' #' @title Bayesian Information Criterion for GKw Regression Models
#' #'
#' #' @description
#' #' Calculates the Bayesian Information Criterion (BIC), also known as Schwarz's
#' #' Bayesian Criterion (SBC), for one or more fitted Generalized Kumaraswamy (GKw)
#' #' regression model objects (class \code{"gkwreg"}). BIC is used for model selection
#' #' and tends to penalize model complexity more heavily than AIC, especially for
#' #' larger datasets.
#' #'
#' #' @param object An object of class \code{"gkwreg"}, typically the result of a
#' #'   call to \code{\link{gkwreg}}.
#' #' @param ... Optionally, one or more additional fitted model objects of class
#' #'   \code{"gkwreg"}, for which BIC should also be calculated.
#' #'
#' #' @details
#' #' The BIC is calculated based on the maximized log-likelihood (\eqn{L}), the
#' #' number of estimated parameters (\eqn{p}) in the model, and the number of
#' #' observations (\eqn{n}):
#' #' \deqn{BIC = -2 \log(L) + p \times \log(n)}
#' #' This function retrieves the log-likelihood, the number of parameters (\code{df}),
#' #' and the number of observations (\code{nobs}) using the \code{\link{logLik.gkwreg}}
#' #' method for the fitted \code{gkwreg} object(s).
#' #'
#' #' Models with lower BIC values are generally preferred. The penalty term \eqn{p \log(n)}
#' #' increases more rapidly with sample size \eqn{n} compared to AIC's penalty \eqn{2p},
#' #' meaning BIC favors simpler models more strongly in larger samples. BIC can be
#' #' motivated from a Bayesian perspective as an approximation related to Bayes factors.
#' #'
#' #' When comparing multiple models passed via \code{...}, the function relies on
#' #' \code{\link[stats]{BIC}}'s default method for creating a comparison table,
#' #' which in turn calls \code{logLik} for each provided object.
#' #'
#' #' @return If just one \code{object} is provided, returns a single numeric BIC value.
#' #'   If multiple objects are provided via \code{...}, returns a \code{data.frame}
#' #'   with rows corresponding to the models and columns for the degrees of freedom
#' #'   (\code{df}) and the BIC values, sorted by BIC.
#' #'
#' #' @author Lopes, J. E.
#' #'
#' #' @references
#' #' Schwarz, G. (1978). Estimating the dimension of a model.
#' #' \emph{The Annals of Statistics}, \strong{6}(2), 461-464.
#' #'
#' #'
#' #' @seealso \code{\link{gkwreg}}, \code{\link{logLik.gkwreg}}, \code{\link{AIC.gkwreg}},
#' #'   \code{\link[stats]{BIC}}
#' #'
#' #' @keywords BIC models likelihood model selection
#' #'
#' #' @examples
#' #' \donttest{
#' #' # Assume 'df' exists with response 'y' and predictors 'x1', 'x2', 'x3'
#' #' # and that rkw() is available and data is appropriate (0 < y < 1).
#' #' set.seed(123)
#' #' n <- 100
#' #' x1 <- runif(n)
#' #' x2 <- rnorm(n)
#' #' x3 <- factor(rbinom(n, 1, 0.4))
#' #' alpha <- exp(0.5 + 0.2 * x1)
#' #' beta <- exp(1.0 - 0.1 * x2 + 0.3 * (x3 == "1"))
#' #' y <- rkw(n, alpha = alpha, beta = beta) # Placeholder if rkw not available
#' #' y <- pmax(pmin(y, 1 - 1e-7), 1e-7)
#' #' df <- data.frame(y = y, x1 = x1, x2 = x2, x3 = x3)
#' #'
#' #' # Fit two competing models
#' #' kw_reg1 <- gkwreg(y ~ x1 | x2, data = df, family = "kw")
#' #' kw_reg2 <- gkwreg(y ~ x1 | x2 + x3, data = df, family = "kw") # More complex beta model
#' #' kw_reg3 <- gkwreg(y ~ 1 | x2 + x3, data = df, family = "kw") # Simpler alpha model
#' #'
#' #' # Calculate BIC for a single model
#' #' bic1 <- BIC(kw_reg1)
#' #' print(bic1)
#' #'
#' #' # Compare models using BIC (lower is better)
#' #' model_comparison_bic <- c(BIC(kw_reg1), BIC(kw_reg2), BIC(kw_reg3))
#' #' print(model_comparison_bic)
#' #' }
#' #'
#' #' @importFrom stats BIC
#' #' @method BIC gkwreg
#' #' @export
#' BIC.gkwreg <- function(object, ...) {
#'   # Check if the object is of class gkwreg
#'   if (!inherits(object, "gkwreg")) {
#'     stop("'object' must be a fitted model object of class 'gkwreg'")
#'   }
#'
#'   # Capture additional model objects
#'   dot_objects <- list(...)
#'
#'   # Handle case with multiple models for comparison using stats::BIC generic logic
#'   if (length(dot_objects) > 0) {
#'     # Ensure all additional objects are also gkwreg models (or compatible)
#'     obj_list <- c(list(object), dot_objects)
#'     classes <- vapply(obj_list, function(o) class(o)[1], character(1))
#'     if (!all(classes == "gkwreg")) {
#'       # Allow comparison if logLik methods exist for other object types
#'       # stats::BIC handles this, just pass them through
#'       warning("Comparing objects of different classes.")
#'     }
#'     # Call the default stats::BIC logic which handles multiple objects
#'     # It uses logLik() on each object.
#'     return(stats::BIC(object = object, ...))
#'   }
#'
#'   # --- Handle single object case ---
#'
#'   # Check if BIC is already computed and stored in the object
#'   # Avoid recalculating if already present
#'   if (!is.null(object$bic)) {
#'     return(object$bic)
#'   }
#'
#'   # Calculate BIC from log-likelihood, number of parameters, and number of observations
#'   ll <- stats::logLik(object) # Use stats::logLik generic, dispatches to logLik.gkwreg
#'
#'   # Extract number of parameters (df) from the logLik object
#'   df <- attr(ll, "df")
#'
#'   # Check if df is valid
#'   if (is.null(df) || is.na(df) || !is.numeric(df) || df < 0) {
#'     warning("Could not extract a valid number of parameters (df) from the logLik object. BIC calculation might be incorrect.")
#'     # Attempt fallback: count coefficients if df is invalid
#'     if (!is.null(object$coefficients)) {
#'       df <- length(object$coefficients)
#'       warning("Using the count of coefficients (", df, ") as the number of parameters.")
#'     } else {
#'       df <- NA_real_ # Cannot determine df
#'       warning("Setting number of parameters to NA.")
#'     }
#'   }
#'
#'   # Extract number of observations (nobs) from the logLik object
#'   n <- attr(ll, "nobs")
#'
#'   # Check if nobs is valid
#'   if (is.null(n) || is.na(n) || !is.numeric(n) || n <= 0) {
#'     warning("Could not extract a valid number of observations (nobs) from the logLik object. BIC calculation might be incorrect.")
#'     # Attempt fallback: try various sources from the object
#'     if (!is.null(object$residuals)) {
#'       n <- length(object$residuals)
#'     } else if (!is.null(object$fitted.values)) {
#'       n <- length(object$fitted.values)
#'     } else if (!is.null(object$y)) {
#'       n <- length(object$y)
#'     } else {
#'       n <- NA_real_ # Cannot determine nobs
#'       warning("Setting number of observations to NA.")
#'     }
#'     if (!is.na(n)) warning("Using length of residuals/fitted/y (", n, ") as the number of observations.")
#'   }
#'
#'   # Check if logLik value is valid
#'   ll_val <- as.numeric(ll)
#'   if (is.null(ll_val) || is.na(ll_val) || !is.finite(ll_val)) {
#'     warning("Invalid log-likelihood value extracted. Cannot compute BIC.")
#'     return(NA_real_)
#'   }
#'
#'   # Check if df and n are valid for calculation
#'   if (is.na(df) || is.na(n) || n <= 0) {
#'     warning("Cannot compute BIC due to missing or invalid 'df' or 'nobs'.")
#'     return(NA_real_)
#'   }
#'
#'   # Calculate and return BIC using the formula -2*logLik + df*log(n)
#'   bic_val <- -2 * ll_val + df * log(n)
#'
#'   return(bic_val)
#' }
