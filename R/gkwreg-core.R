#' Fit Generalized Kumaraswamy Regression Models
#'
#' @description
#' Fits regression models for response variables strictly bounded in the (0, 1)
#' interval using the Generalized Kumaraswamy (GKw) family of distributions.
#' It provides a unified interface for fitting the GKw and six nested submodels,
#' allowing flexible regression on multiple distribution parameters.
#'
#' @details
#' Estimation is performed via Maximum Likelihood using automatic differentiation
#' through the TMB (Template Model Builder) package for efficiency and accuracy.
#' The modeling interface uses \code{\link[Formula]{Formula}} syntax, similar to
#' \code{\link[stats]{glm}} and \code{\link[betareg]{betareg}}, allowing different
#' linear predictors for each distributional parameter.
#'
#' \subsection{Distribution Families}{
#' The \code{family} argument selects a distribution from the GKw hierarchy.
#' Simpler families (\code{"beta"}, \code{"kw"}) are often sufficient and
#' computationally faster. More complex families (\code{"bkw"}, \code{"gkw"})
#' add flexibility but risk overfitting. Use \code{\link[stats]{AIC}},
#' \code{\link[stats]{BIC}}, or \code{\link{anova.gkwreg}} to compare nested models.
#' }
#'
#' \subsection{Formula Specification}{
#' The extended formula syntax \code{y ~ model_alpha | model_beta | ...} allows
#' each parameter's linear predictor to be specified. Parameters are in the
#' order \eqn{\alpha, \beta, \gamma, \delta, \lambda}.
#'
#' If a part is omitted (e.g., \code{y ~ model_alpha}), the remaining parameters
#' are modeled as intercept-only. Parts corresponding to fixed parameters
#' (defined by \code{family}) are ignored.
#'
#' \preformatted{
#' # Kw family ("kw") parameters: alpha, beta
#'
#' # alpha ~ x1 + x2, beta ~ x3
#' y ~ x1 + x2 | x3
#'
#' # alpha ~ x1, beta ~ 1 (intercept only)
#' y ~ x1 | 1
#'
#' # alpha ~ x1, beta ~ x1 (same formula for both)
#' y ~ x1
#' }
#' }
#'
#' \subsection{Optimization and Convergence}{
#' Fitting is performed by \code{\link[stats]{optim}} or \code{\link[stats]{nlminb}}.
#' If convergence fails, consider:
#' \itemize{
#'   \item Checking data for separation, outliers, or collinearity.
#'   \item Rescaling predictors.
#'   \item Trying a different optimizer, e.g.,
#'     \code{control = gkw_control(method = "BFGS")}.
#'   \item Simplifying the model (fewer predictors or a simpler \code{family}).
#'   \item Providing starting values via \code{control = gkw_control(start = ...)}.
#' }
#' }
#'
#' @param formula An object of class \code{\link[Formula]{Formula}}. The formula
#'   uses extended syntax to specify different linear predictors for each
#'   parameter:
#'   \code{y ~ model_alpha | model_beta | model_gamma | model_delta | model_lambda}
#'
#'   \describe{
#'     \item{\code{y}}{The response variable, which must be in the (0, 1) interval.}
#'     \item{\code{model_alpha}}{Predictors for the \eqn{\alpha} parameter.}
#'     \item{\code{model_beta}}{Predictors for the \eqn{\beta} parameter.}
#'     \item{\code{model_gamma}}{Predictors for the \eqn{\gamma} parameter.}
#'     \item{\code{model_delta}}{Predictors for the \eqn{\delta} parameter.}
#'     \item{\code{model_lambda}}{Predictors for the \eqn{\lambda} parameter.}
#'   }
#'   See Details for examples.
#'
#' @param data A data frame containing the variables specified in \code{formula}.
#'
#' @param family A character string specifying the distribution family. Must be one of:
#'   \describe{
#'     \item{\code{"gkw"}}{Generalized Kumaraswamy (default). Parameters:
#'       \eqn{\alpha, \beta, \gamma, \delta, \lambda}.}
#'     \item{\code{"bkw"}}{Beta-Kumaraswamy. Parameters:
#'       \eqn{\alpha, \beta, \gamma, \delta} (fixes \eqn{\lambda = 1}).}
#'     \item{\code{"kkw"}}{Kumaraswamy-Kumaraswamy. Parameters:
#'       \eqn{\alpha, \beta, \delta, \lambda} (fixes \eqn{\gamma = 1}).}
#'     \item{\code{"ekw"}}{Exponentiated Kumaraswamy. Parameters:
#'       \eqn{\alpha, \beta, \lambda} (fixes \eqn{\gamma = 1, \delta = 0}).}
#'     \item{\code{"mc"}}{McDonald (Beta Power). Parameters:
#'       \eqn{\gamma, \delta, \lambda} (fixes \eqn{\alpha = 1, \beta = 1}).}
#'     \item{\code{"kw"}}{Kumaraswamy. Parameters: \eqn{\alpha, \beta}
#'       (fixes \eqn{\gamma = 1, \delta = 0, \lambda = 1}).}
#'     \item{\code{"beta"}}{Beta distribution. Parameters: \eqn{\gamma, \delta}
#'       (fixes \eqn{\alpha = 1, \beta = 1, \lambda = 1}). Corresponds to
#'       \code{shape1 = \eqn{\gamma}}, \code{shape2 = \eqn{\delta}}.}
#'   }
#'
#' @param link Specifies the link function(s) for the distributional parameters.
#'   Can be a single character string (applied to all parameters) or a
#'   named list for parameter-specific links (e.g.,
#'   \code{link = list(alpha = "log", delta = "logit")}).
#'
#'   If \code{NULL} (default), links are: \code{"log"} for \eqn{\alpha, \beta,
#'   \gamma, \lambda} and \code{"logit"} for \eqn{\delta}.
#'
#'   Available links:
#'   \describe{
#'     \item{\code{"log"}}{Log link, maps \eqn{(0, \infty) \to R}. Ensures positivity.}
#'     \item{\code{"logit"}}{Logit link, maps \eqn{(0, 1) \to R}.}
#'     \item{\code{"probit"}}{Probit link, maps \eqn{(0, 1) \to R}.}
#'     \item{\code{"cloglog"}}{Complementary log-log link, maps \eqn{(0, 1) \to R}.}
#'     \item{\code{"cauchy"}}{Cauchy link, maps \eqn{(0, 1) \to R}.}
#'     \item{\code{"identity"}}{Identity link (no transformation). Use with
#'       caution as it does not enforce parameter constraints.}
#'     \item{\code{"sqrt"}}{Square root link, maps \eqn{x \to \sqrt{x}}.}
#'     \item{\code{"inverse"}}{Inverse link, maps \eqn{x \to 1/x}.}
#'     \item{\code{"inverse-square"}}{Inverse squared link, maps \eqn{x \to 1/x^2}.}
#'   }
#'
#' @param link_scale Numeric scale factor(s) for link functions. Can be a
#'   single numeric value (applied to all) or a named list.
#'   Default is \code{NULL}, which uses \code{10} for \eqn{\alpha, \beta,
#'   \gamma, \lambda} and \code{1} for \eqn{\delta}. Smaller scales produce
#'   steeper response curves; larger scales produce more gradual ones.
#'
#' @param subset Optional vector specifying a subset of observations.
#'
#' @param weights Optional numeric vector of prior weights. Currently experimental.
#'
#' @param offset Optional vector or matrix specifying an \emph{a priori} known
#'   component to be included in the linear predictor(s). Offsets are added
#'   \emph{before} the link function is applied.
#'
#' @param na.action A function specifying how to handle \code{NA}s.
#'   Defaults to \code{getOption("na.action")}. See
#'   \code{\link[stats]{na.omit}} and \code{\link[stats]{na.exclude}}.
#'
#' @param contrasts Optional list specifying contrasts for factor variables.
#'   See \code{\link[stats]{contrasts}} and \code{contrasts.arg} in
#'   \code{\link[stats]{model.matrix}}.
#'
#' @param control A list of control parameters from \code{\link{gkw_control}}
#'   specifying optimization details (e.g., \code{method}, \code{start},
#'   \code{hessian}, \code{maxit}). See \code{\link{gkw_control}} for all
#'   options and defaults.
#'
#' @param model Logical (default \code{TRUE}). If \code{TRUE}, the model frame
#'   is returned as component \code{model}.
#'
#' @param x Logical (default \code{FALSE}). If \code{TRUE}, the list of
#'   model matrices (one for each modeled parameter) is returned as component \code{x}.
#'
#' @param y Logical (default \code{TRUE}). If \code{TRUE}, the response
#'   vector is returned as component \code{y}.
#'
#' @param ... Arguments for backward compatibility (e.g., \code{method},
#'   \code{start}, \code{hessian}, \code{silent}). These are deprecated and
#'   should be passed via the \code{control} argument. Using them will
#'   trigger a warning.
#'
#' @return
#' An object of class \code{"gkwreg"}, a list containing:
#'
#' \strong{Model Specification:}
#' \describe{
#'   \item{\code{call}}{The matched function call.}
#'   \item{\code{formula}}{The \code{Formula} object used.}
#'   \item{\code{family}}{The distribution family used.}
#'   \item{\code{link}}{Named list of link functions used.}
#'   \item{\code{link_scale}}{Named list of link scale values.}
#'   \item{\code{control}}{The \code{gkw_control} object used.}
#' }
#'
#' \strong{Parameter Estimates:}
#' \describe{
#'   \item{\code{coefficients}}{Named vector of estimated regression
#'     coefficients (on the link scale).}
#'   \item{\code{fitted_parameters}}{Named list of mean parameter values
#'     (e.g., \eqn{\alpha, \beta}) averaged across observations.}
#'   \item{\code{parameter_vectors}}{Named list of observation-specific
#'     parameter vectors (e.g., \code{alphaVec}, \code{betaVec}).}
#' }
#'
#' \strong{Fitted Values and Residuals:}
#' \describe{
#'   \item{\code{fitted.values}}{Vector of fitted mean values \eqn{E[Y|X]}.}
#'   \item{\code{residuals}}{Vector of response residuals (observed - fitted).}
#' }
#'
#' \strong{Inference:} (Only if \code{control$hessian = TRUE})
#' \describe{
#'   \item{\code{vcov}}{Variance-covariance matrix of coefficients.}
#'   \item{\code{se}}{Vector of standard errors of coefficients.}
#' }
#'
#' \strong{Model Fit Statistics:}
#' \describe{
#'   \item{\code{loglik}}{Maximized log-likelihood.}
#'   \item{\code{aic}}{Akaike Information Criterion.}
#'   \item{\code{bic}}{Bayesian Information Criterion.}
#'   \item{\code{deviance}}{Deviance (-2 * loglik).}
#'   \item{\code{df.residual}}{Residual degrees of freedom (nobs - npar).}
#'   \item{\code{nobs}}{Number of observations.}
#'   \item{\code{npar}}{Total number of estimated parameters.}
#' }
#'
#' \strong{Diagnostic Statistics:}
#' \describe{
#'   \item{\code{rmse}}{Root Mean Squared Error.}
#'   \item{\code{efron_r2}}{Efron's pseudo R-squared.}
#'   \item{\code{mean_absolute_error}}{Mean Absolute Error.}
#' }
#'
#' \strong{Optimization Details:}
#' \describe{
#'   \item{\code{convergence}}{Logical, \code{TRUE} if converged.}
#'   \item{\code{message}}{Convergence message from optimizer.}
#'   \item{\code{iterations}}{Number of optimizer iterations.}
#' }
#'
#' \strong{Optional Components:}
#' \describe{
#'   \item{\code{model}}{The model frame (if \code{model = TRUE}).}
#'   \item{\code{x}}{List of model matrices (if \code{x = TRUE}).}
#'   \item{\code{y}}{The response vector (if \code{y = TRUE}).}
#' }
#'
#' \strong{Internal:}
#' \describe{
#'   \item{\code{tmb_object}}{The raw object from \code{\link[TMB]{MakeADFun}}.}
#' }
#'
#' @references
#' \strong{Generalized Kumaraswamy Distribution:}
#'
#' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized distributions.
#' \emph{Journal of Statistical Computation and Simulation}, \strong{81}(7), 883-898.
#' \doi{10.1080/00949650903530745}
#'
#' \strong{Kumaraswamy Distribution:}
#'
#' Kumaraswamy, P. (1980). A generalized probability density function for
#' double-bounded random processes. \emph{Journal of Hydrology}, \strong{46}(1-2), 79-88.
#' \doi{10.1016/0022-1694(80)90036-0}
#'
#' Jones, M. C. (2009). Kumaraswamy's distribution: A beta-type distribution with
#' some tractability advantages. \emph{Statistical Methodology}, \strong{6}(1), 70-81.
#' \doi{10.1016/j.stamet.2008.04.001}
#'
#' \strong{Beta Regression:}
#'
#' Ferrari, S. L. P., & Cribari-Neto, F. (2004). Beta regression for modelling
#' rates and proportions. \emph{Journal of Applied Statistics}, \strong{31}(7), 799-815.
#' \doi{10.1080/0266476042000214501}
#'
#' \strong{Template Model Builder (TMB):}
#'
#' Kristensen, K., Nielsen, A., Berg, C. W., Skaug, H., & Bell, B. M. (2016).
#' TMB: Automatic Differentiation and Laplace Approximation. \emph{Journal of
#' Statistical Software}, \strong{70}(5), 1-21.
#' \doi{10.18637/jss.v070.i05}
#'
#' @author Lopes, J. E.
#'
#' Maintainer: Lopes, J. E.
#'
#' @seealso
#' \strong{Core Functions:}
#' \code{\link{gkw_control}} for fitting options.
#'
#' \strong{S3 Methods:}
#' \code{\link{summary.gkwreg}}, \code{\link{print.gkwreg}},
#' \code{\link{plot.gkwreg}}, \code{\link{predict.gkwreg}},
#' \code{\link{residuals.gkwreg}}, \code{\link{coef.gkwreg}},
#' \code{\link{vcov.gkwreg}}, \code{\link{logLik.gkwreg}},
#' \code{\link{confint.gkwreg}}, \code{\link{anova.gkwreg}}
#'
#' \strong{Distributions:}
#' \code{\link[gkwdist]{dgkw}}, \code{\link[gkwdist]{pgkw}},
#' \code{\link[gkwdist]{qgkw}}, \code{\link[gkwdist]{rgkw}}
#'
#' \strong{Related Packages:}
#' \code{\link[betareg]{betareg}}, \code{\link[Formula]{Formula}},
#' \code{\link[TMB]{TMB}}
#'
#' @examples
#' \donttest{
#' # SECTION 1: Basic Usage - Getting Started
#' # Load packages and data
#' library(gkwreg)
#' library(gkwdist)
#' data(GasolineYield)
#' 
#' # Example 1.1: Simplest possible model (intercept-only, all defaults)
#' fit_basic <- gkwreg(yield ~ 1, data = GasolineYield, family = "kw")
#' summary(fit_basic)
#' 
#' # Example 1.2: Model with predictors (uses all defaults)
#' # Default: family = "gkw", method = "nlminb", hessian = TRUE
#' fit_default <- gkwreg(yield ~ batch + temp, data = GasolineYield)
#' summary(fit_default)
#' 
#' # Example 1.3: Kumaraswamy model (two-parameter family)
#' # Default link functions: log for both alpha and beta
#' fit_kw <- gkwreg(yield ~ batch + temp, data = GasolineYield, family = "kw")
#' summary(fit_kw)
#' plot(fit_kw, ask = FALSE, use_ggplot = TRUE, arrange_plots = TRUE)
#' 
#' # Example 1.4: Beta model for comparison
#' # Default links: log for gamma and delta
#' fit_beta <- gkwreg(yield ~ batch + temp, data = GasolineYield, family = "beta")
#' 
#' # Compare models using AIC/BIC
#' AIC(fit_kw, fit_beta)
#' BIC(fit_kw, fit_beta)
#' 
#' # SECTION 2: Using gkw_control() for Customization
#' # Example 2.1: Change optimization method to BFGS
#' fit_bfgs <- gkwreg(
#'   yield ~ batch + temp,
#'   data = GasolineYield,
#'   family = "kw",
#'   control = gkw_control(method = "BFGS")
#' )
#' summary(fit_bfgs)
#' 
#' # Example 2.3: Fast fitting without standard errors
#' # Useful for model exploration or large datasets
#' fit_fast <- gkwreg(
#'   yield ~ batch + temp,
#'   data = GasolineYield,
#'   family = "kw",
#'   control = gkw_control(hessian = FALSE)
#' )
#' # Note: Cannot compute confint() without hessian
#' coef(fit_fast) # Point estimates still available
#' 
#' # SECTION 3: Advanced Formula Specifications
#' # Example 3.1: Different predictors for different parameters
#' # alpha depends on batch, beta depends on temp
#' fit_diff <- gkwreg(
#'   yield ~ batch | temp,
#'   data = GasolineYield,
#'   family = "kw"
#' )
#' summary(fit_diff)
#' 
#' # Example 3.2: Intercept-only for one parameter
#' # alpha varies with predictors, beta is constant
#' fit_partial <- gkwreg(
#'   yield ~ batch + temp | 1,
#'   data = GasolineYield,
#'   family = "kw"
#' )
#' 
#' # Example 3.3: Complex model with interactions
#' fit_interact <- gkwreg(
#'   yield ~ batch * temp | temp + I(temp^2),
#'   data = GasolineYield,
#'   family = "kw"
#' )
#' 
#' # SECTION 4: Working with Different Families
#' # Example 4.1: Fit multiple families and compare
#' families <- c("beta", "kw", "ekw", "bkw", "gkw")
#' fits <- lapply(families, function(fam) {
#'   gkwreg(yield ~ batch + temp, data = GasolineYield, family = fam)
#' })
#' names(fits) <- families
#' 
#' # Compare via information criteria
#' comparison <- data.frame(
#'   Family = families,
#'   LogLik = sapply(fits, logLik),
#'   AIC = sapply(fits, AIC),
#'   BIC = sapply(fits, BIC),
#'   npar = sapply(fits, function(x) x$npar)
#' )
#' print(comparison)
#' 
#' # Example 4.2: Formal nested model testing
#' fit_kw <- gkwreg(yield ~ batch + temp, GasolineYield, family = "kw")
#' fit_ekw <- gkwreg(yield ~ batch + temp, GasolineYield, family = "ekw")
#' fit_gkw <- gkwreg(yield ~ batch + temp, GasolineYield, family = "gkw")
#' anova(fit_kw, fit_ekw, fit_gkw)
#' 
#' # SECTION 5: Link Functions and Scales
#' # Example 5.1: Custom link functions
#' fit_links <- gkwreg(
#'   yield ~ batch + temp,
#'   data = GasolineYield,
#'   family = "kw",
#'   link = list(alpha = "sqrt", beta = "log")
#' )
#' 
#' # Example 5.2: Custom link scales
#' # Smaller scale = steeper response curve
#' fit_scale <- gkwreg(
#'   yield ~ batch + temp,
#'   data = GasolineYield,
#'   family = "kw",
#'   link_scale = list(alpha = 5, beta = 15)
#' )
#' 
#' # SECTION 6: Prediction and Inference
#' # Fit model for prediction examples
#' fit <- gkwreg(yield ~ batch + temp, GasolineYield, family = "kw")
#' 
#' # Example 6.1: Confidence intervals
#' confint(fit, level = 0.95) # 95% CI
#' 
#' # SECTION 7: Diagnostic Plots and Model Checking
#' fit <- gkwreg(yield ~ batch + temp, GasolineYield, family = "kw")
#' 
#' # Example 7.1: All diagnostic plots (default)
#' plot(fit, ask = FALSE, use_ggplot = TRUE, arrange_plots = TRUE)
#' 
#' # Example 7.2: Select specific plots
#' plot(fit, which = c(2, 4, 5), use_ggplot = TRUE, arrange_plots = TRUE)
#' 
#' # SECTION 8: Real Data Example - Food Expenditure
#' # Load and prepare data
#' data(FoodExpenditure, package = "betareg")
#' food_data <- FoodExpenditure
#' food_data$prop <- food_data$food / food_data$income
#' 
#' # Example 8.1: Basic model
#' fit_food <- gkwreg(
#'   prop ~ persons | income,
#'   data = food_data,
#'   family = "kw"
#' )
#' summary(fit_food)
#' 
#' # Example 8.2: Compare with Beta regression
#' fit_food_beta <- gkwreg(
#'   prop ~ persons | income,
#'   data = food_data,
#'   family = "beta"
#' )
#' 
#' # Which fits better?
#' AIC(fit_food, fit_food_beta)
#' 
#' # Example 8.3: Interpretation via effects
#' income_seq <- seq(min(food_data$income), max(food_data$income), length = 50)
#' pred_data <- data.frame(
#'   persons = median(food_data$persons),
#'   income = income_seq
#' )
#' pred_food <- predict(fit_food, newdata = pred_data, type = "response")
#' 
#' plot(food_data$income, food_data$prop,
#'      xlab = "Income", ylab = "Proportion Spent on Food",
#'      main = "Food Expenditure Pattern"
#' )
#' lines(income_seq, pred_food, col = "red", lwd = 2)
#' 
#' # SECTION 9: Simulation Studies
#' # Example 9.1: Simple Kumaraswamy simulation
#' set.seed(123)
#' n <- 500
#' x1 <- runif(n, -2, 2)
#' x2 <- rnorm(n)
#' 
#' # True model: log(alpha) = 0.8 + 0.3*x1, log(beta) = 1.2 - 0.2*x2
#' eta_alpha <- 0.8 + 0.3 * x1
#' eta_beta <- 1.2 - 0.2 * x2
#' alpha_true <- exp(eta_alpha)
#' beta_true <- exp(eta_beta)
#' 
#' # Generate response
#' y <- rkw(n, alpha = alpha_true, beta = beta_true)
#' sim_data <- data.frame(y = y, x1 = x1, x2 = x2)
#' 
#' # Fit and check parameter recovery
#' fit_sim <- gkwreg(y ~ x1 | x2, data = sim_data, family = "kw")
#' 
#' # Compare estimated vs true coefficients
#' cbind(
#'   True = c(0.8, 0.3, 1.2, -0.2),
#'   Estimated = coef(fit_sim),
#'   SE = fit_sim$se
#' )
#' 
#' # SECTION 10: Memory and Performance Optimization
#' 
#' # Example 10.1: Minimal object for large datasets
#' fit_minimal <- gkwreg(
#'   yield ~ batch + temp,
#'   data = GasolineYield,
#'   family = "kw",
#'   model = FALSE, # Don't store model frame
#'   x = FALSE,     # Don't store design matrices
#'   y = FALSE,     # Don't store response
#'   control = gkw_control(hessian = FALSE) # Skip Hessian
#' )
#' 
#' # Much smaller object
#' object.size(fit_minimal)
#' 
#' # Trade-off: Limited post-fitting capabilities
#' # Can still use: coef(), logLik(), AIC(), BIC()
#' # Cannot use: predict(), some diagnostics
#' 
#' # SECTION 11: Model Selection and Comparison
#' 
#' # Example 11.1: Nested model testing
#' fit1 <- gkwreg(yield ~ 1, GasolineYield, family = "kw")
#' fit2 <- gkwreg(yield ~ batch, GasolineYield, family = "kw")
#' fit3 <- gkwreg(yield ~ batch + temp, GasolineYield, family = "kw")
#' 
#' # Likelihood ratio tests
#' anova(fit1, fit2, fit3)
#' }
#' @keywords regression models
#' @export
gkwreg <- function(formula,
                   data,
                   family = c("gkw", "bkw", "kkw", "ekw", "mc", "kw", "beta"),
                   link = NULL,
                   link_scale = NULL,
                   subset = NULL,
                   weights = NULL,
                   offset = NULL,
                   na.action = getOption("na.action"),
                   contrasts = NULL,
                   control = gkw_control(),
                   model = TRUE,
                   x = FALSE,
                   y = TRUE,
                   ...) {
  # BACKWARD COMPATIBILITY LAYER
  # Handle deprecated arguments with informative warnings


  dots <- list(...)

  # Check for completely removed arguments (violated Separation of Concerns)
  if ("plot" %in% names(dots)) {
    warning(
      "Argument 'plot' has been removed from gkwreg().\n",
      "  Reason: This argument had no effect and was misleading.\n",
      "  Solution: Use plot(fitted_model) after fitting to generate diagnostics.\n",
      "  Example: fit <- gkwreg(...); plot(fit)",
      call. = FALSE
    )
    dots$plot <- NULL
  }

  if ("conf.level" %in% names(dots)) {
    warning(
      "Argument 'conf.level' has been removed from gkwreg().\n",
      "  Reason: Confidence intervals can be computed at any level without refitting.\n",
      "  Solution: Use confint(fitted_model, level = your_level) after fitting.\n",
      "  Example: fit <- gkwreg(...); confint(fit, level = 0.90)",
      call. = FALSE
    )
    dots$conf.level <- NULL
  }

  # Check for arguments that moved to control object
  moved_to_control <- c(
    "method", "start", "fixed", "hessian",
    "silent", "optimizer.control"
  )
  found_moved <- intersect(names(dots), moved_to_control)

  if (length(found_moved) > 0) {
    warning(
      "Arguments ", paste0("'", found_moved, "'", collapse = ", "),
      " should now be passed via 'control' argument.\n",
      "  New syntax: gkwreg(..., control = gkw_control(",
      paste(found_moved, "= ...", collapse = ", "), "))\n",
      "  Attempting to use provided values, but please update your code.",
      call. = FALSE, immediate. = TRUE
    )

    # Build control from deprecated arguments
    deprecated_control_args <- list()
    for (arg in found_moved) {
      if (arg == "optimizer.control") {
        # Special handling: merge optimizer.control into control
        opt_ctrl <- dots[[arg]]
        if (is.list(opt_ctrl)) {
          deprecated_control_args <- utils::modifyList(deprecated_control_args, opt_ctrl)
        }
      } else {
        deprecated_control_args[[arg]] <- dots[[arg]]
      }
      dots[[arg]] <- NULL
    }

    # Merge deprecated args with user-provided control (user takes precedence)
    if (inherits(control, "gkw_control")) {
      # Convert control object to list, merge, then recreate
      control_list <- as.list(control)
      control_list <- utils::modifyList(deprecated_control_args, control_list)
      control <- do.call(gkw_control, control_list[names(control_list) %in%
        names(formals(gkw_control))])
    } else {
      # User provided a list directly, merge it
      deprecated_control_args <- utils::modifyList(
        deprecated_control_args,
        as.list(control)
      )
      control <- do.call(gkw_control, deprecated_control_args)
    }
  }

  # Warn about completely unused arguments
  if (length(dots) > 0) {
    warning(
      "Unused arguments: ", paste(names(dots), collapse = ", "),
      call. = FALSE
    )
  }

  # VALIDATION & SETUP
  # Store the call for the result
  call <- match.call()

  # Match family argument
  family <- match.arg(family, choices = c("gkw", "bkw", "kkw", "ekw", "mc", "kw", "beta"))

  # Validate control object
  if (!inherits(control, "gkw_control")) {
    stop(
      "'control' must be an object from gkw_control().\n",
      "  Use: control = gkw_control(...)",
      call. = FALSE
    )
  }

  # Check required packages
  if (!requireNamespace("Formula", quietly = TRUE)) {
    stop(
      "Package 'Formula' is required but not installed.\n",
      "  Install with: install.packages('Formula')",
      call. = FALSE
    )
  }

  if (!requireNamespace("TMB", quietly = TRUE)) {
    stop(
      "Package 'TMB' is required but not installed.\n",
      "  Install with: install.packages('TMB')",
      call. = FALSE
    )
  }

  # Extract control parameters for easier access
  method <- control$method
  use_nlminb <- (method == "nlminb")
  silent <- control$silent

  # PARAMETER & FORMULA PROCESSING
  # Using existing internal functions

  # Get parameter information for the specified family
  param_info <- .get_family_param_info(family)
  param_names <- param_info$names
  fixed_params <- param_info$fixed
  param_positions <- param_info$positions

  # Convert to Formula object
  formula_obj <- Formula::as.Formula(formula)

  # Process formula parts for each parameter
  formula_list <- .process_formula_parts(
    formula_obj, param_names,
    fixed_params, data
  )

  # Process link functions
  link_list <- .process_link(link, param_names, fixed_params)

  # Process link scales
  link_scale_list <- .process_link_scale(
    link_scale, link_list,
    param_names, fixed_params
  )

  # Convert link strings to integers for TMB
  link_ints <- .convert_links_to_int(link_list)

  # DATA EXTRACTION & VALIDATION
  # Using existing internal functions

  # Extract model frames, responses, and model matrices
  model_data <- .extract_model_data(
    formula_list = formula_list,
    data = data,
    subset = subset,
    weights = weights,
    na.action = na.action,
    offset = offset,
    contrasts = contrasts,
    original_call = call
  )

  # Get response variable
  y_var <- model_data$y

  # Validate response is in (0, 1)
  invisible(.validate_data(y_var, length(param_names)))

  # FIXED PARAMETERS PROCESSING
  # Using existing internal function

  # Process fixed parameters (from control or family definition)
  fixed_processed <- .process_fixed(control$fixed, param_names, fixed_params)

  # TMB PREPARATION
  # Using existing internal functions

  # Prepare TMB data with correct matrices based on family
  tmb_data <- .prepare_tmb_data(
    model_data, family, param_names, fixed_processed,
    link_ints, link_scale_list, y_var, param_positions
  )

  # Prepare TMB parameters in the correct structure for the family
  tmb_params <- .prepare_tmb_params(
    model_data, family, param_names, fixed_processed,
    param_positions
  )

  # Override with user-provided starting values if available
  if (!is.null(control$start)) {
    for (param_name in names(control$start)) {
      if (param_name %in% param_names) {
        pos <- param_positions[[param_name]]
        tmb_param_name <- paste0("beta", pos)
        if (tmb_param_name %in% names(tmb_params)) {
          user_start <- control$start[[param_name]]
          expected_length <- length(tmb_params[[tmb_param_name]])
          if (length(user_start) == expected_length) {
            tmb_params[[tmb_param_name]] <- user_start
          } else {
            warning(
              "Starting value for '", param_name, "' has length ",
              length(user_start), " but expected ", expected_length,
              ". Ignoring user-provided starting value.",
              call. = FALSE
            )
          }
        }
      }
    }
  }

  # TMB MODEL COMPILATION
  # Using existing internal function

  # Determine TMB model name based on family
  if (family == "beta") {
    dll_name <- "gkwbetareg"
  } else {
    dll_name <- paste0(family, "reg")
  }

  if (!silent) {
    message("Using TMB model: ", dll_name)
  }

  # Check and compile TMB model if needed
  .check_and_compile_TMB_code(dll_name, verbose = !silent)

  # Create TMB object
  obj <- TMB::MakeADFun(
    data = tmb_data,
    parameters = tmb_params,
    DLL = dll_name,
    silent = silent
  )

  # OPTIMIZATION

  if (use_nlminb) {
    # Use nlminb optimizer
    if (!silent) message("Optimizing with nlminb...")

    opt <- tryCatch(
      stats::nlminb(
        start = obj$par,
        objective = obj$fn,
        gradient = obj$gr,
        control = control$nlminb_control
      ),
      error = function(e) {
        stop("Optimization with nlminb failed: ", e$message, call. = FALSE)
      }
    )

    # Extract results
    fit_result <- list(
      coefficients = obj$env$last.par,
      loglik = -opt$objective,
      convergence = (opt$convergence == 0),
      message = opt$message,
      iterations = opt$iterations
    )
  } else {
    # Use optim with specified method
    if (!silent) {
      message("Optimizing with optim (method = ", method, ")...")
    }

    opt <- tryCatch(
      stats::optim(
        par = obj$par,
        fn = obj$fn,
        gr = obj$gr,
        method = method,
        control = control$optim_control
      ),
      error = function(e) {
        stop("Optimization with optim (", method, ") failed: ",
          e$message,
          call. = FALSE
        )
      }
    )

    # Extract results
    fit_result <- list(
      coefficients = opt$par,
      loglik = -opt$value,
      convergence = (opt$convergence == 0),
      message = if (opt$convergence == 0) {
        "Successful convergence"
      } else {
        paste("Optimization failed (code ", opt$convergence, ")")
      },
      iterations = opt$counts[1]
    )
  }

  # COEFFICIENT NAMING
  # Using existing internal function

  # Format coefficient names with parameter mappings
  coef_names <- .format_coefficient_names(param_names, model_data, param_positions)
  names(fit_result$coefficients) <- coef_names

  # HESSIAN & STANDARD ERRORS


  if (control$hessian) {
    if (!silent) message("Computing standard errors...")

    tryCatch(
      {
        sd_report <- TMB::sdreport(obj, getJointPrecision = TRUE)
        fit_result$se <- as.vector(sd_report$sd)
        fit_result$vcov <- as.matrix(sd_report$cov)

        # Add names
        names(fit_result$se) <- coef_names
        rownames(fit_result$vcov) <- coef_names
        colnames(fit_result$vcov) <- coef_names
      },
      error = function(e) {
        warning(
          "Could not compute standard errors: ", e$message, "\n",
          "  Model results are still valid for point estimates.",
          call. = FALSE
        )
        fit_result$se <- rep(NA_real_, length(fit_result$coefficients))
        fit_result$vcov <- matrix(
          NA_real_,
          length(fit_result$coefficients),
          length(fit_result$coefficients)
        )
        names(fit_result$se) <- coef_names
        rownames(fit_result$vcov) <- coef_names
        colnames(fit_result$vcov) <- coef_names
      }
    )
  } else {
    # Hessian not requested
    fit_result$se <- NULL
    fit_result$vcov <- NULL
  }


  # EXTRACT TMB REPORT


  tmb_report <- obj$report()

  # Extract parameter means
  alpha_mean <- tmb_report$alpha_mean
  beta_mean <- tmb_report$beta_mean
  gamma_mean <- tmb_report$gamma_mean
  delta_mean <- tmb_report$delta_mean
  lambda_mean <- tmb_report$lambda_mean

  # Extract parameter vectors for each observation
  alphaVec <- if ("alphaVec" %in% names(tmb_report)) {
    tmb_report$alphaVec
  } else {
    rep(alpha_mean, length(y_var))
  }

  betaVec <- if ("betaVec" %in% names(tmb_report)) {
    tmb_report$betaVec
  } else {
    rep(beta_mean, length(y_var))
  }

  gammaVec <- if ("gammaVec" %in% names(tmb_report)) {
    tmb_report$gammaVec
  } else {
    rep(gamma_mean, length(y_var))
  }

  deltaVec <- if ("deltaVec" %in% names(tmb_report)) {
    tmb_report$deltaVec
  } else {
    rep(delta_mean, length(y_var))
  }

  lambdaVec <- if ("lambdaVec" %in% names(tmb_report)) {
    tmb_report$lambdaVec
  } else {
    rep(lambda_mean, length(y_var))
  }

  # Extract fitted values
  fitted_values <- if ("fitted" %in% names(tmb_report)) {
    tmb_report$fitted
  } else {
    rep(NA_real_, length(y_var))
  }


  # FIT STATISTICS


  # Calculate residuals
  response_residuals <- y_var - fitted_values

  # Sample size and parameter count
  n_obs <- length(y_var)
  npar <- length(fit_result$coefficients)

  # Information criteria
  nll <- -fit_result$loglik
  deviance <- 2.0 * nll
  aic <- deviance + 2.0 * npar
  bic <- deviance + log(n_obs) * npar

  # Degrees of freedom
  df.residual <- n_obs - npar

  # Additional fit statistics
  rmse <- sqrt(mean(response_residuals^2, na.rm = TRUE))
  sse <- sum((y_var - fitted_values)^2, na.rm = TRUE)
  sst <- sum((y_var - mean(y_var, na.rm = TRUE))^2, na.rm = TRUE)
  efron_r2 <- if (sst > 0) 1 - sse / sst else NA_real_
  mean_absolute_error <- mean(abs(response_residuals), na.rm = TRUE)


  # BUILD RESULT OBJECT
  # CRITICAL: Maintain exact same structure and names as original


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
    link_scale = link_scale_list,
    param_names = param_names,
    fixed_params = fixed_params,
    loglik = fit_result$loglik,
    aic = aic,
    bic = bic,
    deviance = deviance,
    df.residual = df.residual,
    nobs = n_obs,
    npar = npar,
    vcov = fit_result$vcov,
    se = fit_result$se,
    convergence = fit_result$convergence,
    message = fit_result$message,
    iterations = fit_result$iterations,
    rmse = rmse,
    efron_r2 = efron_r2,
    mean_absolute_error = mean_absolute_error,
    method = method,
    control = control, # Store control for reference
    tmb_object = obj
  )

  # Add optional components
  if (x) result$x <- model_data$matrices
  if (y) result$y <- y_var
  if (model) result$model <- model_data$model

  # Set class
  class(result) <- "gkwreg"

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


#' Prepare TMB Data for GKw Regression
#'
#' @param model_data List of model data.
#' @param family Family name.
#' @param param_names Names of parameters.
#' @param fixed List of fixed parameters and coefficients.
#' @param link_ints List of link function integers.
#' @param link_scale_list List of link scale values.
#' @param y Response variable.
#' @param param_positions Parameter position mapping for the family.
#' @return A list with TMB data.
#' @keywords internal
.prepare_tmb_data <- function(model_data, family, param_names, fixed, link_ints, link_scale_list, y, param_positions) {
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

      # If link_scale exists for this parameter, use it
      if (param %in% names(link_scale_list)) {
        tmb_data[[scale_name]] <- link_scale_list[[param]]
      }
    }
  }

  return(tmb_data)
}


#' #' Fit Generalized Kumaraswamy Regression Models
#' #'
#' #' @description
#' #' Fits regression models using the Generalized Kumaraswamy (GKw) family of
#' #' distributions for response variables strictly bounded in the interval (0, 1).
#' #' The function allows modeling parameters from all seven submodels of the GKw
#' #' family as functions of predictors using appropriate link functions. Estimation
#' #' is performed using Maximum Likelihood via the TMB (Template Model Builder) package.
#' #' Requires the \code{Formula} and \code{TMB} packages.
#' #'
#' #' @param formula An object of class \code{\link[Formula]{Formula}} (or one that
#' #'   can be coerced to that class). It should be structured as
#' #'   \code{y ~ model_alpha | model_beta | model_gamma | model_delta | model_lambda},
#' #'   where \code{y} is the response variable and each \code{model_*} part specifies
#' #'   the linear predictor for the corresponding parameter (\eqn{\alpha}, \eqn{\beta},
#' #'   \eqn{\gamma}, \eqn{\delta}, \eqn{\lambda}). If a part is omitted or specified
#' #'   as \code{~ 1} or \code{.}, an intercept-only model is used for that parameter.
#' #'   See Details for parameter correspondence in subfamilies.
#' #' @param data A data frame containing the variables specified in the \code{formula}.
#' #' @param family A character string specifying the desired distribution family.
#' #'   Defaults to \code{"gkw"}. Supported families are:
#' #'   \itemize{
#' #'     \item \code{"gkw"}: Generalized Kumaraswamy (5 parameters: \eqn{\alpha, \beta, \gamma, \delta, \lambda})
#' #'     \item \code{"bkw"}: Beta-Kumaraswamy (4 parameters: \eqn{\alpha, \beta, \gamma, \delta}; \eqn{\lambda = 1} fixed)
#' #'     \item \code{"kkw"}: Kumaraswamy-Kumaraswamy (4 parameters: \eqn{\alpha, \beta, \delta, \lambda}; \eqn{\gamma = 1} fixed)
#' #'     \item \code{"ekw"}: Exponentiated Kumaraswamy (3 parameters: \eqn{\alpha, \beta, \lambda}; \eqn{\gamma = 1, \delta = 0} fixed)
#' #'     \item \code{"mc"}: McDonald / Beta Power (3 parameters: \eqn{\gamma, \delta, \lambda}; \eqn{\alpha = 1, \beta = 1} fixed)
#' #'     \item \code{"kw"}: Kumaraswamy (2 parameters: \eqn{\alpha, \beta}; \eqn{\gamma = 1, \delta = 0, \lambda = 1} fixed)
#' #'     \item \code{"beta"}: Beta distribution (2 parameters: \eqn{\gamma, \delta}; \eqn{\alpha = 1, \beta = 1, \lambda = 1} fixed)
#' #'   }
#' #' @param link Either a single character string specifying the same link function
#' #'   for all relevant parameters, or a named list specifying the link function for each
#' #'   modeled parameter (e.g., \code{list(alpha = "log", beta = "log", delta = "logit")}).
#' #'   Defaults are \code{"log"} for \eqn{\alpha, \beta, \gamma, \lambda} (parameters > 0)
#' #'   and \code{"logit"} for \eqn{\delta} (parameter in (0, 1)).
#' #'   Supported link functions are:
#' #'   \itemize{
#' #'     \item \code{"log"}: logarithmic link, maps \eqn{(0, \infty)} to \eqn{(-\infty, \infty)}
#' #'     \item \code{"identity"}: identity link, no transformation
#' #'     \item \code{"inverse"}: inverse link, maps \eqn{x} to \eqn{1/x}
#' #'     \item \code{"sqrt"}: square root link, maps \eqn{x} to \eqn{\sqrt{x}}
#' #'     \item \code{"inverse-square"}: inverse squared link, maps \eqn{x} to \eqn{1/x^2}
#' #'     \item \code{"logit"}: logistic link, maps \eqn{(0, 1)} to \eqn{(-\infty, \infty)}
#' #'     \item \code{"probit"}: probit link, using normal CDF
#' #'     \item \code{"cloglog"}: complementary log-log
#' #'     \item \code{"cauchy"}: Cauchy link, using Cauchy CDF
#' #'   }
#' #' @param link_scale Either a single numeric value specifying the same scale for all
#' #'   link functions, or a named list specifying the scale for each parameter's link
#' #'   function (e.g., \code{list(alpha = 10, beta = 5, delta = 1)}). The scale affects
#' #'   how the link function transforms the linear predictor. Default is 10 for most
#' #'   parameters and 1 for parameters using probability-type links (such as \code{delta}).
#' #'   For probability-type links (logit, probit, cloglog, cauchy), smaller values
#' #'   produce more extreme transformations.
#' #' @param start An optional named list providing initial values for the regression
#' #'   coefficients. Parameter names should match the distribution parameters (alpha,
#' #'   beta, etc.), and values should be vectors corresponding to the coefficients
#' #'   in the respective linear predictors (including intercept). If \code{NULL}
#' #'   (default), suitable starting values are automatically determined based on
#' #'   global parameter estimates.
#' #' @param fixed An optional named list specifying parameters or coefficients to be
#' #'   held fixed at specific values during estimation. Currently not fully implemented.
#' #' @param method Character string specifying the optimization algorithm to use.
#' #'   Options are \code{"nlminb"} (default, using \code{\link[stats]{nlminb}}),
#' #'   \code{"BFGS"}, \code{"Nelder-Mead"}, \code{"CG"}, \code{"SANN"}, or \code{"L-BFGS-B"}.
#' #'   If \code{"nlminb"} is selected, R's \code{\link[stats]{nlminb}} function is used;
#' #'   otherwise, R's \code{\link[stats]{optim}} function is used with the specified method.
#' #' @param hessian Logical. If \code{TRUE} (default), the Hessian matrix is computed
#' #'   via \code{\link[TMB]{sdreport}} to obtain standard errors and the covariance
#' #'   matrix of the estimated coefficients. Setting to \code{FALSE} speeds up fitting
#' #'   but prevents calculation of standard errors and confidence intervals.
#' #' @param plot Logical. If \code{TRUE} (default), enables the generation of diagnostic
#' #'   plots when calling the generic \code{plot()} function on the fitted object.
#' #'   Actual plotting is handled by the \code{plot.gkwreg} method.
#' #' @param conf.level Numeric. The confidence level (1 - alpha) for constructing
#' #'   confidence intervals for the parameters. Default is 0.95. Used only if
#' #'   \code{hessian = TRUE}.
#' #' @param optimizer.control A list of control parameters passed directly to the
#' #'   chosen optimizer (\code{\link[stats]{nlminb}} or \code{\link[stats]{optim}}).
#' #'   See their respective documentation for details.
#' #' @param subset An optional vector specifying a subset of observations from \code{data}
#' #'   to be used in the fitting process.
#' #' @param weights An optional vector of prior weights (e.g., frequency weights)
#' #'   to be used in the fitting process. Should be \code{NULL} or a numeric vector
#' #'   of non-negative values.
#' #' @param offset An optional numeric vector or matrix specifying an a priori known
#' #'   component to be included *on the scale of the linear predictor* for each parameter.
#' #'   If a vector, it's applied to the predictor of the first parameter in the standard
#' #'   order (\eqn{\alpha}). If a matrix, columns must correspond to parameters in the
#' #'   order \eqn{\alpha, \beta, \gamma, \delta, \lambda}.
#' #' @param na.action A function which indicates what should happen when the data
#' #'   contain \code{NA}s. The default (\code{na.fail}) stops if \code{NA}s are
#' #'   present. Other options include \code{\link[stats]{na.omit}} or
#' #'   \code{\link[stats]{na.exclude}}.
#' #' @param contrasts An optional list specifying the contrasts to be used for factor
#' #'   variables in the model. See the \code{contrasts.arg} of
#' #'   \code{\link[stats]{model.matrix.default}}.
#' #' @param x Logical. If \code{TRUE}, the list of model matrices (one for each modeled
#' #'   parameter) is returned as component \code{x} of the fitted object. Default \code{FALSE}.
#' #' @param y Logical. If \code{TRUE} (default), the response variable (after processing
#' #'   by \code{na.action}, \code{subset}) is returned as component \code{y}.
#' #' @param model Logical. If \code{TRUE} (default), the model frame (containing all
#' #'   variables used from \code{data}) is returned as component \code{model}.
#' #' @param silent Logical. If \code{TRUE} (default), suppresses progress messages
#' #'   from TMB compilation and optimization. Set to \code{FALSE} for verbose output.
#' #' @param ... Additional arguments, currently unused or passed down to internal
#' #'   methods (potentially).
#' #'
#' #' @return An object of class \code{gkwreg}. This is a list containing the
#' #'   following components:
#' #'   \item{call}{The matched function call.}
#' #'   \item{family}{The specified distribution family string.}
#' #'   \item{formula}{The \code{\link[Formula]{Formula}} object used.}
#' #'   \item{coefficients}{A named vector of estimated regression coefficients.}
#' #'   \item{fitted.values}{Vector of estimated means (expected values) of the response.}
#' #'   \item{residuals}{Vector of response residuals (observed - fitted mean).}
#' #'   \item{fitted_parameters}{List containing the estimated mean value for each distribution parameter (\eqn{\alpha, \beta, \gamma, \delta, \lambda}).}
#' #'   \item{parameter_vectors}{List containing vectors of the estimated parameters (\eqn{\alpha, \beta, \gamma, \delta, \lambda}) for each observation, evaluated on the link scale.}
#' #'   \item{link}{List of link functions used for each parameter.}
#' #'   \item{param_names}{Character vector of names of the parameters modeled by the family.}
#' #'   \item{fixed_params}{Named list indicating which parameters are fixed by the family definition.}
#' #'   \item{loglik}{The maximized log-likelihood value.}
#' #'   \item{aic}{Akaike Information Criterion.}
#' #'   \item{bic}{Bayesian Information Criterion.}
#' #'   \item{deviance}{The deviance ( -2 * loglik ).}
#' #'   \item{df.residual}{Residual degrees of freedom (nobs - npar).}
#' #'   \item{nobs}{Number of observations used in the fit.}
#' #'   \item{npar}{Total number of estimated parameters (coefficients).}
#' #'   \item{vcov}{The variance-covariance matrix of the coefficients (if \code{hessian = TRUE}).}
#' #'   \item{se}{Standard errors of the coefficients (if \code{hessian = TRUE}).}
#' #'   \item{convergence}{Convergence code from the optimizer (0 typically indicates success).}
#' #'   \item{message}{Convergence message from the optimizer.}
#' #'   \item{iterations}{Number of iterations used by the optimizer.}
#' #'   \item{rmse}{Root Mean Squared Error of response residuals.}
#' #'   \item{efron_r2}{Efron's pseudo R-squared.}
#' #'   \item{mean_absolute_error}{Mean Absolute Error of response residuals.}
#' #'   \item{x}{List of model matrices (if \code{x = TRUE}).}
#' #'   \item{y}{The response vector (if \code{y = TRUE}).}
#' #'   \item{model}{The model frame (if \code{model = TRUE}).}
#' #'   \item{tmb_object}{The raw object returned by \code{\link[TMB]{MakeADFun}}.}
#' #'
#' #' @details
#' #' The \code{gkwreg} function provides a regression framework for the Generalized
#' #' Kumaraswamy (GKw) family and its submodels, extending density estimation
#' #' to include covariates. The response variable must be strictly bounded in the
#' #' (0, 1) interval.
#' #'
#' #' \strong{Model Specification:}
#' #' The extended \code{\link[Formula]{Formula}} syntax is crucial for specifying
#' #' potentially different linear predictors for each relevant distribution parameter.
#' #' The parameters (\eqn{\alpha, \beta, \gamma, \delta, \lambda}) correspond sequentially
#' #' to the parts of the formula separated by \code{|}. For subfamilies where some
#' #' parameters are fixed by definition (see \code{family} argument), the corresponding
#' #' parts of the formula are automatically ignored. For example, in a \code{family = "kw"}
#' #' model, only the first two parts (for \eqn{\alpha} and \eqn{\beta}) are relevant.
#' #'
#' #' \strong{Parameter Constraints and Link Functions:}
#' #' The parameters \eqn{\alpha, \beta, \gamma, \lambda} are constrained to be positive,
#' #' while \eqn{\delta} is constrained to the interval (0, 1). The default link functions
#' #' (\code{"log"} for positive parameters, \code{"logit"} for \eqn{\delta}) ensure these
#' #' constraints during estimation. Users can specify alternative link functions suitable
#' #' for the parameter's domain via the \code{link} argument.
#' #'
#' #' \strong{Link Scales:}
#' #' The \code{link_scale} parameter allows users to control how aggressively the link
#' #' function transforms the linear predictor. For probability-type links (logit, probit,
#' #' cloglog, cauchy), smaller values (e.g., 1) produce more extreme transformations,
#' #' while larger values (e.g., 10) produce more gradual transformations. For continuous
#' #' parameters, scale values control the sensitivity of the transformation.
#' #'
#' #' \strong{Families and Parameters:}
#' #' The function automatically handles parameter fixing based on the chosen \code{family}:
#' #' \itemize{
#' #'   \item \strong{GKw}: All 5 parameters (\eqn{\alpha, \beta, \gamma, \delta, \lambda}) modeled.
#' #'   \item \strong{BKw}: Models \eqn{\alpha, \beta, \gamma, \delta}; fixes \eqn{\lambda = 1}.
#' #'   \item \strong{KKw}: Models \eqn{\alpha, \beta, \delta, \lambda}; fixes \eqn{\gamma = 1}.
#' #'   \item \strong{EKw}: Models \eqn{\alpha, \beta, \lambda}; fixes \eqn{\gamma = 1, \delta = 0}.
#' #'   \item \strong{Mc} (McDonald): Models \eqn{\gamma, \delta, \lambda}; fixes \eqn{\alpha = 1, \beta = 1}.
#' #'   \item \strong{Kw} (Kumaraswamy): Models \eqn{\alpha, \beta}; fixes \eqn{\gamma = 1, \delta = 0, \lambda = 1}.
#' #'   \item \strong{Beta}: Models \eqn{\gamma, \delta}; fixes \eqn{\alpha = 1, \beta = 1, \lambda = 1}. This parameterization corresponds to the standard Beta distribution with shape1 = \eqn{\gamma} and shape2 = \eqn{\delta}.
#' #' }
#' #'
#' #' \strong{Estimation Engine:}
#' #' Maximum Likelihood Estimation (MLE) is performed using C++ templates via the
#' #' \code{TMB} package, which provides automatic differentiation and efficient
#' #' optimization capabilities. The specific TMB template used depends on the chosen \code{family}.
#' #'
#' #' \strong{Optimizer Method (\code{method} argument):}
#' #' \itemize{
#' #'   \item \code{"nlminb"}: Uses R's built-in \code{stats::nlminb} optimizer. Good for problems with box constraints. Default option.
#' #'   \item \code{"Nelder-Mead"}: Uses R's \code{stats::optim} with the Nelder-Mead simplex algorithm, which doesn't require derivatives.
#' #'   \item \code{"BFGS"}: Uses R's \code{stats::optim} with the BFGS quasi-Newton method for unconstrained optimization.
#' #'   \item \code{"CG"}: Uses R's \code{stats::optim} with conjugate gradients method for unconstrained optimization.
#' #'   \item \code{"SANN"}: Uses R's \code{stats::optim} with simulated annealing, a global optimization method useful for problems with multiple local minima.
#' #'   \item \code{"L-BFGS-B"}: Uses R's \code{stats::optim} with the limited-memory BFGS method with box constraints.
#' #' }
#' #'
#' #' @examples
#' #' \donttest{
#' #'
#' #' ## -------------------------------------------------------------------------
#' #' ## 1. Real-world Case Studies
#' #' ## -------------------------------------------------------------------------
#' #'
#' #' ## Example 1: Food Expenditure Data
#' #' # Load required package
#' #' require(gkwdist)
#' #' require(gkwreg)
#' #'
#' #' # Get FoodExpenditure data and create response variable 'y' as proportion of income spent on food
#' #' data(FoodExpenditure)
#' #' food_data <- FoodExpenditure
#' #' food_data <- within(food_data, {
#' #'   y <- food / income
#' #' })
#' #'
#' #' # Define formula: y depends on 'persons' with 'income' as predictor for second parameter
#' #' formu_fe <- y ~ persons | income
#' #'
#' #' # Fit Kumaraswamy model with log link for both parameters
#' #' kw_model <- gkwreg(formu_fe, food_data,
#' #'   family = "kw",
#' #'   link = rep("log", 2), method = "nlminb"
#' #' )
#' #'
#' #' # Display model summary and diagnostics
#' #' summary(kw_model)
#' #' plot(kw_model, use_ggplot = TRUE, arrange_plots = TRUE, sub.caption = "")
#' #'
#' #' ## Example 2: Gasoline Yield Data
#' #' # Load GasolineYield dataset
#' #' data(GasolineYield)
#' #' gasoline_data <- GasolineYield
#' #'
#' #' # Formula: yield depends on batch and temperature
#' #' # First part (for alpha/gamma) includes batch and temp
#' #' # Second part (for beta/delta/phi) includes only temp
#' #' formu_gy <- yield ~ batch + temp | temp
#' #'
#' #' # Fit Kumaraswamy model with log link and BFGS optimization
#' #' kw_model_gas <- gkwreg(formu_gy, gasoline_data,
#' #'   family = "kw",
#' #'   link = rep("log", 2), method = "BFGS"
#' #' )
#' #'
#' #' # Display results
#' #' summary(kw_model_gas)
#' #' plot(kw_model_gas, use_ggplot = TRUE, arrange_plots = TRUE, sub.caption = "")
#' #'
#' #' ## -------------------------------------------------------------------------
#' #' ## 2. Simulation Studies
#' #' ## -------------------------------------------------------------------------
#' #'
#' #' ## Example 1: Simple Kumaraswamy Regression Model
#' #' # Set seed for reproducibility
#' #' set.seed(123)
#' #' n <- 1000
#' #' x1 <- runif(n, -2, 2)
#' #' x2 <- rnorm(n)
#' #'
#' #' # Define true regression coefficients
#' #' alpha_coef <- c(0.8, 0.3, -0.2) # Intercept, x1, x2
#' #' beta_coef <- c(1.2, -0.4, 0.1) # Intercept, x1, x2
#' #'
#' #' # Generate linear predictors and transform using exponential link
#' #' eta_alpha <- alpha_coef[1] + alpha_coef[2] * x1 + alpha_coef[3] * x2
#' #' eta_beta <- beta_coef[1] + beta_coef[2] * x1 + beta_coef[3] * x2
#' #' alpha_true <- exp(eta_alpha)
#' #' beta_true <- exp(eta_beta)
#' #'
#' #' # Generate responses from Kumaraswamy distribution
#' #' y <- rkw(n, alpha = alpha_true, beta = beta_true)
#' #' df1 <- data.frame(y = y, x1 = x1, x2 = x2)
#' #'
#' #' # Fit Kumaraswamy regression model with formula notation
#' #' # Model: alpha ~ x1 + x2 and beta ~ x1 + x2
#' #' kw_reg <- gkwreg(y ~ x1 + x2 | x1 + x2, data = df1, family = "kw", silent = TRUE)
#' #'
#' #' # Alternative model with custom link scales
#' #' kw_reg2 <- gkwreg(y ~ x1 + x2 | x1 + x2,
#' #'   data = df1, family = "kw",
#' #'   link_scale = list(alpha = 5, beta = 8), silent = TRUE
#' #' )
#' #'
#' #' # Display model summary
#' #' summary(kw_reg)
#' #'
#' #' ## Example 2: Generalized Kumaraswamy Regression
#' #' # Set seed for reproducibility
#' #' set.seed(456)
#' #' n <- 1000
#' #' x1 <- runif(n, -1, 1)
#' #' x2 <- rnorm(n)
#' #' x3 <- factor(rbinom(n, 1, 0.5), labels = c("A", "B")) # Factor variable
#' #'
#' #' # Define true regression coefficients for all parameters
#' #' alpha_coef <- c(0.5, 0.2) # Intercept, x1
#' #' beta_coef <- c(0.8, -0.3, 0.1) # Intercept, x1, x2
#' #' gamma_coef <- c(0.6, 0.4) # Intercept, x3B
#' #' delta_coef <- c(0.0, 0.2) # Intercept, x3B (logit scale)
#' #' lambda_coef <- c(-0.2, 0.1) # Intercept, x2
#' #'
#' #' # Create design matrices
#' #' X_alpha <- model.matrix(~x1, data = data.frame(x1 = x1))
#' #' X_beta <- model.matrix(~ x1 + x2, data = data.frame(x1 = x1, x2 = x2))
#' #' X_gamma <- model.matrix(~x3, data = data.frame(x3 = x3))
#' #' X_delta <- model.matrix(~x3, data = data.frame(x3 = x3))
#' #' X_lambda <- model.matrix(~x2, data = data.frame(x2 = x2))
#' #'
#' #' # Generate parameters through linear predictors and appropriate link functions
#' #' alpha <- exp(X_alpha %*% alpha_coef)
#' #' beta <- exp(X_beta %*% beta_coef)
#' #' gamma <- exp(X_gamma %*% gamma_coef)
#' #' delta <- plogis(X_delta %*% delta_coef) # logit link for delta
#' #' lambda <- exp(X_lambda %*% lambda_coef)
#' #'
#' #' # Generate response from Generalized Kumaraswamy distribution
#' #' y <- rgkw(n, alpha = alpha, beta = beta, gamma = gamma, delta = delta, lambda = lambda)
#' #' df2 <- data.frame(y = y, x1 = x1, x2 = x2, x3 = x3)
#' #'
#' #' # Fit GKw regression with parameter-specific formulas
#' #' gkw_reg <- gkwreg(y ~ x1 | x1 + x2 | x3 | x3 | x2, data = df2, family = "gkw")
#' #'
#' #' # Alternative model with custom link scales
#' #' gkw_reg2 <- gkwreg(y ~ x1 | x1 + x2 | x3 | x3 | x2,
#' #'   data = df2, family = "gkw",
#' #'   link_scale = list(
#' #'     alpha = 12, beta = 12, gamma = 12,
#' #'     delta = 0.8, lambda = 12
#' #'   )
#' #' )
#' #'
#' #' # Compare true vs. estimated coefficients
#' #' print("Estimated Coefficients (GKw):")
#' #' print(coef(gkw_reg))
#' #' print("True Coefficients (approx):")
#' #' print(list(
#' #'   alpha = alpha_coef, beta = beta_coef, gamma = gamma_coef,
#' #'   delta = delta_coef, lambda = lambda_coef
#' #' ))
#' #'
#' #' ## Example 3: Beta Regression for Comparison
#' #' # Set seed for reproducibility
#' #' set.seed(789)
#' #' n <- 1000
#' #' x1 <- runif(n, -1, 1)
#' #'
#' #' # True coefficients for Beta parameters (gamma = shape1, delta = shape2)
#' #' gamma_coef <- c(1.0, 0.5) # Intercept, x1 (log scale)
#' #' delta_coef <- c(1.5, -0.7) # Intercept, x1 (log scale)
#' #'
#' #' # Generate parameters through linear predictors and log link
#' #' X_beta_eg <- model.matrix(~x1, data.frame(x1 = x1))
#' #' gamma_true <- exp(X_beta_eg %*% gamma_coef)
#' #' delta_true <- exp(X_beta_eg %*% delta_coef)
#' #'
#' #' # Generate response from Beta distribution
#' #' y <- rbeta_(n, gamma_true, delta_true)
#' #' df_beta <- data.frame(y = y, x1 = x1)
#' #'
#' #' # Fit Beta regression model using gkwreg
#' #' beta_reg <- gkwreg(y ~ x1 | x1,
#' #'   data = df_beta, family = "beta",
#' #'   link = list(gamma = "log", delta = "log")
#' #' )
#' #'
#' #' ## Example 4: Model Comparison using AIC/BIC
#' #' # Fit an alternative model (Kumaraswamy) to the same beta-generated data
#' #' kw_reg2 <- try(gkwreg(y ~ x1 | x1, data = df_beta, family = "kw"))
#' #'
#' #' # Compare models using information criteria
#' #' print("AIC Comparison (Beta vs Kw):")
#' #' c(AIC(beta_reg), AIC(kw_reg2))
#' #' print("BIC Comparison (Beta vs Kw):")
#' #' c(BIC(beta_reg), BIC(kw_reg2))
#' #'
#' #' ## Example 5: Prediction with Fitted Models
#' #' # Create new data for predictions
#' #' newdata <- data.frame(x1 = seq(-1, 1, length.out = 20))
#' #'
#' #' # Predict expected response (mean of the Beta distribution)
#' #' pred_response <- predict(beta_reg, newdata = newdata, type = "response")
#' #'
#' #' # Predict parameters on the scale of the link function
#' #' pred_link <- predict(beta_reg, newdata = newdata, type = "link")
#' #'
#' #' # Predict parameters on the original scale
#' #' pred_params <- predict(beta_reg, newdata = newdata, type = "parameter")
#' #'
#' #' # Visualize fitted model and data
#' #' plot(df_beta$x1, df_beta$y,
#' #'   pch = 20, col = "grey", xlab = "x1", ylab = "y",
#' #'   main = "Beta Regression Fit (using gkwreg)"
#' #' )
#' #' lines(newdata$x1, pred_response, col = "red", lwd = 2)
#' #' legend("topright", legend = "Predicted Mean", col = "red", lty = 1, lwd = 2)
#' #' }
#' #'
#' #' @references
#' #' Kumaraswamy, P. (1980). A generalized probability density function for double-bounded
#' #' random processes. \emph{Journal of Hydrology}, \strong{46}(1-2), 79-88.
#' #'
#' #' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized distributions.
#' #' \emph{Journal of Statistical Computation and Simulation}, \strong{81}(7), 883-898.
#' #'
#' #' Ferrari, S. L. P., & Cribari-Neto, F. (2004). Beta regression for modelling rates and
#' #' proportions. \emph{Journal of Applied Statistics}, \strong{31}(7), 799-815.
#' #'
#' #' Kristensen, K., Nielsen, A., Berg, C. W., Skaug, H., & Bell, B. M. (2016). TMB:
#' #' Automatic Differentiation and Laplace Approximation. \emph{Journal of Statistical
#' #' Software}, \strong{70}(5), 1-21.
#' #' (Underlying TMB package)
#' #'
#' #' Zeileis, A., Kleiber, C., Jackman, S. (2008). Regression Models for Count Data in R.
#' #' \emph{Journal of Statistical Software}, \strong{27}(8), 1-25.
#' #'
#' #'
#' #' Smithson, M., & Verkuilen, J. (2006). A Better Lemon Squeezer? Maximum-Likelihood
#' #' Regression with Beta-Distributed Dependent Variables. \emph{Psychological Methods},
#' #' \strong{11}(1), 54–71.
#' #'
#' #' @author Lopes, J. E.
#' #'
#' #' @seealso \code{\link{summary.gkwreg}}, \code{\link{predict.gkwreg}},
#' #'   \code{\link{plot.gkwreg}}, \code{\link{coef.gkwreg}}, \code{\link{vcov.gkwreg}},
#' #'   \code{\link[stats]{logLik}}, \code{\link[stats]{AIC}},
#' #'   \code{\link[Formula]{Formula}}, \code{\link[TMB]{MakeADFun}},
#' #'   \code{\link[TMB]{sdreport}}
#' #'
#' #' @keywords regression models
#' #' @author  Lopes, J. E.
#' #' @export
#' gkwreg <- function(formula,
#'                    data,
#'                    family = c("gkw", "bkw", "kkw", "ekw", "mc", "kw", "beta"),
#'                    link = NULL,
#'                    link_scale = NULL,
#'                    start = NULL,
#'                    fixed = NULL,
#'                    method = c("nlminb", "BFGS", "Nelder-Mead", "CG", "SANN", "L-BFGS-B"),
#'                    hessian = TRUE,
#'                    plot = TRUE,
#'                    conf.level = 0.95,
#'                    optimizer.control = list(),
#'                    subset = NULL,
#'                    weights = NULL,
#'                    offset = NULL,
#'                    na.action = getOption("na.action"),
#'                    contrasts = NULL,
#'                    x = FALSE,
#'                    y = TRUE,
#'                    model = TRUE,
#'                    silent = TRUE,
#'                    ...) {
#'   # Match arguments
#'   family <- match.arg(family, choices = c("gkw", "bkw", "kkw", "ekw", "mc", "kw", "beta"))
#'   method <- match.arg(method, choices = c("nlminb", "BFGS", "Nelder-Mead", "CG", "SANN", "L-BFGS-B"))
#'
#'   # Determine if we're using nlminb or optim (with specified method)
#'   use_nlminb <- method == "nlminb"
#'
#'   call <- match.call()
#'
#'   # Load Formula package for multi-part formula
#'   if (!requireNamespace("Formula", quietly = TRUE)) {
#'     stop("The 'Formula' package is required for this function. Please install it.")
#'   }
#'
#'   # Load TMB package for model fitting
#'   if (!requireNamespace("TMB", quietly = TRUE)) {
#'     stop("The 'TMB' package is required for this function. Please install it.")
#'   }
#'
#'   # Get parameter information for the specified family
#'   param_info <- .get_family_param_info(family)
#'   param_names <- param_info$names
#'   fixed_params <- param_info$fixed
#'   param_positions <- param_info$positions
#'
#'   # Convert to Formula object
#'   formula_obj <- Formula::as.Formula(formula)
#'
#'   # Process Formula object to get individual formulas for each parameter
#'   formula_list <- .process_formula_parts(formula_obj, param_names, fixed_params, data)
#'
#'   # Process link functions
#'   link_list <- .process_link(link, param_names, fixed_params)
#'
#'   # Process link scales
#'   link_scale_list <- .process_link_scale(link_scale, link_list, param_names, fixed_params)
#'
#'   # Convert link strings to integers for TMB
#'   link_ints <- .convert_links_to_int(link_list)
#'
#'   # Extract model frames, responses, and model matrices
#'   model_data <- .extract_model_data(
#'     formula_list = formula_list,
#'     data = data,
#'     subset = subset,
#'     weights = weights,
#'     na.action = na.action,
#'     offset = offset,
#'     contrasts = contrasts,
#'     original_call = call # Pass the original call for correct scoping
#'   )
#'
#'   # Validate response variable is in (0, 1)
#'   y_var <- model_data$y
#'   invisible(.validate_data(y_var, length(param_names)))
#'
#'   # Initialize result list
#'   result <- list(
#'     call = call,
#'     family = family,
#'     formula = formula,
#'     link = link_list,
#'     link_scale = link_scale_list,
#'     param_names = param_names,
#'     fixed_params = fixed_params
#'   )
#'
#'   # Process fixed parameters
#'   fixed_processed <- .process_fixed(fixed, param_names, fixed_params)
#'
#'   # Prepare TMB data with correct matrices based on family
#'   tmb_data <- .prepare_tmb_data(
#'     model_data, family, param_names, fixed_processed,
#'     link_ints, link_scale_list, y_var, param_positions
#'   )
#'
#'   # Prepare TMB parameters in the correct structure required by each family
#'   tmb_params <- .prepare_tmb_params(
#'     model_data, family, param_names, fixed_processed,
#'     param_positions
#'   )
#'
#'   # Compile and load the appropriate TMB model based on the family
#'   if (family == "beta") {
#'     dll_name <- "gkwbetareg"
#'   } else {
#'     dll_name <- paste0(family, "reg")
#'   }
#'
#'   if (!silent) {
#'     message("Using TMB model: ", dll_name)
#'   }
#'
#'   # Use the existing function to check and compile the TMB model
#'   .check_and_compile_TMB_code(dll_name, verbose = !silent)
#'
#'   # Create TMB object
#'   obj <- TMB::MakeADFun(
#'     data = tmb_data,
#'     parameters = tmb_params,
#'     DLL = dll_name,
#'     silent = silent
#'   )
#'
#'   # Set up optimizer control parameters
#'   if (use_nlminb) {
#'     default_control <- list(eval.max = 500, iter.max = 300, trace = ifelse(silent, 0, 1))
#'   } else { # optim methods
#'     default_control <- list(maxit = 500, trace = ifelse(silent, 0, 1))
#'   }
#'
#'   opt_control <- utils::modifyList(default_control, optimizer.control)
#'
#'   # Optimize the model
#'   if (use_nlminb) {
#'     if (!silent) message("Optimizing with nlminb...")
#'     opt <- tryCatch(
#'       stats::nlminb(
#'         start = obj$par,
#'         objective = obj$fn,
#'         gradient = obj$gr,
#'         control = opt_control
#'       ),
#'       error = function(e) {
#'         stop("Optimization with nlminb failed: ", e$message)
#'       }
#'     )
#'     fit_result <- list(
#'       coefficients = obj$env$last.par,
#'       loglik = -opt$objective,
#'       convergence = opt$convergence == 0,
#'       message = opt$message,
#'       iterations = opt$iterations
#'     )
#'   } else { # optim methods
#'     if (!silent) message(paste("Optimizing with optim method:", method, "..."))
#'     opt <- tryCatch(
#'       stats::optim(
#'         par = obj$par,
#'         fn = obj$fn,
#'         gr = obj$gr,
#'         method = method,
#'         control = opt_control
#'       ),
#'       error = function(e) {
#'         stop(paste("Optimization with optim method", method, "failed:", e$message))
#'       }
#'     )
#'     fit_result <- list(
#'       coefficients = opt$par,
#'       loglik = -opt$value,
#'       convergence = opt$convergence == 0,
#'       message = if (opt$convergence == 0) "Successful convergence" else "Optimization failed",
#'       iterations = opt$counts[1]
#'     )
#'   }
#'
#'   # Format coefficient names with parameter mappings
#'   coef_names <- .format_coefficient_names(param_names, model_data, param_positions)
#'   names(fit_result$coefficients) <- coef_names
#'
#'   # Calculate standard errors and covariance matrix if requested
#'   if (hessian) {
#'     if (!silent) message("Computing standard errors...")
#'     tryCatch(
#'       {
#'         sd_report <- TMB::sdreport(obj, getJointPrecision = TRUE)
#'         fit_result$se <- as.vector(sd_report$sd)
#'         fit_result$vcov <- as.matrix(sd_report$cov)
#'
#'         # Add parameter names to standard errors
#'         names(fit_result$se) <- coef_names
#'
#'         # Add confidence intervals
#'         alpha <- 1 - conf.level
#'         z_value <- stats::qnorm(1 - alpha / 2)
#'         fit_result$ci_lower <- fit_result$coefficients - z_value * fit_result$se
#'         fit_result$ci_upper <- fit_result$coefficients + z_value * fit_result$se
#'       },
#'       error = function(e) {
#'         warning("Could not compute standard errors: ", e$message)
#'         fit_result$se <- rep(NA, length(fit_result$coefficients))
#'         fit_result$vcov <- matrix(NA, length(fit_result$coefficients), length(fit_result$coefficients))
#'         fit_result$ci_lower <- fit_result$ci_upper <- rep(NA, length(fit_result$coefficients))
#'         names(fit_result$se) <- coef_names
#'       }
#'     )
#'   }
#'
#'   # Extract TMB report
#'   tmb_report <- obj$report()
#'
#'   # Extract parameter means
#'   alpha_mean <- tmb_report$alpha_mean
#'   beta_mean <- tmb_report$beta_mean
#'   gamma_mean <- tmb_report$gamma_mean
#'   delta_mean <- tmb_report$delta_mean
#'   lambda_mean <- tmb_report$lambda_mean
#'
#'   # Extract parameter vectors for each observation
#'   alphaVec <- if ("alphaVec" %in% names(tmb_report)) tmb_report$alphaVec else rep(alpha_mean, length(y_var))
#'   betaVec <- if ("betaVec" %in% names(tmb_report)) tmb_report$betaVec else rep(beta_mean, length(y_var))
#'   gammaVec <- if ("gammaVec" %in% names(tmb_report)) tmb_report$gammaVec else rep(gamma_mean, length(y_var))
#'   deltaVec <- if ("deltaVec" %in% names(tmb_report)) tmb_report$deltaVec else rep(delta_mean, length(y_var))
#'   lambdaVec <- if ("lambdaVec" %in% names(tmb_report)) tmb_report$lambdaVec else rep(lambda_mean, length(y_var))
#'
#'   # Extract fitted values
#'   fitted_values <- if ("fitted" %in% names(tmb_report)) tmb_report$fitted else rep(NA, length(y_var))
#'
#'   # Calculate residuals
#'   response_residuals <- y_var - fitted_values
#'
#'   # Calculate fit statistics
#'   n_obs <- length(y_var)
#'   npar <- length(fit_result$coefficients)
#'
#'   # Deviance, AIC, BIC
#'   nll <- -fit_result$loglik
#'   deviance <- 2.0 * nll
#'   aic <- deviance + 2.0 * npar
#'   bic <- deviance + log(n_obs) * npar
#'
#'   # Residual degrees of freedom
#'   df.residual <- n_obs - npar
#'
#'   # Additional statistics
#'   rmse <- sqrt(mean(response_residuals^2, na.rm = TRUE))
#'   sse <- sum((y_var - fitted_values)^2, na.rm = TRUE)
#'   sst <- sum((y_var - mean(y_var))^2, na.rm = TRUE)
#'   efron_r2 <- if (sst > 0) 1 - sse / sst else NA
#'   mean_absolute_error <- mean(abs(response_residuals), na.rm = TRUE)
#'
#'   # Build final result with all necessary components
#'   result <- list(
#'     call = call,
#'     family = family,
#'     formula = formula,
#'     coefficients = fit_result$coefficients,
#'     fitted.values = as.vector(fitted_values),
#'     residuals = as.vector(response_residuals),
#'     fitted_parameters = list(
#'       alpha = alpha_mean,
#'       beta = beta_mean,
#'       gamma = gamma_mean,
#'       delta = delta_mean,
#'       lambda = lambda_mean
#'     ),
#'     parameter_vectors = list(
#'       alphaVec = alphaVec,
#'       betaVec = betaVec,
#'       gammaVec = gammaVec,
#'       deltaVec = deltaVec,
#'       lambdaVec = lambdaVec
#'     ),
#'     link = link_list,
#'     link_scale = link_scale_list,
#'     param_names = param_names,
#'     fixed_params = fixed_params,
#'     loglik = fit_result$loglik,
#'     aic = aic,
#'     bic = bic,
#'     deviance = deviance,
#'     df.residual = df.residual,
#'     nobs = n_obs,
#'     npar = npar,
#'     vcov = if (hessian) fit_result$vcov else NULL,
#'     se = if (hessian) fit_result$se else NULL,
#'     convergence = fit_result$convergence,
#'     message = fit_result$message,
#'     iterations = fit_result$iterations,
#'     rmse = rmse,
#'     efron_r2 = efron_r2,
#'     mean_absolute_error = mean_absolute_error,
#'     method = method
#'   )
#'
#'   # Add extra information if requested
#'   if (x) result$x <- model_data$matrices
#'   if (y) result$y <- y_var
#'   if (model) result$model <- model_data$model
#'
#'   # Store TMB object
#'   result$tmb_object <- obj
#'
#'   # Set class for S3 methods
#'   class(result) <- "gkwreg"
#'
#'   # Return the final result
#'   return(result)
#' }
