# Tests for distribution fitting functions in the gkwreg package
library(testthat)

# Common parameters for testing
set.seed(12345)
n <- 500
tolerance <- 1e-5

# === gkwfit function tests ===
test_that("gkwfit works correctly for Kumaraswamy distribution", {
  # Generate data from Kw(2, 3)
  alpha <- 2
  beta <- 3
  data <- rkw(n, alpha, beta)

  # Fit Kw model
  fit_kw <- gkwfit(data, family = "kw", silent = TRUE)

  # Check basic properties of the fit
  expect_s3_class(fit_kw, "gkwfit")
  expect_equal(fit_kw$family, "kw")
  expect_equal(fit_kw$convergence, TRUE) # Should converge
  expect_equal(length(fit_kw$coefficients), 2) # alpha and beta

  # Check that estimated parameters are close to true values
  expect_lt(abs(fit_kw$coefficients["alpha"] - alpha), 0.5)
  expect_lt(abs(fit_kw$coefficients["beta"] - beta), 0.5)

  # Check additional fit components
  expect_true(is.numeric(fit_kw$loglik))
  expect_true(is.numeric(fit_kw$nobs))
  expect_equal(fit_kw$nobs, n)
})

test_that("gkwfit works correctly for Beta distribution", {
  # Generate data from Beta(2, 4)
  shape1 <- 2
  shape2 <- 4
  data <- stats::rbeta(n, shape1, shape2)

  # Fit Beta model (using the gamma, delta+1 parameterization)
  fit_beta <- gkwfit(data, family = "beta", silent = TRUE)

  # Check basic properties of the fit
  expect_s3_class(fit_beta, "gkwfit")
  expect_equal(fit_beta$family, "beta")
  expect_equal(fit_beta$convergence, TRUE) # Should converge
  expect_equal(length(fit_beta$coefficients), 2) # gamma and delta

  # Check that estimated parameters are close to true values
  # In the gamma, delta+1 parameterization, gamma = shape1, delta = shape2 - 1
  expect_lt(abs(fit_beta$coefficients["gamma"] - shape1), 0.5)
  expect_lt(abs(fit_beta$coefficients["delta"] - (shape2 - 1)), 0.5)
})

test_that("gkwfit works correctly for Generalized Kumaraswamy distribution", {
  # Generate data from GKw with specified parameters
  alpha <- 2
  beta <- 3
  gamma <- 1.5
  delta <- 0.5
  lambda <- 1.2
  data <- rgkw(n, alpha, beta, gamma, delta, lambda)

  # Fit GKw model
  fit_gkw <- gkwfit(data, family = "gkw", silent = TRUE)

  # Check basic properties of the fit
  expect_s3_class(fit_gkw, "gkwfit")
  expect_equal(fit_gkw$family, "gkw")
  expect_equal(fit_gkw$convergence, TRUE) # Should converge
  expect_equal(length(fit_gkw$coefficients), 5) # alpha, beta, gamma, delta, lambda

  # Check that estimated parameters are reasonably close to true values
  # Note: GKw has more parameters so estimates might be less precise
  expect_lt(abs(fit_gkw$coefficients["alpha"] - alpha), 5.0)
  expect_lt(abs(fit_gkw$coefficients["beta"] - beta), 5.0)
  expect_lt(abs(fit_gkw$coefficients["gamma"] - gamma), 5.0)
  expect_lt(abs(fit_gkw$coefficients["delta"] - delta), 5.0)
  expect_lt(abs(fit_gkw$coefficients["lambda"] - lambda), 5.0)
})

test_that("gkwfit supports fixing parameters", {
  # Generate data from GKw with specified parameters
  alpha <- 2
  beta <- 3
  gamma <- 1.5
  delta <- 0.5
  lambda <- 1.2
  data <- rgkw(n, alpha, beta, gamma, delta, lambda)

  # Fit with fixed lambda
  fixed_lambda <- 1.2
  fit_fixed <- gkwfit(data, family = "gkw", fixed = list(lambda = fixed_lambda), silent = TRUE)

  # Check that lambda is fixed at the specified value
  expect_equal(fit_fixed$coefficients[["lambda"]], fixed_lambda)
  expect_equal(length(fit_fixed$coefficients), 5) # Still reports all parameters
  expect_equal(length(fit_fixed$std.errors), 5) # SEs only for non-fixed parameters
})

test_that("gkwfit handles different optimization methods", {
  # Generate data from Kw
  alpha <- 2
  beta <- 3
  data <- rkw(n, alpha, beta)

  # Test different optimization methods
  fit_tmb <- gkwfit(data, family = "kw", fit = "tmb", silent = TRUE)
  fit_nlminb <- gkwfit(data, family = "kw", fit = "nlminb", silent = TRUE)

  # Check that both methods converge
  expect_equal(fit_tmb$convergence, TRUE)
  expect_equal(fit_nlminb$convergence, TRUE)

  # Check that estimates are similar
  expect_lt(abs(fit_tmb$coefficients["alpha"] - fit_nlminb$coefficients["alpha"]), 0.1)
  expect_lt(abs(fit_tmb$coefficients["beta"] - fit_nlminb$coefficients["beta"]), 0.1)
})

# === nrgkw (Newton-Raphson optimization) tests ===
test_that("nrgkw works for Kumaraswamy distribution", {
  # Generate data from Kw
  alpha <- 2
  beta <- 3
  data <- rkw(n, alpha, beta)

  # Use Newton-Raphson for optimization
  nr_result <- nrgkw(NULL, data, family = "kw", verbose = FALSE, get_num_hess = TRUE)

  # Check basic properties
  expect_true(is.list(nr_result))
  expect_true("parameters" %in% names(nr_result))
  expect_true("loglik" %in% names(nr_result))
  expect_true("aic" %in% names(nr_result))

  # Check parameter estimates
  expect_lt(abs(nr_result$parameters[[1]] - alpha), 0.5)
  expect_lt(abs(nr_result$parameters[[2]] - beta), 0.5)
})

# === logLik, AIC, BIC methods for gkwfit tests ===
test_that("logLik method works correctly for gkwfit", {
  # Generate data and fit model
  data <- rkw(n, alpha = 2, beta = 3)
  fit <- gkwfit(data, family = "kw", silent = TRUE)

  # Test logLik method
  ll <- logLik(fit)

  # Check class and attributes
  expect_s3_class(ll, "logLik")
  expect_true(is.numeric(ll))
  expect_equal(attr(ll, "nobs"), n)
  expect_equal(attr(ll, "df"), 2) # Kw has 2 parameters
})

test_that("AIC method works correctly for gkwfit", {
  # Generate data and fit models
  data <- rkw(n, alpha = 2, beta = 3)
  fit_kw <- gkwfit(data, family = "kw", silent = TRUE)
  fit_beta <- gkwfit(data, family = "beta", silent = TRUE)

  # Test AIC for a single model
  aic_kw <- AIC(fit_kw)
  expect_true(is.numeric(aic_kw))
  expect_length(aic_kw, 1)

  # Test AIC for model comparison
  aic_comparison <- AIC(fit_kw, fit_beta)
  expect_s3_class(aic_comparison, "data.frame")
  expect_equal(nrow(aic_comparison), 2)
  expect_true("AIC" %in% colnames(aic_comparison))
  expect_true("df" %in% colnames(aic_comparison))
})

test_that("BIC method works correctly for gkwfit", {
  # Generate data and fit models
  data <- rkw(n, alpha = 2, beta = 3)
  fit_kw <- gkwfit(data, family = "kw", silent = TRUE)
  fit_beta <- gkwfit(data, family = "beta", silent = TRUE)

  # Test BIC for a single model
  bic_kw <- BIC(fit_kw)
  expect_true(is.numeric(bic_kw))
  expect_length(bic_kw, 1)

  # Test BIC for model comparison
  bic_comparison <- BIC(fit_kw, fit_beta)
  expect_s3_class(bic_comparison, "data.frame")
  expect_equal(nrow(bic_comparison), 2)
  expect_true("BIC" %in% colnames(bic_comparison))
  expect_true("df" %in% colnames(bic_comparison))
})

# === coef, confint, vcov, anova methods for gkwfit tests ===
test_that("coef method works correctly for gkwfit", {
  # Generate data and fit model
  data <- rkw(n, alpha = 2, beta = 3)
  fit <- gkwfit(data, family = "kw", silent = TRUE)

  # Test coef method
  coefficients <- coef(fit)

  # Check properties
  expect_true(is.numeric(coefficients))
  expect_equal(length(coefficients), 2) # Kw has 2 parameters
  expect_true("alpha" %in% names(coefficients))
  expect_true("beta" %in% names(coefficients))
})

test_that("confint method works correctly for gkwfit", {
  # Generate data and fit model
  data <- rkw(n, alpha = 2, beta = 3)
  fit <- gkwfit(data, family = "kw", silent = TRUE, hessian = TRUE)

  # Test confint method
  ci <- confint(fit)

  # Check properties
  expect_true(is.matrix(ci))
  expect_equal(nrow(ci), 2) # Kw has 2 parameters
  expect_equal(ncol(ci), 2) # Lower and upper bounds
  expect_equal(rownames(ci), c("alpha", "beta"))

  # Check that CIs are reasonable
  expect_true(all(ci[, 1] < ci[, 2])) # Lower bound < upper bound
  expect_true(fit$coefficients["alpha"] > ci["alpha", 1] && fit$coefficients["alpha"] < ci["alpha", 2])
  expect_true(fit$coefficients["beta"] > ci["beta", 1] && fit$coefficients["beta"] < ci["beta", 2])
})

test_that("vcov method works correctly for gkwfit", {
  # Generate data and fit model
  data <- rkw(n, alpha = 2, beta = 3)
  fit <- gkwfit(data, family = "kw", silent = TRUE, hessian = TRUE)

  # Test vcov method
  cov_matrix <- vcov(fit)

  # Check properties
  expect_true(is.matrix(cov_matrix))
  expect_equal(dim(cov_matrix), c(2, 2)) # Kw has 2 parameters
  expect_equal(rownames(cov_matrix), c("alpha", "beta"))
  expect_equal(colnames(cov_matrix), c("alpha", "beta"))

  # Check that the matrix is symmetric and positive definite
  expect_true(isSymmetric(cov_matrix))
  eigenvalues <- eigen(cov_matrix)$values
  expect_true(all(eigenvalues > 0))
})

test_that("anova method works correctly for gkwfit", {
  # Generate data
  data <- rgkw(n, alpha = 2, beta = 3, gamma = 1.5, delta = 0.5, lambda = 1.2)

  # Fit nested models
  fit_kw <- gkwfit(data, family = "kw", silent = TRUE)
  fit_gkw <- gkwfit(data, family = "gkw", silent = TRUE)

  # Test anova method
  anova_result <- anova(fit_kw, fit_gkw)

  # Check properties
  expect_s3_class(anova_result, "anova")
  expect_equal(nrow(anova_result), 2)
  expect_true(all(c("N.Par", "AIC", "BIC", "LogLik", "Test", "LR stat.", "Pr(>Chi)") %in% colnames(anova_result)))

  # Check that nested models have increasing numbers of parameters
  expect_lt(anova_result$N.Par[1], anova_result$N.Par[[2]])
})

# === print and summary methods for gkwfit tests ===
test_that("print method works for gkwfit", {
  # Generate data and fit model
  data <- rkw(n, alpha = 2, beta = 3)
  fit <- gkwfit(data, family = "kw", silent = TRUE)

  # Capture print output
  output <- capture.output(print(fit))

  # Basic checks on output
  expect_true(length(output) > 0)
  expect_true(any(grepl("Kumaraswamy", output) | grepl("kw", output)))
  expect_true(any(grepl("alpha", output)))
  expect_true(any(grepl("beta", output)))
})

test_that("summary method works for gkwfit", {
  # Generate data and fit model
  data <- rkw(n, alpha = 2, beta = 3)
  fit <- gkwfit(data, family = "kw", silent = TRUE, hessian = TRUE)

  # Get summary
  summary_obj <- summary(fit)

  # Capture summary output
  output <- capture.output(print(summary_obj))

  # Basic checks on output
  expect_true(length(output) > 0)
  expect_true(any(grepl("Family", output)))
  expect_true(any(grepl("Coefficients", output)))
  expect_true(any(grepl("alpha", output)))
  expect_true(any(grepl("beta", output)))
  expect_true(any(grepl("Std. Error", output)))
})

# === gkwgof (goodness-of-fit) tests ===
test_that("gkwgof works correctly", {
  # Generate data
  data <- rkw(n, alpha = 2, beta = 3)
  fit <- gkwfit(data, family = "kw", silent = TRUE)

  # Run goodness-of-fit analysis
  gof <- gkwgof(fit, plot = FALSE, print_summary = FALSE)

  # Check basic properties
  expect_true(is.list(gof))
  expect_true("information_criteria" %in% names(gof))
  expect_true("distance_tests" %in% names(gof))
  expect_true("probability_plots" %in% names(gof))

  # Check information criteria
  expect_true(all(c("AIC", "BIC", "AICc") %in% names(gof$information_criteria)))

  # Check GOF tests
  expect_true(all(c("ks", "cvm", "ad") %in% names(gof$distance_tests)))
})

test_that("extract_gof_stats works correctly", {
  # Generate data and fit multiple models
  data <- rkw(n, alpha = 2, beta = 3)
  fit_kw <- gkwfit(data, family = "kw", silent = TRUE)
  fit_beta <- gkwfit(data, family = "beta", silent = TRUE)

  # Calculate GOF statistics
  gof_kw <- gkwgof(fit_kw, plot = FALSE, print_summary = FALSE)
  gof_beta <- gkwgof(fit_beta, plot = FALSE, print_summary = FALSE)

  # Extract and compare statistics
  stats <- extract_gof_stats(gof_kw, gof_beta)

  # Check basic properties
  expect_s3_class(stats, "data.frame")
  expect_equal(nrow(stats), 2)
  expect_true(all(c("AIC", "BIC", "KS", "CvM", "AD") %in% colnames(stats)))
})

# test_that("compare_gof_plot works correctly", {
#   # Skip this test in non-interactive sessions or when ggplot2 isn't available
#   skip_if_not_installed("ggplot2")
#   skip_if(!interactive())
#
#   # Generate data and fit multiple models
#   data <- rkw(n, alpha = 2, beta = 3)
#   fit_kw <- gkwfit(data, family = "kw", silent = TRUE)
#   fit_beta <- gkwfit(data, family = "beta", silent = TRUE)
#
#   # Calculate GOF statistics
#   gof_kw <- gkwgof(fit_kw, plot = FALSE, print_summary = FALSE)
#   gof_beta <- gkwgof(fit_beta, plot = FALSE, print_summary = FALSE)
# })
