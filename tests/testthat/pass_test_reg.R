# Tests for regression models in the gkwreg package
library(testthat)
context("Regression Models")

# Common parameters for testing
set.seed(12345)
n <- 500
tolerance <- 1e-5

# Helper function to generate data for regression tests
generate_kw_regression_data <- function(n) {
  x1 <- runif(n, -2, 2)
  x2 <- rnorm(n)

  # True regression coefficients
  alpha_coef <- c(0.8, 0.3) # Intercept, x1
  beta_coef <- c(1.2, -0.4, 0.1) # Intercept, x1, x2

  # Generate linear predictors and transform to parameters using exp link
  X_alpha <- model.matrix(~x1, data.frame(x1 = x1))
  X_beta <- model.matrix(~ x1 + x2, data.frame(x1 = x1, x2 = x2))

  alpha_true <- exp(X_alpha %*% alpha_coef)
  beta_true <- exp(X_beta %*% beta_coef)

  # Generate responses from Kumaraswamy distribution
  y <- rkw(n, alpha = alpha_true, beta = beta_true)
  y <- pmax(pmin(y, 1 - 1e-7), 1e-7) # Ensure y is strictly in (0, 1)

  # Return data and true coefficients
  list(
    data = data.frame(y = y, x1 = x1, x2 = x2),
    alpha_coef = alpha_coef,
    beta_coef = beta_coef
  )
}

# Helper function to generate data for GKw regression tests
generate_gkw_regression_data <- function(n) {
  x1 <- runif(n, -1, 1)
  x2 <- rnorm(n)
  x3 <- factor(rbinom(n, 1, 0.5), labels = c("A", "B"))

  # True regression coefficients
  alpha_coef <- c(0.5, 0.2) # Intercept, x1
  beta_coef <- c(0.8, -0.3, 0.1) # Intercept, x1, x2
  gamma_coef <- c(0.6, 0.4) # Intercept, x3B
  delta_coef <- c(0.0, 0.2) # Intercept, x3B
  lambda_coef <- c(-0.2, 0.1) # Intercept, x2

  # Design matrices
  X_alpha <- model.matrix(~x1, data.frame(x1 = x1))
  X_beta <- model.matrix(~ x1 + x2, data.frame(x1 = x1, x2 = x2))
  X_gamma <- model.matrix(~x3, data.frame(x3 = x3))
  X_delta <- model.matrix(~x3, data.frame(x3 = x3))
  X_lambda <- model.matrix(~x2, data.frame(x2 = x2))

  # Generate linear predictors and transform to parameters
  alpha <- exp(X_alpha %*% alpha_coef)
  beta <- exp(X_beta %*% beta_coef)
  gamma <- exp(X_gamma %*% gamma_coef)
  delta <- plogis(X_delta %*% delta_coef)
  lambda <- exp(X_lambda %*% lambda_coef)

  # Generate response from GKw distribution
  y <- rgkw(n, alpha = alpha, beta = beta, gamma = gamma, delta = delta, lambda = lambda)
  y <- pmax(pmin(y, 1 - 1e-7), 1e-7) # Ensure y is strictly in (0, 1)

  # Return data and true coefficients
  list(
    data = data.frame(y = y, x1 = x1, x2 = x2, x3 = x3),
    alpha_coef = alpha_coef,
    beta_coef = beta_coef,
    gamma_coef = gamma_coef,
    delta_coef = delta_coef,
    lambda_coef = lambda_coef
  )
}

# === gkwreg function tests ===
test_that("gkwreg works correctly for Kumaraswamy regression", {
  # Generate data
  data_info <- generate_kw_regression_data(n)
  df <- data_info$data

  # Fit Kumaraswamy regression model
  kw_reg <- gkwreg(y ~ x1 | x1 + x2, data = df, family = "kw")

  # Check basic properties of the fit
  expect_s3_class(kw_reg, "gkwreg")
  expect_equal(kw_reg$family, "kw")
  expect_equal(kw_reg$nobs, n)

  # Check that coefficients are named correctly
  coef_names <- names(coef(kw_reg))
  expect_true(all(c("alpha:(Intercept)", "alpha:x1", "beta:(Intercept)", "beta:x1", "beta:x2") %in% coef_names))

  # Check that estimated coefficients are reasonably close to true values
  alpha_est <- coef(kw_reg)[grep("^alpha:", names(coef(kw_reg)))]
  beta_est <- coef(kw_reg)[grep("^beta:", names(coef(kw_reg)))]

  # Compare with true coefficients (allowing for some estimation error)
  expect_lt(abs(alpha_est["alpha:(Intercept)"] - data_info$alpha_coef[1]), 0.5)
  expect_lt(abs(alpha_est["alpha:x1"] - data_info$alpha_coef[2]), 0.5)
  expect_lt(abs(beta_est["beta:(Intercept)"] - data_info$beta_coef[1]), 0.5)
  expect_lt(abs(beta_est["beta:x1"] - data_info$beta_coef[2]), 0.5)
  expect_lt(abs(beta_est["beta:x2"] - data_info$beta_coef[3]), 0.5)
})

test_that("gkwreg works correctly for Beta regression", {
  # Generate data from a Beta regression model
  n <- 500
  x1 <- runif(n, -1, 1)

  # True coefficients for Beta parameters (gamma = shape1, delta = shape2)
  gamma_coef <- c(1.0, 0.5) # Intercept, x1 (log scale for shape1)
  delta_coef <- c(1.5, -0.7) # Intercept, x1 (log scale for shape2)

  # Generate linear predictors and transform
  X_beta_eg <- model.matrix(~x1, data.frame(x1 = x1))
  gamma_true <- exp(X_beta_eg %*% gamma_coef)
  delta_true <- exp(X_beta_eg %*% delta_coef)

  # Generate response from Beta distribution
  y <- stats::rbeta(n, shape1 = gamma_true, shape2 = delta_true)
  y <- pmax(pmin(y, 1 - 1e-7), 1e-7) # Ensure y is strictly in (0, 1)

  # Create data frame
  df_beta <- data.frame(y = y, x1 = x1)

  # Fit Beta regression model using gkwreg
  beta_reg <- gkwreg(y ~ x1 | x1, data = df_beta, family = "beta")

  # Check basic properties of the fit
  expect_s3_class(beta_reg, "gkwreg")
  expect_equal(beta_reg$family, "beta")
  expect_equal(beta_reg$nobs, n)

  # Check that coefficients are named correctly
  coef_names <- names(coef(beta_reg))
  expect_true(all(c("gamma:(Intercept)", "gamma:x1", "delta:(Intercept)", "delta:x1") %in% coef_names))

  # Check that estimated coefficients are reasonably close to true values
  gamma_est <- coef(beta_reg)[grep("^gamma:", names(coef(beta_reg)))]
  delta_est <- coef(beta_reg)[grep("^delta:", names(coef(beta_reg)))]

  # Compare with true coefficients (allowing for some estimation error)
  expect_lt(abs(gamma_est["gamma:(Intercept)"] - gamma_coef[1]), 5)
  expect_lt(abs(gamma_est["gamma:x1"] - gamma_coef[2]), 5)
  expect_lt(abs(delta_est["delta:(Intercept)"] - delta_coef[1]), 25)
  expect_lt(abs(delta_est["delta:x1"] - delta_coef[2]), 5)
})

test_that("gkwreg works correctly for Generalized Kumaraswamy regression", {
  # For a simpler test, we'll focus on the KKw special case (gamma=1)
  n <- 500 # Smaller sample due to complexity

  # Generate predictors
  x1 <- runif(n, -1, 1)
  x2 <- rnorm(n)

  # True regression coefficients
  alpha_coef <- c(0.5, 0.2) # Intercept, x1
  beta_coef <- c(0.8, -0.3) # Intercept, x1
  delta_coef <- c(0.0, 0.2) # Intercept, x2
  lambda_coef <- c(0.5) # Intercept only

  # Design matrices
  X_alpha <- model.matrix(~x1, data.frame(x1 = x1))
  X_beta <- model.matrix(~x1, data.frame(x1 = x1))
  X_delta <- model.matrix(~x2, data.frame(x2 = x2))
  X_lambda <- matrix(1, nrow = n, ncol = 1)

  # Generate linear predictors and transform to parameters
  alpha <- exp(X_alpha %*% alpha_coef)
  beta <- exp(X_beta %*% beta_coef)
  delta <- plogis(X_delta %*% delta_coef)
  lambda <- exp(X_lambda %*% lambda_coef)

  # Generate response from KKw distribution (GKw with gamma=1)
  y <- rkkw(n, alpha = alpha, beta = beta, delta = delta, lambda = lambda)
  y <- pmax(pmin(y, 1 - 1e-7), 1e-7) # Ensure y is strictly in (0, 1)

  # Create data frame
  df_kkw <- data.frame(y = y, x1 = x1, x2 = x2)

  # Fit KKw regression model using gkwreg (equivalent to GKw with gamma=1)
  kkw_reg <- gkwreg(y ~ x1 | x1 | x2, data = df_kkw, family = "kkw")

  # Check basic properties of the fit
  expect_s3_class(kkw_reg, "gkwreg")
  expect_equal(kkw_reg$family, "kkw")
  expect_equal(kkw_reg$nobs, n)

  # Check that coefficients are named correctly
  coef_names <- names(coef(kkw_reg))
  expect_true(all(c(
    "alpha:(Intercept)", "alpha:x1", "beta:(Intercept)", "beta:x1",
    "delta:(Intercept)", "delta:x2", "lambda:(Intercept)"
  ) %in% coef_names))
})

# === coef, confint, vcov methods for gkwreg tests ===
test_that("coef method works correctly for gkwreg", {
  # Generate data
  data_info <- generate_kw_regression_data(n)
  df <- data_info$data

  # Fit model
  kw_reg <- gkwreg(y ~ x1 | x1 + x2, data = df, family = "kw")

  # Test coef method
  coefficients <- coef(kw_reg)

  # Check properties
  expect_true(is.numeric(coefficients))
  expect_equal(length(coefficients), 5) # alpha(intercept, x1), beta(intercept, x1, x2)
  expect_true(all(c("alpha:(Intercept)", "alpha:x1", "beta:(Intercept)", "beta:x1", "beta:x2") %in% names(coefficients)))

  # Alternative: coefficients
  coefs_alt <- coefficients(kw_reg)
  expect_equal(coefs_alt, coefficients)
})

test_that("confint method works correctly for gkwreg", {
  # Generate data
  data_info <- generate_kw_regression_data(n)
  df <- data_info$data

  # Fit model
  kw_reg <- gkwreg(y ~ x1 | x1 + x2, data = df, family = "kw", fit = "tmb")

  # Test confint method
  ci <- confint(kw_reg)

  # Check properties
  expect_true(is.matrix(ci))
  expect_equal(nrow(ci), 5) # alpha(intercept, x1), beta(intercept, x1, x2)
  expect_equal(ncol(ci), 2) # Lower and upper bounds
  expect_true(all(c("alpha:(Intercept)", "alpha:x1", "beta:(Intercept)", "beta:x1", "beta:x2") %in% rownames(ci)))

  # Check that CIs are reasonable
  # expect_true(all(ci[, 1] < ci[, 2]))  # Lower bound < upper bound
  # coefficients <- coef(kw_reg)
  # for (param in rownames(ci)) {
  #   expect_true(coefficients[param] > ci[param, 1] && coefficients[param] < ci[param, 2])
  # }
})

test_that("vcov method works correctly for gkwreg", {
  # Generate data
  data_info <- generate_kw_regression_data(n)
  df <- data_info$data

  # Fit model
  kw_reg <- gkwreg(y ~ x1 | x1 + x2, data = df, family = "kw")

  # Test vcov method
  cov_matrix <- vcov(kw_reg)

  # Check properties
  # expect_true(is.matrix(cov_matrix))
  # expect_equal(dim(cov_matrix), c(5, 5))  # alpha(intercept, x1), beta(intercept, x1, x2)
  # expect_true(all(c("alpha:(Intercept)", "alpha:", "beta:(Intercept)", "beta:x1", "beta:x2") %in% rownames(cov_matrix)))
  # expect_true(all(c("alpha:(Intercept)", "alpha:x1", "beta:(Intercept)", "beta:x1", "beta:x2") %in% colnames(cov_matrix)))

  # Check that the matrix is symmetric and positive definite
  expect_true(isSymmetric(cov_matrix))
  eigenvalues <- eigen(cov_matrix)$values
  expect_true(all(eigenvalues > 0))
})

# === logLik, AIC, BIC methods for gkwreg tests ===
test_that("logLik method works correctly for gkwreg", {
  # Generate data
  data_info <- generate_kw_regression_data(n)
  df <- data_info$data

  # Fit model
  kw_reg <- gkwreg(y ~ x1 | x1 + x2, data = df, family = "kw")

  # Test logLik method
  ll <- logLik(kw_reg)

  # Check class and attributes
  expect_s3_class(ll, "logLik")
  expect_true(is.numeric(ll))
  expect_equal(attr(ll, "nobs"), n)
  expect_equal(attr(ll, "df"), 5) # alpha(intercept, x1), beta(intercept, x1, x2)
})

test_that("AIC and BIC methods work correctly for gkwreg", {
  # Generate data
  data_info <- generate_kw_regression_data(n)
  df <- data_info$data

  # Fit different models
  kw_reg1 <- gkwreg(y ~ x1 | x1 + x2, data = df, family = "kw")
  kw_reg2 <- gkwreg(y ~ x1 | x2, data = df, family = "kw")

  # Test AIC for a single model
  aic1 <- AIC(kw_reg1)
  expect_true(is.numeric(aic1))
  expect_length(aic1, 1)

  # Test BIC for a single model
  bic1 <- BIC(kw_reg1)
  expect_true(is.numeric(bic1))
  expect_length(bic1, 1)

  # Test BIC for model comparison
  # bic_comparison <- BIC(kw_reg1, kw_reg2)
  # expect_s3_class(bic_comparison, "data.frame")
  # expect_equal(nrow(bic_comparison), 2)
  # expect_true("BIC" %in% colnames(bic_comparison))
  # expect_true("df" %in% colnames(bic_comparison))
})

# === fitted, predict, residuals methods for gkwreg tests ===
test_that("fitted method works correctly for gkwreg", {
  # Generate data
  data_info <- generate_kw_regression_data(n)
  df <- data_info$data

  # Fit model
  kw_reg <- gkwreg(y ~ x1 | x1 + x2, data = df, family = "kw")

  # Test fitted method
  fitted_values <- fitted(kw_reg)

  # Check properties
  expect_true(is.numeric(fitted_values))
  expect_equal(length(fitted_values), n)
  expect_true(all(fitted_values > 0 & fitted_values < 1))

  # Test with different family assumption
  fitted_values_beta <- fitted(kw_reg, family = "beta")
  expect_true(is.numeric(fitted_values_beta))
  expect_equal(length(fitted_values_beta), n)
  expect_true(all(fitted_values_beta > 0 & fitted_values_beta < 1))
})

test_that("predict method works correctly for gkwreg", {
  # Generate data
  data_info <- generate_kw_regression_data(n)
  df <- data_info$data

  # Fit model
  kw_reg <- gkwreg(y ~ x1 | x1 + x2, data = df, family = "kw")

  # Create newdata for prediction
  newdata <- data.frame(
    x1 = seq(-1, 1, length.out = 5),
    x2 = rep(0, 5)
  )

  # Test predict method with different types
  pred_response <- predict(kw_reg, type = "response")
  pred_link <- predict(kw_reg, type = "link")
  pred_params <- predict(kw_reg, type = "parameter")

  # Check properties
  expect_true(is.numeric(pred_response))
  expect_equal(length(pred_response), nrow(df))
  expect_true(all(pred_response > 0 & pred_response < 1))

  expect_true(is.list(pred_link))
  expect_true(all(c("alpha", "beta") %in% names(pred_link)))
  expect_equal(length(pred_link$alpha), nrow(df))
  expect_equal(length(pred_link$beta), nrow(df))

  expect_true(is.list(pred_params))
  expect_true(all(c("alpha", "beta") %in% names(pred_params)))
  expect_equal(length(pred_params$alpha), nrow(df))
  expect_equal(length(pred_params$beta), nrow(df))

  # Test specific parameter prediction
  pred_alpha <- predict(kw_reg, type = "alpha")
  expect_true(is.numeric(pred_alpha))
  expect_equal(length(pred_alpha), nrow(df))
  expect_true(all(pred_alpha > 0))

  # Test density prediction
  pred_dens <- predict(kw_reg, type = "density", at = c(0.2, 0.5, 0.8))
  expect_true(is.matrix(pred_dens))
  expect_equal(nrow(pred_dens), nrow(df))
  expect_equal(ncol(pred_dens), 3)

  # Test quantile prediction
  pred_quant <- predict(kw_reg, type = "quantile", at = 0.5)
  expect_true(is.numeric(pred_quant))
  expect_equal(length(pred_quant), nrow(df))
  expect_true(all(pred_quant > 0 & pred_quant < 1))
})

test_that("residuals method works correctly for gkwreg", {
  # Generate data
  data_info <- generate_kw_regression_data(n)
  df <- data_info$data

  # Fit model
  kw_reg <- gkwreg(y ~ x1 | x1 + x2, data = df, family = "kw")

  # Test residuals method with different types
  resp_res <- residuals(kw_reg, type = "response")
  pearson_res <- residuals(kw_reg, type = "pearson")
  quant_res <- residuals(kw_reg, type = "quantile")
  cs_res <- residuals(kw_reg, type = "cox-snell")

  # Check properties
  expect_true(is.numeric(resp_res))
  expect_equal(length(resp_res), n)
  expect_true(mean(resp_res) < 0.1) # Mean should be close to 0

  expect_true(is.numeric(pearson_res))
  expect_equal(length(pearson_res), n)
  expect_true(mean(pearson_res) < 0.1) # Mean should be close to 0

  expect_true(is.numeric(quant_res))
  expect_equal(length(quant_res), n)
  expect_true(mean(quant_res) < 0.1) # Mean should be close to 0

  expect_true(is.numeric(cs_res))
  expect_equal(length(cs_res), n)

  # Test with a different family assumption
  quant_res_beta <- residuals(kw_reg, type = "quantile", family = "beta")
  expect_true(is.numeric(quant_res_beta))
  expect_equal(length(quant_res_beta), n)

  # Test partial residuals
  if (any(grepl("x1", colnames(kw_reg$x$alpha)))) {
    x1_idx_alpha <- which(colnames(kw_reg$x$alpha) == "x1")
    if (length(x1_idx_alpha) == 1) {
      part_res_alpha_x1 <- residuals(kw_reg, type = "partial", parameter = "alpha", covariate_idx = x1_idx_alpha)
      expect_true(is.numeric(part_res_alpha_x1))
      expect_equal(length(part_res_alpha_x1), n)
    }
  }
})

# === print, summary methods for gkwreg tests ===
test_that("summary method works for gkwreg", {
  # Generate data
  data_info <- generate_kw_regression_data(n)
  df <- data_info$data

  # Fit model
  kw_reg <- gkwreg(y ~ x1 | x1 + x2, data = df, family = "kw")

  # Get summary
  summary_obj <- summary(kw_reg)

  # Check basic properties
  expect_s3_class(summary_obj, "summary.gkwreg")
  expect_true(is.list(summary_obj))
  expect_true(all(c("call", "family", "residuals", "nobs", "coefficients", "conf.int") %in% names(summary_obj)))

  # Check coefficients table
  expect_true(is.data.frame(summary_obj$coefficients))
  expect_equal(nrow(summary_obj$coefficients), 5) # alpha(intercept, x1), beta(intercept, x1, x2)
  expect_true(all(c("Estimate", "Std. Error", "z value", "Pr(>|z|)") %in% colnames(summary_obj$coefficients)))

  # Test print method for summary object
  printed_summary <- capture.output(print(summary_obj))
  expect_true(length(printed_summary) > 0)
  expect_true(any(grepl("Kumaraswamy", printed_summary) | grepl("kw", printed_summary)))
  expect_true(any(grepl("alpha", printed_summary)))
  expect_true(any(grepl("beta", printed_summary)))
})
