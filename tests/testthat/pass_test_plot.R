# Tests for plot methods and profile likelihoods in the gkwreg package
library(testthat)

# Common parameters for testing
set.seed(12345)
n <- 500
tolerance <- 1e-5

# === plot.gkwfit tests ===
test_that("plot.gkwfit basic functionality works", {
  # Skip tests in non-interactive sessions
  skip_if(!interactive())

  # Generate data
  alpha <- 2
  beta <- 3
  data <- rkw(n, alpha, beta)

  # Fit model without plotting initially
  fit <- gkwfit(data, family = "kw", plot = FALSE, silent = TRUE)

  # Test basic plot functionality
  expect_no_error(plot(fit))

  # Test plot with specific panels
  expect_no_error(plot(fit, which = c(1, 3)))

  # Test plot with custom parameters
  expect_no_error(plot(fit, which = 1, main = "Custom Title", ylim = c(0, 3)))
})

test_that("plot.gkwfit profile likelihoods work", {
  # Skip tests in non-interactive sessions
  skip_if(!interactive())

  set.seed(100)
  # Generate data
  alpha <- 2
  beta <- 3
  data <- rkw(n, alpha, beta)

  # Fit model with profile likelihood
  fit <- gkwfit(data, family = "kw", profile = TRUE, npoints = 10, plot = FALSE, silent = TRUE)

  # Test profile likelihood plot
  expect_no_error(plot(fit, profile = TRUE))

  # Check that profile likelihood data exists
  expect_true(!is.null(fit$profile))
  expect_true(is.list(fit$profile))
  expect_true(all(c("alpha", "beta") %in% names(fit$profile)))

  # Check structure of profile data
  expect_true(is.list(fit$profile))
  expect_equal(ncol(fit$profile$alpha), 3)
})

test_that("plot.gkwreg diagnostic plots work", {
  # Skip tests in non-interactive sessions
  skip_if(!interactive())

  set.seed(100)
  # Generate data
  x1 <- runif(n, -2, 2)
  x2 <- rnorm(n)
  alpha_coef <- c(0.8, 0.3)
  beta_coef <- c(1.2, -0.4, 0.1)

  # Generate linear predictors and transform to parameters
  X_alpha <- model.matrix(~x1, data.frame(x1 = x1))
  X_beta <- model.matrix(~ x1 + x2, data.frame(x1 = x1, x2 = x2))
  alpha_true <- exp(X_alpha %*% alpha_coef)
  beta_true <- exp(X_beta %*% beta_coef)

  # Generate responses
  y <- rkw(n, alpha = alpha_true, beta = beta_true)
  y <- pmax(pmin(y, 1 - 1e-7), 1e-7) # Ensure y is strictly in (0, 1)
  df <- data.frame(y = y, x1 = x1, x2 = x2)

  # Fit model
  kw_reg <- gkwreg(y ~ x1 | x1 + x2, data = df, family = "kw")

  # Test individual plot types
  expect_no_error(plot(kw_reg, which = 1, type = "quantile"))
  expect_no_error(plot(kw_reg, which = 2, type = "pearson")) # Residuals vs Index
  expect_no_error(plot(kw_reg, which = 3, type = "quantile")) # QQ Plot
  expect_no_error(plot(kw_reg, which = 4)) # Scale-Location
  expect_no_error(plot(kw_reg, which = 5)) # Cook's Distance
  expect_no_error(plot(kw_reg, which = 6)) # Predicted vs Observed

  # Test ggplot2-based plots if available
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    expect_no_error(plot(kw_reg, which = 1, use_ggplot = TRUE))
  }
})

# === gkwgof tests with plots ===
test_that("gkwgof plotting functionality works", {
  # Skip tests in non-interactive sessions
  skip_if(!interactive())

  set.seed(100)

  # Generate data
  alpha <- 2
  beta <- 3
  data <- rkw(n, alpha, beta)

  # Fit model
  fit <- gkwfit(data, family = "kw", silent = TRUE)

  # Run goodness-of-fit analysis with plots
  expect_no_error(gkwgof(fit, plot = TRUE, print_summary = FALSE))

  # Test with specific plot panels
  expect_no_error(gkwgof(fit, plot = TRUE, which = c("density", "cdf"), print_summary = FALSE))

  # Test with simulated p-values
  expect_no_error(gkwgof(fit, plot = TRUE, simulate_p_values = TRUE, n_bootstrap = 50, print_summary = FALSE))
})

test_that("compare_gof_plot functionality works", {
  # Skip tests in non-interactive sessions
  skip_if(!interactive())
  set.seed(100)
  # Generate data
  data <- rkw(n, alpha = 2, beta = 3)

  # Fit different models
  fit_kw <- gkwfit(data, family = "kw", silent = TRUE)
  fit_beta <- gkwfit(data, family = "beta", silent = TRUE)
  fit_gkw <- gkwfit(data, family = "gkw", silent = TRUE)

  # Calculate GOF statistics
  gof_kw <- gkwgof(fit_kw, plot = FALSE, print_summary = FALSE)
  gof_beta <- gkwgof(fit_beta, plot = FALSE, print_summary = FALSE)
  gof_gkw <- gkwgof(fit_gkw, plot = FALSE, print_summary = FALSE)

  # Test comparison plot functionality
  expect_no_error(compare_gof_plot(list(KW = gof_kw, Beta = gof_beta, GKW = gof_gkw), plot_type = "all"))

  # Test individual plot types
  for (plot_type in c("radar", "bar", "table")) {
    expect_no_error(compare_gof_plot(list(KW = gof_kw, Beta = gof_beta), plot_type = plot_type))
  }

  # Check the return value
  plots <- compare_gof_plot(list(KW = gof_kw, Beta = gof_beta),
    plot_type = "all", return_plots = TRUE
  )
  expect_true(is.list(plots))
})

# === Likelihood function tests in cases not covered before ===
test_that("Negative log-likelihood functions handle edge cases", {
  # Generate data
  set.seed(100)
  alpha <- 2
  beta <- 3
  data <- rkw(50, alpha, beta)

  # Test data with values very close to boundaries
  edge_data <- c(1e-10, 0.5, 1 - 1e-10)
  expect_no_error(llkw(c(alpha = 2, beta = 3), edge_data))
  expect_no_error(llgkw(c(alpha = 2, beta = 3, gamma = 1, delta = 0, lambda = 1), edge_data))
})

# === nrgkw (Newton-Raphson) more detailed tests ===
test_that("nrgkw handles different optimization methods", {
  # Generate data
  set.seed(100)
  alpha <- 2
  beta <- 3
  data <- rkw(100, alpha, beta) # Smaller sample for faster computation

  # Test different optimization methods
  methods <- c("trust-region", "newton-raphson", "hybrid")
  for (method in methods) {
    result <- nrgkw(NULL, data, family = "kw", optimization_method = method, verbose = FALSE)
    expect_true(is.list(result))
    expect_true("parameters" %in% names(result))
    expect_true("loglik" %in% names(result))
    expect_true("aic" %in% names(result))
  }

  # Test with diagnostic information
  result_diag <- nrgkw(NULL, data, family = "kw", get_num_hess = TRUE, verbose = FALSE)
  expect_true("condition_number" %in% names(result_diag))
  expect_true("std_errors" %in% names(result_diag))
  expect_true("z_values" %in% names(result_diag))
  expect_true("p_values" %in% names(result_diag))
})

# === gkwreg with different link functions ===
test_that("gkwreg handles different link functions", {
  # Generate data
  set.seed(100)
  n <- 200 # Smaller sample for faster computation
  x1 <- runif(n, -1, 1)

  # Generate response with different link for beta parameter
  alpha_true <- exp(0.5 + 0.2 * x1) # log link
  beta_true <- 1 + 0.5 * x1 # identity link
  y <- rkw(n, alpha = alpha_true, beta = beta_true)
  y <- pmax(pmin(y, 1 - 1e-7), 1e-7) # Ensure y is strictly in (0, 1)
  df <- data.frame(y = y, x1 = x1)

  # Fit model with custom links
  kw_reg_custom <- gkwreg(y ~ x1 | x1,
    data = df, family = "kw",
    link = list(alpha = "log", beta = "identity")
  )

  # Check basic properties
  expect_s3_class(kw_reg_custom, "gkwreg")
  expect_equal(kw_reg_custom$link$alpha, "log")
  expect_equal(kw_reg_custom$link$beta, "identity")

  pred <- predict(kw_reg_custom, type = "response")
  expect_length(pred, nrow(df))
  expect_true(all(pred > 0 & pred < 1))
})


# === comprehensive family comparison case ===
test_that("All families can be fitted to the same data", {
  # Generate data from a generic GKw distribution
  set.seed(100)
  alpha <- 2
  beta <- 3
  gamma <- 1.5
  delta <- 0.5
  lambda <- 1.2
  data <- rgkw(n, alpha, beta, gamma, delta, lambda)

  # Fit all possible family models
  fit_gkw <- gkwfit(data, family = "gkw", silent = TRUE)
  fit_kw <- gkwfit(data, family = "kw", silent = TRUE)
  fit_beta <- gkwfit(data, family = "beta", silent = TRUE)
  fit_bkw <- gkwfit(data, family = "bkw", silent = TRUE)
  fit_ekw <- gkwfit(data, family = "ekw", silent = TRUE)
  fit_kkw <- gkwfit(data, family = "kkw", silent = TRUE)
  fit_mc <- gkwfit(data, family = "mc", silent = TRUE)

  # Check that all fits converged
  expect_equal(fit_gkw$convergence, TRUE)
  expect_equal(fit_kw$convergence, TRUE)
  expect_equal(fit_beta$convergence, TRUE)
  expect_equal(fit_bkw$convergence, TRUE)
  expect_equal(fit_ekw$convergence, TRUE)
  expect_equal(fit_kkw$convergence, TRUE)
  expect_equal(fit_mc$convergence, TRUE)

  # Compare likelihoods and information criteria
  models <- list(
    gkw = fit_gkw, kw = fit_kw, beta = fit_beta,
    bkw = fit_bkw, ekw = fit_ekw, kkw = fit_kkw, mc = fit_mc
  )

  likelihoods <- sapply(models, function(m) logLik(m)[1])
  expect_true(likelihoods["gkw"] >= max(likelihoods[names(likelihoods) != "gkw"]))

  aics <- sapply(models, AIC)
  expect_true(is.numeric(aics))

  bics <- sapply(models, BIC)
  expect_true(is.numeric(bics))
})

# === tests for more complex special cases ===
test_that("Newton-Raphson with different starting points converges to same solution", {
  # Generate data
  alpha <- 2
  beta <- 3
  data <- rkw(100, alpha, beta) # Smaller sample for faster computation

  # Different starting points
  start1 <- c(alpha = 1, beta = 1)
  start2 <- c(alpha = 3, beta = 4)

  # Fit models with different starting points
  fit1 <- nrgkw(start1, data, family = "kw", verbose = FALSE)
  fit2 <- nrgkw(start2, data, family = "kw", verbose = FALSE)

  # Check that both converge to similar solutions
  expect_equal(fit1$parameters["alpha"], fit2$parameters["alpha"], tolerance = 0.1)
  expect_equal(fit1$parameters["beta"], fit2$parameters["beta"], tolerance = 0.1)
})


# === gkwfit profile likelihood tests ===
test_that("gkwfit with profile=TRUE creates valid profile likelihoods", {
  # Generate data
  set.seed(100)
  alpha <- 2
  beta <- 3
  data <- rkw(100, alpha, beta) # Smaller sample for faster computation

  # Fit with profile likelihood
  fit <- gkwfit(data, family = "kw", profile = TRUE, npoints = 10, plot = FALSE, silent = TRUE)

  # Check profile likelihood structure
  expect_true(is.list(fit$profile))
  expect_true(all(c("alpha", "beta") %in% names(fit$profile)))

  # Check alpha profile
  expect_true(is.data.frame(fit$profile$alpha))
  expect_equal(ncol(fit$profile$alpha), 3)


  # Check that minimum negative log-likelihood is at the MLE
  alpha_mle_idx <- which(fit$profile$alpha[, 1] == alpha_mle)
  alpha_mle_nll <- fit$profile$alpha[alpha_mle_idx, 2]
  expect_true(all(fit$profile$alpha[, 2] >= alpha_mle_nll))
})
