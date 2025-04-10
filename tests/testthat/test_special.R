# Tests for package integration and special cases
library(testthat)

# Common parameters for testing
set.seed(12345)
n <- 500
tolerance <- 1e-5

# === Integration tests across multiple components ===
test_that("Workflow combining fitting and regression works end-to-end", {
  # 1. Generate data from a specific distribution
  alpha <- 2
  beta <- 3
  raw_data <- rkw(n, alpha, beta)

  # 2. Fit univariate distribution
  fit <- gkwfit(raw_data, family = "kw", silent = TRUE)

  # 3. Generate data for regression model using fitted parameters
  x1 <- runif(n, -1, 1)
  x2 <- rnorm(n)

  # Use fitted parameters as starting points for regression
  alpha_true <- exp(fit$coefficients["alpha"] * (1 + 0.2 * x1))
  beta_true <- exp(fit$coefficients["beta"] * (1 - 0.1 * x2))
  y <- rkw(n, alpha = alpha_true, beta = beta_true)
  y <- pmax(pmin(y, 1 - 1e-7), 1e-7)
  df <- data.frame(y = y, x1 = x1, x2 = x2)

  # 4. Fit regression model
  kw_reg <- gkwreg(y ~ x1 | x2, data = df, family = "kw")

  # 5. Calculate goodness-of-fit statistics
  gof <- gkwgof(fit, plot = FALSE, print_summary = FALSE)

  pred <- predict(kw_reg, newdata = NULL, type = "response")

  # Check that all components work together
  expect_s3_class(fit, "gkwfit")
  expect_s3_class(kw_reg, "gkwreg")
  expect_true(is.list(gof))
  expect_true(is.numeric(pred))
})

test_that("Model comparison workflow across different families works", {
  # Generate data from a more complex distribution (GKw)
  alpha <- 2
  beta <- 3
  gamma <- 1.5
  delta <- 0.5
  lambda <- 1.2
  data <- rgkw(n, alpha, beta, gamma, delta, lambda)

  # Fit various families
  fit_gkw <- gkwfit(data, family = "gkw", silent = TRUE)
  fit_kw <- gkwfit(data, family = "kw", silent = TRUE)
  fit_beta <- gkwfit(data, family = "beta", silent = TRUE)

  # Calculate goodness-of-fit for each model
  gof_gkw <- gkwgof(fit_gkw, plot = FALSE, print_summary = FALSE)
  gof_kw <- gkwgof(fit_kw, plot = FALSE, print_summary = FALSE)
  gof_beta <- gkwgof(fit_beta, plot = FALSE, print_summary = FALSE)

  # Extract GOF statistics
  gof_stats <- extract_gof_stats(gof_gkw, gof_kw, gof_beta)

  # Compare models using AIC and BIC
  aic_models <- AIC(fit_gkw)
  bic_models <- BIC(fit_gkw)

  # Perform likelihood ratio tests for nested models
  lrt_kw_gkw <- anova(fit_kw, fit_gkw)

  # Check that all components work together
  expect_s3_class(lrt_kw_gkw, "anova")
})

# === Special case: handling extreme parameter values ===
test_that("gkwfit handles extreme parameter values", {
  # Generate data with extreme parameter values
  n <- 1000
  alpha <- 0.1 # Very small alpha
  beta <- 10 # Large beta
  data <- rkw(n, alpha, beta)

  # Generate data with large gamma and delta
  gamma <- 5 # Large gamma
  delta <- 10 # Large delta
  data_beta <- rbeta_(n, gamma, delta)

  # Fit model
  fit_beta <- gkwfit(data_beta, family = "beta", silent = TRUE)

  # Check convergence and parameter estimates
  expect_equal(fit_beta$convergence, TRUE)
  expect_lt(abs(fit_beta$coefficients["gamma"] - gamma) / gamma, 0.5)
  expect_lt(abs(fit_beta$coefficients["delta"] - delta) / delta, 0.5)
})



# === Special case: data generation consistency ===
test_that("Random generation functions produce consistent results with set.seed", {
  # Test each random generation function
  seed <- 12345
  n_test <- 10

  # Kw distribution
  set.seed(seed)
  kw_samples1 <- rkw(n_test, alpha = 2, beta = 3)
  set.seed(seed)
  kw_samples2 <- rkw(n_test, alpha = 2, beta = 3)
  expect_equal(kw_samples1, kw_samples2)

  # GKw distribution
  set.seed(seed)
  gkw_samples1 <- rgkw(n_test, alpha = 2, beta = 3, gamma = 1.5, delta = 0.5, lambda = 1.2)
  set.seed(seed)
  gkw_samples2 <- rgkw(n_test, alpha = 2, beta = 3, gamma = 1.5, delta = 0.5, lambda = 1.2)
  expect_equal(gkw_samples1, gkw_samples2)

  # Beta_ distribution
  set.seed(seed)
  beta_samples1 <- rbeta_(n_test, gamma = 2, delta = 3)
  set.seed(seed)
  beta_samples2 <- rbeta_(n_test, gamma = 2, delta = 3)
  expect_equal(beta_samples1, beta_samples2)

  # Mc distribution
  set.seed(seed)
  mc_samples1 <- rmc(n_test, gamma = 2, delta = 3, lambda = 1.5)
  set.seed(seed)
  mc_samples2 <- rmc(n_test, gamma = 2, delta = 3, lambda = 1.5)
  expect_equal(mc_samples1, mc_samples2)
})

# === Special case: parameter relationships in nested models ===
test_that("Fitted parameters follow expected relationships in nested models", {
  set.seed(123)
  # Generate data from Kw distribution
  alpha <- 2
  beta <- 3
  data <- rkw(n, alpha, beta)

  # Fit more general models
  fit_kw <- gkwfit(data, family = "kw", silent = TRUE)
  fit_ekw <- gkwfit(data, family = "ekw", silent = TRUE)
  fit_bkw <- gkwfit(data, family = "bkw", silent = TRUE)
  fit_gkw <- gkwfit(data, family = "gkw", silent = TRUE)

  # Kw is EKw with lambda = 1
  expect_lt(abs(fit_ekw$coefficients["lambda"] - 1), 1)

  # Kw is BKw with gamma = 1, delta = 0
  expect_lt(abs(fit_bkw$coefficients["gamma"] - 1), 5)
  expect_lt(abs(fit_bkw$coefficients["delta"]), 5)

  # Kw is GKw with gamma = 1, delta = 0, lambda = 1
  expect_lt(abs(fit_gkw$coefficients["gamma"] - 1), 5)
  expect_lt(abs(fit_gkw$coefficients["delta"]), 5)
  expect_lt(abs(fit_gkw$coefficients["lambda"] - 1), 5)
})

# === Special case: prediction examples ===
test_that("predict.gkwreg provides various useful outputs", {
  # Generate regression data
  n <- 200
  x1 <- runif(n, -1, 1)
  x2 <- rnorm(n)
  alpha_true <- exp(0.5 + 0.2 * x1)
  beta_true <- exp(0.8 - 0.3 * x1 + 0.1 * x2)
  y <- rkw(n, alpha = alpha_true, beta = beta_true)
  y <- pmax(pmin(y, 1 - 1e-7), 1e-7)
  df <- data.frame(y = y, x1 = x1, x2 = x2)

  # Fit model
  kw_reg <- gkwreg(y ~ x1 | x1 + x2, data = df, family = "kw")

  # Create prediction grid
  newdata <- expand.grid(
    x1 = seq(-1, 1, length.out = 5),
    x2 = seq(-2, 2, length.out = 5)
  )

  # Test different prediction types

  # 1. Expected response
  pred_mean <- predict(kw_reg, type = "response")

  expect_length(pred_mean, nrow(df))
  expect_true(all(pred_mean > 0 & pred_mean < 1))

  # 3. Distribution parameters
  pred_params <- predict(kw_reg, NULL, type = "parameter")
  expect_true(is.list(pred_params))
  expect_true(all(c("alpha", "beta") %in% names(pred_params)))
  expect_length(pred_params$alpha, nrow(df))
  expect_length(pred_params$beta, nrow(df))
  expect_true(all(pred_params$alpha > 0))
  expect_true(all(pred_params$beta > 0))

  # 4. Quantiles at various probabilities
  probs <- c(0.1, 0.5, 0.9)
  for (p in probs) {
    pred_quant <- predict(kw_reg, NULL, type = "quantile", at = p)
    expect_length(pred_quant, nrow(df))
    expect_true(all(pred_quant > 0 & pred_quant < 1))
  }

  # 5. Density at various points
  points <- c(0.2, 0.5, 0.8)
  pred_dens <- predict(kw_reg, NULL, type = "density", at = points)
  expect_true(is.matrix(pred_dens))
  expect_equal(nrow(pred_dens), nrow(df))
  expect_equal(ncol(pred_dens), length(points))
  expect_true(all(pred_dens >= 0))
})

# === Special case: using fixed parameters in regression ===
test_that("gkwreg with intercept-only formulas fixes parameters across observations", {
  # Generate regression data
  n <- 200
  x1 <- runif(n, -1, 1)
  x2 <- rnorm(n)

  # Only vary alpha with x1, keep beta constant
  alpha_true <- exp(0.5 + 0.2 * x1)
  beta_true <- rep(exp(0.8), n) # Constant beta

  y <- rkw(n, alpha = alpha_true, beta = beta_true)
  y <- pmax(pmin(y, 1 - 1e-7), 1e-7)
  df <- data.frame(y = y, x1 = x1, x2 = x2)

  # Fit model with intercept-only for beta
  kw_reg <- gkwreg(y ~ x1 | 1, data = df, family = "kw", x = TRUE)

  # Check model matrices
  expect_equal(ncol(kw_reg$x$alpha), 2) # Intercept, x1
  expect_equal(ncol(kw_reg$x$beta), 1) # Intercept only

  # Check coefficient names
  coef_names <- names(coef(kw_reg))
  expect_true("alpha:(Intercept)" %in% coef_names)
  expect_true("alpha:x1" %in% coef_names)
  expect_true("beta:(Intercept)" %in% coef_names)
  expect_equal(length(coef_names), 3)

  # Verify that predictions vary only in alpha, not beta
  newdata <- data.frame(x1 = c(-1, 0, 1), x2 = c(0, 0, 0))
  pred_params <- predict(kw_reg, type = "parameter")

  # Alpha should vary
  expect_true(length(unique(pred_params$alpha)) > 3)
  # Beta should be constant
  expect_true(length(unique(pred_params$beta)) == 1)
})

# === Special case: test against known theoretical values ===
test_that("Moments of distributions match theoretical values", {
  # Generate large samples for accuracy
  n_large <- 10000

  # 1. Kw distribution
  alpha <- 2
  beta <- 3

  kw_samples <- rkw(n_large, alpha, beta)

  # Expected mean for Kw(alpha, beta) ≈ beta / (alpha*beta + 1)
  expected_mean_kw <- beta / (1 + alpha * beta)
  sample_mean_kw <- mean(kw_samples)
  expect_equal(sample_mean_kw, expected_mean_kw, tolerance = 0.10)

  # 2. Beta distribution
  gamma <- 2 # shape1
  delta <- 3 # shape2 - 1

  beta_samples <- rbeta_(n_large, gamma, delta)

  # Expected mean for Beta(gamma, delta+1) = gamma / (gamma + delta + 1)
  expected_mean_beta <- gamma / (gamma + delta + 1)
  sample_mean_beta <- mean(beta_samples)
  expect_equal(sample_mean_beta, expected_mean_beta, tolerance = 0.05)

  # 3. Check variance for Beta distribution
  # Var[Beta(a,b)] = a*b / ((a+b)^2 * (a+b+1))
  expected_var_beta <- (gamma * (delta + 1)) / ((gamma + delta + 1)^2 * (gamma + delta + 2))
  sample_var_beta <- var(beta_samples)
  expect_equal(sample_var_beta, expected_var_beta, tolerance = 0.05)
})

# === Special case: boundary behavior of distributions ===
test_that("Distributions handle values near boundaries correctly", {
  # Test values very close to 0 and 1
  x_near_0 <- 1e-10
  x_near_1 <- 1 - 1e-10

  alpha <- 2
  beta <- 3

  # 1. Kw distribution
  # Near 0, density should approach 0 for alpha > 1
  density_near_0 <- dkw(x_near_0, alpha, beta)
  expect_true(density_near_0 > 0 && density_near_0 < 1e5) # Positive but finite

  # Near 1, density should approach 0 for beta > 1
  density_near_1 <- dkw(x_near_1, alpha, beta)
  expect_true(density_near_1 > 0 && density_near_1 < 1e5) # Positive but finite

  # 2. Quantile function behavior at boundaries
  # q(0) should be 0, q(1) should be 1
  q0 <- qkw(0, alpha, beta)
  q1 <- qkw(1, alpha, beta)
  expect_equal(q0, 0)
  expect_equal(q1, 1)

  # Values very close to 0 and 1
  q_near_0 <- qkw(1e-10, alpha, beta)
  q_near_1 <- qkw(1 - 1e-10, alpha, beta)
  expect_true(q_near_0 > 0 && q_near_0 < 0.1)
  expect_true(q_near_1 < 1 && q_near_1 > 0.9)

  # 3. CDF function behavior at boundaries
  # F(0) should be 0, F(1) should be 1
  p0 <- pkw(0, alpha, beta)
  p1 <- pkw(1, alpha, beta)
  expect_equal(p0, 0)
  expect_equal(p1, 1)
})

# === Special case: multimodal distribution fit ===
test_that("gkwfit can handle bimodal data with appropriate family", {
  # Generate bimodal data by mixing distributions
  n1 <- 250
  n2 <- 250

  set.seed(111)
  # First mode around 0.3
  data1 <- rbeta(n1, shape1 = 3, shape2 = 7)

  # Second mode around 0.7
  data2 <- rbeta(n2, shape1 = 7, shape2 = 3)

  # Combined data
  bimodal_data <- c(data1, data2)

  # Fit with different families
  fit_beta <- gkwfit(bimodal_data, family = "beta", silent = TRUE)
  fit_kw <- gkwfit(bimodal_data, family = "kw", silent = TRUE)
  fit_gkw <- gkwfit(bimodal_data, family = "gkw", silent = TRUE)

  # All models should converge
  expect_equal(fit_beta$convergence, TRUE)
  expect_equal(fit_kw$convergence, TRUE)
  expect_equal(fit_gkw$convergence, TRUE)

  # Compare log-likelihoods
  ll_beta <- logLik(fit_beta)[1]
  ll_kw <- logLik(fit_kw)[1]
  ll_gkw <- logLik(fit_gkw)[1]

  # GKw should have highest likelihood as it's the most flexible
  expect_true(ll_gkw > ll_beta && ll_gkw > ll_kw)

  # Check KS test for GKw fit (might not pass, but should be better than others)
  ks_gkw <- ks.test(bimodal_data, function(x) {
    pgkw(
      x,
      fit_gkw$coefficients["alpha"],
      fit_gkw$coefficients["beta"],
      fit_gkw$coefficients["gamma"],
      fit_gkw$coefficients["delta"],
      fit_gkw$coefficients["lambda"]
    )
  })

  ks_beta <- ks.test(bimodal_data, function(x) {
    stats::pbeta(
      x,
      fit_beta$coefficients["gamma"],
      fit_beta$coefficients["delta"] + 1
    )
  })

  # GKw KS p-value should be higher than Beta
  expect_true(ks_gkw$p.value > ks_beta$p.value)
})


# === Special case: test reproducibility of random seeds ===
test_that("Setting seeds ensures reproducible results in complex workflows", {
  # Create a complex workflow with multiple random components
  run_workflow <- function(seed) {
    set.seed(seed)

    # 1. Generate training data with noise
    n <- 200
    x1 <- runif(n, -1, 1)
    x2 <- rnorm(n)

    alpha_true <- exp(0.5 + 0.2 * x1)
    beta_true <- exp(0.8 - 0.3 * x1 + 0.1 * x2)
    y <- rkw(n, alpha = alpha_true, beta = beta_true)

    # Add some noise
    y <- y * (1 + rnorm(n, 0, 0.01))
    y <- pmax(pmin(y, 1 - 1e-7), 1e-7)
    df <- data.frame(y = y, x1 = x1, x2 = x2)

    # 2. Split into training/testing
    train_idx <- sample(1:n, 0.7 * n)
    df_train <- df[train_idx, ]
    df_test <- df[-train_idx, ]

    # 3. Fit model
    kw_reg <- gkwreg(y ~ x1 | x1 + x2, data = df_train, family = "kw")

    # 4. Make predictions
    pred_test <- predict(kw_reg, type = "response")

    # Return results
    list(coefficients = coef(kw_reg), error = kw_reg$residuals)
  }

  # Run the workflow twice with the same seed
  result1 <- run_workflow(54321)
  result2 <- run_workflow(54321)

  # Results should be identical
  expect_equal(result1$coefficients, result2$coefficients)
  expect_equal(result1$error, result2$error)

  # Run with a different seed
  result3 <- run_workflow(12345)

  # Results should be different
  expect_false(identical(result1$error, result3$error))
})
