# Tests for distribution functions in the gkwreg package
library(testthat)

# Common parameters for testing
set.seed(12345)
n <- 500
tolerance <- 1e-5

# === Kumaraswamy (Kw) distribution tests ===
test_that("Kumaraswamy density function works correctly", {
  # Basic functionality
  x <- 0.5
  alpha <- 2
  beta <- 3
  density <- as.numeric(dkw(x, alpha, beta))
  expected <- alpha * beta * x^(alpha - 1) * (1 - x^alpha)^(beta - 1)
  expect_equal(density, expected, tolerance = tolerance)

  # Vectorized inputs
  x_vec <- c(0.2, 0.5, 0.8)
  densities <- dkw(x_vec, alpha, beta)
  expect_length(densities, length(x_vec))

  # Log density
  log_density <- as.numeric(dkw(x, alpha, beta, log_prob = TRUE))
  expect_equal(log_density, log(density), tolerance = tolerance)

  # Edge cases
  expect_equal(as.numeric(dkw(0, alpha, beta)), 0)
  expect_equal(as.numeric(dkw(1, alpha, beta)), 0)
})

test_that("Kumaraswamy CDF function works correctly", {
  # Basic functionality
  q <- 0.5
  alpha <- 2
  beta <- 3
  cdf <- as.numeric(pkw(q, alpha, beta))
  expected <- 1 - (1 - q^alpha)^beta
  expect_equal(cdf, expected, tolerance = tolerance)

  # Vectorized inputs
  q_vec <- c(0.2, 0.5, 0.8)
  cdfs <- pkw(q_vec, alpha, beta)
  expect_length(cdfs, length(q_vec))

  # Upper tail and log probabilities
  upper_tail <- as.numeric(pkw(q, alpha, beta, lower_tail = FALSE))
  expect_equal(upper_tail, 1 - cdf, tolerance = tolerance)

  log_cdf <- as.numeric(pkw(q, alpha, beta, log_p = TRUE))
  expect_equal(log_cdf, log(cdf), tolerance = tolerance)

  # Edge cases
  expect_equal(as.numeric(pkw(0, alpha, beta)), 0)
  expect_equal(as.numeric(pkw(1, alpha, beta)), 1)
})

test_that("Kumaraswamy quantile function works correctly", {
  # Basic functionality
  p <- 0.5
  alpha <- 2
  beta <- 3
  q <- as.numeric(qkw(p, alpha, beta))
  expect_equal(as.numeric(pkw(q, alpha, beta)), p, tolerance = tolerance)

  # Vectorized inputs
  p_vec <- c(0.2, 0.5, 0.8)
  quantiles <- as.numeric(qkw(p_vec, alpha, beta))
  expect_length(quantiles, length(p_vec))

  # Upper tail and log probabilities
  q_upper <- as.numeric(qkw(p, alpha, beta, lower_tail = FALSE))
  expect_equal(q_upper, as.numeric(qkw(1 - p, alpha, beta)), tolerance = tolerance)

  q_log <- as.numeric(qkw(log(p), alpha, beta, log_p = TRUE))
  expect_equal(q_log, q, tolerance = tolerance)

  # Edge cases
  expect_equal(as.numeric(qkw(0, alpha, beta)), 0)
  expect_equal(as.numeric(qkw(1, alpha, beta)), 1)
})

test_that("Kumaraswamy random generation works correctly", {
  # Basic functionality
  alpha <- 2
  beta <- 3
  samples <- rkw(n, alpha, beta)

  expect_length(samples, n)
  expect_true(all(samples >= 0 & samples <= 1))

  # Test that the sample mean and variance are close to theoretical values
  # For Kw(alpha, beta), mean = alpha*beta*B(1+1/alpha, beta)/(alpha*beta)
  # where B() is the beta function
  # This is an approximation for quick testing
  expect_true(abs(mean(samples) - 0.5) < 0.1)

  # Test reproducibility with set.seed
  set.seed(12345)
  samples1 <- rkw(5, alpha, beta)
  set.seed(12345)
  samples2 <- rkw(5, alpha, beta)
  expect_equal(samples1, samples2)

  # Vectorized parameters
  alphas <- c(2, 3, 4)
  samples_vec <- rkw(length(alphas), alpha = alphas, beta = beta)
  expect_length(samples_vec, length(alphas))
})

test_that("Generalized Kumaraswamy CDF function works correctly", {
  # Basic functionality
  q <- 0.5
  alpha <- 2
  beta <- 3
  gamma <- 1
  delta <- 0
  lambda <- 1

  # For these parameters, GKw reduces to Kw
  cdf_gkw <- pgkw(q, alpha, beta, gamma, delta, lambda)
  cdf_kw <- pkw(q, alpha, beta)
  expect_equal(cdf_gkw, cdf_kw, tolerance = tolerance)

  # Different parameterization
  gamma <- 2
  delta <- 1
  lambda <- 1.5
  cdf <- pgkw(q, alpha, beta, gamma, delta, lambda)
  expect_true(is.numeric(cdf) && cdf >= 0 && cdf <= 1)

  # Upper tail and log probabilities
  upper_tail <- pgkw(q, alpha, beta, gamma, delta, lambda, lower_tail = FALSE)
  expect_equal(upper_tail, 1 - cdf, tolerance = tolerance)

  log_cdf <- pgkw(q, alpha, beta, gamma, delta, lambda, log_p = TRUE)
  expect_equal(log_cdf, log(cdf), tolerance = tolerance)
})

test_that("Generalized Kumaraswamy quantile function works correctly", {
  # Basic functionality
  p <- 0.5
  alpha <- 2
  beta <- 3
  gamma <- 1
  delta <- 0
  lambda <- 1

  q <- as.numeric(qgkw(p, alpha, beta, gamma, delta, lambda))
  expect_equal(as.numeric(pgkw(q, alpha, beta, gamma, delta, lambda)), p, tolerance = tolerance)

  # Different parameterization
  gamma <- 2
  delta <- 1
  lambda <- 1.5
  q_diff <- as.numeric(qgkw(p, alpha, beta, gamma, delta, lambda))
  expect_equal(as.numeric(pgkw(q_diff, alpha, beta, gamma, delta, lambda)), p, tolerance = tolerance)
})

test_that("Generalized Kumaraswamy random generation works correctly", {
  # Basic functionality
  alpha <- 2
  beta <- 3
  gamma <- 1
  delta <- 0
  lambda <- 1

  samples <- rgkw(n, alpha, beta, gamma, delta, lambda)
  expect_length(samples, n)
  expect_true(all(samples >= 0 & samples <= 1))

  # Test that the generated samples follow the distribution
  # For these parameters, GKw reduces to Kw
  ks_test <- suppressWarnings(ks.test(samples, function(x) pkw(x, alpha, beta)))
  expect_true(ks_test$p.value > 0.01)
})

# === Beta (gamma, delta+1 parameterization) tests ===
test_that("Beta_ density function works correctly", {
  # Basic functionality
  x <- 0.5
  gamma <- 2
  delta <- 3 # This means shape2 = delta + 1 = 4 in standard parameterization

  density_beta_ <- as.numeric(dbeta_(x, gamma, delta))
  density_beta <- stats::dbeta(x, shape1 = gamma, shape2 = delta + 1)
  expect_equal(density_beta_, density_beta, tolerance = tolerance)

  # Log density
  log_density_beta_ <- as.numeric(dbeta_(x, gamma, delta, log_prob = TRUE))
  log_density_beta <- stats::dbeta(x, shape1 = gamma, shape2 = delta + 1, log = TRUE)
  expect_equal(log_density_beta_, log_density_beta, tolerance = tolerance)
})

test_that("Beta_ CDF function works correctly", {
  # Basic functionality
  q <- 0.5
  gamma <- 2
  delta <- 3 # This means shape2 = delta + 1 = 4 in standard parameterization

  cdf_beta_ <- as.numeric(pbeta_(q, gamma, delta))
  cdf_beta <- stats::pbeta(q, shape1 = gamma, shape2 = delta + 1)
  expect_equal(cdf_beta_, cdf_beta, tolerance = tolerance)

  # Upper tail and log probabilities
  upper_tail_beta_ <- as.numeric(pbeta_(q, gamma, delta, lower_tail = FALSE))
  upper_tail_beta <- stats::pbeta(q, shape1 = gamma, shape2 = delta + 1, lower.tail = FALSE)
  expect_equal(upper_tail_beta_, upper_tail_beta, tolerance = tolerance)

  log_cdf_beta_ <- as.numeric(pbeta_(q, gamma, delta, log_p = TRUE))
  log_cdf_beta <- stats::pbeta(q, shape1 = gamma, shape2 = delta + 1, log.p = TRUE)
  expect_equal(log_cdf_beta_, log_cdf_beta, tolerance = tolerance)
})

test_that("Beta_ quantile function works correctly", {
  # Basic functionality
  p <- 0.5
  gamma <- 2
  delta <- 3 # This means shape2 = delta + 1 = 4 in standard parameterization

  q_beta_ <- as.numeric(qbeta_(p, gamma, delta))
  q_beta <- stats::qbeta(p, shape1 = gamma, shape2 = delta + 1)
  expect_equal(q_beta_, q_beta, tolerance = tolerance)

  # Upper tail and log probabilities
  q_upper_beta_ <- as.numeric(qbeta_(p, gamma, delta, lower_tail = FALSE))
  q_upper_beta <- stats::qbeta(p, shape1 = gamma, shape2 = delta + 1, lower.tail = FALSE)
  expect_equal(q_upper_beta_, q_upper_beta, tolerance = tolerance)

  q_log_beta_ <- as.numeric(qbeta_(log(p), gamma, delta, log_p = TRUE))
  q_log_beta <- stats::qbeta(log(p), shape1 = gamma, shape2 = delta + 1, log.p = TRUE)
  expect_equal(q_log_beta_, q_log_beta, tolerance = tolerance)
})

test_that("Beta_ random generation works correctly", {
  # Basic functionality
  gamma <- 2
  delta <- 3 # This means shape2 = delta + 1 = 4 in standard parameterization

  samples_beta_ <- as.numeric(rbeta_(n, gamma, delta))

  expect_length(samples_beta_, n)
  expect_true(all(samples_beta_ >= 0 & samples_beta_ <= 1))

  # Test that the generated samples follow the distribution
  ks_test <- suppressWarnings(ks.test(samples_beta_, function(x) as.numeric(pbeta_(x, gamma, delta))))
  expect_true(ks_test$p.value > 0.01)
})

# === McDonald (Mc) distribution tests ===
test_that("McDonald (Mc) density function works correctly", {
  # Basic functionality
  x <- 0.5
  gamma <- 2
  delta <- 3
  lambda <- 1 # With lambda=1, Mc is equivalent to Beta(gamma, delta+1)

  density_mc <- as.numeric(dmc(x, gamma, delta, lambda))
  density_beta <- stats::dbeta(x, shape1 = gamma, shape2 = delta + 1)
  expect_equal(density_mc, density_beta, tolerance = tolerance)

  # Different lambda
  lambda <- 2
  density_mc_lambda <- as.numeric(dmc(x, gamma, delta, lambda))
  expect_true(is.numeric(density_mc_lambda) && density_mc_lambda >= 0)

  # Log density
  log_density_mc <- as.numeric(dmc(x, gamma, delta, lambda, log_prob = TRUE))
  expect_equal(log_density_mc, log(density_mc_lambda), tolerance = tolerance)
})

test_that("McDonald (Mc) CDF function works correctly", {
  # Basic functionality
  q <- 0.5
  gamma <- 2
  delta <- 3
  lambda <- 1 # With lambda=1, Mc is equivalent to Beta(gamma, delta+1)

  cdf_mc <- as.numeric(pmc(q, gamma, delta, lambda))
  cdf_beta <- stats::pbeta(q, shape1 = gamma, shape2 = delta + 1)
  expect_equal(cdf_mc, cdf_beta, tolerance = tolerance)

  # Upper tail and log probabilities
  upper_tail_mc <- as.numeric(pmc(q, gamma, delta, lambda, lower_tail = FALSE))
  expect_equal(upper_tail_mc, 1 - cdf_mc, tolerance = tolerance)

  log_cdf_mc <- as.numeric(pmc(q, gamma, delta, lambda, log_p = TRUE))
  expect_equal(log_cdf_mc, log(cdf_mc), tolerance = tolerance)
})

test_that("McDonald (Mc) quantile function works correctly", {
  # Basic functionality
  p <- 0.5
  gamma <- 2
  delta <- 3
  lambda <- 1 # With lambda=1, Mc is equivalent to Beta(gamma, delta+1)

  q_mc <- as.numeric(qmc(p, gamma, delta, lambda))
  q_beta <- stats::qbeta(p, shape1 = gamma, shape2 = delta + 1)
  expect_equal(q_mc, q_beta, tolerance = tolerance)

  # Different lambda
  lambda <- 2
  q_mc_lambda <- as.numeric(qmc(p, gamma, delta, lambda))
  expect_equal(as.numeric(pmc(q_mc_lambda, gamma, delta, lambda)), p, tolerance = tolerance)
})

test_that("McDonald (Mc) random generation works correctly", {
  # Basic functionality
  gamma <- 2
  delta <- 3
  lambda <- 1 # With lambda=1, Mc is equivalent to Beta(gamma, delta+1)

  samples_mc <- rmc(n, gamma, delta, lambda)

  expect_length(samples_mc, n)
  expect_true(all(samples_mc >= 0 & samples_mc <= 1))

  # Test that the generated samples follow the distribution
  ks_test <- suppressWarnings(ks.test(samples_mc, function(x) pmc(x, gamma, delta, lambda)))
  expect_true(ks_test$p.value > 0.01)
})

test_that("Beta-Kumaraswamy (BKw) CDF function works correctly", {
  # Basic functionality
  q <- 0.5
  alpha <- 2
  beta <- 3
  gamma <- 1
  delta <- 0.5

  cdf_bkw <- pbkw(q, alpha, beta, gamma, delta)
  # When gamma=1, BKw is equivalent to GKw with lambda=1
  cdf_gkw <- pgkw(q, alpha, beta, gamma, delta, lambda = 1)
  expect_equal(cdf_bkw, cdf_gkw, tolerance = tolerance)

  # Upper tail and log probabilities
  upper_tail_bkw <- pbkw(q, alpha, beta, gamma, delta, lower_tail = FALSE)
  expect_equal(upper_tail_bkw, 1 - cdf_bkw, tolerance = tolerance)

  log_cdf_bkw <- pbkw(q, alpha, beta, gamma, delta, log_p = TRUE)
  expect_equal(log_cdf_bkw, log(cdf_bkw), tolerance = tolerance)
})

test_that("Beta-Kumaraswamy (BKw) quantile function works correctly", {
  # Basic functionality
  p <- 0.5
  alpha <- 2
  beta <- 3
  gamma <- 1
  delta <- 0.5

  q_bkw <- as.numeric(qbkw(p, alpha, beta, gamma, delta))
  # Verify that quantile and CDF are inverses
  expect_equal(as.numeric(pbkw(q_bkw, alpha, beta, gamma, delta)), p, tolerance = tolerance)
})

test_that("Beta-Kumaraswamy (BKw) random generation works correctly", {
  # Basic functionality
  alpha <- 2
  beta <- 3
  gamma <- 1
  delta <- 0.5

  samples_bkw <- rbkw(n, alpha, beta, gamma, delta)

  expect_length(samples_bkw, n)
  expect_true(all(samples_bkw >= 0 & samples_bkw <= 1))

  # Test that the generated samples follow the distribution
  ks_test <- suppressWarnings(ks.test(samples_bkw, function(x) pbkw(x, alpha, beta, gamma, delta)))
  expect_true(ks_test$p.value > 0.01)
})

test_that("Exponentiated Kumaraswamy (EKw) CDF function works correctly", {
  # Basic functionality
  q <- 0.5
  alpha <- 2
  beta <- 3
  lambda <- 1 # With lambda=1, EKw is equivalent to Kw

  cdf_ekw <- pekw(q, alpha, beta, lambda)
  cdf_kw <- pkw(q, alpha, beta)
  expect_equal(cdf_ekw, cdf_kw, tolerance = tolerance)

  # Upper tail and log probabilities
  upper_tail_ekw <- pekw(q, alpha, beta, lambda, lower_tail = FALSE)
  expect_equal(upper_tail_ekw, 1 - cdf_ekw, tolerance = tolerance)

  log_cdf_ekw <- pekw(q, alpha, beta, lambda, log_p = TRUE)
  expect_equal(log_cdf_ekw, log(cdf_ekw), tolerance = tolerance)
})

test_that("Exponentiated Kumaraswamy (EKw) quantile function works correctly", {
  # Basic functionality
  p <- 0.5
  alpha <- 2
  beta <- 3
  lambda <- 1 # With lambda=1, EKw is equivalent to Kw

  q_ekw <- as.numeric(qekw(p, alpha, beta, lambda))
  q_kw <- as.numeric(qkw(p, alpha, beta))
  expect_equal(q_ekw, q_kw, tolerance = tolerance)

  # Different lambda
  lambda <- 2
  q_ekw_lambda <- as.numeric(qekw(p, alpha, beta, lambda))
  # Verify that quantile and CDF are inverses
  expect_equal(as.numeric(pekw(q_ekw_lambda, alpha, beta, lambda)), p, tolerance = tolerance)
})

test_that("Exponentiated Kumaraswamy (EKw) random generation works correctly", {
  # Basic functionality
  alpha <- 2
  beta <- 3
  lambda <- 2

  samples_ekw <- rekw(n, alpha, beta, lambda)

  expect_length(samples_ekw, n)
  expect_true(all(samples_ekw >= 0 & samples_ekw <= 1))

  # Test that the generated samples follow the distribution
  ks_test <- suppressWarnings(ks.test(samples_ekw, function(x) pekw(x, alpha, beta, lambda)))
  expect_true(ks_test$p.value > 0.01)
})


test_that("Kumaraswamy-Kumaraswamy (KKw) CDF function works correctly", {
  # Basic functionality
  q <- 0.5
  alpha <- 2
  beta <- 3
  delta <- 0.5
  lambda <- 1.5

  cdf_kkw <- pkkw(q, alpha, beta, delta, lambda)
  # KKw is equivalent to GKw with gamma=1
  cdf_gkw <- pgkw(q, alpha, beta, gamma = 1, delta, lambda)
  expect_equal(cdf_kkw, cdf_gkw, tolerance = tolerance)

  # Upper tail and log probabilities
  upper_tail_kkw <- pkkw(q, alpha, beta, delta, lambda, lower_tail = FALSE)
  expect_equal(upper_tail_kkw, 1 - cdf_kkw, tolerance = tolerance)

  log_cdf_kkw <- pkkw(q, alpha, beta, delta, lambda, log_p = TRUE)
  expect_equal(log_cdf_kkw, log(cdf_kkw), tolerance = tolerance)
})

test_that("Kumaraswamy-Kumaraswamy (KKw) quantile function works correctly", {
  # Basic functionality
  p <- 0.5
  alpha <- 2
  beta <- 3
  delta <- 0.5
  lambda <- 1.5

  q_kkw <- as.numeric(qkkw(p, alpha, beta, delta, lambda))
  # Verify that quantile and CDF are inverses
  expect_equal(as.numeric(pkkw(q_kkw, alpha, beta, delta, lambda)), p, tolerance = tolerance)
})

test_that("Kumaraswamy-Kumaraswamy (KKw) random generation works correctly", {
  # Basic functionality
  alpha <- 2
  beta <- 3
  delta <- 0.5
  lambda <- 1.5

  samples_kkw <- rkkw(n, alpha, beta, delta, lambda)

  expect_length(samples_kkw, n)
  expect_true(all(samples_kkw >= 0 & samples_kkw <= 1))

  # Test that the generated samples follow the distribution
  ks_test <- suppressWarnings(ks.test(samples_kkw, function(x) pkkw(x, alpha, beta, delta, lambda)))
  expect_true(ks_test$p.value > 0.01)
})
