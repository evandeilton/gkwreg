
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gkwreg: Generalized Kumaraswamy Regression Models

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/gkwreg)](https://cran.r-project.org/package=gkwreg)
<!-- [![R-CMD-check](https://github.com/evandeilton/gkwreg/workflows/R-CMD-check/badge.svg)](https://github.com/evandeilton/gkwreg/actions) -->
[![R-CMD-check](https://github.com/evandeilton/gkwreg/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/evandeilton/gkwreg/actions/workflows/R-CMD-check.yaml)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Downloads](https://cranlogs.r-pkg.org/badges/gkwreg)](https://cran.r-project.org/package=gkwreg)

## Overview

`gkwreg` is a comprehensive R package for statistical modeling of
bounded continuous data on the unit interval (0,1). The package
implements both the parameter estimation and regression modeling for the
five-parameter Generalized Kumaraswamy (GKw) distribution family using
high-performance C++ code via RcppArmadillo and Template Model Builder
(TMB).

The package offers two main components:

1.  **Distribution fitting** (`gkwfit`): Efficient maximum likelihood
    estimation for all seven distribution families in the GKw hierarchy,
    including special cases like Kumaraswamy and Beta distributions.

2.  **Regression modeling** (`gkwreg`): A powerful regression framework
    where all five distribution parameters (α, β, γ, δ, λ) can be
    independently modeled as functions of covariates through flexible
    link functions.

This implementation significantly extends existing limited-parameter
models for bounded data, providing more flexibility for modeling complex
bounded response variables in fields like environmental science,
economics, biology, and social sciences.

## Mathematical Foundation

### The Generalized Kumaraswamy Distribution

The GKw distribution is a highly flexible five-parameter family
applicable to bounded continuous random variables. For $y \in (0,1)$,
the probability density function is:

$$f(y; \alpha, \beta, \gamma, \delta, \lambda) = \frac{\lambda \alpha \beta y^{\alpha-1} (1-y^\alpha)^{\beta-1}}{B(\gamma, \delta+1)} \left[1-(1-y^\alpha)^\beta\right]^{\gamma\lambda-1} \left\{1-\left[1-(1-y^\alpha)^\beta\right]^\lambda\right\}^{\delta}$$

where $\alpha, \beta, \gamma, \delta, \lambda > 0$ and
$B(\gamma, \delta+1)$ is the complete Beta function.

The associated cumulative distribution function is:

$$F(y; \alpha, \beta, \gamma, \delta, \lambda) = I_{[1-(1-y^\alpha)^\beta]^\lambda}(\gamma, \delta+1)$$

where $I_z(a,b)$ is the regularized incomplete Beta function.

### Distribution Hierarchy

The GKw family includes seven nested distributions:

| Distribution | Parameters    | Constraint          | Flexibility |
|--------------|---------------|---------------------|-------------|
| GKw          | α, β, γ, δ, λ | None                | Highest     |
| BKw          | α, β, γ, δ    | λ = 1               | High        |
| KKw          | α, β, δ, λ    | γ = 1               | High        |
| EKw          | α, β, λ       | γ = 1, δ = 0        | Moderate    |
| Mc           | γ, δ, λ       | α = 1, β = 1        | Moderate    |
| Kw           | α, β          | γ = 1, δ = 0, λ = 1 | Basic       |
| Beta         | γ, δ          | α = 1, β = 1, λ = 1 | Basic       |

This hierarchy allows model simplification through likelihood ratio
tests, facilitating parsimony while maintaining adequate fit.

### Regression Framework

The regression model extends the GKw distribution by modeling each
parameter as a function of covariates:

$$\theta_i = g_\theta^{-1}(\mathbf{x}_{\theta i}^T \boldsymbol{\beta}_\theta), \quad \theta \in \{\alpha, \beta, \gamma, \delta, \lambda\}$$

where: - $\mathbf{x}_{\theta i}$ is the covariate vector for observation
$i$ and parameter $\theta$ - $\boldsymbol{\beta}_\theta$ is the
corresponding coefficient vector - $g_\theta^{-1}$ is the inverse link
function for parameter $\theta$

The package supports nine inverse link functions:

1.  **Log**: $g^{-1}(\eta) = \exp(\eta)$
2.  **Logit**:
    $g^{-1}(\eta) = \text{scale} \cdot \frac{\exp(\eta)}{1+\exp(\eta)}$
3.  **Probit**: $g^{-1}(\eta) = \text{scale} \cdot \Phi(\eta)$
4.  **Cauchy**:
    $g^{-1}(\eta) = \text{scale} \cdot \Bigl(\frac{1}{\pi}\arctan(\eta) + \frac{1}{2}\Bigr)$
5.  **Cloglog**:
    $g^{-1}(\eta) = \text{scale} \cdot \bigl(1 - \exp(-\exp(\eta))\bigr)$
6.  **Identity**: $g^{-1}(\eta) = \eta$
7.  **Square Root**: $g^{-1}(\eta) = \eta^2$
8.  **Inverse**: $g^{-1}(\eta) = \frac{1}{\eta}$
9.  **Inverse-Square**: $g^{-1}(\eta) = \frac{1}{\sqrt{\eta}}$

The `scale` parameter allows rescaling bounded link functions (logit,
probit, cauchy, cloglog) to appropriate parameter ranges.

## Installation

``` r
# Install from CRAN
install.packages("gkwreg")

# Or development version from GitHub
# install.packages("devtools")
devtools::install_github("evandeilton/gkwreg")
```

### System Requirements

- R (≥ 4.0.0)
- TMB (≥ 1.8.0)
- C++ compiler with C++11 support
- RcppArmadillo (≥ 0.11.0)

## Core Functionality

### Distribution Fitting with gkwfit()

The `gkwfit()` function fits any of the seven distributions in the GKw
family using maximum likelihood estimation:

``` r
# Generate sample data from a Kumaraswamy distribution
set.seed(123)
data <- rbeta(500, 2, 3)^(1/2) # Simple transform to approximate Kw(2,3)

# Fit the full 5-parameter GKw model
fit_gkw <- gkwfit(data, family = "gkw", fit = "tmb", hessian = TRUE)
summary(fit_gkw)

# Fit a Kumaraswamy model (2 parameters)
fit_kw <- gkwfit(data, family = "kw", fit = "tmb")
summary(fit_kw)

# Compare models with likelihood ratio test
# This is automatic if submodels = TRUE in the full model fit
fit_full <- gkwfit(data, family = "gkw", submodels = TRUE)
fit_full$lrt$kw  # Shows LRT for GKw vs Kw

# Plot diagnostics
plot(fit_kw)
```

Key features of `gkwfit()`:

- Choice of estimation methods: TMB (default) or Newton-Raphson
- Profile likelihood computation for uncertainty quantification
- Comprehensive diagnostics including goodness-of-fit tests
- Automatic fitting of submodels for nested model comparison

### Regression Modeling with gkwreg()

The `gkwreg()` function implements the full GKw regression framework:

``` r
# Simulate regression data
set.seed(456)
n <- 1000
x1 <- rnorm(n)
x2 <- rbinom(n, 1, 0.5)
x3 <- runif(n, -1, 1)

# True parameters as functions of covariates
alpha <- exp(0.5 + 0.3 * x1)
beta <- exp(0.8 - 0.4 * x2)
gamma <- exp(0)  # Constant
delta <- exp(-0.2 - 0.5 * x3)
lambda <- exp(0.1)  # Constant

# Generate response using probability integral transform
u <- runif(n)
y <- qgkw(u, alpha, beta, gamma, delta, lambda)  # Requires gkwreg package functions

# Create data frame
dat <- data.frame(y = y, x1 = x1, x2 = x2, x3 = x3)

# Fit GKw regression model
# Formula notation: y ~ alpha_covariates | beta_covariates | gamma_covariates | delta_covariates | lambda_covariates
model <- gkwreg(
  formula = y ~ x1 | x2 | 1 | x3 | 1,
  data = dat,
  link = c(1, 1, 1, 1, 1),  # All log links
  scale = c(10, 10, 10, 10, 10)
)

# Model summary
summary(model)

# Diagnostic plots
plot(model, which = 1:6, use_ggplot = TRUE, arrange_plots = TRUE)

# Predictions
new_data <- data.frame(x1 = c(0, 1), x2 = c(0, 1), x3 = c(-0.5, 0.5))
predict(model, newdata = new_data, type = "response")
predict(model, newdata = new_data, type = "parameter")  # Get all parameters
```

Key features of `gkwreg()`:

- Formula interface for specifying different covariates for each
  parameter
- Nine link functions for flexible parameter modeling
- High-performance TMB implementation with automatic differentiation
- Comprehensive diagnostics with multiple residual types
- Extensive prediction capabilities

## Advanced Usage

### Custom Starting Values and Optimization Control

For complex models or challenging datasets, customizing the starting
values and optimization parameters can improve convergence:

``` r
# Custom starting values for gkwreg()
model_custom <- gkwreg(
  formula = y ~ x1 | x2 | 1 | x3 | 1,
  data = dat,
  link = c(1, 1, 1, 1, 1),
  start = list(
    beta1 = c(0.4, 0.3),  # Intercept and coefficient for x1
    beta2 = c(0.8, -0.4), # Intercept and coefficient for x2
    beta3 = c(0),         # Intercept only for gamma
    beta4 = c(-0.2, -0.5),# Intercept and coefficient for x3
    beta5 = c(0.1)        # Intercept only for lambda
  ),
  control = list(
    optimizer = "nlminb",
    eval.max = 10000,
    iter.max = 10000,
    rel.tol = 1e-12
  )
)

# Custom starting values for gkwfit()
fit_custom <- gkwfit(
  data = data,
  family = "gkw",
  start = list(alpha = 2, beta = 3, gamma = 1, delta = 0.5, lambda = 1),
  optimizer.control = list(tol = 1e-8, max_iter = 200, step_size = 0.5)
)
```

### Model Comparison and Selection

The package provides several tools for model comparison:

``` r
# Fit nested models
model_full <- gkwreg(y ~ x1 + x2 | x1 + x2 | x1 + x2 | x1 + x2 | x1 + x2, data = dat)
model_reduced <- gkwreg(y ~ x1 + x2 | x1 | 1 | x1 | 1, data = dat)

# Compare with AIC/BIC
AIC(model_full, model_reduced)
BIC(model_full, model_reduced)

# Likelihood ratio test (requires manual calculation for gkwreg)
lr_stat <- 2 * (logLik(model_full) - logLik(model_reduced))
p_value <- 1 - pchisq(lr_stat, df = model_full$npar - model_reduced$npar)
cat("LR statistic:", lr_stat, "p-value:", p_value, "\n")
```

### Visualization and Diagnostics

The package provides extensive visualization capabilities:

``` r
# Basic diagnostic plots for distribution fit
plot(fit_gkw, which = 1:4)

# Advanced regression diagnostics
plot(model, which = 1:6, type = "quantile", use_ggplot = TRUE, arrange_plots = TRUE)

# Custom residual plots
pearson_resid <- residuals(model, type = "pearson")
plot(model$fitted.values, pearson_resid,
     xlab = "Fitted Values", ylab = "Pearson Residuals",
     main = "Pearson Residuals vs Fitted Values")
abline(h = 0, lty = 2, col = "red")
```

## Comparison with Other Packages

The `gkwreg` package extends several existing approaches for modeling
bounded data:

| Package | Distribution | Parameters | Regression | Link Functions |
|----|----|----|----|----|
| gkwreg | GKw family (7 models) | 2-5 | All parameters | 9 types |
| betareg | Beta | 2 | Mean and precision | Several |
| gamlss | Many | Varies | All parameters | Several |
| zoib | Beta with inflation | 2 + inflation | Mean and precision | Limited |
| brms | Beta | 2 | Mean and precision | Several |

Key advantages of `gkwreg`: - More flexible shape modeling through 5
parameters - Unified framework for 7 distribution families - Enhanced
numerical stability through TMB implementation - Extended link function
repertoire - Comprehensive diagnostic toolkit

## Implementation Details

The `gkwreg` package is implemented using a combination of R, Rcpp,
RcppArmadillo, and TMB:

1.  **High-Performance Computing**: TMB provides automatic
    differentiation for efficient gradient and Hessian computation,
    accelerating convergence and improving numerical stability.

2.  **Memory Efficiency**: The implementation uses sparse matrix
    operations where appropriate to handle large datasets with many
    covariates.

3.  **Numerical Safeguards**: Careful handling of boundary cases
    prevents common numerical issues in bounded distributions.

4.  **Object-Oriented Design**: Consistent S3 class implementation
    ensures compatibility with R’s modeling ecosystem.

The log-likelihood function used for maximum likelihood estimation is:

$$\ell(\boldsymbol{\Theta}) = \sum_{i=1}^{n} \Big[ \log(\lambda_i) + \log(\alpha_i) + \log(\beta_i) - \log\{B(\gamma_i,\delta_i+1)\} + (\alpha_i-1)\log(y_i) + (\beta_i-1)\log(1-y_i^{\alpha_i}) + (\gamma_i\lambda_i-1)\log\{1-(1-y_i^{\alpha_i})^{\beta_i}\} + \delta_i\log\{1-[1-(1-y_i^{\alpha_i})^{\beta_i}]^{\lambda_i}\} \Big]$$

## References

Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized
distributions. *Journal of Statistical Computation and Simulation*,
81(7), 883-898.

Kumaraswamy, P. (1980). A generalized probability density function for
double-bounded random processes. *Journal of Hydrology*, 46(1-2), 79-88.

Jones, M. C. (2009). Kumaraswamy’s distribution: A beta-type
distribution with some tractability advantages. *Statistical
Methodology*, 6(1), 70-81.

Mitnik, P. A., & Baek, S. (2013). The Kumaraswamy distribution:
median-dispersion re-parameterizations for regression modeling and
simulation-based estimation. *Statistical Papers*, 54(1), 177-192.

Kristensen, K., Nielsen, A., Berg, C. W., Skaug, H., & Bell, B. M.
(2016). TMB: Automatic differentiation and Laplace approximation.
*Journal of Statistical Software*, 70(5), 1-21.

## License

This package is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License, version 3.
