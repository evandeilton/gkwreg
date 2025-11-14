# Summary Method for Generalized Kumaraswamy Regression Models

Computes and returns a detailed statistical summary for a fitted
Generalized Kumaraswamy (GKw) regression model object of class
`"gkwreg"`.

## Usage

``` r
# S3 method for class 'gkwreg'
summary(object, conf.level = 0.95, ...)
```

## Arguments

- object:

  An object of class `"gkwreg"`, typically the result of a call to
  [`gkwreg`](https://evandeilton.github.io/gkwreg/reference/gkwreg.md).

- conf.level:

  Numeric. The desired confidence level for constructing confidence
  intervals for the regression coefficients. Default is 0.95.

- ...:

  Additional arguments, currently ignored by this method.

## Value

An object of class `"summary.gkwreg"`, which is a list containing the
following components:

- call:

  The original function call that created the `object`.

- family:

  Character string specifying the distribution family.

- coefficients:

  A data frame (matrix) containing the coefficient estimates, standard
  errors, z-values, and p-values.

- conf.int:

  A matrix containing the lower and upper bounds of the confidence
  intervals for the coefficients (if standard errors are available).

- link:

  A list of character strings specifying the link functions used.

- fitted_parameters:

  A list containing the mean values of the estimated distribution
  parameters.

- residuals:

  A named numeric vector containing summary statistics for the response
  residuals.

- nobs:

  Number of observations used in the fit.

- npar:

  Total number of estimated regression coefficients.

- df.residual:

  Residual degrees of freedom.

- loglik:

  The maximized log-likelihood value.

- aic:

  Akaike Information Criterion.

- bic:

  Bayesian Information Criterion.

- rmse:

  Root Mean Squared Error of the residuals.

- efron_r2:

  Efron's pseudo-R-squared value.

- mean_absolute_error:

  Mean Absolute Error of the residuals.

- convergence:

  Convergence code from the optimizer.

- iterations:

  Number of iterations reported by the optimizer.

- conf.level:

  The confidence level used for calculating intervals.

## Details

This method provides a comprehensive summary of the fitted `gkwreg`
model. It calculates z-values and p-values for the regression
coefficients based on the estimated standard errors (if available) and
computes confidence intervals at the specified `conf.level`. The summary
includes:

- The model call.

- The distribution family used.

- A table of coefficients including estimates, standard errors,
  z-values, and p-values. Note: Significance stars are typically added
  by the corresponding `print.summary.gkwreg` method.

- Confidence intervals for the coefficients.

- Link functions used for each parameter.

- Mean values of the fitted distribution parameters (\\\alpha, \beta,
  \gamma, \delta, \lambda\\).

- A five-number summary (Min, Q1, Median, Q3, Max) plus the mean of the
  response residuals.

- Key model fit statistics (Log-likelihood, AIC, BIC, RMSE, Efron's
  R^2).

- Information about model convergence and optimizer iterations.

If standard errors were not computed (e.g., `hessian = FALSE` in the
original `gkwreg` call), the coefficient table will only contain
estimates, and confidence intervals will not be available.

## See also

[`gkwreg`](https://evandeilton.github.io/gkwreg/reference/gkwreg.md),
[`print.summary.gkwreg`](https://evandeilton.github.io/gkwreg/reference/print.summary.gkwreg.md),
[`coef`](https://rdrr.io/r/stats/coef.html),
[`confint`](https://rdrr.io/r/stats/confint.html)

## Author

Lopes, J. E.

## Examples

``` r
# \donttest{
set.seed(123)
n <- 100
x1 <- runif(n, -2, 2)
x2 <- rnorm(n)
alpha_coef <- c(0.8, 0.3, -0.2)
beta_coef <- c(1.2, -0.4, 0.1)
eta_alpha <- alpha_coef[1] + alpha_coef[2] * x1 + alpha_coef[3] * x2
eta_beta <- beta_coef[1] + beta_coef[2] * x1 + beta_coef[3] * x2
alpha_true <- exp(eta_alpha)
beta_true <- exp(eta_beta)
# Use stats::rbeta as a placeholder if rkw is unavailable
y <- stats::rbeta(n, shape1 = alpha_true, shape2 = beta_true)
y <- pmax(pmin(y, 1 - 1e-7), 1e-7)
df <- data.frame(y = y, x1 = x1, x2 = x2)

# Fit a Kumaraswamy regression model
kw_reg <- gkwreg(y ~ x1 + x2 | x1 + x2, data = df, family = "kw")

# Generate detailed summary using the summary method
summary_kw <- summary(kw_reg)

# Print the summary object (uses print.summary.gkwreg)
print(summary_kw)
#> 
#> Generalized Kumaraswamy Regression Model Summary
#> 
#> Family: kw 
#> 
#> Call:
#> gkwreg(formula = y ~ x1 + x2 | x1 + x2, data = df, family = "kw")
#> 
#> Residuals:
#>     Min  Q1.25%  Median    Mean  Q3.75%     Max 
#> -0.4910 -0.0942  0.0028  0.0042  0.1329  0.5094 
#> 
#> Coefficients:
#>                   Estimate Std. Error z value Pr(>|z|)    
#> alpha:(Intercept)  0.73589    0.10135   7.261 3.85e-13 ***
#> alpha:x1           0.25878    0.08860   2.921  0.00349 ** 
#> alpha:x2          -0.12441    0.08825  -1.410  0.15863    
#> beta:(Intercept)   1.45325    0.17992   8.077 6.64e-16 ***
#> beta:x1           -0.52161    0.17460  -2.987  0.00281 ** 
#> beta:x2            0.03766    0.18694   0.201  0.84034    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Confidence intervals (95%):
#>                        3%     98%
#> alpha:(Intercept)  0.5372  0.9345
#> alpha:x1           0.0851  0.4324
#> alpha:x2          -0.2974  0.0486
#> beta:(Intercept)   1.1006  1.8059
#> beta:x1           -0.8638 -0.1794
#> beta:x2           -0.3287  0.4041
#> 
#> Link functions:
#> alpha: log
#> beta: log
#> 
#> Fitted parameter means:
#> alpha: 2.207
#> beta: 5.072
#> gamma: 1
#> delta: 0
#> lambda: 1
#> 
#> Model fit statistics:
#> Number of observations: 100 
#> Number of parameters: 6 
#> Residual degrees of freedom: 94 
#> Log-likelihood: 46.13 
#> AIC: -80.27 
#> BIC: -64.63 
#> RMSE: 0.1728 
#> Efron's R2: 0.5496 
#> Mean Absolute Error: 0.1365 
#> 
#> Convergence status: Successful 
#> Iterations: 27 
#> 

# Extract coefficient table directly from the summary object
coef_table <- coef(summary_kw) # Equivalent to summary_kw$coefficients
print(coef_table)
#>                      Estimate Std. Error    z value     Pr(>|z|)
#> alpha:(Intercept)  0.73589120 0.10135318  7.2606619 3.852007e-13
#> alpha:x1           0.25877630 0.08860253  2.9206422 3.493107e-03
#> alpha:x2          -0.12440703 0.08825176 -1.4096833 1.586332e-01
#> beta:(Intercept)   1.45324717 0.17992343  8.0770311 6.636262e-16
#> beta:x1           -0.52161306 0.17460352 -2.9874143 2.813482e-03
#> beta:x2            0.03766118 0.18693698  0.2014646 8.403353e-01
# }
```
