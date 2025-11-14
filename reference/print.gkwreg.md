# Print Method for Generalized Kumaraswamy Regression Models

Print method for objects of class `"gkwreg"`. Provides a concise summary
of the fitted model following the style of
[`print.lm`](https://rdrr.io/r/stats/lm.html).

## Usage

``` r
# S3 method for class 'gkwreg'
print(x, digits = max(3, getOption("digits") - 3), ...)
```

## Arguments

- x:

  An object of class `"gkwreg"`, typically obtained from
  [`gkwreg`](https://evandeilton.github.io/gkwreg/reference/gkwreg.md).

- digits:

  Minimum number of significant digits to print. Default is
  `max(3, getOption("digits") - 3)`.

- ...:

  Additional arguments passed to or from other methods.

## Value

The object `x`, invisibly.

## Details

The print method provides a concise overview of the fitted model,
showing: the call, deviance residuals summary, coefficient estimates,
link functions, and basic fit statistics. For more detailed output
including standard errors and significance tests, use
[`summary.gkwreg`](https://evandeilton.github.io/gkwreg/reference/summary.gkwreg.md).

## See also

[`gkwreg`](https://evandeilton.github.io/gkwreg/reference/gkwreg.md),
[`summary.gkwreg`](https://evandeilton.github.io/gkwreg/reference/summary.gkwreg.md)

## Author

Lopes, J. E.

## Examples

``` r
# \donttest{
data(GasolineYield)
fit <- gkwreg(yield ~ batch + temp, data = GasolineYield, family = "kw")
#> Warning: NaNs produced
print(fit)
#> 
#> Call:  gkwreg(formula = yield ~ batch + temp, data = GasolineYield, 
#>     family = "kw")
#> 
#> Deviance Residuals: 
#>       Min        1Q    Median        3Q       Max 
#> -0.028057 -0.002998  0.003669  0.009988  0.017351 
#> 
#> Coefficients:
#> alpha:(Intercept)      alpha:batch1      alpha:batch2      alpha:batch3 
#>          0.669355          0.838130          0.591274          0.710287 
#>      alpha:batch4      alpha:batch5      alpha:batch6      alpha:batch7 
#>          0.466756          0.520260          0.449785          0.224769 
#>      alpha:batch8      alpha:batch9        alpha:temp  beta:(Intercept) 
#>          0.203183          0.141845          0.005282         28.881698 
#> 
#> Link functions:  alpha = log, beta = log
#> 
#> Degrees of Freedom: 31 Total (i.e. Null);  20 Residual
#> Residual Deviance: -193.9    AIC: -169.9
#> Log-Lik: 96.97,  BIC: -152.3 
#> R-squared: 0.9856,  RMSE: 0.01266 
#> 
#> Number of Fisher Scoring iterations: 68
#> 

# With more digits
print(fit, digits = 5)
#> 
#> Call:  gkwreg(formula = yield ~ batch + temp, data = GasolineYield, 
#>     family = "kw")
#> 
#> Deviance Residuals: 
#>        Min         1Q     Median         3Q        Max 
#> -0.0280573 -0.0029985  0.0036691  0.0099881  0.0173508 
#> 
#> Coefficients:
#> alpha:(Intercept)      alpha:batch1      alpha:batch2      alpha:batch3 
#>         0.6693546         0.8381300         0.5912742         0.7102867 
#>      alpha:batch4      alpha:batch5      alpha:batch6      alpha:batch7 
#>         0.4667565         0.5202601         0.4497854         0.2247689 
#>      alpha:batch8      alpha:batch9        alpha:temp  beta:(Intercept) 
#>         0.2031834         0.1418454         0.0052823        28.8816975 
#> 
#> Link functions:  alpha = log, beta = log
#> 
#> Degrees of Freedom: 31 Total (i.e. Null);  20 Residual
#> Residual Deviance: -193.94   AIC: -169.94
#> Log-Lik: 96.969,  BIC: -152.35 
#> R-squared: 0.98561,  RMSE: 0.01266 
#> 
#> Number of Fisher Scoring iterations: 68
#> 
# }
```
