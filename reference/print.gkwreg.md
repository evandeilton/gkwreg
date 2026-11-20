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
#> -0.029770 -0.008759  0.003532  0.008115  0.016937 
#> 
#> Coefficients:
#> alpha:(Intercept)      alpha:batch1      alpha:batch2      alpha:batch3 
#>          0.669891          0.838120          0.591222          0.710176 
#>      alpha:batch4      alpha:batch5      alpha:batch6      alpha:batch7 
#>          0.466629          0.520248          0.449796          0.224739 
#>      alpha:batch8      alpha:batch9        alpha:temp  beta:(Intercept) 
#>          0.203245          0.141747          0.005282         28.892380 
#> 
#> Link functions:  alpha = log, beta = log
#> 
#> Degrees of Freedom: 31 Total (i.e. Null);  20 Residual
#> Residual Deviance: -193.9    AIC: -169.9
#> Log-Lik: 96.97,  BIC: -152.4 
#> R-squared: 0.9853,  RMSE: 0.01278 
#> 
#> Number of Fisher Scoring iterations: 72
#> 

# With more digits
print(fit, digits = 5)
#> 
#> Call:  gkwreg(formula = yield ~ batch + temp, data = GasolineYield, 
#>     family = "kw")
#> 
#> Deviance Residuals: 
#>        Min         1Q     Median         3Q        Max 
#> -0.0297697 -0.0087587  0.0035316  0.0081155  0.0169374 
#> 
#> Coefficients:
#> alpha:(Intercept)      alpha:batch1      alpha:batch2      alpha:batch3 
#>          0.669891          0.838120          0.591222          0.710176 
#>      alpha:batch4      alpha:batch5      alpha:batch6      alpha:batch7 
#>          0.466629          0.520248          0.449796          0.224739 
#>      alpha:batch8      alpha:batch9        alpha:temp  beta:(Intercept) 
#>          0.203245          0.141747          0.005282         28.892380 
#> 
#> Link functions:  alpha = log, beta = log
#> 
#> Degrees of Freedom: 31 Total (i.e. Null);  20 Residual
#> Residual Deviance: -193.94   AIC: -169.94
#> Log-Lik: 96.969,  BIC: -152.35 
#> R-squared: 0.98533,  RMSE: 0.012784 
#> 
#> Number of Fisher Scoring iterations: 72
#> 
# }
```
