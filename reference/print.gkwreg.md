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
#> -0.030691 -0.007432  0.003587  0.008642  0.017537 
#> 
#> Coefficients:
#> alpha:(Intercept)      alpha:batch1      alpha:batch2      alpha:batch3 
#>          0.669333          0.838133          0.591261          0.710279 
#>      alpha:batch4      alpha:batch5      alpha:batch6      alpha:batch7 
#>          0.466760          0.520279          0.449768          0.224772 
#>      alpha:batch8      alpha:batch9        alpha:temp  beta:(Intercept) 
#>          0.203167          0.141873          0.005282         28.881137 
#> 
#> Link functions:  alpha = log, beta = log
#> 
#> Degrees of Freedom: 31 Total (i.e. Null);  20 Residual
#> Residual Deviance: -193.9    AIC: -169.9
#> Log-Lik: 96.97,  BIC: -152.3 
#> R-squared: 0.985,  RMSE: 0.01292 
#> 
#> Number of Fisher Scoring iterations: 77
#> 

# With more digits
print(fit, digits = 5)
#> 
#> Call:  gkwreg(formula = yield ~ batch + temp, data = GasolineYield, 
#>     family = "kw")
#> 
#> Deviance Residuals: 
#>        Min         1Q     Median         3Q        Max 
#> -0.0306910 -0.0074323  0.0035873  0.0086417  0.0175366 
#> 
#> Coefficients:
#> alpha:(Intercept)      alpha:batch1      alpha:batch2      alpha:batch3 
#>         0.6693333         0.8381330         0.5912612         0.7102789 
#>      alpha:batch4      alpha:batch5      alpha:batch6      alpha:batch7 
#>         0.4667599         0.5202794         0.4497678         0.2247716 
#>      alpha:batch8      alpha:batch9        alpha:temp  beta:(Intercept) 
#>         0.2031671         0.1418731         0.0052825        28.8811373 
#> 
#> Link functions:  alpha = log, beta = log
#> 
#> Degrees of Freedom: 31 Total (i.e. Null);  20 Residual
#> Residual Deviance: -193.94   AIC: -169.94
#> Log-Lik: 96.969,  BIC: -152.35 
#> R-squared: 0.98502,  RMSE: 0.012917 
#> 
#> Number of Fisher Scoring iterations: 77
#> 
# }
```
