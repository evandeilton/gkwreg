# Confidence Intervals for Generalized Kumaraswamy Regression Parameters

Computes confidence intervals for model parameters in fitted gkwreg
objects using Wald (normal approximation) method based on asymptotic
theory.

## Usage

``` r
# S3 method for class 'gkwreg'
confint(object, parm, level = 0.95, ...)
```

## Arguments

- object:

  An object of class `"gkwreg"` from
  [`gkwreg`](https://evandeilton.github.io/gkwreg/reference/gkwreg.md).

- parm:

  A specification of which parameters are to be given confidence
  intervals, either a vector of numbers or a vector of names. If
  missing, all parameters are considered.

- level:

  The confidence level required. Default is 0.95.

- ...:

  Additional arguments (currently unused).

## Value

A matrix (or vector) with columns giving lower and upper confidence
limits for each parameter. These will be labeled as (1-level)/2 and 1 -
(1-level)/2 in percent (by default 2.5 percent and 97.5 percent).

## Details

The confidence intervals are computed using the Wald method based on
asymptotic normality of maximum likelihood estimators: \$\$CI =
\hat{\theta} \pm z\_{\alpha/2} \times SE(\hat{\theta})\$\$ where
\\z\_{\alpha/2}\\ is the appropriate normal quantile and
\\SE(\hat{\theta})\\ is the standard error from the Hessian matrix.

The model must have been fitted with `hessian = TRUE` (the default) in
[`gkw_control`](https://evandeilton.github.io/gkwreg/reference/gkw_control.md).
If standard errors are not available, an error is raised.

## See also

[`gkwreg`](https://evandeilton.github.io/gkwreg/reference/gkwreg.md),
[`summary.gkwreg`](https://evandeilton.github.io/gkwreg/reference/summary.gkwreg.md),
[`confint`](https://rdrr.io/r/stats/confint.html)

## Author

Lopes, J. E.

## Examples

``` r
# \donttest{
data(GasolineYield)
fit <- gkwreg(yield ~ batch + temp, data = GasolineYield, family = "kw")
#> Warning: NaNs produced

# 95 percent confidence intervals
confint(fit)
#> Warning: some standard errors are NA or infinite; intervals may be unreliable
#>                         2.5 %     97.5 %
#> alpha:(Intercept)  0.61360892  0.7261721
#> alpha:batch1       0.78583358  0.8904055
#> alpha:batch2       0.53532330  0.6471204
#> alpha:batch3       0.65443679  0.7659146
#> alpha:batch4       0.41441364  0.5188438
#> alpha:batch5       0.46463674  0.5758602
#> alpha:batch6       0.39394438  0.5056469
#> alpha:batch7       0.17260900  0.2768690
#> alpha:batch8       0.14783914  0.2586511
#> alpha:batch9       0.07933947  0.2041538
#> alpha:temp                NaN        NaN
#> beta:(Intercept)  28.04744330 29.7373164

# 90 percent confidence intervals
confint(fit, level = 0.90)
#> Warning: some standard errors are NA or infinite; intervals may be unreliable
#>                           5 %       95 %
#> alpha:(Intercept)  0.62265751  0.7171235
#> alpha:batch1       0.79423978  0.8819993
#> alpha:batch2       0.54431031  0.6381334
#> alpha:batch3       0.66339813  0.7569532
#> alpha:batch4       0.42280844  0.5104490
#> alpha:batch5       0.47357764  0.5669193
#> alpha:batch6       0.40292379  0.4966675
#> alpha:batch7       0.18099013  0.2684879
#> alpha:batch8       0.15674695  0.2497433
#> alpha:batch9       0.08937289  0.1941204
#> alpha:temp                NaN        NaN
#> beta:(Intercept)  28.18328674 29.6014730

# Specific parameters
confint(fit, parm = "alpha:(Intercept)")
#> Warning: some standard errors are NA or infinite; intervals may be unreliable
#>                       2.5 %    97.5 %
#> alpha:(Intercept) 0.6136089 0.7261721
confint(fit, parm = 1:3)
#> Warning: some standard errors are NA or infinite; intervals may be unreliable
#>                       2.5 %    97.5 %
#> alpha:(Intercept) 0.6136089 0.7261721
#> alpha:batch1      0.7858336 0.8904055
#> alpha:batch2      0.5353233 0.6471204
# }
```
