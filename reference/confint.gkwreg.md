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
#> alpha:(Intercept)  0.61295476  0.7257544
#> alpha:batch1       0.78585088  0.8904090
#> alpha:batch2       0.53536141  0.6471869
#> alpha:batch3       0.65450034  0.7660730
#> alpha:batch4       0.41450496  0.5190080
#> alpha:batch5       0.46463077  0.5758895
#> alpha:batch6       0.39393615  0.5056346
#> alpha:batch7       0.17262525  0.2769125
#> alpha:batch8       0.14778262  0.2585841
#> alpha:batch9       0.07938588  0.2043050
#> alpha:temp                NaN        NaN
#> beta:(Intercept)  28.03403986 29.7293552

# 90 percent confidence intervals
confint(fit, level = 0.90)
#> Warning: some standard errors are NA or infinite; intervals may be unreliable
#>                           5 %       95 %
#> alpha:(Intercept)  0.62202236  0.7166868
#> alpha:batch1       0.79425597  0.8820039
#> alpha:batch2       0.54435070  0.6381976
#> alpha:batch3       0.66346930  0.7571040
#> alpha:batch4       0.42290563  0.5106073
#> alpha:batch5       0.47357450  0.5669458
#> alpha:batch6       0.40291523  0.4966556
#> alpha:batch7       0.18100857  0.2685292
#> alpha:batch8       0.15668960  0.2496771
#> alpha:batch9       0.08942772  0.1942632
#> alpha:temp                NaN        NaN
#> beta:(Intercept)  28.17032079 29.5930743

# Specific parameters
confint(fit, parm = "alpha:(Intercept)")
#> Warning: some standard errors are NA or infinite; intervals may be unreliable
#>                       2.5 %    97.5 %
#> alpha:(Intercept) 0.6129548 0.7257544
confint(fit, parm = 1:3)
#> Warning: some standard errors are NA or infinite; intervals may be unreliable
#>                       2.5 %    97.5 %
#> alpha:(Intercept) 0.6129548 0.7257544
#> alpha:batch1      0.7858509 0.8904090
#> alpha:batch2      0.5353614 0.6471869
# }
```
