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
#> alpha:(Intercept)  0.61309049  0.7255762
#> alpha:batch1       0.78581980  0.8904463
#> alpha:batch2       0.53531608  0.6472063
#> alpha:batch3       0.65446611  0.7660917
#> alpha:batch4       0.41447097  0.5190488
#> alpha:batch5       0.46459924  0.5759595
#> alpha:batch6       0.39388604  0.5056496
#> alpha:batch7       0.17258723  0.2769559
#> alpha:batch8       0.14773407  0.2586001
#> alpha:batch9       0.07935147  0.2043948
#> alpha:temp                NaN        NaN
#> beta:(Intercept)  28.03819235 29.7240822

# 90 percent confidence intervals
confint(fit, level = 0.90)
#> Warning: some standard errors are NA or infinite; intervals may be unreliable
#>                          5 %       95 %
#> alpha:(Intercept)  0.6221328  0.7165338
#> alpha:batch1       0.7942304  0.8820357
#> alpha:batch2       0.5443106  0.6382118
#> alpha:batch3       0.6634393  0.7571185
#> alpha:batch4       0.4228776  0.5106421
#> alpha:batch5       0.4735511  0.5670076
#> alpha:batch6       0.4028704  0.4966653
#> alpha:batch7       0.1809771  0.2685661
#> alpha:batch8       0.1566462  0.2496880
#> alpha:batch9       0.0894033  0.1943430
#> alpha:temp               NaN        NaN
#> beta:(Intercept)  28.1737156 29.5885589

# Specific parameters
confint(fit, parm = "alpha:(Intercept)")
#> Warning: some standard errors are NA or infinite; intervals may be unreliable
#>                       2.5 %    97.5 %
#> alpha:(Intercept) 0.6130905 0.7255762
confint(fit, parm = 1:3)
#> Warning: some standard errors are NA or infinite; intervals may be unreliable
#>                       2.5 %    97.5 %
#> alpha:(Intercept) 0.6130905 0.7255762
#> alpha:batch1      0.7858198 0.8904463
#> alpha:batch2      0.5353161 0.6472063
# }
```
