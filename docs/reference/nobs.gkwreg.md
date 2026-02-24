# Number of Observations for GKw Regression Models

Extracts the number of observations from a fitted Generalized
Kumaraswamy regression model.

## Usage

``` r
# S3 method for class 'gkwreg'
nobs(object, ...)
```

## Arguments

- object:

  An object of class `"gkwreg"`, typically obtained from
  [`gkwreg`](https://evandeilton.github.io/gkwreg/reference/gkwreg.md).

- ...:

  Currently not used.

## Value

Integer representing the number of observations used in model fitting.

## See also

[`gkwreg`](https://evandeilton.github.io/gkwreg/reference/gkwreg.md)

## Author

Lopes, J. E.

## Examples

``` r
# \donttest{
data(GasolineYield)
fit <- gkwreg(yield ~ batch + temp, data = GasolineYield, family = "kw")
#> Warning: NaNs produced
nobs(fit)
#> [1] 32
# }
```
