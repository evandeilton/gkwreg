# Extract Family from GKw Regression Model

Extracts the family specification from a fitted Generalized Kumaraswamy
regression model object.

## Usage

``` r
# S3 method for class 'gkwreg'
family(object, ...)
```

## Arguments

- object:

  An object of class `"gkwreg"`.

- ...:

  Currently not used.

## Value

A character string indicating the family used in the model.

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
family(fit)
#> [1] "kw"
# }
```
