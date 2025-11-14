# Get Call from GKw Regression Model

Extracts the call that was used to fit a Generalized Kumaraswamy
regression model.

## Usage

``` r
# S3 method for class 'gkwreg'
getCall(x, ...)
```

## Arguments

- x:

  An object of class `"gkwreg"`.

- ...:

  Currently not used.

## Value

The matched call.

## See also

[`gkwreg`](https://evandeilton.github.io/gkwreg/reference/gkwreg.md),
[`update.gkwreg`](https://evandeilton.github.io/gkwreg/reference/update.gkwreg.md)

## Author

Lopes, J. E.

## Examples

``` r
# \donttest{
data(GasolineYield)
fit <- gkwreg(yield ~ batch + temp, data = GasolineYield, family = "kw")
#> Warning: NaNs produced
getCall(fit)
#> gkwreg(formula = yield ~ batch + temp, data = GasolineYield, 
#>     family = "kw")
# }
```
