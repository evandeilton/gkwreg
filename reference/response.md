# Extract Response Variable from GKw Regression Model

Extracts the response variable from a fitted Generalized Kumaraswamy
regression model object.

## Usage

``` r
response(object, ...)

# S3 method for class 'gkwreg'
response(object, ...)
```

## Arguments

- object:

  An object of class `"gkwreg"`.

- ...:

  Currently not used.

## Value

A numeric vector containing the response variable values.

## See also

[`gkwreg`](https://evandeilton.github.io/gkwreg/reference/gkwreg.md),
[`fitted.gkwreg`](https://evandeilton.github.io/gkwreg/reference/fitted.gkwreg.md),
[`residuals.gkwreg`](https://evandeilton.github.io/gkwreg/reference/residuals.gkwreg.md)

## Author

Lopes, J. E.

## Examples

``` r
# \donttest{
data(GasolineYield)
fit <- gkwreg(yield ~ batch + temp, data = GasolineYield, family = "kw")
#> Warning: NaNs produced
y <- response(fit)
head(y)
#>     1     2     3     4     5     6 
#> 0.122 0.223 0.347 0.457 0.080 0.131 
# }
```
