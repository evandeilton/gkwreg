# Extract Model Frame from GKw Regression Model

Extracts the model frame from a fitted Generalized Kumaraswamy
regression model object.

## Usage

``` r
# S3 method for class 'gkwreg'
model.frame(formula, ...)
```

## Arguments

- formula:

  An object of class `"gkwreg"`.

- ...:

  Currently not used.

## Value

A data frame containing the variables used in fitting the model.

## See also

[`gkwreg`](https://evandeilton.github.io/gkwreg/reference/gkwreg.md),
[`model.matrix.gkwreg`](https://evandeilton.github.io/gkwreg/reference/model.matrix.gkwreg.md)

## Author

Lopes, J. E.

## Examples

``` r
# \donttest{
data(GasolineYield)
fit <- gkwreg(yield ~ batch + temp, data = GasolineYield, family = "kw")
#> Warning: NaNs produced
head(model.frame(fit))
#>   yield batch temp
#> 1 0.122     1  205
#> 2 0.223     1  275
#> 3 0.347     1  345
#> 4 0.457     1  407
#> 5 0.080     2  218
#> 6 0.131     2  273
# }
```
