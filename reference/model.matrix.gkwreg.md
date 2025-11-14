# Extract Model Matrix from GKw Regression Model

Extracts the model matrix (design matrix) from a fitted Generalized
Kumaraswamy regression model object.

## Usage

``` r
# S3 method for class 'gkwreg'
model.matrix(object, ...)
```

## Arguments

- object:

  An object of class `"gkwreg"`.

- ...:

  Currently not used.

## Value

A design matrix.

## See also

[`gkwreg`](https://evandeilton.github.io/gkwreg/reference/gkwreg.md),
[`model.frame.gkwreg`](https://evandeilton.github.io/gkwreg/reference/model.frame.gkwreg.md)

## Author

Lopes, J. E.

## Examples

``` r
# \donttest{
data(GasolineYield)
fit <- gkwreg(yield ~ batch + temp, data = GasolineYield, family = "kw")
#> Warning: NaNs produced
head(model.matrix(fit))
#>   (Intercept) batch1 batch2 batch3 batch4 batch5 batch6 batch7 batch8 batch9
#> 1           1      1      0      0      0      0      0      0      0      0
#> 2           1      1      0      0      0      0      0      0      0      0
#> 3           1      1      0      0      0      0      0      0      0      0
#> 4           1      1      0      0      0      0      0      0      0      0
#> 5           1      0      1      0      0      0      0      0      0      0
#> 6           1      0      1      0      0      0      0      0      0      0
#>   temp
#> 1  205
#> 2  275
#> 3  345
#> 4  407
#> 5  218
#> 6  273
# }
```
