# Extract Terms from GKw Regression Model

Extracts the terms object from a fitted Generalized Kumaraswamy
regression model.

## Usage

``` r
# S3 method for class 'gkwreg'
terms(x, ...)
```

## Arguments

- x:

  An object of class `"gkwreg"`.

- ...:

  Currently not used.

## Value

A terms object.

## See also

[`gkwreg`](https://evandeilton.github.io/gkwreg/reference/gkwreg.md),
[`formula.gkwreg`](https://evandeilton.github.io/gkwreg/reference/formula.gkwreg.md)

## Author

Lopes, J. E.

## Examples

``` r
# \donttest{
data(GasolineYield)
fit <- gkwreg(yield ~ batch + temp, data = GasolineYield, family = "kw")
#> Warning: NaNs produced
terms(fit)
#> yield ~ batch + temp
#> attr(,"variables")
#> list(yield, batch, temp)
#> attr(,"factors")
#>       batch temp
#> yield     0    0
#> batch     1    0
#> temp      0    1
#> attr(,"term.labels")
#> [1] "batch" "temp" 
#> attr(,"order")
#> [1] 1 1
#> attr(,"intercept")
#> [1] 1
#> attr(,"response")
#> [1] 1
#> attr(,".Environment")
#> <environment: 0x55e56d880640>
#> attr(,"predvars")
#> list(yield, batch, temp)
#> attr(,"dataClasses")
#>     yield     batch      temp 
#> "numeric"  "factor" "numeric" 
# }
```
