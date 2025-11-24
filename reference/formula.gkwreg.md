# Extract Formula from GKw Regression Model

Extracts the model formula from a fitted Generalized Kumaraswamy
regression model object. Properly handles formulas with up to 5 parts.

## Usage

``` r
# S3 method for class 'gkwreg'
formula(x, ...)
```

## Arguments

- x:

  An object of class `"gkwreg"`.

- ...:

  Currently not used.

## Value

The formula used to fit the model. For multi-part formulas, returns an
object of class `"Formula"`.

## See also

[`gkwreg`](https://evandeilton.github.io/gkwreg/reference/gkwreg.md),
[`update.gkwreg`](https://evandeilton.github.io/gkwreg/reference/update.gkwreg.md)

## Author

Lopes, J. E.

## Examples

``` r
# \donttest{
data(GasolineYield)

# Simple formula
fit1 <- gkwreg(yield ~ batch + temp, data = GasolineYield, family = "kw")
#> Warning: NaNs produced
formula(fit1)
#> yield ~ batch + temp
#> <environment: 0x56033592e068>

# Two-part formula
fit2 <- gkwreg(yield ~ temp | batch, data = GasolineYield, family = "kw")
formula(fit2)
#> yield ~ temp | batch
#> <environment: 0x56033592e068>

# Five-part formula
fit3 <- gkwreg(yield ~ temp | batch | temp | 1 | 1,
  data = GasolineYield, family = "gkw"
)
#> using C++ compiler: ‘g++ (Ubuntu 13.3.0-6ubuntu2~24.04) 13.3.0’
formula(fit3)
#> yield ~ temp | batch | temp | 1 | 1
#> <environment: 0x56033592e068>
# }
```
