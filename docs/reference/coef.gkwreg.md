# Extract Coefficients from a Fitted GKw Regression Model

Extracts the estimated regression coefficients from a fitted Generalized
Kumaraswamy (GKw) regression model object of class `"gkwreg"`. This is
an S3 method for the generic [`coef`](https://rdrr.io/r/stats/coef.html)
function.

## Usage

``` r
# S3 method for class 'gkwreg'
coef(object, ...)
```

## Arguments

- object:

  An object of class `"gkwreg"`, typically the result of a call to
  [`gkwreg`](https://evandeilton.github.io/gkwreg/reference/gkwreg.md).

- ...:

  Additional arguments, currently ignored by this method.

## Value

A named numeric vector containing the estimated regression coefficients
for all modeled parameters. The names indicate the parameter (e.g.,
`alpha`, `beta`) and the corresponding predictor variable (e.g.,
`(Intercept)`, `x1`).

## Details

This function provides the standard way to access the estimated
regression coefficients from a model fitted with
[`gkwreg`](https://evandeilton.github.io/gkwreg/reference/gkwreg.md). It
simply extracts the `coefficients` component from the fitted model
object. The function [`coefficients`](https://rdrr.io/r/stats/coef.html)
is an alias for this function.

## See also

[`gkwreg`](https://evandeilton.github.io/gkwreg/reference/gkwreg.md),
[`summary.gkwreg`](https://evandeilton.github.io/gkwreg/reference/summary.gkwreg.md),
[`coef`](https://rdrr.io/r/stats/coef.html),
[`confint`](https://rdrr.io/r/stats/confint.html)

## Author

Lopes, J. E.
