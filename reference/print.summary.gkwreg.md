# Print Method for Generalized Kumaraswamy Regression Summaries

Formats and prints the summary output of a fitted Generalized
Kumaraswamy (GKw) regression model (objects of class
`"summary.gkwreg"`).

## Usage

``` r
# S3 method for class 'summary.gkwreg'
print(
  x,
  digits = max(3, getOption("digits") - 3),
  signif.stars = getOption("show.signif.stars"),
  ...
)
```

## Arguments

- x:

  An object of class `"summary.gkwreg"`, typically the result of a call
  to
  [`summary.gkwreg`](https://evandeilton.github.io/gkwreg/reference/summary.gkwreg.md).

- digits:

  Integer, controlling the number of significant digits to print.
  Defaults to `max(3, getOption("digits") - 3)`.

- signif.stars:

  Logical. If `TRUE`, significance stars are printed next to the
  p-values in the coefficient table. Defaults to the value of
  `getOption("show.signif.stars")`.

- ...:

  Additional arguments, currently ignored by this method.

## Value

Invisibly returns the original input object `x`. This allows the output
of [`print()`](https://rdrr.io/r/base/print.html) to be assigned, but it
primarily prints the formatted summary to the console.

## Details

This is the print method for objects created by
[`summary.gkwreg`](https://evandeilton.github.io/gkwreg/reference/summary.gkwreg.md).
It formats the summary information for display in the console. It is
typically invoked automatically when
[`print()`](https://rdrr.io/r/base/print.html) is called on a
`summary.gkwreg` object, or simply by typing the name of the summary
object in an interactive R session.

The output includes:

- Model family and the original function call.

- Summary statistics for residuals.

- A coefficient table with estimates, standard errors, z-values, and
  p-values, optionally marked with significance stars (using
  [`printCoefmat`](https://rdrr.io/r/stats/printCoefmat.html)).

- Confidence intervals for coefficients (if available).

- Link functions used for each parameter.

- Mean values of the fitted distribution parameters.

- Key model fit statistics (LogLik, AIC, BIC, RMSE, R^2, etc.).

- Convergence status and number of iterations.

## See also

[`summary.gkwreg`](https://evandeilton.github.io/gkwreg/reference/summary.gkwreg.md),
[`gkwreg`](https://evandeilton.github.io/gkwreg/reference/gkwreg.md),
[`printCoefmat`](https://rdrr.io/r/stats/printCoefmat.html)

## Author

Lopes, J. E.
