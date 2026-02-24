# Print Method for ANOVA of GKw Models

Print method for analysis of deviance tables produced by
[`anova.gkwreg`](https://evandeilton.github.io/gkwreg/reference/anova.gkwreg.md).

## Usage

``` r
# S3 method for class 'anova.gkwreg'
print(
  x,
  digits = max(getOption("digits") - 2L, 3L),
  signif.stars = getOption("show.signif.stars", TRUE),
  dig.tst = digits,
  ...
)
```

## Arguments

- x:

  An object of class `"anova.gkwreg"` from
  [`anova.gkwreg`](https://evandeilton.github.io/gkwreg/reference/anova.gkwreg.md).

- digits:

  Minimum number of significant digits to print. Default is
  `max(getOption("digits") - 2, 3)`.

- signif.stars:

  Logical; if `TRUE` (default), significance stars are printed alongside
  p-values. Can be controlled globally via
  `options(show.signif.stars = FALSE)`.

- dig.tst:

  Number of digits for test statistics. Default is `digits`.

- ...:

  Additional arguments (currently ignored).

## Value

The object `x`, invisibly.

## See also

[`anova.gkwreg`](https://evandeilton.github.io/gkwreg/reference/anova.gkwreg.md)

## Author

Lopes, J. E.
