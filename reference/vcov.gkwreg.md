# Extract Variance-Covariance Matrix from a Generalized Kumaraswamy Regression Model

This function extracts the variance-covariance matrix of the estimated
parameters from a fitted Generalized Kumaraswamy regression model. The
variance-covariance matrix is essential for statistical inference,
including hypothesis testing and confidence interval calculation.

## Usage

``` r
# S3 method for class 'gkwreg'
vcov(object, complete = TRUE, ...)
```

## Arguments

- object:

  An object of class `"gkwreg"`, typically the result of a call to
  [`gkwreg`](https://evandeilton.github.io/gkwreg/reference/gkwreg.md).

- complete:

  Logical indicating whether the complete variance-covariance matrix
  should be returned in case some coefficients were omitted from the
  original fit. Currently ignored for `gkwreg` objects.

- ...:

  Additional arguments (currently not used).

## Value

A square matrix with row and column names corresponding to the
coefficients in the model. If the variance-covariance matrix is not
available (for example, if the model was fitted with `hessian = FALSE`),
the function returns `NULL` with a warning.

## Details

The variance-covariance matrix is estimated based on the observed
information matrix, which is derived from the second derivatives of the
log-likelihood function with respect to the model parameters. For
`gkwreg` objects, this matrix is typically computed using the TMB
(Template Model Builder) automatic differentiation framework during
model fitting.

The diagonal elements of the variance-covariance matrix correspond to
the squared standard errors of the parameter estimates, while the
off-diagonal elements represent the covariances between pairs of
parameters.

## See also

[`gkwreg`](https://evandeilton.github.io/gkwreg/reference/gkwreg.md),
[`confint`](https://rdrr.io/r/stats/confint.html),
[`summary.gkwreg`](https://evandeilton.github.io/gkwreg/reference/summary.gkwreg.md)
