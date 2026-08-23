# Extract Log-Likelihood from Generalized Kumaraswamy Regression Models

Extracts the log-likelihood value from a fitted Generalized Kumaraswamy
(GKw) regression model object.

## Usage

``` r
# S3 method for class 'gkwreg'
logLik(object, ...)
```

## Arguments

- object:

  An object of class `"gkwreg"`, typically obtained from
  [`gkwreg`](https://evandeilton.github.io/gkwreg/reference/gkwreg.md).

- ...:

  Currently not used.

## Value

An object of class `"logLik"` containing the log-likelihood value with
the following attributes:

- `df`:

  Number of estimated parameters

- `nobs`:

  Number of observations

## Details

The log-likelihood is extracted from the fitted model object and
returned as an object of class `"logLik"` with appropriate attributes
for the number of parameters (`df`) and observations (`nobs`). These
attributes are required for information criteria calculations.

For a GKw regression model with parameter vector \\\theta\\, the
log-likelihood is defined as: \$\$\ell(\theta \mid y) = \sum\_{i=1}^n
\log f(y_i; \alpha_i, \beta_i, \gamma_i, \delta_i, \lambda_i)\$\$ where
\\f(\cdot)\\ is the probability density function of the specified GKw
family distribution, and the parameters may depend on covariates through
link functions.

## See also

[`gkwreg`](https://evandeilton.github.io/gkwreg/reference/gkwreg.md),
[`AIC.gkwreg`](https://evandeilton.github.io/gkwreg/reference/AIC.gkwreg.md),
[`BIC.gkwreg`](https://evandeilton.github.io/gkwreg/reference/BIC.gkwreg.md)

## Author

Lopes, J. E.

## Examples

``` r
# \donttest{
# Load example data
data(GasolineYield)

# Fit a Kumaraswamy regression model
fit <- gkwreg(yield ~ batch + temp, data = GasolineYield, family = "kw")
#> Warning: NaNs produced

# Extract log-likelihood
ll <- logLik(fit)
print(ll)
#> 'log Lik.' 96.96932 (df=12)

# Access attributes
cat("Log-likelihood:", as.numeric(ll), "\n")
#> Log-likelihood: 96.96932 
cat("Parameters:", attr(ll, "df"), "\n")
#> Parameters: 12 
cat("Observations:", attr(ll, "nobs"), "\n")
#> Observations: 32 
# }
```
