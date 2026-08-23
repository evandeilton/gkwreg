# Bayesian Information Criterion for GKw Regression Models

Calculates the Bayesian Information Criterion (BIC), also known as the
Schwarz Information Criterion (SIC), for fitted Generalized Kumaraswamy
regression models.

## Usage

``` r
# S3 method for class 'gkwreg'
BIC(object, ...)
```

## Arguments

- object:

  An object of class `"gkwreg"`, typically obtained from
  [`gkwreg`](https://evandeilton.github.io/gkwreg/reference/gkwreg.md).

- ...:

  Optionally more fitted model objects.

## Value

If only one object is provided, returns a numeric value with the BIC. If
multiple objects are provided, returns a data frame with columns `df`
and `BIC`, with rows named according to the object names in the call.

## Details

The BIC is computed as: \$\$BIC = -2\ell(\hat{\theta}) + p \cdot
\log(n)\$\$ where \\\ell(\hat{\theta})\\ is the maximized
log-likelihood, \\p\\ is the number of estimated parameters, and \\n\\
is the sample size.

When multiple objects are provided, a data frame comparing all models is
returned. Lower BIC values indicate better models. BIC penalizes model
complexity more heavily than AIC, particularly for large samples, and
tends to favor more parsimonious models.

The BIC can be derived from a Bayesian perspective as an approximation
to the logarithm of the Bayes factor, under certain regularity
conditions and assuming uniform priors.

## References

Schwarz, G. (1978). Estimating the dimension of a model. *The Annals of
Statistics*, **6**(2), 461–464.
[doi:10.1214/aos/1176344136](https://doi.org/10.1214/aos/1176344136)

## See also

[`gkwreg`](https://evandeilton.github.io/gkwreg/reference/gkwreg.md),
[`logLik.gkwreg`](https://evandeilton.github.io/gkwreg/reference/logLik.gkwreg.md),
[`AIC.gkwreg`](https://evandeilton.github.io/gkwreg/reference/AIC.gkwreg.md)

## Author

Lopes, J. E.

## Examples

``` r
# \donttest{
# Load example data
data(GasolineYield)

# Fit competing models
fit1 <- gkwreg(yield ~ batch, data = GasolineYield, family = "kw")
fit2 <- gkwreg(yield ~ batch + temp, data = GasolineYield, family = "kw")
#> Warning: NaNs produced
fit3 <- gkwreg(yield ~ temp, data = GasolineYield, family = "kw")

# Calculate BIC for single model
BIC(fit1)
#> [1] -26.81126

# Compare multiple models (with proper names)
BIC(fit1, fit2, fit3)
#>      df        BIC
#> fit1 11  -26.81126
#> fit2 12 -152.34901
#> fit3  3  -70.50538
# }
```
