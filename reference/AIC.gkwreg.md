# Akaike Information Criterion for GKw Regression Models

Calculates the Akaike Information Criterion (AIC) for fitted Generalized
Kumaraswamy regression models.

## Usage

``` r
# S3 method for class 'gkwreg'
AIC(object, ..., k = 2)
```

## Arguments

- object:

  An object of class `"gkwreg"`, typically obtained from
  [`gkwreg`](https://evandeilton.github.io/gkwreg/reference/gkwreg.md).

- ...:

  Optionally more fitted model objects.

- k:

  Numeric, the penalty per parameter. Default is `k = 2` for classical
  AIC. Setting `k = log(n)` gives BIC-equivalent penalty.

## Value

If only one object is provided, returns a numeric value with the AIC. If
multiple objects are provided, returns a data frame with columns `df`
and `AIC`, with rows named according to the object names in the call.

## Details

The AIC is computed as: \$\$AIC = -2\ell(\hat{\theta}) + k \cdot p\$\$
where \\\ell(\hat{\theta})\\ is the maximized log-likelihood and \\p\\
is the number of estimated parameters.

When multiple objects are provided, a data frame comparing all models is
returned. Lower AIC values indicate better models, balancing
goodness-of-fit against model complexity.

For small sample sizes, consider the corrected AIC (AICc): \$\$AICc =
AIC + \frac{2p(p+1)}{n-p-1}\$\$ where \\n\\ is the sample size. This
correction is not automatically applied but can be calculated manually.

## References

Akaike, H. (1974). A new look at the statistical model identification.
*IEEE Transactions on Automatic Control*, **19**(6), 716–723.
[doi:10.1109/TAC.1974.1100705](https://doi.org/10.1109/TAC.1974.1100705)

Burnham, K. P., & Anderson, D. R. (2004). Multimodel inference:
Understanding AIC and BIC in model selection. *Sociological Methods &
Research*, **33**(2), 261–304.
[doi:10.1177/0049124104268644](https://doi.org/10.1177/0049124104268644)

## See also

[`gkwreg`](https://evandeilton.github.io/gkwreg/reference/gkwreg.md),
[`logLik.gkwreg`](https://evandeilton.github.io/gkwreg/reference/logLik.gkwreg.md),
[`BIC.gkwreg`](https://evandeilton.github.io/gkwreg/reference/BIC.gkwreg.md)

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

# Calculate AIC for single model
AIC(fit1)
#> [1] -42.93436

# Compare multiple models (with proper names)
AIC(fit1, fit2, fit3)
#>      df        AIC
#> fit1 11  -42.93436
#> fit2 12 -169.93784
#> fit3  3  -74.90259

# Use different penalty
AIC(fit1, k = 4)
#> [1] -20.93436
# }
```
