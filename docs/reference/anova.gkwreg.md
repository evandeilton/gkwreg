# Analysis of Deviance for GKw Regression Models

Computes an analysis of deviance table for one or more fitted
Generalized Kumaraswamy (GKw) regression model objects. When multiple
models are provided, likelihood ratio tests are performed to compare
nested models.

## Usage

``` r
# S3 method for class 'gkwreg'
anova(object, ..., test = c("Chisq", "none"))
```

## Arguments

- object:

  An object of class `"gkwreg"`, typically obtained from
  [`gkwreg`](https://evandeilton.github.io/gkwreg/reference/gkwreg.md).

- ...:

  Additional objects of class `"gkwreg"` for model comparison. Models
  must be nested and fitted to the same dataset.

- test:

  A character string specifying the test statistic to use. Currently
  only `"Chisq"` (default) is supported, which performs likelihood ratio
  tests. Can also be `"none"` for no tests.

## Value

An object of class `c("anova.gkwreg", "anova", "data.frame")`, with the
following columns:

- `Resid. Df`:

  Residual degrees of freedom

- `Resid. Dev`:

  Residual deviance (-2 \<U+00D7\> log-likelihood)

- `Df`:

  Change in degrees of freedom (for model comparisons)

- `Deviance`:

  Change in deviance (for model comparisons)

- `Pr(>Chi)`:

  P-value from the chi-squared test (if `test = "Chisq"`)

## Details

When a single model is provided, the function returns a table showing
the residual degrees of freedom and deviance.

When multiple models are provided, the function compares them using
likelihood ratio tests (LRT). Models are automatically ordered by their
complexity (degrees of freedom). The LRT statistic is computed as:
\$\$LRT = 2(\ell_1 - \ell_0)\$\$ where \\\ell_1\\ is the log-likelihood
of the more complex model and \\\ell_0\\ is the log-likelihood of the
simpler (nested) model. Under the null hypothesis that the simpler model
is adequate, the LRT statistic follows a chi-squared distribution with
degrees of freedom equal to the difference in the number of parameters
between the models.

**Important**: This method assumes that the models being compared are
nested (i.e., one model is a special case of the other) and fitted to
the same data. Comparing non-nested models or models fitted to different
datasets will produce unreliable results. Use
[`AIC`](https://rdrr.io/r/stats/AIC.html) or
[`BIC`](https://rdrr.io/r/stats/AIC.html) for comparing non-nested
models.

The deviance is defined as \\-2 \times \text{log-likelihood}\\. For
models fitted by maximum likelihood, smaller (more negative) deviances
indicate better fit. Note that deviance can be negative when the
log-likelihood is positive, which occurs when density values exceed 1
(common in continuous distributions on bounded intervals). What matters
for inference is the *change* in deviance between models, which should
be positive when the more complex model fits better.

## References

Wilks, S. S. (1938). The large-sample distribution of the likelihood
ratio for testing composite hypotheses. *The Annals of Mathematical
Statistics*, **9**(1), 60–62.
[doi:10.1214/aoms/1177732360](https://doi.org/10.1214/aoms/1177732360)

Pawitan, Y. (2001). *In All Likelihood: Statistical Modelling and
Inference Using Likelihood*. Oxford University Press.

## See also

[`gkwreg`](https://evandeilton.github.io/gkwreg/reference/gkwreg.md),
[`logLik.gkwreg`](https://evandeilton.github.io/gkwreg/reference/logLik.gkwreg.md),
[`AIC.gkwreg`](https://evandeilton.github.io/gkwreg/reference/AIC.gkwreg.md),
[`BIC.gkwreg`](https://evandeilton.github.io/gkwreg/reference/BIC.gkwreg.md),
[`lrtest`](https://evandeilton.github.io/gkwreg/reference/lrtest.md)

## Author

Lopes, J. E.

## Examples

``` r
# \donttest{
# Load example data
data(GasolineYield)

# Fit a series of nested models
fit1 <- gkwreg(yield ~ 1, data = GasolineYield, family = "kw")
fit2 <- gkwreg(yield ~ temp, data = GasolineYield, family = "kw")
fit3 <- gkwreg(yield ~ batch + temp, data = GasolineYield, family = "kw")
#> Warning: NaNs produced

# ANOVA table for single model
anova(fit3)
#> Analysis of Deviance Table
#> 
#> Model: gkwreg(formula = yield ~ batch + temp, data = GasolineYield, 
#> Model:     family = "kw")
#> 
#>   Resid. Df Resid. Dev
#> 1  20.00000 -193.93865

# Compare nested models using likelihood ratio tests
anova(fit1, fit2, fit3)
#> Analysis of Deviance Table
#> 
#> Model 1: yield ~ 1
#> Model 2: yield ~ temp
#> Model 3: yield ~ batch + temp
#> 
#>      Resid. Df Resid. Dev Df  Deviance Pr(>Chi)    
#> fit1  30.00000  -57.02258                          
#> fit2  29.00000  -80.90259  1  23.88001  < 1e-04 ***
#> fit3  20.00000 -193.93865  9 113.03606  < 1e-04 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> Model 1 vs 2: Adding temperature is highly significant (p < 0.001)
#> Model 2 vs 3: Adding batch is highly significant (p < 0.001)

# Compare two models
anova(fit2, fit3, test = "Chisq")
#> Analysis of Deviance Table
#> 
#> Model 1: yield ~ temp
#> Model 2: yield ~ batch + temp
#> 
#>      Resid. Df Resid. Dev Df  Deviance Pr(>Chi)    
#> fit2  29.00000  -80.90259                          
#> fit3  20.00000 -193.93865  9 113.03606  < 1e-04 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Suppress test statistics
anova(fit1, fit2, fit3, test = "none")
#> Analysis of Deviance Table
#> 
#> Model 1: yield ~ 1
#> Model 2: yield ~ temp
#> Model 3: yield ~ batch + temp
#> 
#>      Resid. Df Resid. Dev Df  Deviance
#> fit1  30.00000  -57.02258             
#> fit2  29.00000  -80.90259  1  23.88001
#> fit3  20.00000 -193.93865  9 113.03606
# }
```
