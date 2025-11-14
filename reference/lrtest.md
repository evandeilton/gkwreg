# Likelihood Ratio Test for Nested GKw Models

Performs a likelihood ratio test to compare two nested Generalized
Kumaraswamy regression models.

## Usage

``` r
lrtest(object, object2)
```

## Arguments

- object:

  A fitted model object of class `"gkwreg"` (the restricted model).

- object2:

  A fitted model object of class `"gkwreg"` (the full model).

## Value

A list with class `"htest"` containing:

- `statistic`:

  The LRT test statistic

- `parameter`:

  Degrees of freedom for the test

- `p.value`:

  P-value from the chi-squared distribution

- `method`:

  Description of the test

- `data.name`:

  Names of the compared models

## Details

This function performs a likelihood ratio test (LRT) to compare two
nested models. The test statistic is: \$\$LRT = 2(\ell\_{\text{full}} -
\ell\_{\text{restricted}})\$\$ which follows a chi-squared distribution
with degrees of freedom equal to the difference in the number of
parameters.

The models must be nested (one is a special case of the other) and
fitted to the same data for the test to be valid.

## See also

[`anova.gkwreg`](https://evandeilton.github.io/gkwreg/reference/anova.gkwreg.md)

## Author

Lopes, J. E.

## Examples

``` r
# \donttest{
data(GasolineYield)

# Fit nested models
fit_restricted <- gkwreg(yield ~ temp, data = GasolineYield, family = "kw")
fit_full <- gkwreg(yield ~ batch + temp, data = GasolineYield, family = "kw")
#> Warning: NaNs produced

# Likelihood ratio test
lrtest(fit_restricted, fit_full)
#> 
#>  Likelihood Ratio Test for Nested GKw Models
#> 
#> data:  fit_restricted vs fit_full
#> LRT = 113.04, df = 9, p-value < 2.2e-16
#> 
# }
```
