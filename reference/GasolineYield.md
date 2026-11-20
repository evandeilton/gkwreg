# Gasoline Yield from Crude Oil

Operational data on the proportion of crude oil converted to gasoline
after distillation and fractionation processes.

## Usage

``` r
GasolineYield
```

## Format

A data frame with 32 observations on 6 variables:

- yield:

  numeric. Proportion of crude oil converted to gasoline after
  distillation and fractionation (response variable).

- gravity:

  numeric. Crude oil gravity in degrees API (American Petroleum
  Institute scale).

- pressure:

  numeric. Vapor pressure of crude oil in pounds per square inch (psi).

- temp10:

  numeric. Temperature in degrees Fahrenheit at which 10\\ crude oil has
  vaporized.

- temp:

  numeric. Temperature in degrees Fahrenheit at which all gasoline has
  vaporized (end point).

- batch:

  factor. Batch indicator distinguishing the 10 different crude oils
  used in the experiment.

## Source

Taken from Prater (1956).

## Details

This dataset was collected by Prater (1956) to study gasoline yield from
crude oil. The dependent variable is the proportion of crude oil after
distillation and fractionation. Atkinson (1985) analyzed this dataset
using linear regression and noted that there is "indication that the
error distribution is not quite symmetrical, giving rise to some unduly
large and small residuals".

The dataset contains 32 observations. It has been noted (Daniel and
Wood, 1971, Chapter 8) that there are only ten sets of values of the
first three explanatory variables which correspond to ten different
crudes subjected to experimentally controlled distillation conditions.
These conditions are captured in variable `batch` and the data were
ordered according to the ascending order of `temp10`.

## References

Atkinson, A.C. (1985). *Plots, Transformations and Regression: An
Introduction to Graphical Methods of Diagnostic Regression Analysis*.
New York: Oxford University Press.

Cribari-Neto, F., and Zeileis, A. (2010). Beta Regression in R. *Journal
of Statistical Software*, **34**(2), 1–24.
[doi:10.18637/jss.v034.i02](https://doi.org/10.18637/jss.v034.i02)

Daniel, C., and Wood, F.S. (1971). *Fitting Equations to Data*. New
York: John Wiley and Sons.

Ferrari, S.L.P., and Cribari-Neto, F. (2004). Beta Regression for
Modeling Rates and Proportions. *Journal of Applied Statistics*,
**31**(7), 799–815.

Prater, N.H. (1956). Estimate Gasoline Yields from Crudes. *Petroleum
Refiner*, **35**(5), 236–238.

## Examples

``` r
# \donttest{
require(gkwreg)
require(gkwdist)

data(GasolineYield)

# Example 1: Kumaraswamy regression with batch effects
# Model mean yield as function of batch and temperature
# Allow precision to vary with temperature (heteroscedasticity)
fit_kw <- gkwreg(yield ~ batch + temp | temp,
  data = GasolineYield,
  family = "kw"
)
summary(fit_kw)
#> 
#> Generalized Kumaraswamy Regression Model Summary
#> 
#> Family: kw 
#> 
#> Call:
#> gkwreg(formula = yield ~ batch + temp | temp, data = GasolineYield, 
#>     family = "kw")
#> 
#> Residuals:
#>     Min  Q1.25%  Median    Mean  Q3.75%     Max 
#> -0.0456 -0.0070  0.0019  0.0025  0.0074  0.0707 
#> 
#> Coefficients:
#>                    Estimate Std. Error   z value Pr(>|z|)    
#> alpha:(Intercept) 2.517e-01  2.900e-02 8.679e+00  < 2e-16 ***
#> alpha:batch1      8.280e-01  2.912e-02 2.843e+01  < 2e-16 ***
#> alpha:batch2      5.846e-01  3.125e-02 1.870e+01  < 2e-16 ***
#> alpha:batch3      7.122e-01  2.893e-02 2.462e+01  < 2e-16 ***
#> alpha:batch4      4.648e-01  2.606e-02 1.783e+01  < 2e-16 ***
#> alpha:batch5      5.234e-01  2.626e-02 1.993e+01  < 2e-16 ***
#> alpha:batch6      4.497e-01  2.722e-02 1.652e+01  < 2e-16 ***
#> alpha:batch7      2.251e-01  2.497e-02 9.017e+00  < 2e-16 ***
#> alpha:batch8      2.028e-01  2.556e-02 7.932e+00 2.15e-15 ***
#> alpha:batch9      1.426e-01  2.991e-02 4.768e+00 1.86e-06 ***
#> alpha:temp        6.565e-03  2.584e-09 2.541e+06  < 2e-16 ***
#> beta:(Intercept)  1.687e+01  1.449e+00 1.164e+01  < 2e-16 ***
#> beta:temp         3.709e-02  2.917e-03 1.271e+01  < 2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Confidence intervals (95%):
#>                        3%     98%
#> alpha:(Intercept)  0.1949  0.3086
#> alpha:batch1       0.7709  0.8850
#> alpha:batch2       0.5233  0.6458
#> alpha:batch3       0.6555  0.7689
#> alpha:batch4       0.4137  0.5158
#> alpha:batch5       0.4719  0.5748
#> alpha:batch6       0.3963  0.5030
#> alpha:batch7       0.1762  0.2741
#> alpha:batch8       0.1527  0.2529
#> alpha:batch9       0.0840  0.2012
#> alpha:temp         0.0066  0.0066
#> beta:(Intercept)  14.0346 19.7151
#> beta:temp          0.0314  0.0428
#> 
#> Link functions:
#> alpha: log
#> beta: log
#> 
#> Fitted parameter means:
#> alpha: 19.06
#> beta: 2.377e+13
#> gamma: 1
#> delta: 0
#> lambda: 1
#> 
#> Model fit statistics:
#> Number of observations: 32 
#> Number of parameters: 13 
#> Residual degrees of freedom: 19 
#> Log-likelihood: 96.75 
#> AIC: -167.5 
#> BIC: -148.4 
#> RMSE: 0.02249 
#> Efron's R2: 0.9546 
#> Mean Absolute Error: 0.01503 
#> 
#> Convergence status: Failed 
#> Iterations: 76 
#> 

# Interpretation:
# - Alpha (mean): Different batches have different baseline yields
#   Temperature affects yield transformation
# - Beta (precision): Higher temperatures may produce more variable yields

# Example 2: Full model with all physical-chemical properties
fit_kw_full <- gkwreg(
  yield ~ gravity + pressure + temp10 + temp |
    temp10 + temp,
  data = GasolineYield,
  family = "kw"
)
#> Warning: NaNs produced
summary(fit_kw_full)
#> 
#> Generalized Kumaraswamy Regression Model Summary
#> 
#> Family: kw 
#> 
#> Call:
#> gkwreg(formula = yield ~ gravity + pressure + temp10 + temp | 
#>     temp10 + temp, data = GasolineYield, family = "kw")
#> 
#> Residuals:
#>     Min  Q1.25%  Median    Mean  Q3.75%     Max 
#> -0.0516 -0.0169 -0.0034 -0.0053  0.0033  0.0444 
#> 
#> Coefficients:
#>                    Estimate Std. Error z value Pr(>|z|)    
#> alpha:(Intercept)  0.457548   0.166036   2.756  0.00586 ** 
#> alpha:gravity      0.001539   0.003912   0.393  0.69400    
#> alpha:pressure     0.024644   0.012936   1.905  0.05676 .  
#> alpha:temp10      -0.001364        NaN     NaN      NaN    
#> alpha:temp         0.005833        NaN     NaN      NaN    
#> beta:(Intercept)  -1.766739   3.731457  -0.473  0.63588    
#> beta:temp10        0.059930   0.013139   4.561 5.08e-06 ***
#> beta:temp          0.008332   0.001735   4.802 1.57e-06 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Confidence intervals (95%):
#>                        3%    98%
#> alpha:(Intercept)  0.1321 0.7830
#> alpha:gravity     -0.0061 0.0092
#> alpha:pressure    -0.0007 0.0500
#> alpha:temp10          NaN    NaN
#> alpha:temp            NaN    NaN
#> beta:(Intercept)  -9.0803 5.5468
#> beta:temp10        0.0342 0.0857
#> beta:temp          0.0049 0.0117
#> 
#> Link functions:
#> alpha: log
#> beta: log
#> 
#> Fitted parameter means:
#> alpha: 9.935
#> beta: 59340415
#> gamma: 1
#> delta: 0
#> lambda: 1
#> 
#> Model fit statistics:
#> Number of observations: 32 
#> Number of parameters: 8 
#> Residual degrees of freedom: 24 
#> Log-likelihood: 80.45 
#> AIC: -144.9 
#> BIC: -133.2 
#> RMSE: 0.02081 
#> Efron's R2: 0.9611 
#> Mean Absolute Error: 0.01542 
#> 
#> Convergence status: Successful 
#> Iterations: 70 
#> 

# Interpretation:
# - Mean model captures effects of crude oil properties
# - Precision varies with vaporization temperatures

# Example 3: Exponentiated Kumaraswamy for extreme yields
# Some batches may produce unusually high/low yields
fit_ekw <- gkwreg(
  yield ~ batch + temp | # alpha: batch effects
    temp | # beta: temperature precision
    batch, # lambda: batch-specific tail behavior
  data = GasolineYield,
  family = "ekw"
)
summary(fit_ekw)
#> 
#> Generalized Kumaraswamy Regression Model Summary
#> 
#> Family: ekw 
#> 
#> Call:
#> gkwreg(formula = yield ~ batch + temp | temp | batch, data = GasolineYield, 
#>     family = "ekw")
#> 
#> Residuals:
#>     Min  Q1.25%  Median    Mean  Q3.75%     Max 
#> -0.0360 -0.0143 -0.0047  0.0017  0.0049  0.2229 
#> 
#> Coefficients:
#>                     Estimate Std. Error z value Pr(>|z|)
#> alpha:(Intercept)  -0.318429        NaN     NaN      NaN
#> alpha:batch1        0.533730        NaN     NaN      NaN
#> alpha:batch2        0.525889        NaN     NaN      NaN
#> alpha:batch3        0.690016        NaN     NaN      NaN
#> alpha:batch4        0.350257        NaN     NaN      NaN
#> alpha:batch5        0.543610        NaN     NaN      NaN
#> alpha:batch6        0.311936        NaN     NaN      NaN
#> alpha:batch7        0.193041        NaN     NaN      NaN
#> alpha:batch8        0.229574        NaN     NaN      NaN
#> alpha:batch9        0.161059        NaN     NaN      NaN
#> alpha:temp          0.006078        NaN     NaN      NaN
#> beta:(Intercept)   10.512482        NaN     NaN      NaN
#> beta:temp           0.011623        NaN     NaN      NaN
#> lambda:(Intercept)  0.356490        NaN     NaN      NaN
#> lambda:batch1      48.038859        NaN     NaN      NaN
#> lambda:batch2       2.279667        NaN     NaN      NaN
#> lambda:batch3       0.561476        NaN     NaN      NaN
#> lambda:batch4       7.008866        NaN     NaN      NaN
#> lambda:batch5      -0.081997        NaN     NaN      NaN
#> lambda:batch6       8.826415        NaN     NaN      NaN
#> lambda:batch7       1.644509        NaN     NaN      NaN
#> lambda:batch8       0.014281        NaN     NaN      NaN
#> lambda:batch9      -0.021619        NaN     NaN      NaN
#> 
#> Confidence intervals (95%):
#>                     3% 98%
#> alpha:(Intercept)  NaN NaN
#> alpha:batch1       NaN NaN
#> alpha:batch2       NaN NaN
#> alpha:batch3       NaN NaN
#> alpha:batch4       NaN NaN
#> alpha:batch5       NaN NaN
#> alpha:batch6       NaN NaN
#> alpha:batch7       NaN NaN
#> alpha:batch8       NaN NaN
#> alpha:batch9       NaN NaN
#> alpha:temp         NaN NaN
#> beta:(Intercept)   NaN NaN
#> beta:temp          NaN NaN
#> lambda:(Intercept) NaN NaN
#> lambda:batch1      NaN NaN
#> lambda:batch2      NaN NaN
#> lambda:batch3      NaN NaN
#> lambda:batch4      NaN NaN
#> lambda:batch5      NaN NaN
#> lambda:batch6      NaN NaN
#> lambda:batch7      NaN NaN
#> lambda:batch8      NaN NaN
#> lambda:batch9      NaN NaN
#> 
#> Link functions:
#> alpha: log
#> beta: log
#> lambda: log
#> 
#> Fitted parameter means:
#> alpha: 8.419
#> beta: 2313893
#> gamma: 1
#> delta: 0
#> lambda: 1.302e+20
#> 
#> Model fit statistics:
#> Number of observations: 32 
#> Number of parameters: 23 
#> Residual degrees of freedom: 9 
#> Log-likelihood: 125.5 
#> AIC: -205 
#> BIC: -171.3 
#> RMSE: 0.04243 
#> Efron's R2: 0.8384 
#> Mean Absolute Error: 0.01919 
#> 
#> Convergence status: Failed 
#> Iterations: 105 
#> 

# Interpretation:
# - Lambda varies by batch: Some crude oils have more extreme
#   yield distributions (heavy tails for very high/low yields)

# Model comparison: Does tail flexibility improve fit?
anova(fit_kw, fit_ekw)
#> Analysis of Deviance Table
#> 
#> Model 1: yield ~ batch + temp | temp
#> Model 2: yield ~ batch + temp | temp | batch
#> 
#>         Resid. Df Resid. Dev Df Deviance Pr(>Chi)    
#> fit_kw   19.00000 -193.49997                         
#> fit_ekw   9.00000 -250.98137 10 57.48139  < 1e-04 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Diagnostic plots
par(mfrow = c(2, 2))
plot(fit_kw, which = c(1, 2, 4, 5))
#> Simulating envelope ( 100 iterations): .......... Done!

par(mfrow = c(1, 1))
# }
```
