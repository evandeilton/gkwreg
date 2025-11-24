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
#> Warning: NaNs produced
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
#> -0.4914 -0.0175 -0.0034 -0.0249  0.0065  0.0745 
#> 
#> Coefficients:
#>                    Estimate Std. Error z value Pr(>|z|)    
#> alpha:(Intercept)  0.319842   0.026078  12.265  < 2e-16 ***
#> alpha:batch1       0.826719   0.028879  28.627  < 2e-16 ***
#> alpha:batch2       0.582466   0.030836  18.889  < 2e-16 ***
#> alpha:batch3       0.710004   0.028636  24.794  < 2e-16 ***
#> alpha:batch4       0.461921   0.025812  17.895  < 2e-16 ***
#> alpha:batch5       0.520953   0.026011  20.028  < 2e-16 ***
#> alpha:batch6       0.447006   0.026950  16.586  < 2e-16 ***
#> alpha:batch7       0.222963   0.024710   9.023  < 2e-16 ***
#> alpha:batch8       0.200138   0.024892   8.040 8.96e-16 ***
#> alpha:batch9       0.138999   0.029271   4.749 2.05e-06 ***
#> alpha:temp         0.006445        NaN     NaN      NaN    
#> beta:(Intercept)  18.458156   1.405043  13.137  < 2e-16 ***
#> beta:temp          0.034533   0.002942  11.737  < 2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Confidence intervals (95%):
#>                        3%     98%
#> alpha:(Intercept)  0.2687  0.3710
#> alpha:batch1       0.7701  0.8833
#> alpha:batch2       0.5220  0.6429
#> alpha:batch3       0.6539  0.7661
#> alpha:batch4       0.4113  0.5125
#> alpha:batch5       0.4700  0.5719
#> alpha:batch6       0.3942  0.4998
#> alpha:batch7       0.1745  0.2714
#> alpha:batch8       0.1514  0.2489
#> alpha:batch9       0.0816  0.1964
#> alpha:temp            NaN     NaN
#> beta:(Intercept)  15.7043 21.2120
#> beta:temp          0.0288  0.0403
#> 
#> Link functions:
#> alpha: log
#> beta: log
#> 
#> Fitted parameter means:
#> alpha: 19.51
#> beta: 4.025e+13
#> gamma: 1
#> delta: 0
#> lambda: 1
#> 
#> Model fit statistics:
#> Number of observations: 32 
#> Number of parameters: 13 
#> Residual degrees of freedom: 19 
#> Log-likelihood: 96.81 
#> AIC: -167.6 
#> BIC: -148.6 
#> RMSE: 0.09951 
#> Efron's R2: 0.111 
#> Mean Absolute Error: 0.04167 
#> 
#> Convergence status: Failed 
#> Iterations: 93 
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
#> beta:(Intercept)  -1.766735   3.731464  -0.473  0.63588    
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
#> beta: 59340318
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
#> -0.0356 -0.0133 -0.0020  0.0048  0.0082  0.2229 
#> 
#> Coefficients:
#>                    Estimate Std. Error z value Pr(>|z|)
#> alpha:(Intercept)  -0.36374        NaN     NaN      NaN
#> alpha:batch1        0.52791        NaN     NaN      NaN
#> alpha:batch2        0.49225        NaN     NaN      NaN
#> alpha:batch3        0.67943        NaN     NaN      NaN
#> alpha:batch4        0.31973        NaN     NaN      NaN
#> alpha:batch5        0.50498        NaN     NaN      NaN
#> alpha:batch6        0.30367        NaN     NaN      NaN
#> alpha:batch7        0.16380        NaN     NaN      NaN
#> alpha:batch8        0.19265        NaN     NaN      NaN
#> alpha:batch9        0.09724        NaN     NaN      NaN
#> alpha:temp          0.00610        NaN     NaN      NaN
#> beta:(Intercept)    9.81853        NaN     NaN      NaN
#> beta:temp           0.01178        NaN     NaN      NaN
#> lambda:(Intercept) -0.14664        NaN     NaN      NaN
#> lambda:batch1      36.09593        NaN     NaN      NaN
#> lambda:batch2       3.21155        NaN     NaN      NaN
#> lambda:batch3       1.08510        NaN     NaN      NaN
#> lambda:batch4       7.26510        NaN     NaN      NaN
#> lambda:batch5       0.63205        NaN     NaN      NaN
#> lambda:batch6       7.47889        NaN     NaN      NaN
#> lambda:batch7       2.23889        NaN     NaN      NaN
#> lambda:batch8       0.72835        NaN     NaN      NaN
#> lambda:batch9       1.47656        NaN     NaN      NaN
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
#> alpha: 7.912
#> beta: 1225127
#> gamma: 1
#> delta: 0
#> lambda: 5.123e+14
#> 
#> Model fit statistics:
#> Number of observations: 32 
#> Number of parameters: 23 
#> Residual degrees of freedom: 9 
#> Log-likelihood: 114.3 
#> AIC: -182.5 
#> BIC: -148.8 
#> RMSE: 0.04235 
#> Efron's R2: 0.839 
#> Mean Absolute Error: 0.01908 
#> 
#> Convergence status: Failed 
#> Iterations: 99 
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
#>         Resid. Df Resid. Dev Df Deviance   Pr(>Chi)    
#> fit_kw   19.00000 -193.61594                           
#> fit_ekw   9.00000 -228.53221 10 34.91627 0.00012904 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Diagnostic plots
par(mfrow = c(2, 2))
plot(fit_kw, which = c(1, 2, 4, 5))
#> Simulating envelope ( 100 iterations): .......... Done!

par(mfrow = c(1, 1))
# }
```
