# Intraocular Gas Decay in Retinal Surgery

Longitudinal data on the recorded decay of intraocular gas
(perfluoropropane) in complex retinal surgeries. The dataset tracks the
proportion of gas remaining over time following vitrectomy procedures.

## Usage

``` r
retinal
```

## Format

A data frame with 40 observations on 7 variables:

- ID:

  integer. Patient identification number for longitudinal tracking.

- Gas:

  numeric. Proportion of intraocular gas remaining (0-1 scale). Response
  variable measuring the fraction of perfluoropropane gas still present
  in the vitreous cavity.

- Time:

  numeric. Time point of measurement (days or weeks post-surgery).

- LogT:

  numeric. Logarithm of time, log(Time). Used to linearize the
  exponential decay pattern.

- LogT2:

  numeric. Squared logarithm of time, (log(Time))^2. Captures nonlinear
  decay patterns.

- Level:

  factor. Initial gas concentration level at the time of injection.
  Different starting concentrations affect decay kinetics.

## Source

Based on clinical data from vitreoretinal surgery patients. Originally
analyzed in Song and Tan (2000).

## Details

This longitudinal dataset comes from a study of gas decay following
vitreoretinal surgery. Perfluoropropane (C3F8) is commonly used as a
temporary tamponade agent in retinal detachment repair and other complex
vitreoretinal procedures.

**Clinical background:** During vitrectomy for retinal detachment, gas
bubbles are injected into the vitreous cavity to help reattach the
retina by providing internal tamponade. The gas gradually absorbs and
dissipates over time. Understanding the decay rate is important for:

- Predicting when patients can resume normal activities (esp. air
  travel)

- Assessing treatment efficacy

- Planning follow-up examinations

**Decay kinetics:** Gas decay typically follows a nonlinear pattern that
can be approximated by exponential or power-law functions. The log
transformation (LogT, LogT2) helps linearize these relationships for
regression modeling.

**Data structure:** This is a longitudinal/panel dataset with repeated
measurements on the same patients over time. Correlation structures
(exchangeable, AR(1), etc.) should be considered when modeling.

The proportional nature of the gas variable (bounded between 0 and 1)
makes this dataset ideal for:

- Simplex marginal models (original application by Song & Tan 2000)

- Beta regression with longitudinal correlation structures

- Kumaraswamy regression with heteroscedastic errors

## References

Meyers, S.M., Ambler, J.S., Tan, M., Werner, J.C., and Huang, S.S.
(1992). Variation of Perfluoropropane Disappearance After Vitrectomy.
*Retina*, **12**, 359–363.

Song, P.X.-K., and Tan, M. (2000). Marginal Models for Longitudinal
Continuous Proportional Data. *Biometrics*, **56**, 496–502.
[doi:10.1111/j.0006-341x.2000.00496.x](https://doi.org/10.1111/j.0006-341x.2000.00496.x)

Song, P.X.-K., Qiu, Z., and Tan, M. (2004). Modelling Heterogeneous
Dispersion in Marginal Models for Longitudinal Proportional Data.
*Biometrical Journal*, **46**, 540–553.

Qiu, Z., Song, P.X.-K., and Tan, M. (2008). Simplex Mixed-Effects Models
for Longitudinal Proportional Data. *Scandinavian Journal of
Statistics*, **35**, 577–596.
[doi:10.1111/j.1467-9469.2008.00603.x](https://doi.org/10.1111/j.1467-9469.2008.00603.x)

Zhang, P., Qiu, Z., and Shi, C. (2016). simplexreg: An R Package for
Regression Analysis of Proportional Data Using the Simplex Distribution.
*Journal of Statistical Software*, **71**(11), 1–21.
[doi:10.18637/jss.v071.i11](https://doi.org/10.18637/jss.v071.i11)

## Examples

``` r
# \donttest{
require(gkwreg)
require(gkwdist)

data(retinal)

# Example 1: Nonlinear time decay with level effects
# Model gas decay as quadratic function of log-time
# Allow precision to vary by initial gas concentration
fit_kw <- gkwreg(
  Gas ~ LogT + LogT2 + Level |
    Level,
  data = retinal,
  family = "kw"
)
summary(fit_kw)
#> 
#> Generalized Kumaraswamy Regression Model Summary
#> 
#> Family: kw 
#> 
#> Call:
#> gkwreg(formula = Gas ~ LogT + LogT2 + Level | Level, data = retinal, 
#>     family = "kw")
#> 
#> Residuals:
#>     Min  Q1.25%  Median    Mean  Q3.75%     Max 
#> -0.5665 -0.2076  0.0087 -0.0436  0.0971  0.4073 
#> 
#> Coefficients:
#>                   Estimate Std. Error z value Pr(>|z|)    
#> alpha:(Intercept)  1.69138    0.22256   7.600 2.97e-14 ***
#> alpha:LogT        -0.25598    0.25054  -1.022  0.30693    
#> alpha:LogT2       -0.12019    0.06061  -1.983  0.04737 *  
#> alpha:Level        0.51892    0.15828   3.279  0.00104 ** 
#> beta:(Intercept)  -0.40918    0.09481  -4.316 1.59e-05 ***
#> beta:Level         0.08475    0.13199   0.642  0.52084    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Confidence intervals (95%):
#>                        3%     98%
#> alpha:(Intercept)  1.2552  2.1276
#> alpha:LogT        -0.7470  0.2351
#> alpha:LogT2       -0.2390 -0.0014
#> alpha:Level        0.2087  0.8291
#> beta:(Intercept)  -0.5950 -0.2233
#> beta:Level        -0.1740  0.3434
#> 
#> Link functions:
#> alpha: log
#> beta: log
#> 
#> Fitted parameter means:
#> alpha: 2.672
#> beta: 0.6776
#> gamma: 1
#> delta: 0
#> lambda: 1
#> 
#> Model fit statistics:
#> Number of observations: 181 
#> Number of parameters: 6 
#> Residual degrees of freedom: 175 
#> Log-likelihood: 103.7 
#> AIC: -195.4 
#> BIC: -176.2 
#> RMSE: 0.2194 
#> Efron's R2: 0.4994 
#> Mean Absolute Error: 0.177 
#> 
#> Convergence status: Successful 
#> Iterations: 21 
#> 

# Interpretation:
# - Alpha: Decay curve shape varies by initial gas concentration
#   LogT + LogT2 capture nonlinear exponential-like decay
# - Beta: Precision differs by concentration level
#   Higher concentration may produce more/less variable decay

# Example 2: Heteroscedastic model
# Variability in gas proportion may change over time
fit_kw_hetero <- gkwreg(
  Gas ~ LogT + LogT2 + Level |
    Level + LogT,
  data = retinal,
  family = "kw"
)
summary(fit_kw_hetero)
#> 
#> Generalized Kumaraswamy Regression Model Summary
#> 
#> Family: kw 
#> 
#> Call:
#> gkwreg(formula = Gas ~ LogT + LogT2 + Level | Level + LogT, data = retinal, 
#>     family = "kw")
#> 
#> Residuals:
#>     Min  Q1.25%  Median    Mean  Q3.75%     Max 
#> -0.5261 -0.1223  0.0192 -0.0081  0.1092  0.5323 
#> 
#> Coefficients:
#>                   Estimate Std. Error z value Pr(>|z|)    
#> alpha:(Intercept)  0.90707    0.34706   2.614  0.00896 ** 
#> alpha:LogT         0.40385    0.29920   1.350  0.17710    
#> alpha:LogT2       -0.19493    0.05992  -3.253  0.00114 ** 
#> alpha:Level        0.41066    0.13956   2.943  0.00326 ** 
#> beta:(Intercept)  -1.25839    0.18008  -6.988 2.79e-12 ***
#> beta:Level         0.10830    0.13838   0.783  0.43384    
#> beta:LogT          0.51739    0.08071   6.410 1.45e-10 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Confidence intervals (95%):
#>                        3%     98%
#> alpha:(Intercept)  0.2268  1.5873
#> alpha:LogT        -0.1826  0.9903
#> alpha:LogT2       -0.3124 -0.0775
#> alpha:Level        0.1371  0.6842
#> beta:(Intercept)  -1.6113 -0.9054
#> beta:Level        -0.1629  0.3795
#> beta:LogT          0.3592  0.6756
#> 
#> Link functions:
#> alpha: log
#> beta: log
#> 
#> Fitted parameter means:
#> alpha: 2.25
#> beta: 1.06
#> gamma: 1
#> delta: 0
#> lambda: 1
#> 
#> Model fit statistics:
#> Number of observations: 181 
#> Number of parameters: 7 
#> Residual degrees of freedom: 174 
#> Log-likelihood: 125.4 
#> AIC: -236.8 
#> BIC: -214.4 
#> RMSE: 0.2059 
#> Efron's R2: 0.5592 
#> Mean Absolute Error: 0.1566 
#> 
#> Convergence status: Successful 
#> Iterations: 29 
#> 

# Interpretation:
# - Beta: Precision varies with both level and time
#   Early measurements may be more variable than late measurements

# Test heteroscedasticity
anova(fit_kw, fit_kw_hetero)
#> Analysis of Deviance Table
#> 
#> Model 1: Gas ~ LogT + LogT2 + Level | Level
#> Model 2: Gas ~ LogT + LogT2 + Level | Level + LogT
#> 
#>               Resid. Df Resid. Dev Df Deviance Pr(>Chi)    
#> fit_kw        175.00000 -207.36034                         
#> fit_kw_hetero 174.00000 -250.75874  1 43.39840  < 1e-04 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Example 3: Exponentiated Kumaraswamy for decay tails
# Gas decay may show different tail behavior at extreme time points
# (very fast initial decay or very slow residual decay)
fit_ekw <- gkwreg(
  Gas ~ LogT + LogT2 + Level | # alpha: decay curve
    Level + LogT | # beta: heteroscedasticity
    Level, # lambda: tail heaviness by level
  data = retinal,
  family = "ekw"
)
summary(fit_ekw)
#> 
#> Generalized Kumaraswamy Regression Model Summary
#> 
#> Family: ekw 
#> 
#> Call:
#> gkwreg(formula = Gas ~ LogT + LogT2 + Level | Level + LogT | 
#>     Level, data = retinal, family = "ekw")
#> 
#> Residuals:
#>     Min  Q1.25%  Median    Mean  Q3.75%     Max 
#> -0.4958 -0.0799  0.0651  0.0446  0.1693  0.5107 
#> 
#> Coefficients:
#>                    Estimate Std. Error z value Pr(>|z|)    
#> alpha:(Intercept)   1.36711    1.11962   1.221  0.22207    
#> alpha:LogT          0.27552    0.32922   0.837  0.40267    
#> alpha:LogT2        -0.16545    0.05963  -2.775  0.00553 ** 
#> alpha:Level        -0.38899    0.39384  -0.988  0.32330    
#> beta:(Intercept)   -1.36470    0.34796  -3.922 8.78e-05 ***
#> beta:Level          0.25158    0.19003   1.324  0.18555    
#> beta:LogT           0.54369    0.12922   4.207 2.58e-05 ***
#> lambda:(Intercept) -0.32755    0.76219  -0.430  0.66738    
#> lambda:Level        0.92335    0.39316   2.349  0.01885 *  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Confidence intervals (95%):
#>                         3%     98%
#> alpha:(Intercept)  -0.8273  3.5615
#> alpha:LogT         -0.3698  0.9208
#> alpha:LogT2        -0.2823 -0.0486
#> alpha:Level        -1.1609  0.3829
#> beta:(Intercept)   -2.0467 -0.6827
#> beta:Level         -0.1209  0.6240
#> beta:LogT           0.2904  0.7970
#> lambda:(Intercept) -1.8214  1.1663
#> lambda:Level        0.1528  1.6939
#> 
#> Link functions:
#> alpha: log
#> beta: log
#> lambda: log
#> 
#> Fitted parameter means:
#> alpha: 2.771
#> beta: 1.082
#> gamma: 1
#> delta: 0
#> lambda: 1.08
#> 
#> Model fit statistics:
#> Number of observations: 181 
#> Number of parameters: 9 
#> Residual degrees of freedom: 172 
#> Log-likelihood: 129.2 
#> AIC: -240.3 
#> BIC: -211.5 
#> RMSE: 0.209 
#> Efron's R2: 0.5458 
#> Mean Absolute Error: 0.1644 
#> 
#> Convergence status: Successful 
#> Iterations: 37 
#> 

# Interpretation:
# - Lambda varies by level: Different initial concentrations may have
#   different rates of extreme decay (very fast or very slow residual gas)
# - Important for predicting complete absorption time

# Example 4: McDonald distribution for asymmetric decay
# Alternative parameterization for skewed decay patterns
fit_mc <- gkwreg(
  Gas ~ LogT + LogT2 + Level | # gamma
    LogT + Level | # delta
    Level, # lambda
  data = retinal,
  family = "mc",
  control = gkw_control(
    method = "BFGS",
    maxit = 1500,
    reltol = 1e-8
  )
)
#> Warning: NaNs produced
summary(fit_mc)
#> 
#> Generalized Kumaraswamy Regression Model Summary
#> 
#> Family: mc 
#> 
#> Call:
#> gkwreg(formula = Gas ~ LogT + LogT2 + Level | LogT + Level | 
#>     Level, data = retinal, family = "mc", control = gkw_control(method = "BFGS", 
#>     maxit = 1500, reltol = 1e-08))
#> 
#> Residuals:
#>     Min  Q1.25%  Median    Mean  Q3.75%     Max 
#> -0.5652 -0.1159  0.0034 -0.0098  0.1214  0.5185 
#> 
#> Coefficients:
#>                    Estimate Std. Error z value Pr(>|z|)    
#> gamma:(Intercept)   0.99232    1.85227   0.536    0.592    
#> gamma:LogT          1.11127    0.23830   4.663 3.11e-06 ***
#> gamma:LogT2        -0.44693    0.05657  -7.900 2.79e-15 ***
#> gamma:Level         0.28447    1.38833   0.205    0.838    
#> delta:(Intercept)  -6.27067    0.70158  -8.938  < 2e-16 ***
#> delta:LogT          0.17905        NaN     NaN      NaN    
#> delta:Level        -0.38972        NaN     NaN      NaN    
#> lambda:(Intercept)  0.29057    1.83059   0.159    0.874    
#> lambda:Level       -0.07583    1.42998  -0.053    0.958    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Confidence intervals (95%):
#>                         3%     98%
#> gamma:(Intercept)  -2.6381  4.6227
#> gamma:LogT          0.6442  1.5783
#> gamma:LogT2        -0.5578 -0.3360
#> gamma:Level        -2.4366  3.0055
#> delta:(Intercept)  -7.6458 -4.8956
#> delta:LogT             NaN     NaN
#> delta:Level            NaN     NaN
#> lambda:(Intercept) -3.2973  3.8785
#> lambda:Level       -2.8785  2.7269
#> 
#> Link functions:
#> gamma: log
#> delta: logit
#> lambda: log
#> 
#> Fitted parameter means:
#> alpha: 1
#> beta: 1
#> gamma: 2.848
#> delta: 0.0269
#> lambda: 1.317
#> 
#> Model fit statistics:
#> Number of observations: 181 
#> Number of parameters: 9 
#> Residual degrees of freedom: 172 
#> Log-likelihood: 66.48 
#> AIC: -115 
#> BIC: -86.18 
#> RMSE: 0.2145 
#> Efron's R2: 0.5216 
#> Mean Absolute Error: 0.1631 
#> 
#> Convergence status: Successful 
#> Iterations: 26 
#> 

# Model comparison
AIC(fit_kw, fit_kw_hetero, fit_ekw, fit_mc)
#>               df       AIC
#> fit_kw         6 -195.3603
#> fit_kw_hetero  7 -236.7587
#> fit_ekw        9 -240.3186
#> fit_mc         9 -114.9668
# }
```
