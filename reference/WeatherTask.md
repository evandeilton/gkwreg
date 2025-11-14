# Weather Task with Priming and Precise and Imprecise Probabilities

Data from a cognitive psychology experiment on probabilistic learning
and probability judgments. Participants estimated probabilities for
weather events under different priming and precision conditions.

## Usage

``` r
WeatherTask
```

## Format

A data frame with 345 observations on 4 variables:

- agreement:

  numeric. Probability indicated by participants, or the average between
  minimum and maximum estimates in the imprecise condition. Response
  variable scaled to (0, 1).

- priming:

  factor with levels `two-fold` (case prime) and `seven-fold` (class
  prime). Indicates the partition priming condition.

- eliciting:

  factor with levels `precise` and `imprecise` (lower and upper limit).
  Indicates whether participants gave point estimates or interval
  estimates.

## Source

Taken from Smithson et al. (2011) supplements.

## Details

All participants in the study were either first- or second-year
undergraduate students in psychology, none of whom had a strong
background in probability or were familiar with imprecise probability
theories.

**Task description:** Participants were asked: "What is the probability
that the temperature at Canberra airport on Sunday will be higher than
'specified temperature'?"

**Experimental manipulations:**

- **Priming:** Two-fold (simple binary: above/below) vs. seven-fold
  (multiple temperature categories)

- **Eliciting:** Precise (single probability estimate) vs. imprecise
  (lower and upper bounds)

The study examines how partition priming (number of response categories)
and elicitation format affect probability judgments. Classical findings
suggest that more categories (seven-fold) lead to different probability
assessments than binary categories (two-fold).

## References

Smithson, M., Merkle, E.C., and Verkuilen, J. (2011). Beta Regression
Finite Mixture Models of Polarization and Priming. *Journal of
Educational and Behavioral Statistics*, **36**(6), 804–831.
[doi:10.3102/1076998610396893](https://doi.org/10.3102/1076998610396893)

Smithson, M., and Segale, C. (2009). Partition Priming in Judgments of
Imprecise Probabilities. *Journal of Statistical Theory and Practice*,
**3**(1), 169–181.

## Examples

``` r
# \donttest{
require(gkwreg)
require(gkwdist)

data(WeatherTask)

# Example 1: Main effects model
# Probability judgments affected by priming and elicitation format
fit_kw <- gkwreg(
  agreement ~ priming + eliciting,
  data = WeatherTask,
  family = "kw"
)
summary(fit_kw)
#> 
#> Generalized Kumaraswamy Regression Model Summary
#> 
#> Family: kw 
#> 
#> Call:
#> gkwreg(formula = agreement ~ priming + eliciting, data = WeatherTask, 
#>     family = "kw")
#> 
#> Residuals:
#>     Min  Q1.25%  Median    Mean  Q3.75%     Max 
#> -0.2941 -0.1053 -0.0494 -0.0063  0.0899  0.7053 
#> 
#> Coefficients:
#>                          Estimate Std. Error z value Pr(>|z|)    
#> alpha:(Intercept)         0.39534    0.06067   6.516 7.24e-11 ***
#> alpha:primingseven-fold  -0.18903    0.05241  -3.606  0.00031 ***
#> alpha:elicitingimprecise  0.21471    0.05233   4.103 4.07e-05 ***
#> beta:(Intercept)          1.81969    0.09799  18.571  < 2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Confidence intervals (95%):
#>                               3%     98%
#> alpha:(Intercept)         0.2764  0.5143
#> alpha:primingseven-fold  -0.2918 -0.0863
#> alpha:elicitingimprecise  0.1122  0.3173
#> beta:(Intercept)          1.6276  2.0117
#> 
#> Link functions:
#> alpha: log
#> beta: log
#> 
#> Fitted parameter means:
#> alpha: 1.524
#> beta: 6.164
#> gamma: 1
#> delta: 0
#> lambda: 1
#> 
#> Model fit statistics:
#> Number of observations: 345 
#> Number of parameters: 4 
#> Residual degrees of freedom: 341 
#> Log-likelihood: 201.6 
#> AIC: -395.2 
#> BIC: -379.9 
#> RMSE: 0.1541 
#> Efron's R2: 0.08752 
#> Mean Absolute Error: 0.1253 
#> 
#> Convergence status: Successful 
#> Iterations: 25 
#> 

# Interpretation:
# - Alpha: Seven-fold priming may shift probability estimates
#   Imprecise elicitation may produce different mean estimates

# Example 2: Interaction model with heteroscedasticity
# Priming effects may differ by elicitation format
# Variability may also depend on conditions
fit_kw_interact <- gkwreg(
  agreement ~ priming * eliciting |
    priming + eliciting,
  data = WeatherTask,
  family = "kw"
)
summary(fit_kw_interact)
#> 
#> Generalized Kumaraswamy Regression Model Summary
#> 
#> Family: kw 
#> 
#> Call:
#> gkwreg(formula = agreement ~ priming * eliciting | priming + 
#>     eliciting, data = WeatherTask, family = "kw")
#> 
#> Residuals:
#>     Min  Q1.25%  Median    Mean  Q3.75%     Max 
#> -0.2930 -0.1159 -0.0320 -0.0052  0.0910  0.6891 
#> 
#> Coefficients:
#>                                            Estimate Std. Error z value Pr(>|z|)
#> alpha:(Intercept)                           0.26198    0.08991   2.914  0.00357
#> alpha:primingseven-fold                     0.10566    0.10814   0.977  0.32850
#> alpha:elicitingimprecise                    0.21199    0.11980   1.769  0.07681
#> alpha:primingseven-fold:elicitingimprecise  0.08181    0.10438   0.784  0.43317
#> beta:(Intercept)                            1.46021    0.15558   9.385  < 2e-16
#> beta:primingseven-fold                      0.85229    0.20297   4.199 2.68e-05
#> beta:elicitingimprecise                     0.09174    0.19953   0.460  0.64566
#>                                               
#> alpha:(Intercept)                          ** 
#> alpha:primingseven-fold                       
#> alpha:elicitingimprecise                   .  
#> alpha:primingseven-fold:elicitingimprecise    
#> beta:(Intercept)                           ***
#> beta:primingseven-fold                     ***
#> beta:elicitingimprecise                       
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Confidence intervals (95%):
#>                                                 3%    98%
#> alpha:(Intercept)                           0.0858 0.4382
#> alpha:primingseven-fold                    -0.1063 0.3176
#> alpha:elicitingimprecise                   -0.0228 0.4468
#> alpha:primingseven-fold:elicitingimprecise -0.1228 0.2864
#> beta:(Intercept)                            1.1553 1.7651
#> beta:primingseven-fold                      0.4545 1.2501
#> beta:elicitingimprecise                    -0.2993 0.4828
#> 
#> Link functions:
#> alpha: log
#> beta: log
#> 
#> Fitted parameter means:
#> alpha: 1.566
#> beta: 7.428
#> gamma: 1
#> delta: 0
#> lambda: 1
#> 
#> Model fit statistics:
#> Number of observations: 345 
#> Number of parameters: 7 
#> Residual degrees of freedom: 338 
#> Log-likelihood: 210.9 
#> AIC: -407.8 
#> BIC: -380.9 
#> RMSE: 0.1537 
#> Efron's R2: 0.09164 
#> Mean Absolute Error: 0.124 
#> 
#> Convergence status: Successful 
#> Iterations: 31 
#> 

# Interpretation:
# - Alpha: Interaction tests if partition priming works differently
#   for precise vs. imprecise probability judgments
# - Beta: Precision varies by experimental condition

# Test interaction
anova(fit_kw, fit_kw_interact)
#> Analysis of Deviance Table
#> 
#> Model 1: agreement ~ priming + eliciting
#> Model 2: agreement ~ priming * eliciting | priming + eliciting
#> 
#>                 Resid. Df Resid. Dev Df Deviance   Pr(>Chi)    
#> fit_kw          341.00000 -403.23406                           
#> fit_kw_interact 338.00000 -421.81480  3 18.58074 0.00033376 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Example 3: McDonald distribution for polarized responses
# Probability judgments often show polarization (clustering at extremes)
# particularly under certain priming conditions
fit_mc <- gkwreg(
  agreement ~ priming * eliciting | # gamma
    priming * eliciting | # delta
    priming, # lambda: priming affects polarization
  data = WeatherTask,
  family = "mc",
  control = gkw_control(method = "BFGS", maxit = 1500)
)
summary(fit_mc)
#> 
#> Generalized Kumaraswamy Regression Model Summary
#> 
#> Family: mc 
#> 
#> Call:
#> gkwreg(formula = agreement ~ priming * eliciting | priming * 
#>     eliciting | priming, data = WeatherTask, family = "mc", control = gkw_control(method = "BFGS", 
#>     maxit = 1500))
#> 
#> Residuals:
#>     Min  Q1.25%  Median    Mean  Q3.75%     Max 
#> -0.2027 -0.0727 -0.0416  0.0351  0.1084  0.7564 
#> 
#> Coefficients:
#>                                             Estimate Std. Error z value
#> gamma:(Intercept)                           0.090319   0.222111   0.407
#> gamma:primingseven-fold                     0.036067   0.426063   0.085
#> gamma:elicitingimprecise                    0.054185   0.522357   0.104
#> gamma:primingseven-fold:elicitingimprecise  0.023760   0.566762   0.042
#> delta:(Intercept)                          -0.039825   0.534925  -0.074
#> delta:primingseven-fold                    -0.010793   0.644386  -0.017
#> delta:elicitingimprecise                   -0.025605   1.096937  -0.023
#> delta:primingseven-fold:elicitingimprecise -0.008725   1.274826  -0.007
#> lambda:(Intercept)                          0.165563   0.026202   6.319
#> lambda:primingseven-fold                    0.062143   0.216997   0.286
#>                                            Pr(>|z|)    
#> gamma:(Intercept)                             0.684    
#> gamma:primingseven-fold                       0.933    
#> gamma:elicitingimprecise                      0.917    
#> gamma:primingseven-fold:elicitingimprecise    0.967    
#> delta:(Intercept)                             0.941    
#> delta:primingseven-fold                       0.987    
#> delta:elicitingimprecise                      0.981    
#> delta:primingseven-fold:elicitingimprecise    0.995    
#> lambda:(Intercept)                         2.64e-10 ***
#> lambda:primingseven-fold                      0.775    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Confidence intervals (95%):
#>                                                 3%    98%
#> gamma:(Intercept)                          -0.3450 0.5256
#> gamma:primingseven-fold                    -0.7990 0.8711
#> gamma:elicitingimprecise                   -0.9696 1.0780
#> gamma:primingseven-fold:elicitingimprecise -1.0871 1.1346
#> delta:(Intercept)                          -1.0883 1.0086
#> delta:primingseven-fold                    -1.2738 1.2522
#> delta:elicitingimprecise                   -2.1756 2.1244
#> delta:primingseven-fold:elicitingimprecise -2.5073 2.4899
#> lambda:(Intercept)                          0.1142 0.2169
#> lambda:primingseven-fold                   -0.3632 0.4874
#> 
#> Link functions:
#> gamma: log
#> delta: logit
#> lambda: log
#> 
#> Fitted parameter means:
#> alpha: 1
#> beta: 1
#> gamma: 1.151
#> delta: 4.851
#> lambda: 1.216
#> 
#> Model fit statistics:
#> Number of observations: 345 
#> Number of parameters: 10 
#> Residual degrees of freedom: 335 
#> Log-likelihood: 175.9 
#> AIC: -331.8 
#> BIC: -293.4 
#> RMSE: 0.1669 
#> Efron's R2: -0.07105 
#> Mean Absolute Error: 0.1242 
#> 
#> Convergence status: Successful 
#> Iterations: 7 
#> 

# Interpretation:
# - Lambda varies by priming: Seven-fold priming may produce more
#   extreme/polarized probability judgments
# }
```
