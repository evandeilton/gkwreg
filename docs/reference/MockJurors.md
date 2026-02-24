# Confidence of Mock Jurors in Their Verdicts

Data from a study examining factors that influence mock juror confidence
in verdicts for criminal trials. The experiment manipulates verdict
options (two-option vs. three-option) and presence of conflicting
testimonial evidence.

## Usage

``` r
MockJurors
```

## Format

A data frame with 104 observations on 3 variables:

- confidence:

  numeric. Juror confidence in their verdict, scaled to the open unit
  interval (0, 1). Original scale was 0-100.

- verdict:

  factor indicating whether a two-option verdict (guilty vs. acquittal)
  or three-option verdict (with Scottish 'not proven' alternative) was
  requested. Sum contrast coding is employed.

- conflict:

  factor. Is there conflicting testimonial evidence? Values are `no` or
  `yes`. Sum contrast coding is employed.

## Source

Data collected by Deady (2004), analyzed by Smithson and Verkuilen
(2006).

## Details

The data were collected by Deady (2004) among first-year psychology
students at Australian National University. The experiment examined how
the availability of a third verdict option ('not proven') and
conflicting evidence affect juror confidence.

Smithson and Verkuilen (2006) employed the data, scaling the original
confidence (on a scale 0-100) to the open unit interval using the
transformation: `((original_confidence/100) * 103 - 0.5) / 104`.

**Important note:** The original coding of `conflict` in the data
provided from Smithson's homepage is -1/1 which Smithson and Verkuilen
(2006) describe to mean no/yes. However, all their results (sample
statistics, histograms, etc.) suggest that it actually means yes/no,
which was employed in the corrected `MockJurors` dataset.

## References

Deady, S. (2004). *The Psychological Third Verdict: 'Not Proven' or 'Not
Willing to Make a Decision'?* Unpublished honors thesis, The Australian
National University, Canberra.

Smithson, M., and Verkuilen, J. (2006). A Better Lemon Squeezer?
Maximum-Likelihood Regression with Beta-Distributed Dependent Variables.
*Psychological Methods*, **11**(1), 54–71.

## Examples

``` r
# \donttest{
require(gkwreg)
require(gkwdist)

data(MockJurors)

# Example 1: Main effects model with heteroscedasticity
# Confidence depends on verdict options and conflicting evidence
# Variability may also depend on these factors
fit_kw <- gkwreg(
  confidence ~ verdict + conflict |
    verdict * conflict,
  data = MockJurors,
  family = "kw"
)
summary(fit_kw)
#> 
#> Generalized Kumaraswamy Regression Model Summary
#> 
#> Family: kw 
#> 
#> Call:
#> gkwreg(formula = confidence ~ verdict + conflict | verdict * 
#>     conflict, data = MockJurors, family = "kw")
#> 
#> Residuals:
#>     Min  Q1.25%  Median    Mean  Q3.75%     Max 
#> -0.6484 -0.1143  0.0515  0.0152  0.1827  0.3416 
#> 
#> Coefficients:
#>                       Estimate Std. Error z value Pr(>|z|)    
#> alpha:(Intercept)      0.70572    0.13722   5.143 2.71e-07 ***
#> alpha:verdict         -0.31540    0.15835  -1.992  0.04640 *  
#> alpha:conflict         0.21806    0.15492   1.408  0.15925    
#> beta:(Intercept)      -0.14787    0.12674  -1.167  0.24332    
#> beta:verdict          -0.35402    0.13337  -2.654  0.00794 ** 
#> beta:conflict          0.10940    0.13351   0.819  0.41258    
#> beta:verdict:conflict -0.12932    0.09992  -1.294  0.19558    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Confidence intervals (95%):
#>                            3%     98%
#> alpha:(Intercept)      0.4368  0.9747
#> alpha:verdict         -0.6258 -0.0050
#> alpha:conflict        -0.0856  0.5217
#> beta:(Intercept)      -0.3963  0.1005
#> beta:verdict          -0.6154 -0.0926
#> beta:conflict         -0.1523  0.3711
#> beta:verdict:conflict -0.3252  0.0665
#> 
#> Link functions:
#> alpha: log
#> beta: log
#> 
#> Fitted parameter means:
#> alpha: 2.184
#> beta: 0.9347
#> gamma: 1
#> delta: 0
#> lambda: 1
#> 
#> Model fit statistics:
#> Number of observations: 104 
#> Number of parameters: 7 
#> Residual degrees of freedom: 97 
#> Log-likelihood: 35.56 
#> AIC: -57.12 
#> BIC: -38.61 
#> RMSE: 0.2084 
#> Efron's R2: 0.0305 
#> Mean Absolute Error: 0.1657 
#> 
#> Convergence status: Successful 
#> Iterations: 18 
#> 

# Interpretation:
# - Alpha (mean): Additive effects of verdict type and conflict
#   Three-option verdicts may reduce confidence
#   Conflicting evidence reduces confidence
# - Beta (precision): Interaction suggests confidence variability
#   depends on combination of verdict options and evidence type

# Example 2: Full interaction in mean model
fit_kw_interact <- gkwreg(
  confidence ~ verdict * conflict |
    verdict * conflict,
  data = MockJurors,
  family = "kw"
)
summary(fit_kw_interact)
#> 
#> Generalized Kumaraswamy Regression Model Summary
#> 
#> Family: kw 
#> 
#> Call:
#> gkwreg(formula = confidence ~ verdict * conflict | verdict * 
#>     conflict, data = MockJurors, family = "kw")
#> 
#> Residuals:
#>     Min  Q1.25%  Median    Mean  Q3.75%     Max 
#> -0.6078 -0.1026  0.0131  0.0076  0.1653  0.3822 
#> 
#> Coefficients:
#>                        Estimate Std. Error z value Pr(>|z|)    
#> alpha:(Intercept)       0.78602    0.13396   5.867 4.43e-09 ***
#> alpha:verdict          -0.31171    0.13396  -2.327  0.01997 *  
#> alpha:conflict          0.30523    0.13396   2.278  0.02270 *  
#> alpha:verdict:conflict  0.41032    0.13396   3.063  0.00219 ** 
#> beta:(Intercept)       -0.09854    0.12599  -0.782  0.43415    
#> beta:verdict           -0.34600    0.12599  -2.746  0.00603 ** 
#> beta:conflict           0.11003    0.12599   0.873  0.38250    
#> beta:verdict:conflict   0.10701    0.12599   0.849  0.39569    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Confidence intervals (95%):
#>                             3%     98%
#> alpha:(Intercept)       0.5235  1.0486
#> alpha:verdict          -0.5743 -0.0491
#> alpha:conflict          0.0427  0.5678
#> alpha:verdict:conflict  0.1478  0.6729
#> beta:(Intercept)       -0.3455  0.1484
#> beta:verdict           -0.5929 -0.0991
#> beta:conflict          -0.1369  0.3570
#> beta:verdict:conflict  -0.1399  0.3539
#> 
#> Link functions:
#> alpha: log
#> beta: log
#> 
#> Fitted parameter means:
#> alpha: 2.549
#> beta: 0.9711
#> gamma: 1
#> delta: 0
#> lambda: 1
#> 
#> Model fit statistics:
#> Number of observations: 104 
#> Number of parameters: 8 
#> Residual degrees of freedom: 96 
#> Log-likelihood: 40.12 
#> AIC: -64.24 
#> BIC: -43.09 
#> RMSE: 0.2078 
#> Efron's R2: 0.03573 
#> Mean Absolute Error: 0.1631 
#> 
#> Convergence status: Successful 
#> Iterations: 18 
#> 

# Interpretation:
# - Full interaction: Third verdict option may have different effects
#   depending on whether evidence is conflicting

# Test interaction significance
anova(fit_kw, fit_kw_interact)
#> Analysis of Deviance Table
#> 
#> Model 1: confidence ~ verdict + conflict | verdict * conflict
#> Model 2: confidence ~ verdict * conflict | verdict * conflict
#> 
#>                 Resid. Df Resid. Dev Df Deviance  Pr(>Chi)   
#> fit_kw           97.00000  -71.12225                         
#> fit_kw_interact  96.00000  -80.24264  1  9.12039 0.0025278 **
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Example 3: McDonald distribution for extreme confidence patterns
# Jurors may show very high confidence (ceiling effects) or very low
# confidence depending on conditions
fit_mc <- gkwreg(
  confidence ~ verdict * conflict | # gamma: full interaction
    verdict * conflict | # delta: full interaction
    verdict + conflict, # lambda: additive extremity effects
  data = MockJurors,
  family = "mc",
  control = gkw_control(
    method = "BFGS",
    maxit = 1500,
    reltol = 1e-8
  )
)
summary(fit_mc)
#> 
#> Generalized Kumaraswamy Regression Model Summary
#> 
#> Family: mc 
#> 
#> Call:
#> gkwreg(formula = confidence ~ verdict * conflict | verdict * 
#>     conflict | verdict + conflict, data = MockJurors, family = "mc", 
#>     control = gkw_control(method = "BFGS", maxit = 1500, reltol = 1e-08))
#> 
#> Residuals:
#>     Min  Q1.25%  Median    Mean  Q3.75%     Max 
#> -0.4395 -0.0471  0.0520  0.1064  0.2500  0.7032 
#> 
#> Coefficients:
#>                        Estimate Std. Error z value Pr(>|z|)    
#> gamma:(Intercept)        3.2186     5.7452   0.560   0.5753    
#> gamma:verdict            2.0199     5.5102   0.367   0.7139    
#> gamma:conflict          -0.7418     0.6489  -1.143   0.2530    
#> gamma:verdict:conflict   0.9145     0.1266   7.224 5.04e-13 ***
#> delta:(Intercept)       -3.6223     0.3943  -9.186  < 2e-16 ***
#> delta:verdict           -0.8930     0.3948  -2.262   0.0237 *  
#> delta:conflict           0.2426     0.3944   0.615   0.5386    
#> delta:verdict:conflict   0.7301     0.3960   1.843   0.0653 .  
#> lambda:(Intercept)      -2.2131     5.7326  -0.386   0.6995    
#> lambda:verdict          -2.3000     5.5063  -0.418   0.6762    
#> lambda:conflict          1.0677     0.6352   1.681   0.0928 .  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Confidence intervals (95%):
#>                              3%     98%
#> gamma:(Intercept)       -8.0418 14.4790
#> gamma:verdict           -8.7798 12.8196
#> gamma:conflict          -2.0136  0.5300
#> gamma:verdict:conflict   0.6664  1.1627
#> delta:(Intercept)       -4.3952 -2.8494
#> delta:verdict           -1.6667 -0.1193
#> delta:conflict          -0.5305  1.0156
#> delta:verdict:conflict  -0.0461  1.5063
#> lambda:(Intercept)     -13.4488  9.0225
#> lambda:verdict         -13.0922  8.4921
#> lambda:conflict         -0.1773  2.3127
#> 
#> Link functions:
#> gamma: log
#> delta: logit
#> lambda: log
#> 
#> Fitted parameter means:
#> alpha: 1
#> beta: 1
#> gamma: 100.7
#> delta: 0.4194
#> lambda: 0.8955
#> 
#> Model fit statistics:
#> Number of observations: 104 
#> Number of parameters: 11 
#> Residual degrees of freedom: 93 
#> Log-likelihood: 19.98 
#> AIC: -17.96 
#> BIC: 11.13 
#> RMSE: 0.2857 
#> Efron's R2: -0.8222 
#> Mean Absolute Error: 0.2204 
#> 
#> Convergence status: Successful 
#> Iterations: 38 
#> 

# Interpretation:
# - Lambda: Models asymmetry and extreme confidence
#   Some conditions produce more polarized confidence (very high or very low)

# Example 4: Exponentiated Kumaraswamy alternative
fit_ekw <- gkwreg(
  confidence ~ verdict * conflict | # alpha
    verdict + conflict | # beta
    conflict, # lambda: conflict affects extremity
  data = MockJurors,
  family = "ekw",
  control = gkw_control(
    method = "BFGS",
    maxit = 1500
  )
)
#> Warning: NaNs produced
summary(fit_ekw)
#> 
#> Generalized Kumaraswamy Regression Model Summary
#> 
#> Family: ekw 
#> 
#> Call:
#> gkwreg(formula = confidence ~ verdict * conflict | verdict + 
#>     conflict | conflict, data = MockJurors, family = "ekw", control = gkw_control(method = "BFGS", 
#>     maxit = 1500))
#> 
#> Residuals:
#>     Min  Q1.25%  Median    Mean  Q3.75%     Max 
#> -0.7253 -0.2066 -0.0576 -0.0722  0.0744  0.2647 
#> 
#> Coefficients:
#>                        Estimate Std. Error z value Pr(>|z|)
#> alpha:(Intercept)       0.45685        NaN     NaN      NaN
#> alpha:verdict          -0.05041    0.12348  -0.408    0.683
#> alpha:conflict          0.09766        NaN     NaN      NaN
#> alpha:verdict:conflict  0.13324    0.10132   1.315    0.189
#> beta:(Intercept)       -0.64393        NaN     NaN      NaN
#> beta:verdict           -0.17246    0.12703  -1.358    0.175
#> beta:conflict          -0.07422        NaN     NaN      NaN
#> lambda:(Intercept)      0.45685        NaN     NaN      NaN
#> lambda:conflict         0.09766        NaN     NaN      NaN
#> 
#> Confidence intervals (95%):
#>                             3%    98%
#> alpha:(Intercept)          NaN    NaN
#> alpha:verdict          -0.2924 0.1916
#> alpha:conflict             NaN    NaN
#> alpha:verdict:conflict -0.0654 0.3318
#> beta:(Intercept)           NaN    NaN
#> beta:verdict           -0.4214 0.0765
#> beta:conflict              NaN    NaN
#> lambda:(Intercept)         NaN    NaN
#> lambda:conflict            NaN    NaN
#> 
#> Link functions:
#> alpha: log
#> beta: log
#> lambda: log
#> 
#> Fitted parameter means:
#> alpha: 1.608
#> beta: 0.5339
#> gamma: 1
#> delta: 0
#> lambda: 1.589
#> 
#> Model fit statistics:
#> Number of observations: 104 
#> Number of parameters: 9 
#> Residual degrees of freedom: 95 
#> Log-likelihood: 9.75 
#> AIC: -1.499 
#> BIC: 22.3 
#> RMSE: 0.2213 
#> Efron's R2: -0.09353 
#> Mean Absolute Error: 0.1677 
#> 
#> Convergence status: Successful 
#> Iterations: 5 
#> 

# Compare 3-parameter models
AIC(fit_ekw, fit_mc)
#>         df        AIC
#> fit_ekw  9  -1.499307
#> fit_mc  11 -17.957086
# }
```
