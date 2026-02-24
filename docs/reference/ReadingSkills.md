# Dyslexia and IQ Predicting Reading Accuracy

Data for assessing the contribution of non-verbal IQ to children's
reading skills in dyslexic and non-dyslexic children. This is a classic
dataset demonstrating beta regression with interaction effects and
heteroscedasticity.

## Usage

``` r
ReadingSkills
```

## Format

A data frame with 44 observations on 4 variables:

- accuracy:

  numeric. Reading accuracy score scaled to the open unit interval (0,
  1). Perfect scores of 1 were replaced with 0.99.

- accuracy1:

  numeric. Unrestricted reading accuracy score in (0, 1), including
  boundary observations.

- dyslexia:

  factor. Is the child dyslexic? Levels: `no` (control group) and `yes`
  (dyslexic group). Sum contrast coding is employed.

- iq:

  numeric. Non-verbal intelligence quotient transformed to z-scores
  (mean = 0, SD = 1).

## Source

Data collected by Pammer and Kevan (2004).

## Details

The data were collected by Pammer and Kevan (2004) and employed by
Smithson and Verkuilen (2006) in their seminal beta regression paper.
The sample includes 19 dyslexic children and 25 controls recruited from
primary schools in the Australian Capital Territory. Children's ages
ranged from 8 years 5 months to 12 years 3 months.

Mean reading accuracy was 0.606 for dyslexic readers and 0.900 for
controls. The study investigates whether dyslexia contributes to reading
accuracy even when controlling for IQ (which is on average lower for
dyslexics).

**Transformation details:** The original reading accuracy score was
transformed by Smithson and Verkuilen (2006) to fit beta regression
requirements:

1.  First, the original accuracy was scaled using the minimal and
    maximal scores (a and b) that can be obtained in the test:
    `accuracy1 = (original - a)/(b - a)` (a and b values are not
    provided).

2.  Subsequently, `accuracy` was obtained from `accuracy1` by replacing
    all observations with a value of 1 with 0.99 to fit the open
    interval (0, 1).

The data clearly show asymmetry and heteroscedasticity (especially in
the control group), making beta regression more appropriate than
standard linear regression.

## References

Cribari-Neto, F., and Zeileis, A. (2010). Beta Regression in R. *Journal
of Statistical Software*, **34**(2), 1–24.
[doi:10.18637/jss.v034.i02](https://doi.org/10.18637/jss.v034.i02)

Gr\<U+00FC\>n, B., Kosmidis, I., and Zeileis, A. (2012). Extended Beta
Regression in R: Shaken, Stirred, Mixed, and Partitioned. *Journal of
Statistical Software*, **48**(11), 1–25.
[doi:10.18637/jss.v048.i11](https://doi.org/10.18637/jss.v048.i11)

Kosmidis, I., and Zeileis, A. (2024). Extended-Support Beta Regression
for (0, 1) Responses. *arXiv:2409.07233*.
[doi:10.48550/arXiv.2409.07233](https://doi.org/10.48550/arXiv.2409.07233)

Pammer, K., and Kevan, A. (2004). *The Contribution of Visual
Sensitivity, Phonological Processing and Nonverbal IQ to Children's
Reading*. Unpublished manuscript, The Australian National University,
Canberra.

Smithson, M., and Verkuilen, J. (2006). A Better Lemon Squeezer?
Maximum-Likelihood Regression with Beta-Distributed Dependent Variables.
*Psychological Methods*, **11**(1), 54–71.

## Examples

``` r
# \donttest{
require(gkwreg)
require(gkwdist)

data(ReadingSkills)

# Example 1: Standard Kumaraswamy with interaction and heteroscedasticity
# Mean: Dyslexia <U+00D7> IQ interaction (do groups differ in IQ effect?)
# Precision: Main effects (variability differs by group and IQ level)
fit_kw <- gkwreg(
  accuracy ~ dyslexia * iq |
    dyslexia + iq,
  data = ReadingSkills,
  family = "kw",
  control = gkw_control(method = "L-BFGS-B", maxit = 2000)
)
summary(fit_kw)
#> 
#> Generalized Kumaraswamy Regression Model Summary
#> 
#> Family: kw 
#> 
#> Call:
#> gkwreg(formula = accuracy ~ dyslexia * iq | dyslexia + iq, data = ReadingSkills, 
#>     family = "kw", control = gkw_control(method = "L-BFGS-B", 
#>         maxit = 2000))
#> 
#> Residuals:
#>     Min  Q1.25%  Median    Mean  Q3.75%     Max 
#> -0.2641 -0.0305  0.0093  0.0159  0.0601  0.2986 
#> 
#> Coefficients:
#>                   Estimate Std. Error z value Pr(>|z|)    
#> alpha:(Intercept)   1.7158     0.2387   7.187 6.61e-13 ***
#> alpha:dyslexia      0.8632     0.2522   3.422 0.000621 ***
#> alpha:iq            1.2727     0.2244   5.671 1.42e-08 ***
#> alpha:dyslexia:iq  -1.1788     0.1920  -6.139 8.29e-10 ***
#> beta:(Intercept)    2.8774     0.5252   5.479 4.29e-08 ***
#> beta:dyslexia       3.5541     0.5621   6.323 2.56e-10 ***
#> beta:iq             1.0895     0.3788   2.876 0.004029 ** 
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Confidence intervals (95%):
#>                        3%     98%
#> alpha:(Intercept)  1.2479  2.1837
#> alpha:dyslexia     0.3689  1.3576
#> alpha:iq           0.8329  1.7125
#> alpha:dyslexia:iq -1.5551 -0.8024
#> beta:(Intercept)   1.8480  3.9068
#> beta:dyslexia      2.4525  4.6558
#> beta:iq            0.3470  1.8320
#> 
#> Link functions:
#> alpha: log
#> beta: log
#> 
#> Fitted parameter means:
#> alpha: 26.55
#> beta: 204.8
#> gamma: 1
#> delta: 0
#> lambda: 1
#> 
#> Model fit statistics:
#> Number of observations: 44 
#> Number of parameters: 7 
#> Residual degrees of freedom: 37 
#> Log-likelihood: 68.45 
#> AIC: -122.9 
#> BIC: -110.4 
#> RMSE: 0.1094 
#> Efron's R2: 0.6179 
#> Mean Absolute Error: 0.0753 
#> 
#> Convergence status: Successful 
#> Iterations: 82 
#> 

# Interpretation:
# - Alpha (mean): Interaction shows dyslexic children benefit less from
#   higher IQ compared to controls
# - Beta (precision): Controls show more variable accuracy (higher precision)
#   IQ increases consistency of performance

# Example 2: Simpler model without interaction
fit_kw_simple <- gkwreg(
  accuracy ~ dyslexia + iq |
    dyslexia + iq,
  data = ReadingSkills,
  family = "kw",
  control = gkw_control(method = "L-BFGS-B", maxit = 2000)
)

# Test if interaction is significant
anova(fit_kw_simple, fit_kw)
#> Analysis of Deviance Table
#> 
#> Model 1: accuracy ~ dyslexia + iq | dyslexia + iq
#> Model 2: accuracy ~ dyslexia * iq | dyslexia + iq
#> 
#>               Resid. Df Resid. Dev Df Deviance   Pr(>Chi)    
#> fit_kw_simple  38.00000 -124.28762                           
#> fit_kw         37.00000 -136.89181  1 12.60419 0.00038488 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Example 3: Exponentiated Kumaraswamy for ceiling effects
# Reading accuracy often shows ceiling effects (many perfect/near-perfect scores)
# Lambda parameter can model this right-skewed asymmetry
fit_ekw <- gkwreg(
  accuracy ~ dyslexia * iq | # alpha
    dyslexia + iq | # beta
    dyslexia, # lambda: ceiling effect by group
  data = ReadingSkills,
  family = "ekw",
  control = gkw_control(method = "L-BFGS-B", maxit = 2000)
)
summary(fit_ekw)
#> 
#> Generalized Kumaraswamy Regression Model Summary
#> 
#> Family: ekw 
#> 
#> Call:
#> gkwreg(formula = accuracy ~ dyslexia * iq | dyslexia + iq | dyslexia, 
#>     data = ReadingSkills, family = "ekw", control = gkw_control(method = "L-BFGS-B", 
#>         maxit = 2000))
#> 
#> Residuals:
#>     Min  Q1.25%  Median    Mean  Q3.75%     Max 
#> -0.2447 -0.0243  0.0120  0.0152  0.0686  0.2096 
#> 
#> Coefficients:
#>                    Estimate Std. Error z value Pr(>|z|)  
#> alpha:(Intercept)   -1.9487     4.6708  -0.417   0.6765  
#> alpha:dyslexia       5.3238     5.0741   1.049   0.2941  
#> alpha:iq             2.3183     1.0622   2.182   0.0291 *
#> alpha:dyslexia:iq   -2.3858     1.1040  -2.161   0.0307 *
#> beta:(Intercept)     6.0559    16.7532   0.361   0.7177  
#> beta:dyslexia        6.3640    16.7322   0.380   0.7037  
#> beta:iq              0.5124     0.2990   1.714   0.0866 .
#> lambda:(Intercept)   2.6537     4.0631   0.653   0.5137  
#> lambda:dyslexia     -3.9242     4.5629  -0.860   0.3898  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Confidence intervals (95%):
#>                          3%     98%
#> alpha:(Intercept)  -11.1033  7.2059
#> alpha:dyslexia      -4.6212 15.2689
#> alpha:iq             0.2364  4.4003
#> alpha:dyslexia:iq   -4.5496 -0.2219
#> beta:(Intercept)   -26.7797 38.8915
#> beta:dyslexia      -26.4305 39.1585
#> beta:iq             -0.0736  1.0984
#> lambda:(Intercept)  -5.3098 10.6172
#> lambda:dyslexia    -12.8673  5.0188
#> 
#> Link functions:
#> alpha: log
#> beta: log
#> lambda: log
#> 
#> Fitted parameter means:
#> alpha: 13.48
#> beta: 84001
#> gamma: 1
#> delta: 0
#> lambda: 409.1
#> 
#> Model fit statistics:
#> Number of observations: 44 
#> Number of parameters: 9 
#> Residual degrees of freedom: 35 
#> Log-likelihood: 70.62 
#> AIC: -123.2 
#> BIC: -107.2 
#> RMSE: 0.09653 
#> Efron's R2: 0.7024 
#> Mean Absolute Error: 0.07239 
#> 
#> Convergence status: Failed 
#> Iterations: 118 
#> 

# Interpretation:
# - Lambda varies by dyslexia status: Controls have stronger ceiling effect
#   (more compression at high accuracy) than dyslexic children

# Test if ceiling effect modeling improves fit
anova(fit_kw, fit_ekw)
#> Analysis of Deviance Table
#> 
#> Model 1: accuracy ~ dyslexia * iq | dyslexia + iq
#> Model 2: accuracy ~ dyslexia * iq | dyslexia + iq | dyslexia
#> 
#>         Resid. Df Resid. Dev Df Deviance Pr(>Chi)  
#> fit_kw   37.00000 -136.89181                       
#> fit_ekw  35.00000 -141.23391  2  4.34210  0.11406  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Example 4: McDonald distribution alternative
# Provides different parameterization for extreme values
fit_mc <- gkwreg(
  accuracy ~ dyslexia * iq | # gamma
    dyslexia + iq | # delta
    dyslexia * iq, # lambda: interaction affects tails
  data = ReadingSkills,
  family = "mc",
  control = gkw_control(method = "L-BFGS-B", maxit = 2000)
)
summary(fit_mc)
#> 
#> Generalized Kumaraswamy Regression Model Summary
#> 
#> Family: mc 
#> 
#> Call:
#> gkwreg(formula = accuracy ~ dyslexia * iq | dyslexia + iq | dyslexia * 
#>     iq, data = ReadingSkills, family = "mc", control = gkw_control(method = "L-BFGS-B", 
#>     maxit = 2000))
#> 
#> Residuals:
#>     Min  Q1.25%  Median    Mean  Q3.75%     Max 
#> -0.1259  0.0350  0.7058  0.5210  0.9899  0.9899 
#> 
#> Coefficients:
#>                    Estimate Std. Error z value Pr(>|z|)
#> gamma:(Intercept)    7.8675        NaN     NaN      NaN
#> gamma:dyslexia       9.1494        NaN     NaN      NaN
#> gamma:iq             0.9036        NaN     NaN      NaN
#> gamma:dyslexia:iq   -2.9325        NaN     NaN      NaN
#> delta:(Intercept)  -23.0435        NaN     NaN      NaN
#> delta:dyslexia      63.6331        NaN     NaN      NaN
#> delta:iq            -1.2728        NaN     NaN      NaN
#> lambda:(Intercept)  -5.4452        NaN     NaN      NaN
#> lambda:dyslexia     -8.5214        NaN     NaN      NaN
#> lambda:iq           -0.3551        NaN     NaN      NaN
#> lambda:dyslexia:iq   2.3403        NaN     NaN      NaN
#> 
#> Confidence intervals (95%):
#>                     3% 98%
#> gamma:(Intercept)  NaN NaN
#> gamma:dyslexia     NaN NaN
#> gamma:iq           NaN NaN
#> gamma:dyslexia:iq  NaN NaN
#> delta:(Intercept)  NaN NaN
#> delta:dyslexia     NaN NaN
#> delta:iq           NaN NaN
#> lambda:(Intercept) NaN NaN
#> lambda:dyslexia    NaN NaN
#> lambda:iq          NaN NaN
#> lambda:dyslexia:iq NaN NaN
#> 
#> Link functions:
#> gamma: log
#> delta: logit
#> lambda: log
#> 
#> Fitted parameter means:
#> alpha: 1
#> beta: 1
#> gamma: 104076816
#> delta: 4.318
#> lambda: 19.5
#> 
#> Model fit statistics:
#> Number of observations: 44 
#> Number of parameters: 11 
#> Residual degrees of freedom: 33 
#> Log-likelihood: 60.78 
#> AIC: -99.55 
#> BIC: -79.93 
#> RMSE: 0.6856 
#> Efron's R2: -14.01 
#> Mean Absolute Error: 0.5337 
#> 
#> Convergence status: Failed 
#> Iterations: 142 
#> 

# Compare 3-parameter models
AIC(fit_ekw, fit_mc)
#>         df        AIC
#> fit_ekw  9 -123.23391
#> fit_mc  11  -99.55191
# }
```
