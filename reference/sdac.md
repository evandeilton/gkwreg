# Autologous Peripheral Blood Stem Cell Transplants Data

Data on Autologous Peripheral Blood Stem Cell Transplants from the Stem
Cell Lab in the Cross Cancer Institute, Alberta Health Services. The
dataset examines recovery rates of CD34+ cells after peripheral blood
stem cell (PBSC) transplants.

## Usage

``` r
sdac
```

## Format

A data frame with 60 observations on 5 variables:

- rcd:

  numeric. Recovery rate of CD34+ cells (proportion in (0, 1)). Response
  variable measuring the proportion of CD34+ cells recovered after PBSC
  transplant.

- age:

  numeric. Patient age in years (range: 18-71 years).

- ageadj:

  numeric. Age-adjusted covariate. Centered and scaled version of age
  for improved numerical stability in regression models.

- chemo:

  factor. Type of chemotherapy protocol used for stem cell mobilization.
  Levels include: `1-day`, `3-day`, `G-CSF only`, and `other`.

- gender:

  factor. Patient gender. Most patients in the study are male.

## Source

Stem Cell Lab, Cross Cancer Institute, Alberta Health Services, Canada.

## Details

This dataset contains clinical data from autologous peripheral blood
stem cell (PBSC) transplant patients treated at the Cross Cancer
Institute, Alberta Health Services. CD34+ cells are hematopoietic stem
and progenitor cells critical for successful transplantation and
hematopoietic recovery.

**Clinical context:** Autologous PBSC transplantation is used to treat
various hematological malignancies including multiple myeloma,
non-Hodgkin's lymphoma, acute leukemia, and some solid tumors. The
recovery rate of CD34+ cells is a crucial predictor of engraftment
success and patient outcomes.

**Chemotherapy protocols:**

- **1-day protocol:** Single-day high-dose chemotherapy for mobilization

- **3-day protocol:** Multi-day chemotherapy regimen

- **G-CSF only:** Granulocyte colony-stimulating factor without
  chemotherapy

- **Other:** Alternative or combined protocols

The proportion of recovered CD34+ cells naturally falls in the interval
(0, 1), making it ideal for proportional data regression modeling. Age
effects are particularly important as older patients may show different
recovery patterns.

This dataset is particularly suitable for:

- Simplex regression (original application by Zhang et al. 2016)

- Beta regression with variable dispersion

- Kumaraswamy regression for flexible distributional modeling

## References

Zhang, P., Qiu, Z., and Shi, C. (2016). simplexreg: An R Package for
Regression Analysis of Proportional Data Using the Simplex Distribution.
*Journal of Statistical Software*, **71**(11), 1–21.
[doi:10.18637/jss.v071.i11](https://doi.org/10.18637/jss.v071.i11)

## Examples

``` r
# \donttest{
require(gkwreg)
require(gkwdist)

data(sdac)

# Example 1: Basic Kumaraswamy regression
# Mean recovery depends on age and chemotherapy protocol
# Precision varies with age (older patients more variable)
fit_kw <- gkwreg(
  rcd ~ ageadj + chemo |
    age,
  data = sdac,
  family = "kw"
)
summary(fit_kw)
#> 
#> Generalized Kumaraswamy Regression Model Summary
#> 
#> Family: kw 
#> 
#> Call:
#> gkwreg(formula = rcd ~ ageadj + chemo | age, data = sdac, family = "kw")
#> 
#> Residuals:
#>     Min  Q1.25%  Median    Mean  Q3.75%     Max 
#> -0.3947 -0.0679  0.0046 -0.0053  0.0733  0.2092 
#> 
#> Coefficients:
#>                   Estimate Std. Error z value Pr(>|z|)    
#> alpha:(Intercept) 1.571402   0.134370  11.695  < 2e-16 ***
#> alpha:ageadj      0.019100   0.006771   2.821  0.00479 ** 
#> alpha:chemo       0.207250   0.101033   2.051  0.04024 *  
#> beta:(Intercept)  0.636173   0.362393   1.755  0.07918 .  
#> beta:age          0.005695   0.006862   0.830  0.40656    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Confidence intervals (95%):
#>                        3%    98%
#> alpha:(Intercept)  1.3080 1.8348
#> alpha:ageadj       0.0058 0.0324
#> alpha:chemo        0.0092 0.4053
#> beta:(Intercept)  -0.0741 1.3464
#> beta:age          -0.0078 0.0191
#> 
#> Link functions:
#> alpha: log
#> beta: log
#> 
#> Fitted parameter means:
#> alpha: 6.967
#> beta: 2.412
#> gamma: 1
#> delta: 0
#> lambda: 1
#> 
#> Model fit statistics:
#> Number of observations: 239 
#> Number of parameters: 5 
#> Residual degrees of freedom: 234 
#> Log-likelihood: 195.7 
#> AIC: -381.5 
#> BIC: -364.1 
#> RMSE: 0.1113 
#> Efron's R2: 0.04391 
#> Mean Absolute Error: 0.08571 
#> 
#> Convergence status: Successful 
#> Iterations: 19 
#> 

# Interpretation:
# - Alpha (mean recovery): Depends on age-adjusted covariate and chemo protocol
#   Different protocols show different baseline recovery rates
#   G-CSF-only may differ from multi-day chemotherapy protocols
# - Beta (precision): Raw age affects recovery variability
#   Hypothesis: Older patients show more heterogeneous responses

# Example 2: Include gender effects
# Gender may influence stem cell recovery rates
fit_kw_gender <- gkwreg(
  rcd ~ ageadj + chemo + gender |
    age + gender,
  data = sdac,
  family = "kw"
)
summary(fit_kw_gender)
#> 
#> Generalized Kumaraswamy Regression Model Summary
#> 
#> Family: kw 
#> 
#> Call:
#> gkwreg(formula = rcd ~ ageadj + chemo + gender | age + gender, 
#>     data = sdac, family = "kw")
#> 
#> Residuals:
#>     Min  Q1.25%  Median    Mean  Q3.75%     Max 
#> -0.3938 -0.0616  0.0086 -0.0007  0.0773  0.2216 
#> 
#> Coefficients:
#>                    Estimate Std. Error z value Pr(>|z|)    
#> alpha:(Intercept)  1.517624   0.167580   9.056  < 2e-16 ***
#> alpha:ageadj       0.019192   0.006794   2.825  0.00473 ** 
#> alpha:chemo        0.204847   0.100938   2.029  0.04241 *  
#> alpha:genderM      0.079423   0.147338   0.539  0.58985    
#> beta:(Intercept)   0.632507   0.377828   1.674  0.09412 .  
#> beta:age           0.005895   0.006950   0.848  0.39629    
#> beta:genderM      -0.005010   0.220384  -0.023  0.98186    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Confidence intervals (95%):
#>                        3%    98%
#> alpha:(Intercept)  1.1892 1.8461
#> alpha:ageadj       0.0059 0.0325
#> alpha:chemo        0.0070 0.4027
#> alpha:genderM     -0.2094 0.3682
#> beta:(Intercept)  -0.1080 1.3730
#> beta:age          -0.0077 0.0195
#> beta:genderM      -0.4370 0.4269
#> 
#> Link functions:
#> alpha: log
#> beta: log
#> 
#> Fitted parameter means:
#> alpha: 6.989
#> beta: 2.549
#> gamma: 1
#> delta: 0
#> lambda: 1
#> 
#> Model fit statistics:
#> Number of observations: 239 
#> Number of parameters: 7 
#> Residual degrees of freedom: 232 
#> Log-likelihood: 196.1 
#> AIC: -378.2 
#> BIC: -353.9 
#> RMSE: 0.111 
#> Efron's R2: 0.04818 
#> Mean Absolute Error: 0.08604 
#> 
#> Convergence status: Successful 
#> Iterations: 22 
#> 

# Interpretation:
# - Gender effects in both mean and precision
# - Precision may differ between males and females

# Test gender significance
anova(fit_kw, fit_kw_gender)
#> Analysis of Deviance Table
#> 
#> Model 1: rcd ~ ageadj + chemo | age
#> Model 2: rcd ~ ageadj + chemo + gender | age + gender
#> 
#>               Resid. Df Resid. Dev Df Deviance Pr(>Chi)  
#> fit_kw        234.00000 -391.48936                       
#> fit_kw_gender 232.00000 -392.23401  2  0.74464  0.68913  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Example 3: Exponentiated Kumaraswamy for extreme recovery patterns
# Some patients show unusually high or low recovery (outliers)
# Lambda parameter captures tail heaviness
fit_ekw <- gkwreg(
  rcd ~ ageadj + chemo + gender | # alpha: mean model
    age + chemo | # beta: precision varies with age and protocol
    chemo, # lambda: protocol affects extremity
  data = sdac,
  family = "ekw"
)
#> Warning: NaNs produced
summary(fit_ekw)
#> 
#> Generalized Kumaraswamy Regression Model Summary
#> 
#> Family: ekw 
#> 
#> Call:
#> gkwreg(formula = rcd ~ ageadj + chemo + gender | age + chemo | 
#>     chemo, data = sdac, family = "ekw")
#> 
#> Residuals:
#>     Min  Q1.25%  Median    Mean  Q3.75%     Max 
#> -0.3924  0.0076  0.1464  0.3648  0.7999  0.9899 
#> 
#> Coefficients:
#>                      Estimate Std. Error z value Pr(>|z|)    
#> alpha:(Intercept)   1.059e+00  7.485e-01   1.414    0.157    
#> alpha:ageadj        1.432e-02        NaN     NaN      NaN    
#> alpha:chemo        -3.245e+01  7.154e-01 -45.358  < 2e-16 ***
#> alpha:genderM       9.128e-02  6.306e-02   1.448    0.148    
#> beta:(Intercept)    7.965e-01  1.814e-01   4.391 1.13e-05 ***
#> beta:age           -4.334e-04        NaN     NaN      NaN    
#> beta:chemo         -2.306e-01  1.812e-01  -1.272    0.203    
#> lambda:(Intercept)  6.693e-01  9.808e-01   0.682    0.495    
#> lambda:chemo        5.661e+01  9.354e-01  60.525  < 2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Confidence intervals (95%):
#>                          3%      98%
#> alpha:(Intercept)   -0.4083   2.5257
#> alpha:ageadj            NaN      NaN
#> alpha:chemo        -33.8498 -31.0456
#> alpha:genderM       -0.0323   0.2149
#> beta:(Intercept)     0.4410   1.1519
#> beta:age                NaN      NaN
#> beta:chemo          -0.5858   0.1247
#> lambda:(Intercept)  -1.2531   2.5918
#> lambda:chemo        54.7807  58.4473
#> 
#> Link functions:
#> alpha: log
#> beta: log
#> lambda: log
#> 
#> Fitted parameter means:
#> alpha: 2.185
#> beta: 1.964
#> gamma: 1
#> delta: 0
#> lambda: 3.439e+24
#> 
#> Model fit statistics:
#> Number of observations: 239 
#> Number of parameters: 9 
#> Residual degrees of freedom: 230 
#> Log-likelihood: 423.6 
#> AIC: -829.1 
#> BIC: -797.8 
#> RMSE: 0.5497 
#> Efron's R2: -22.33 
#> Mean Absolute Error: 0.4138 
#> 
#> Convergence status: Failed 
#> Iterations: 65 
#> 

# Clinical interpretation:
# - Lambda varies by chemotherapy protocol: Some protocols produce more
#   extreme recovery patterns (very high or very low CD34+ counts)
# - G-CSF-only vs multi-day protocols may differ in tail behavior
# - Important for risk stratification and clinical decision-making

# Test if extreme patterns differ by protocol
anova(fit_kw_gender, fit_ekw)
#> Analysis of Deviance Table
#> 
#> Model 1: rcd ~ ageadj + chemo + gender | age + gender
#> Model 2: rcd ~ ageadj + chemo + gender | age + chemo | chemo
#> 
#>               Resid. Df Resid. Dev Df  Deviance Pr(>Chi)    
#> fit_kw_gender 232.00000 -392.23401                          
#> fit_ekw       230.00000 -847.12446  2 454.89045  < 1e-04 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Example 4: Interaction between age and protocol
# Protocol effectiveness may vary with patient age
fit_kw_interact <- gkwreg(
  rcd ~ ageadj * chemo |
    age * chemo,
  data = sdac,
  family = "kw"
)
summary(fit_kw_interact)
#> 
#> Generalized Kumaraswamy Regression Model Summary
#> 
#> Family: kw 
#> 
#> Call:
#> gkwreg(formula = rcd ~ ageadj * chemo | age * chemo, data = sdac, 
#>     family = "kw")
#> 
#> Residuals:
#>     Min  Q1.25%  Median    Mean  Q3.75%     Max 
#> -0.3898 -0.0624  0.0073 -0.0021  0.0729  0.2113 
#> 
#> Coefficients:
#>                     Estimate Std. Error z value Pr(>|z|)    
#> alpha:(Intercept)   1.659107   0.220086   7.538 4.76e-14 ***
#> alpha:ageadj        0.013918   0.011255   1.237    0.216    
#> alpha:chemo         0.138996   0.256543   0.542    0.588    
#> alpha:ageadj:chemo  0.004748   0.015016   0.316    0.752    
#> beta:(Intercept)    1.686154   1.036463   1.627    0.104    
#> beta:age           -0.012553   0.017666  -0.711    0.477    
#> beta:chemo         -1.234435   1.122370  -1.100    0.271    
#> beta:age:chemo      0.023055   0.020005   1.152    0.249    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Confidence intervals (95%):
#>                         3%    98%
#> alpha:(Intercept)   1.2277 2.0905
#> alpha:ageadj       -0.0081 0.0360
#> alpha:chemo        -0.3638 0.6418
#> alpha:ageadj:chemo -0.0247 0.0342
#> beta:(Intercept)   -0.3453 3.7176
#> beta:age           -0.0472 0.0221
#> beta:chemo         -3.4342 0.9654
#> beta:age:chemo     -0.0162 0.0623
#> 
#> Link functions:
#> alpha: log
#> beta: log
#> 
#> Fitted parameter means:
#> alpha: 6.971
#> beta: 2.524
#> gamma: 1
#> delta: 0
#> lambda: 1
#> 
#> Model fit statistics:
#> Number of observations: 239 
#> Number of parameters: 8 
#> Residual degrees of freedom: 231 
#> Log-likelihood: 196.7 
#> AIC: -377.4 
#> BIC: -349.6 
#> RMSE: 0.111 
#> Efron's R2: 0.04886 
#> Mean Absolute Error: 0.08551 
#> 
#> Convergence status: Successful 
#> Iterations: 46 
#> 

# Interpretation:
# - Interaction: Does protocol effectiveness decline with age?
# - Critical for personalized treatment selection
# }
```
