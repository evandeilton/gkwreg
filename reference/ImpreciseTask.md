# Imprecise Probabilities for Sunday Weather and Boeing Stock Task

Data from a cognitive psychology experiment where participants estimated
upper and lower probabilities for events to occur and not to occur. The
study examines judgment under uncertainty with imprecise probability
assessments.

## Usage

``` r
ImpreciseTask
```

## Format

A data frame with 242 observations on 3 variables:

- task:

  factor with levels `Boeing stock` and `Sunday weather`. Indicates
  which task the participant performed.

- location:

  numeric. Average of the lower estimate for the event not to occur and
  the upper estimate for the event to occur (proportion).

- difference:

  numeric. Difference between upper and lower probability estimates,
  measuring imprecision or uncertainty.

## Source

Taken from Smithson et al. (2011) supplements.

## Details

All participants in the study were either first- or second-year
undergraduate students in psychology at Australian universities, none of
whom had a strong background in probability theory or were familiar with
imprecise probability theories.

For the Sunday weather task, participants were asked to estimate the
probability that the temperature at Canberra airport on Sunday would be
higher than a specified value.

For the Boeing stock task, participants were asked to estimate the
probability that Boeing's stock would rise more than those in a list of
30 companies.

For each task, participants were asked to provide lower and upper
estimates for the event to occur and not to occur.

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

data(ImpreciseTask)

# Example 1: Basic model with task effects
# Probability location varies by task type and uncertainty level
fit_kw <- gkwreg(location ~ task * difference,
  data = ImpreciseTask,
  family = "kw"
)
summary(fit_kw)
#> 
#> Generalized Kumaraswamy Regression Model Summary
#> 
#> Family: kw 
#> 
#> Call:
#> gkwreg(formula = location ~ task * difference, data = ImpreciseTask, 
#>     family = "kw")
#> 
#> Residuals:
#>     Min  Q1.25%  Median    Mean  Q3.75%     Max 
#> -0.3556 -0.0788  0.0274  0.0019  0.0944  0.3499 
#> 
#> Coefficients:
#>                                     Estimate Std. Error z value Pr(>|z|)    
#> alpha:(Intercept)                    1.25134    0.11593  10.794   <2e-16 ***
#> alpha:taskSunday weather            -0.09543    0.11108  -0.859   0.3903    
#> alpha:difference                     0.28623    0.21709   1.318   0.1874    
#> alpha:taskSunday weather:difference  0.43791    0.25964   1.687   0.0917 .  
#> beta:(Intercept)                     2.74022    0.15340  17.864   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Confidence intervals (95%):
#>                                          3%    98%
#> alpha:(Intercept)                    1.0241 1.4786
#> alpha:taskSunday weather            -0.3131 0.1223
#> alpha:difference                    -0.1393 0.7117
#> alpha:taskSunday weather:difference -0.0710 0.9468
#> beta:(Intercept)                     2.4396 3.0409
#> 
#> Link functions:
#> alpha: log
#> beta: log
#> 
#> Fitted parameter means:
#> alpha: 3.805
#> beta: 15.47
#> gamma: 1
#> delta: 0
#> lambda: 1
#> 
#> Model fit statistics:
#> Number of observations: 242 
#> Number of parameters: 5 
#> Residual degrees of freedom: 237 
#> Log-likelihood: 162.2 
#> AIC: -314.3 
#> BIC: -296.9 
#> RMSE: 0.1198 
#> Efron's R2: 0.101 
#> Mean Absolute Error: 0.0981 
#> 
#> Convergence status: Successful 
#> Iterations: 21 
#> 

# Interpretation:
# - Alpha: Task type and uncertainty (difference) interact to affect
#   probability estimates
# - Different tasks may have different baseline probability assessments

# Example 2: Heteroscedastic model
# Precision of estimates may vary by task and uncertainty
fit_kw_hetero <- gkwreg(
  location ~ task * difference |
    task + difference,
  data = ImpreciseTask,
  family = "kw"
)
summary(fit_kw_hetero)
#> 
#> Generalized Kumaraswamy Regression Model Summary
#> 
#> Family: kw 
#> 
#> Call:
#> gkwreg(formula = location ~ task * difference | task + difference, 
#>     data = ImpreciseTask, family = "kw")
#> 
#> Residuals:
#>     Min  Q1.25%  Median    Mean  Q3.75%     Max 
#> -0.3607 -0.0787  0.0333  0.0010  0.0907  0.3574 
#> 
#> Coefficients:
#>                                     Estimate Std. Error z value Pr(>|z|)    
#> alpha:(Intercept)                     1.7403     0.1329  13.092  < 2e-16 ***
#> alpha:taskSunday weather             -0.4913     0.1365  -3.598 0.000320 ***
#> alpha:difference                     -0.1068     0.2145  -0.498 0.618539    
#> alpha:taskSunday weather:difference   0.3001     0.2410   1.245 0.213134    
#> beta:(Intercept)                      4.4688     0.4707   9.493  < 2e-16 ***
#> beta:taskSunday weather              -1.4832     0.4354  -3.407 0.000658 ***
#> beta:difference                      -1.3788     0.5970  -2.309 0.020925 *  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Confidence intervals (95%):
#>                                          3%     98%
#> alpha:(Intercept)                    1.4798  2.0009
#> alpha:taskSunday weather            -0.7589 -0.2237
#> alpha:difference                    -0.5272  0.3136
#> alpha:taskSunday weather:difference -0.1723  0.7725
#> beta:(Intercept)                     3.5462  5.3914
#> beta:taskSunday weather             -2.3366 -0.6298
#> beta:difference                     -2.5489 -0.2086
#> 
#> Link functions:
#> alpha: log
#> beta: log
#> 
#> Fitted parameter means:
#> alpha: 4.164
#> beta: 25.37
#> gamma: 1
#> delta: 0
#> lambda: 1
#> 
#> Model fit statistics:
#> Number of observations: 242 
#> Number of parameters: 7 
#> Residual degrees of freedom: 235 
#> Log-likelihood: 170.6 
#> AIC: -327.3 
#> BIC: -302.9 
#> RMSE: 0.1197 
#> Efron's R2: 0.102 
#> Mean Absolute Error: 0.09725 
#> 
#> Convergence status: Successful 
#> Iterations: 26 
#> 

# Interpretation:
# - Beta: Variability in estimates differs between tasks
#   Higher uncertainty (difference) may lead to less precise estimates

# Example 3: McDonald distribution for extreme uncertainty
# Some participants may show very extreme probability assessments
fit_mc <- gkwreg(
  location ~ task * difference | # gamma: full interaction
    task * difference | # delta: full interaction
    task, # lambda: task affects extremity
  data = ImpreciseTask,
  family = "mc",
  control = gkw_control(
    method = "BFGS",
    maxit = 1500
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
#> gkwreg(formula = location ~ task * difference | task * difference | 
#>     task, data = ImpreciseTask, family = "mc", control = gkw_control(method = "BFGS", 
#>     maxit = 1500))
#> 
#> Residuals:
#>     Min  Q1.25%  Median    Mean  Q3.75%     Max 
#> -0.2631  0.1058  0.2147  0.1669  0.2255  0.6379 
#> 
#> Coefficients:
#>                                     Estimate Std. Error z value Pr(>|z|)
#> gamma:(Intercept)                    0.12000        NaN     NaN      NaN
#> gamma:taskSunday weather             0.08321    0.37727   0.221    0.825
#> gamma:difference                     0.03513        NaN     NaN      NaN
#> gamma:taskSunday weather:difference  0.01876        NaN     NaN      NaN
#> delta:(Intercept)                   -0.08309        NaN     NaN      NaN
#> delta:taskSunday weather            -0.05747    0.58059  -0.099    0.921
#> delta:difference                    -0.02637        NaN     NaN      NaN
#> delta:taskSunday weather:difference -0.01482    1.06142  -0.014    0.989
#> lambda:(Intercept)                   0.25123        NaN     NaN      NaN
#> lambda:taskSunday weather            0.17376        NaN     NaN      NaN
#> 
#> Confidence intervals (95%):
#>                                          3%    98%
#> gamma:(Intercept)                       NaN    NaN
#> gamma:taskSunday weather            -0.6562 0.8226
#> gamma:difference                        NaN    NaN
#> gamma:taskSunday weather:difference     NaN    NaN
#> delta:(Intercept)                       NaN    NaN
#> delta:taskSunday weather            -1.1954 1.0805
#> delta:difference                        NaN    NaN
#> delta:taskSunday weather:difference -2.0952 2.0655
#> lambda:(Intercept)                      NaN    NaN
#> lambda:taskSunday weather               NaN    NaN
#> 
#> Link functions:
#> gamma: log
#> delta: logit
#> lambda: log
#> 
#> Fitted parameter means:
#> alpha: 1
#> beta: 1
#> gamma: 1.212
#> delta: 4.668
#> lambda: 1.457
#> 
#> Model fit statistics:
#> Number of observations: 242 
#> Number of parameters: 10 
#> Residual degrees of freedom: 232 
#> Log-likelihood: 20.49 
#> AIC: -20.97 
#> BIC: 13.92 
#> RMSE: 0.2123 
#> Efron's R2: -1.824 
#> Mean Absolute Error: 0.1885 
#> 
#> Convergence status: Successful 
#> Iterations: 7 
#> 

# Interpretation:
# - Lambda varies by task: Weather vs. stock may produce
#   different patterns of extreme probability assessments
# }
```
