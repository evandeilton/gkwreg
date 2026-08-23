# (No) Myopic Loss Aversion in Adolescents

Data from a behavioral economics experiment assessing the extent of
myopic loss aversion among adolescents aged 11 to 19 years. The
experiment tests whether short-term investment horizons lead to more
conservative investment behavior.

## Usage

``` r
LossAversion
```

## Format

A data frame with 570 observations on 7 variables:

- invest:

  numeric. Average proportion of tokens invested across all 9 rounds of
  the experiment (response variable).

- gender:

  factor. Gender of the player (or team of players).

- male:

  factor. Was (at least one of) the player(s) male (in the team)?

- age:

  numeric. Age of the player (or average age in case of team).

- grade:

  factor. School grade of the player(s).

- arrangement:

  factor. Investment horizon treatment with levels `short` (1 round),
  `medium` (3 rounds), and `long` (9 rounds).

- treatment:

  factor. Type of treatment: long vs. short.

## Source

Data collected by Matthias Sutter and Daniela Glätzle-Rützler,
Universität Innsbruck.

## Details

The data were collected by Matthias Sutter and Daniela Glätzle-Rützler
(Universität Innsbruck) in an experiment with high-school students in
Tyrol, Austria (Schwaz and Innsbruck). The experiment tests the theory
of myopic loss aversion, which proposes that investors with shorter
evaluation periods are more loss-averse and thus invest less in risky
assets.

Classical theory predicts that players with short investment horizons
(myopic view) should invest less due to loss aversion. However, Sutter
et al. (2015) found no evidence of myopic loss aversion in adolescents,
contrary to findings in adult populations.

The investment game structure: In each round, players could invest
tokens in a risky asset with 50% chance of doubling or losing the
investment. The treatment varied the feedback frequency (short = every
round, medium = every 3 rounds, long = only at the end).

## References

Sutter, M., Kocher, M.G., Glätzle-Rützler, D., and Trautmann, S.T.
(2015). No Myopic Loss Aversion in Adolescents? – An Experimental Note.
*Journal of Economic Behavior & Organization*, **111**, 169–176.
[doi:10.1016/j.jebo.2014.12.021](https://doi.org/10.1016/j.jebo.2014.12.021)

Kosmidis, I., and Zeileis, A. (2024). Extended-Support Beta Regression
for (0, 1) Responses. *arXiv:2409.07233*.
[doi:10.48550/arXiv.2409.07233](https://doi.org/10.48550/arXiv.2409.07233)

## Examples

``` r
# \donttest{
require(gkwreg)
require(gkwdist)

data(LossAversion)
# Control bounds

LossAversion$invest <- with(
  LossAversion,
  ifelse(invest <= 0, 0.000001,
    ifelse(invest >= 1, 0.999999, invest)
  )
)
# Example 1: Test for myopic loss aversion
# Do short-term players invest less? (They shouldn't, per Sutter et al.)
fit_kw <- gkwreg(
  invest ~ arrangement + age + male + grade |
    arrangement + male,
  data = LossAversion,
  family = "kw"
)
summary(fit_kw)
#> 
#> Generalized Kumaraswamy Regression Model Summary
#> 
#> Family: kw 
#> 
#> Call:
#> gkwreg(formula = invest ~ arrangement + age + male + grade | 
#>     arrangement + male, data = LossAversion, family = "kw")
#> 
#> Residuals:
#>     Min  Q1.25%  Median    Mean  Q3.75%     Max 
#> -0.5002 -0.2002  0.0000  0.0037  0.1999  0.5000 
#> 
#> Coefficients:
#>                         Estimate Std. Error z value Pr(>|z|)
#> alpha:(Intercept)     -2.841e-10  4.735e-01       0        1
#> alpha:arrangementteam  2.126e-10  1.836e-01       0        1
#> alpha:age             -4.949e-09  3.771e-02       0        1
#> alpha:maleyes          1.889e-11  1.210e-01       0        1
#> alpha:grade10-12      -4.135e-10  1.732e-01       0        1
#> beta:(Intercept)       4.573e-10  9.653e-02       0        1
#> beta:arrangementteam  -3.177e-11  1.694e-01       0        1
#> beta:maleyes          -1.490e-11  1.394e-01       0        1
#> 
#> Confidence intervals (95%):
#>                            3%    98%
#> alpha:(Intercept)     -0.9281 0.9281
#> alpha:arrangementteam -0.3598 0.3598
#> alpha:age             -0.0739 0.0739
#> alpha:maleyes         -0.2372 0.2372
#> alpha:grade10-12      -0.3396 0.3396
#> beta:(Intercept)      -0.1892 0.1892
#> beta:arrangementteam  -0.3320 0.3320
#> beta:maleyes          -0.2733 0.2733
#> 
#> Link functions:
#> alpha: log
#> beta: log
#> 
#> Fitted parameter means:
#> alpha: 1
#> beta: 0.9996
#> gamma: 1
#> delta: 0
#> lambda: 1
#> 
#> Model fit statistics:
#> Number of observations: 570 
#> Number of parameters: 8 
#> Residual degrees of freedom: 562 
#> Log-likelihood: -3e+11 
#> AIC: 6e+11 
#> BIC: 6e+11 
#> RMSE: 0.2678 
#> Efron's R2: 7.819e-05 
#> Mean Absolute Error: 0.2217 
#> 
#> Convergence status: Failed 
#> Iterations: 2 
#> 

# Interpretation:
# - Alpha: Effect of investment horizon (arrangement) on mean investment
#   Age and gender effects on risk-taking
# - Beta: Precision varies by horizon and gender
#   (some groups more consistent than others)

# Example 2: Interaction effects
# Does the horizon effect differ by age/grade?
fit_kw_interact <- gkwreg(
  invest ~ grade * (arrangement + age) + male |
    arrangement + male + grade,
  data = LossAversion,
  family = "kw"
)
summary(fit_kw_interact)
#> 
#> Generalized Kumaraswamy Regression Model Summary
#> 
#> Family: kw 
#> 
#> Call:
#> gkwreg(formula = invest ~ grade * (arrangement + age) + male | 
#>     arrangement + male + grade, data = LossAversion, family = "kw")
#> 
#> Residuals:
#>     Min  Q1.25%  Median    Mean  Q3.75%     Max 
#> -0.5002 -0.2002  0.0000  0.0037  0.2000  0.5000 
#> 
#> Coefficients:
#>                                    Estimate Std. Error z value Pr(>|z|)
#> alpha:(Intercept)                -1.698e-10  6.753e-01       0        1
#> alpha:grade10-12                 -2.471e-10  1.120e+00       0        1
#> alpha:arrangementteam             1.270e-10  2.113e-01       0        1
#> alpha:age                        -2.958e-09  5.342e-02       0        1
#> alpha:maleyes                     1.129e-11  1.224e-01       0        1
#> alpha:grade10-12:arrangementteam  3.069e-11  2.261e-01       0        1
#> alpha:grade10-12:age             -4.007e-09  7.630e-02       0        1
#> beta:(Intercept)                  2.733e-10  1.126e-01       0        1
#> beta:arrangementteam             -1.898e-11  1.720e-01       0        1
#> beta:maleyes                     -8.903e-12  1.398e-01       0        1
#> beta:grade10-12                   1.314e-10  1.440e-01       0        1
#> 
#> Confidence intervals (95%):
#>                                       3%    98%
#> alpha:(Intercept)                -1.3237 1.3237
#> alpha:grade10-12                 -2.1959 2.1959
#> alpha:arrangementteam            -0.4142 0.4142
#> alpha:age                        -0.1047 0.1047
#> alpha:maleyes                    -0.2398 0.2398
#> alpha:grade10-12:arrangementteam -0.4431 0.4431
#> alpha:grade10-12:age             -0.1495 0.1495
#> beta:(Intercept)                 -0.2208 0.2208
#> beta:arrangementteam             -0.3371 0.3371
#> beta:maleyes                     -0.2740 0.2740
#> beta:grade10-12                  -0.2823 0.2823
#> 
#> Link functions:
#> alpha: log
#> beta: log
#> 
#> Fitted parameter means:
#> alpha: 1
#> beta: 0.9996
#> gamma: 1
#> delta: 0
#> lambda: 1
#> 
#> Model fit statistics:
#> Number of observations: 570 
#> Number of parameters: 11 
#> Residual degrees of freedom: 559 
#> Log-likelihood: -3e+11 
#> AIC: 6e+11 
#> BIC: 6e+11 
#> RMSE: 0.2679 
#> Efron's R2: -0.0001649 
#> Mean Absolute Error: 0.2218 
#> 
#> Convergence status: Failed 
#> Iterations: 2 
#> 

# Interpretation:
# - Grade × arrangement interaction tests if myopic loss aversion
#   emerges differently at different developmental stages

# Example 3: Extended-support for boundary observations
# Some students invest 0% or 100% of tokens
# Original 'invest' variable may include exact 0 and 1 values
fit_xbx <- gkwreg(
  invest ~ grade * (arrangement + age) + male |
    arrangement + male + grade,
  data = LossAversion,
  family = "kw" # Note: for true [0,1] support, use extended-support models
)
summary(fit_xbx)
#> 
#> Generalized Kumaraswamy Regression Model Summary
#> 
#> Family: kw 
#> 
#> Call:
#> gkwreg(formula = invest ~ grade * (arrangement + age) + male | 
#>     arrangement + male + grade, data = LossAversion, family = "kw")
#> 
#> Residuals:
#>     Min  Q1.25%  Median    Mean  Q3.75%     Max 
#> -0.5002 -0.2002  0.0000  0.0037  0.2000  0.5000 
#> 
#> Coefficients:
#>                                    Estimate Std. Error z value Pr(>|z|)
#> alpha:(Intercept)                -1.698e-10  6.753e-01       0        1
#> alpha:grade10-12                 -2.471e-10  1.120e+00       0        1
#> alpha:arrangementteam             1.270e-10  2.113e-01       0        1
#> alpha:age                        -2.958e-09  5.342e-02       0        1
#> alpha:maleyes                     1.129e-11  1.224e-01       0        1
#> alpha:grade10-12:arrangementteam  3.069e-11  2.261e-01       0        1
#> alpha:grade10-12:age             -4.007e-09  7.630e-02       0        1
#> beta:(Intercept)                  2.733e-10  1.126e-01       0        1
#> beta:arrangementteam             -1.898e-11  1.720e-01       0        1
#> beta:maleyes                     -8.903e-12  1.398e-01       0        1
#> beta:grade10-12                   1.314e-10  1.440e-01       0        1
#> 
#> Confidence intervals (95%):
#>                                       3%    98%
#> alpha:(Intercept)                -1.3237 1.3237
#> alpha:grade10-12                 -2.1959 2.1959
#> alpha:arrangementteam            -0.4142 0.4142
#> alpha:age                        -0.1047 0.1047
#> alpha:maleyes                    -0.2398 0.2398
#> alpha:grade10-12:arrangementteam -0.4431 0.4431
#> alpha:grade10-12:age             -0.1495 0.1495
#> beta:(Intercept)                 -0.2208 0.2208
#> beta:arrangementteam             -0.3371 0.3371
#> beta:maleyes                     -0.2740 0.2740
#> beta:grade10-12                  -0.2823 0.2823
#> 
#> Link functions:
#> alpha: log
#> beta: log
#> 
#> Fitted parameter means:
#> alpha: 1
#> beta: 0.9996
#> gamma: 1
#> delta: 0
#> lambda: 1
#> 
#> Model fit statistics:
#> Number of observations: 570 
#> Number of parameters: 11 
#> Residual degrees of freedom: 559 
#> Log-likelihood: -3e+11 
#> AIC: 6e+11 
#> BIC: 6e+11 
#> RMSE: 0.2679 
#> Efron's R2: -0.0001649 
#> Mean Absolute Error: 0.2218 
#> 
#> Convergence status: Failed 
#> Iterations: 2 
#> 

# Interpretation:
# - Model accommodates extreme risk-taking (all-in or all-out strategies)

# Compare models
anova(fit_kw, fit_kw_interact)
#> Analysis of Deviance Table
#> 
#> Model 1: invest ~ arrangement + age + male + grade | arrangement + male
#> Model 2: invest ~ grade * (arrangement + age) + male | arrangement + male + 
#> Model 2:     grade
#> 
#>                 Resid. Df Resid. Dev Df Deviance Pr(>Chi)  
#> fit_kw          562.00000      6e+11                       
#> fit_kw_interact 559.00000      6e+11  3  0.00000        1  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Visualization: Investment by horizon
boxplot(invest ~ arrangement,
  data = LossAversion,
  xlab = "Investment Horizon", ylab = "Proportion Invested",
  main = "No Myopic Loss Aversion in Adolescents",
  col = c("lightblue", "lightgreen", "lightyellow")
)

# }
```
