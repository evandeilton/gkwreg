# Update and Re-fit a GKw Regression Model

Updates and (by default) re-fits a Generalized Kumaraswamy regression
model. This method allows modification of the model formula, data, or
other arguments without having to completely re-specify the model call.
Supports formulas with up to 5 parts (alpha, beta, gamma, delta, lambda)
using the Formula package.

## Usage

``` r
# S3 method for class 'gkwreg'
update(object, formula., ..., data. = NULL, evaluate = TRUE)
```

## Arguments

- object:

  An object of class `"gkwreg"`, typically obtained from
  [`gkwreg`](https://evandeilton.github.io/gkwreg/reference/gkwreg.md).

- formula.:

  Changes to the formula. This is a formula where `.` refers to the
  corresponding part of the old formula. For multi-part formulas (e.g.,
  `y ~ x1 | x2 | x3`), you can update each part separately using the `|`
  separator.

- ...:

  Additional arguments to the call, or arguments with changed values.
  Use `name = NULL` to remove an argument.

- data.:

  Optional. A new data frame in which to evaluate the updated model. If
  omitted, the original data is used.

- evaluate:

  Logical. If `TRUE` (default), the updated model is fitted. If `FALSE`,
  the updated call is returned without fitting.

## Value

If `evaluate = TRUE`, a new fitted model object of class `"gkwreg"`. If
`evaluate = FALSE`, an updated call.

## Details

The `update` method allows you to modify a fitted model and re-fit it
with the changes. The GKw regression model supports formulas with up to
5 parts:
`y ~ model_alpha | model_beta | model_gamma | model_delta | model_lambda`

Each part can be updated independently using `.` to refer to the current
specification:

- `. ~ . + x | . | . | . | .` - Add `x` to alpha only

- `. ~ . | . + x | . | . | .` - Add `x` to beta only

- `. ~ . | . | . + x | . | .` - Add `x` to gamma only

- `. ~ . + x | . + x | . | . | .` - Add `x` to alpha and beta

- `. ~ . - x | . | . | . | .` - Remove `x` from alpha

Omitting parts at the end is allowed (they default to `.`):

- `. ~ . + x | .` is equivalent to `. ~ . + x | . | . | . | .`

- `. ~ . | . + x` is equivalent to `. ~ . | . + x | . | . | .`

## See also

[`gkwreg`](https://evandeilton.github.io/gkwreg/reference/gkwreg.md),
[`update`](https://rdrr.io/r/stats/update.html),
[`Formula`](https://rdrr.io/pkg/Formula/man/Formula.html)

## Author

Lopes, J. E.

## Examples

``` r
# \donttest{
# Load example data
require(gkwreg)

data(GasolineYield)

# EXAMPLE 1: Simple formulas (1 part - alpha only)

m1_0 <- gkwreg(yield ~ 1, data = GasolineYield, family = "kw")
m1_1 <- update(m1_0, . ~ . + temp)
m1_2 <- update(m1_1, . ~ . + batch)
#> Warning: NaNs produced
m1_3 <- update(m1_2, . ~ . - temp)

anova(m1_0, m1_1, m1_2)
#> Analysis of Deviance Table
#> 
#> Model 1: yield ~ 1
#> Model 2: yield ~ temp
#> Model 3: yield ~ temp + batch
#> 
#>      Resid. Df Resid. Dev Df  Deviance Pr(>Chi)    
#> m1_0  30.00000  -57.02258                          
#> m1_1  29.00000  -80.90259  1  23.88001  < 1e-04 ***
#> m1_2  20.00000 -193.93742  9 113.03483  < 1e-04 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
AIC(m1_0, m1_1, m1_2, m1_3)
#>      df        AIC
#> m1_0  2  -53.02258
#> m1_1  3  -74.90259
#> m1_2 12 -169.93742
#> m1_3 11  -42.93436
BIC(m1_0, m1_1, m1_2, m1_3)
#>      df        BIC
#> m1_0  2  -50.09111
#> m1_1  3  -70.50538
#> m1_2 12 -152.34859
#> m1_3 11  -26.81126

# EXAMPLE 2: Two-part formulas (alpha | beta)

# Start with intercept-only for both
m2_0 <- gkwreg(yield ~ 1 | 1, data = GasolineYield, family = "kw")

# Add temp to alpha
m2_1 <- update(m2_0, . ~ . + temp | .)

# Add batch to beta
m2_2 <- update(m2_1, . ~ . | . + batch)
#> Warning: NaNs produced

# Add batch to alpha too
m2_3 <- update(m2_2, . ~ . + batch | .)

anova(m2_0, m2_1, m2_2, m2_3)
#> Analysis of Deviance Table
#> 
#> Model 1: yield ~ 1 | 1
#> Model 2: yield ~ temp
#> Model 3: yield ~ temp | batch
#> Model 4: yield ~ temp + batch | batch
#> 
#>      Resid. Df Resid. Dev Df  Deviance   Pr(>Chi)    
#> m2_0  30.00000  -57.02258                            
#> m2_1  29.00000  -80.90259  1  23.88001    < 1e-04 ***
#> m2_2  20.00000 -183.27483  9 102.37224    < 1e-04 ***
#> m2_3  11.00000 -215.08643  9  31.81161 0.00021462 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
AIC(m2_0, m2_1, m2_2, m2_3)
#>      df        AIC
#> m2_0  2  -53.02258
#> m2_1  3  -74.90259
#> m2_2 12 -159.27483
#> m2_3 21 -173.08643

# EXAMPLE 3: Three-part formulas (alpha | beta | gamma)

m3_0 <- gkwreg(yield ~ 1,
  data = GasolineYield,
  family = "gkw",
  control = gkw_control(method = "BFGS", maxit = 2000)
)
#> Warning: NaNs produced

m3_1 <- update(m3_0, . ~ . + temp | . | .)
#> Warning: NaNs produced
m3_2 <- update(m3_1, . ~ . | . + batch | .)
#> Warning: NaNs produced
m3_3 <- update(m3_2, . ~ . | . | . + temp)
#> Warning: NaNs produced

anova(m3_0, m3_1, m3_2, m3_3)
#> Warning: negative deviance change detected; models may not be nested
#> Warning: negative deviance change detected; models may not be nested
#> Analysis of Deviance Table
#> 
#> Model 1: yield ~ 1
#> Model 2: yield ~ temp
#> Model 3: yield ~ temp | batch
#> Model 4: yield ~ temp | batch | temp
#> 
#>      Resid. Df Resid. Dev Df Deviance Pr(>Chi)  
#> m3_0  27.00000  -48.98863                       
#> m3_1  26.00000  -52.38707  1  3.39844 0.065258 .
#> m3_2  17.00000  -52.38706  9  0.00000           
#> m3_3  16.00000  -43.00466  1 -9.38240           
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# EXAMPLE 4: Practical nested model comparison

# Null model
fit0 <- gkwreg(yield ~ 1,
  data = GasolineYield,
  family = "kw",
  control = gkw_control(method = "BFGS", maxit = 2000)
)

# Add main effects to alpha
fit1 <- update(fit0, . ~ . + temp)
fit2 <- update(fit1, . ~ . + batch)

# Model beta parameter
fit3 <- update(fit2, . ~ . | temp)
fit4 <- update(fit3, . ~ . | . + batch)

# Full comparison
anova(fit0, fit1, fit2, fit3, fit4)
#> Warning: negative deviance change detected; models may not be nested
#> Analysis of Deviance Table
#> 
#> Model 1: yield ~ 1
#> Model 2: yield ~ temp
#> Model 3: yield ~ temp + batch
#> Model 4: yield ~ temp + batch | temp
#> Model 5: yield ~ temp + batch | temp + batch
#> 
#>      Resid. Df Resid. Dev Df Deviance  Pr(>Chi)   
#> fit0  30.00000  -16.60035                         
#> fit1  29.00000   -8.36052  1 -8.23983             
#> fit2  20.00000   -8.36052  9    1e-05         1   
#> fit3  19.00000  -18.17449  1  9.81397 0.0017319 **
#> fit4  10.00000  -18.17450  9    1e-05         1   
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
AIC(fit0, fit1, fit2, fit3, fit4)
#>      df        AIC
#> fit0  2 -12.600346
#> fit1  3  -2.360515
#> fit2 12  15.639476
#> fit3 13   7.825509
#> fit4 22  25.825501
BIC(fit0, fit1, fit2, fit3, fit4)
#>      df       BIC
#> fit0  2 -9.668874
#> fit1  3  2.036692
#> fit2 12 33.228307
#> fit3 13 26.880076
#> fit4 22 58.071691

# EXAMPLE 5: Changing other parameters

# Change family
fit_gkw <- update(fit2, family = "gkw")
#> Warning: NaNs produced

# Change link function
fit_logit <- update(fit2, link = list(alpha = "logit"))

# View call without fitting
update(fit2, . ~ . | . + temp, evaluate = FALSE)
#> gkwreg(formula = yield ~ temp + batch | temp, data = GasolineYield, 
#>     family = "kw", control = gkw_control(method = "BFGS", maxit = 2000))
# }
```
