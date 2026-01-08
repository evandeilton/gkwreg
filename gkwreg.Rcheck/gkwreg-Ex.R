pkgname <- "gkwreg"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "gkwreg-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('gkwreg')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("AIC.gkwreg")
### * AIC.gkwreg

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: AIC.gkwreg
### Title: Akaike Information Criterion for GKw Regression Models
### Aliases: AIC.gkwreg

### ** Examples

## No test: 
# Load example data
data(GasolineYield)

# Fit competing models
fit1 <- gkwreg(yield ~ batch, data = GasolineYield, family = "kw")
fit2 <- gkwreg(yield ~ batch + temp, data = GasolineYield, family = "kw")
fit3 <- gkwreg(yield ~ temp, data = GasolineYield, family = "kw")

# Calculate AIC for single model
AIC(fit1)

# Compare multiple models (with proper names)
AIC(fit1, fit2, fit3)

# Use different penalty
AIC(fit1, k = 4)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("AIC.gkwreg", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("BIC.gkwreg")
### * BIC.gkwreg

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: BIC.gkwreg
### Title: Bayesian Information Criterion for GKw Regression Models
### Aliases: BIC.gkwreg

### ** Examples

## No test: 
# Load example data
data(GasolineYield)

# Fit competing models
fit1 <- gkwreg(yield ~ batch, data = GasolineYield, family = "kw")
fit2 <- gkwreg(yield ~ batch + temp, data = GasolineYield, family = "kw")
fit3 <- gkwreg(yield ~ temp, data = GasolineYield, family = "kw")

# Calculate BIC for single model
BIC(fit1)

# Compare multiple models (with proper names)
BIC(fit1, fit2, fit3)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("BIC.gkwreg", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("CarTask")
### * CarTask

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: CarTask
### Title: Partition-Primed Probability Judgement Task for Car Dealership
### Aliases: CarTask
### Keywords: datasets

### ** Examples

## No test: 
require(gkwreg)
require(gkwdist)

data(CarTask)

# Example 1: Task effects on probability judgments
# Do people judge probabilities differently for car vs. salesperson?
fit_kw <- gkwreg(
  probability ~ task,
  data = CarTask,
  family = "kw"
)
summary(fit_kw)

# Interpretation:
# - Alpha: Task type affects mean probability estimate
#   Salesperson condition (1/4 = 0.25) vs. car type (unclear baseline)

# Example 2: Individual differences model
# Need for Closure/Certainty may moderate probability judgments
fit_kw_nfcc <- gkwreg(
  probability ~ task * NFCCscale |
    task,
  data = CarTask,
  family = "kw"
)
summary(fit_kw_nfcc)

# Interpretation:
# - Interaction: NFCC may have different effects depending on task
#   People high in need for certainty may respond differently to
#   explicit partitions (4 salespersons) vs. implicit partitions (car types)
# - Beta: Precision varies by task type

# Example 3: Exponentiated Kumaraswamy for extreme estimates
# Some participants may give very extreme probability estimates
fit_ekw <- gkwreg(
  probability ~ task * NFCCscale | # alpha
    task | # beta
    task, # lambda: extremity differs by task
  data = CarTask,
  family = "ekw"
)
summary(fit_ekw)

# Interpretation:
# - Lambda varies by task: Salesperson condition (explicit partition)
#   may produce more extreme estimates (closer to 0 or 1)

# Visualization: Probability by task and NFCC
plot(probability ~ NFCCscale,
  data = CarTask,
  col = c("blue", "red")[task], pch = 19,
  xlab = "Need for Closure/Certainty", ylab = "Probability Estimate",
  main = "Car Task: Individual Differences in Probability Judgment"
)
legend("topright",
  legend = levels(CarTask$task),
  col = c("blue", "red"), pch = 19
)

# Distribution comparison
boxplot(probability ~ task,
  data = CarTask,
  xlab = "Task Condition", ylab = "Probability Estimate",
  main = "Partition Priming Effects",
  col = c("lightblue", "lightcoral")
)
abline(h = 0.25, lty = 2, col = "gray")
text(1.5, 0.27, "Uniform (1/4)", col = "gray")
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("CarTask", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("FoodExpenditure")
### * FoodExpenditure

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: FoodExpenditure
### Title: Proportion of Household Income Spent on Food
### Aliases: FoodExpenditure
### Keywords: datasets

### ** Examples

## No test: 
require(gkwreg)
require(gkwdist)

data(FoodExpenditure)
FoodExpenditure$prop <- FoodExpenditure$food / FoodExpenditure$income

# Example 1: Basic Kumaraswamy regression
# Proportion spent on food decreases with income (Engel's law)
# Larger households spend more on food
fit_kw <- gkwreg(prop ~ income + persons,
  data = FoodExpenditure,
  family = "kw"
)
summary(fit_kw)

# Interpretation:
# - Alpha: Negative income effect (Engel's law)
#   Positive household size effect
# - Beta: Constant precision (homoscedastic model)

# Example 2: Heteroscedastic model
# Variability in food proportion may differ by income and household size
fit_kw_hetero <- gkwreg(
  prop ~ income + persons |
    income + persons,
  data = FoodExpenditure,
  family = "kw"
)
summary(fit_kw_hetero)

# Interpretation:
# - Beta: Precision varies with both income and household size
#   Wealthier or larger households may show different spending variability

# Test for heteroscedasticity
anova(fit_kw, fit_kw_hetero)

# Example 3: Exponentiated Kumaraswamy for extreme spending patterns
# Some households may have unusual food spending (very frugal or lavish)
fit_ekw <- gkwreg(
  prop ~ income + persons | # alpha
    persons | # beta: household size affects precision
    income, # lambda: income affects extremity
  data = FoodExpenditure,
  family = "ekw"
)
summary(fit_ekw)

# Interpretation:
# - Lambda: Income level affects tail behavior
#   Rich households may show more extreme (unusual) spending patterns

# Visualization: Engel curve
plot(prop ~ income,
  data = FoodExpenditure,
  xlab = "Annual Income ($)", ylab = "Proportion Spent on Food",
  main = "Engel Curve for Food Expenditure"
)
# Add fitted values
FoodExpenditure$fitted_kw <- fitted(fit_kw)
points(FoodExpenditure$income, FoodExpenditure$fitted_kw,
  col = "blue", pch = 19, cex = 0.8
)
legend("topright",
  legend = c("Observed", "Fitted"),
  col = c("black", "blue"), pch = c(1, 19)
)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("FoodExpenditure", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("GasolineYield")
### * GasolineYield

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: GasolineYield
### Title: Gasoline Yield from Crude Oil
### Aliases: GasolineYield
### Keywords: datasets

### ** Examples

## No test: 
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
summary(fit_kw_full)

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

# Interpretation:
# - Lambda varies by batch: Some crude oils have more extreme
#   yield distributions (heavy tails for very high/low yields)

# Model comparison: Does tail flexibility improve fit?
anova(fit_kw, fit_ekw)

# Diagnostic plots
par(mfrow = c(2, 2))
plot(fit_kw, which = c(1, 2, 4, 5))
par(mfrow = c(1, 1))
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("GasolineYield", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("ImpreciseTask")
### * ImpreciseTask

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ImpreciseTask
### Title: Imprecise Probabilities for Sunday Weather and Boeing Stock Task
### Aliases: ImpreciseTask
### Keywords: datasets

### ** Examples

## No test: 
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
summary(fit_mc)

# Interpretation:
# - Lambda varies by task: Weather vs. stock may produce
#   different patterns of extreme probability assessments
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ImpreciseTask", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("LossAversion")
### * LossAversion

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: LossAversion
### Title: (No) Myopic Loss Aversion in Adolescents
### Aliases: LossAversion
### Keywords: datasets

### ** Examples

## No test: 
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

# Interpretation:
# - Model accommodates extreme risk-taking (all-in or all-out strategies)

# Compare models
anova(fit_kw, fit_kw_interact)

# Visualization: Investment by horizon
boxplot(invest ~ arrangement,
  data = LossAversion,
  xlab = "Investment Horizon", ylab = "Proportion Invested",
  main = "No Myopic Loss Aversion in Adolescents",
  col = c("lightblue", "lightgreen", "lightyellow")
)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("LossAversion", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("MockJurors")
### * MockJurors

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: MockJurors
### Title: Confidence of Mock Jurors in Their Verdicts
### Aliases: MockJurors
### Keywords: datasets

### ** Examples

## No test: 
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

# Interpretation:
# - Full interaction: Third verdict option may have different effects
#   depending on whether evidence is conflicting

# Test interaction significance
anova(fit_kw, fit_kw_interact)

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
summary(fit_ekw)

# Compare 3-parameter models
AIC(fit_ekw, fit_mc)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("MockJurors", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ReadingSkills")
### * ReadingSkills

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ReadingSkills
### Title: Dyslexia and IQ Predicting Reading Accuracy
### Aliases: ReadingSkills
### Keywords: datasets

### ** Examples

## No test: 
require(gkwreg)
require(gkwdist)

data(ReadingSkills)

# Example 1: Standard Kumaraswamy with interaction and heteroscedasticity
# Mean: Dyslexia × IQ interaction (do groups differ in IQ effect?)
# Precision: Main effects (variability differs by group and IQ level)
fit_kw <- gkwreg(
  accuracy ~ dyslexia * iq |
    dyslexia + iq,
  data = ReadingSkills,
  family = "kw",
  control = gkw_control(method = "L-BFGS-B", maxit = 2000)
)
summary(fit_kw)

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

# Interpretation:
# - Lambda varies by dyslexia status: Controls have stronger ceiling effect
#   (more compression at high accuracy) than dyslexic children

# Test if ceiling effect modeling improves fit
anova(fit_kw, fit_ekw)

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

# Compare 3-parameter models
AIC(fit_ekw, fit_mc)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ReadingSkills", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("StressAnxiety")
### * StressAnxiety

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: StressAnxiety
### Title: Dependency of Anxiety on Stress
### Aliases: StressAnxiety
### Keywords: datasets

### ** Examples

## No test: 
require(gkwreg)
require(gkwdist)

data(StressAnxiety)

# Example 1: Basic heteroscedastic relationship
# Mean anxiety increases with stress
# Variability in anxiety also changes with stress
fit_kw <- gkwreg(
  anxiety ~ stress |
    stress,
  data = StressAnxiety,
  family = "kw"
)
summary(fit_kw)

# Interpretation:
# - Alpha: Positive relationship between stress and mean anxiety
# - Beta: Precision changes with stress level
#   (anxiety becomes more/less variable at different stress levels)

# Compare to homoscedastic model
fit_kw_homo <- gkwreg(anxiety ~ stress,
  data = StressAnxiety, family = "kw"
)
anova(fit_kw_homo, fit_kw)

# Example 2: Nonlinear stress effects via polynomial
# Stress-anxiety relationship often shows threshold or saturation effects
fit_kw_poly <- gkwreg(
  anxiety ~ poly(stress, 2) | # quadratic mean
    poly(stress, 2), # quadratic precision
  data = StressAnxiety,
  family = "kw"
)
summary(fit_kw_poly)

# Interpretation:
# - Quadratic terms allow for:
#   * Threshold effects (anxiety accelerates at high stress)
#   * Saturation effects (anxiety plateaus at extreme stress)

# Test nonlinearity
anova(fit_kw, fit_kw_poly)

# Example 3: Exponentiated Kumaraswamy for extreme anxiety patterns
# Some individuals may show very extreme anxiety responses to stress
fit_ekw <- gkwreg(
  anxiety ~ poly(stress, 2) | # alpha: quadratic mean
    poly(stress, 2) | # beta: quadratic precision
    stress, # lambda: linear tail effect
  data = StressAnxiety,
  family = "ekw"
)
summary(fit_ekw)

# Interpretation:
# - Lambda: Linear component captures asymmetry at extreme stress levels
#   (very high stress may produce different tail behavior)

# Example 4: McDonald distribution for highly skewed anxiety
# Anxiety distributions are often right-skewed (ceiling effects)
fit_mc <- gkwreg(
  anxiety ~ poly(stress, 2) | # gamma
    poly(stress, 2) | # delta
    stress, # lambda: extremity
  data = StressAnxiety,
  family = "mc",
  control = gkw_control(method = "BFGS", maxit = 1500)
)
summary(fit_mc)

# Compare models
AIC(fit_kw, fit_kw_poly, fit_ekw, fit_mc)

# Visualization: Stress-Anxiety relationship
plot(anxiety ~ stress,
  data = StressAnxiety,
  xlab = "Stress Level", ylab = "Anxiety Level",
  main = "Stress-Anxiety Relationship with Heteroscedasticity",
  pch = 19, col = rgb(0, 0, 1, 0.3)
)

# Add fitted curve
stress_seq <- seq(min(StressAnxiety$stress), max(StressAnxiety$stress),
  length.out = 100
)
pred_mean <- predict(fit_kw, newdata = data.frame(stress = stress_seq))
lines(stress_seq, pred_mean, col = "red", lwd = 2)

# Add lowess smooth for comparison
lines(lowess(StressAnxiety$stress, StressAnxiety$anxiety),
  col = "blue", lwd = 2, lty = 2
)
legend("topleft",
  legend = c("Kumaraswamy fit", "Lowess smooth"),
  col = c("red", "blue"), lwd = 2, lty = c(1, 2)
)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("StressAnxiety", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("WeatherTask")
### * WeatherTask

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: WeatherTask
### Title: Weather Task with Priming and Precise and Imprecise
###   Probabilities
### Aliases: WeatherTask
### Keywords: datasets

### ** Examples

## No test: 
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

# Interpretation:
# - Alpha: Interaction tests if partition priming works differently
#   for precise vs. imprecise probability judgments
# - Beta: Precision varies by experimental condition

# Test interaction
anova(fit_kw, fit_kw_interact)

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

# Interpretation:
# - Lambda varies by priming: Seven-fold priming may produce more
#   extreme/polarized probability judgments
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("WeatherTask", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("anova.gkwreg")
### * anova.gkwreg

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: anova.gkwreg
### Title: Analysis of Deviance for GKw Regression Models
### Aliases: anova.gkwreg

### ** Examples

## No test: 
# Load example data
data(GasolineYield)

# Fit a series of nested models
fit1 <- gkwreg(yield ~ 1, data = GasolineYield, family = "kw")
fit2 <- gkwreg(yield ~ temp, data = GasolineYield, family = "kw")
fit3 <- gkwreg(yield ~ batch + temp, data = GasolineYield, family = "kw")

# ANOVA table for single model
anova(fit3)

# Compare nested models using likelihood ratio tests
anova(fit1, fit2, fit3)
#> Model 1 vs 2: Adding temperature is highly significant (p < 0.001)
#> Model 2 vs 3: Adding batch is highly significant (p < 0.001)

# Compare two models
anova(fit2, fit3, test = "Chisq")

# Suppress test statistics
anova(fit1, fit2, fit3, test = "none")
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("anova.gkwreg", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("confint.gkwreg")
### * confint.gkwreg

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: confint.gkwreg
### Title: Confidence Intervals for Generalized Kumaraswamy Regression
###   Parameters
### Aliases: confint.gkwreg

### ** Examples

## No test: 
data(GasolineYield)
fit <- gkwreg(yield ~ batch + temp, data = GasolineYield, family = "kw")

# 95 percent confidence intervals
confint(fit)

# 90 percent confidence intervals
confint(fit, level = 0.90)

# Specific parameters
confint(fit, parm = "alpha:(Intercept)")
confint(fit, parm = 1:3)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("confint.gkwreg", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("family.gkwreg")
### * family.gkwreg

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: family.gkwreg
### Title: Extract Family from GKw Regression Model
### Aliases: family.gkwreg

### ** Examples

## No test: 
data(GasolineYield)
fit <- gkwreg(yield ~ batch + temp, data = GasolineYield, family = "kw")
family(fit)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("family.gkwreg", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("fitted.gkwreg")
### * fitted.gkwreg

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: fitted.gkwreg
### Title: Extract Fitted Values from a Generalized Kumaraswamy Regression
###   Model
### Aliases: fitted.gkwreg
### Keywords: fitted methods regression

### ** Examples

## No test: 
require(gkwreg)
require(gkwdist)

# Example 1: Basic usage with FoodExpenditure data
data(FoodExpenditure)
FoodExpenditure$prop <- FoodExpenditure$food / FoodExpenditure$income

fit_kw <- gkwreg(prop ~ income + persons | income,
  data = FoodExpenditure,
  family = "kw"
)

# Extract fitted values
fitted_vals <- fitted(fit_kw)

# Visualize fit quality
plot(FoodExpenditure$prop, fitted_vals,
  xlab = "Observed Proportion",
  ylab = "Fitted Values",
  main = "Observed vs Fitted: Food Expenditure",
  pch = 19, col = rgb(0, 0, 1, 0.5)
)
abline(0, 1, col = "red", lwd = 2)

# Calculate R-squared analogue
cor(FoodExpenditure$prop, fitted_vals)^2

# Example 2: Comparing fitted values across families
data(GasolineYield)

fit_ekw <- gkwreg(yield ~ batch + temp | temp | batch,
  data = GasolineYield,
  family = "ekw"
)

# Fitted values under different family assumptions
fitted_ekw <- fitted(fit_ekw)
fitted_kw <- fitted(fit_ekw, family = "kw")
fitted_beta <- fitted(fit_ekw, family = "beta")

# Compare differences
comparison <- data.frame(
  EKW = fitted_ekw,
  KW = fitted_kw,
  Beta = fitted_beta,
  Diff_EKW_KW = fitted_ekw - fitted_kw,
  Diff_EKW_Beta = fitted_ekw - fitted_beta
)
head(comparison)

# Visualize differences
par(mfrow = c(1, 2))
plot(fitted_ekw, fitted_kw,
  xlab = "EKW Fitted", ylab = "KW Fitted",
  main = "EKW vs KW Family Assumptions",
  pch = 19, col = "darkblue"
)
abline(0, 1, col = "red", lty = 2)

plot(fitted_ekw, fitted_beta,
  xlab = "EKW Fitted", ylab = "Beta Fitted",
  main = "EKW vs Beta Family Assumptions",
  pch = 19, col = "darkgreen"
)
abline(0, 1, col = "red", lty = 2)
par(mfrow = c(1, 1))

# Example 3: Diagnostic plot with confidence bands
data(ReadingSkills)

fit_mc <- gkwreg(
  accuracy ~ dyslexia * iq | dyslexia + iq | dyslexia,
  data = ReadingSkills,
  family = "mc"
)

fitted_vals <- fitted(fit_mc)

# Residual plot
residuals_resp <- ReadingSkills$accuracy - fitted_vals

plot(fitted_vals, residuals_resp,
  xlab = "Fitted Values",
  ylab = "Raw Residuals",
  main = "Residual Plot: Reading Accuracy",
  pch = 19, col = ReadingSkills$dyslexia,
  ylim = range(residuals_resp) * 1.2
)
abline(h = 0, col = "red", lwd = 2, lty = 2)
lowess_fit <- lowess(fitted_vals, residuals_resp)
lines(lowess_fit, col = "blue", lwd = 2)
legend("topright",
  legend = c("Control", "Dyslexic", "Zero Line", "Lowess"),
  col = c("black", "red", "red", "blue"),
  pch = c(19, 19, NA, NA),
  lty = c(NA, NA, 2, 1),
  lwd = c(NA, NA, 2, 2)
)

# Example 4: Large dataset efficiency check
set.seed(2024)
n <- 5000
x1 <- rnorm(n)
x2 <- runif(n, -2, 2)
alpha <- exp(0.3 + 0.5 * x1)
beta <- exp(1.2 - 0.4 * x2)
y <- rkw(n, alpha, beta)
large_data <- data.frame(y = y, x1 = x1, x2 = x2)

fit_large <- gkwreg(y ~ x1 | x2,
  data = large_data,
  family = "kw"
)

# Time the extraction
system.time({
  fitted_large <- fitted(fit_large)
})

# Verify extraction
length(fitted_large)
summary(fitted_large)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("fitted.gkwreg", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("formula.gkwreg")
### * formula.gkwreg

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: formula.gkwreg
### Title: Extract Formula from GKw Regression Model
### Aliases: formula.gkwreg

### ** Examples

## No test: 
data(GasolineYield)

# Simple formula
fit1 <- gkwreg(yield ~ batch + temp, data = GasolineYield, family = "kw")
formula(fit1)

# Two-part formula
fit2 <- gkwreg(yield ~ temp | batch, data = GasolineYield, family = "kw")
formula(fit2)

# Five-part formula
fit3 <- gkwreg(yield ~ temp | batch | temp | 1 | 1,
  data = GasolineYield, family = "gkw"
)
formula(fit3)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("formula.gkwreg", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("getCall.gkwreg")
### * getCall.gkwreg

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: getCall.gkwreg
### Title: Get Call from GKw Regression Model
### Aliases: getCall.gkwreg

### ** Examples

## No test: 
data(GasolineYield)
fit <- gkwreg(yield ~ batch + temp, data = GasolineYield, family = "kw")
getCall(fit)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("getCall.gkwreg", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("gkw_control")
### * gkw_control

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: gkw_control
### Title: Control Parameters for Generalized Kumaraswamy Regression
### Aliases: gkw_control print.gkw_control

### ** Examples

## No test: 
# Default control (used automatically if not specified)
ctrl <- gkw_control()
print(ctrl)

# Increase iterations for difficult problem
ctrl_robust <- gkw_control(maxit = 1000, trace = 1)

# Try alternative optimizer
ctrl_bfgs <- gkw_control(method = "BFGS")

# Fast fitting without standard errors
ctrl_fast <- gkw_control(hessian = FALSE)

# Verbose debugging
ctrl_debug <- gkw_control(silent = FALSE, trace = 2)

# Custom starting values
ctrl_start <- gkw_control(
  start = list(
    alpha = c(0.5, 0.2),
    beta = c(1.0, -0.3)
  )
)

# Configure Nelder-Mead with custom reflection/contraction
ctrl_nm <- gkw_control(
  method = "Nelder-Mead",
  alpha = 1.5,
  beta = 0.75
)

# Configure L-BFGS-B for bounded optimization
ctrl_lbfgsb <- gkw_control(
  method = "L-BFGS-B",
  factr = 1e6,
  lmm = 10
)

# Configure SANN for rough surfaces
ctrl_sann <- gkw_control(
  method = "SANN",
  temp = 20,
  tmax = 20,
  maxit = 20000
)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("gkw_control", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("gkwreg")
### * gkwreg

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: gkwreg
### Title: Fit Generalized Kumaraswamy Regression Models
### Aliases: gkwreg
### Keywords: models regression

### ** Examples

## No test: 
# SECTION 1: Basic Usage - Getting Started
# Load packages and data
library(gkwreg)
library(gkwdist)
data(GasolineYield)

# Example 1.1: Simplest possible model (intercept-only, all defaults)
fit_basic <- gkwreg(yield ~ 1, data = GasolineYield, family = "kw")
summary(fit_basic)

# Example 1.2: Model with predictors (uses all defaults)
# Default: family = "gkw", method = "nlminb", hessian = TRUE
fit_default <- gkwreg(yield ~ batch + temp, data = GasolineYield)
summary(fit_default)

# Example 1.3: Kumaraswamy model (two-parameter family)
# Default link functions: log for both alpha and beta
fit_kw <- gkwreg(yield ~ batch + temp, data = GasolineYield, family = "kw")
summary(fit_kw)

par(mfrow = c(3, 2))
plot(fit_kw, ask = FALSE)

# Example 1.4: Beta model for comparison
# Default links: log for gamma and delta
fit_beta <- gkwreg(yield ~ batch + temp, data = GasolineYield, family = "beta")

# Compare models using AIC/BIC
AIC(fit_kw, fit_beta)
BIC(fit_kw, fit_beta)

# SECTION 2: Using gkw_control() for Customization

# Example 2.1: Change optimization method to BFGS
fit_bfgs <- gkwreg(
  yield ~ batch + temp,
  data = GasolineYield,
  family = "kw",
  control = gkw_control(method = "BFGS")
)
summary(fit_bfgs)

# Example 2.2: Increase iterations and enable verbose output
fit_verbose <- gkwreg(
  yield ~ batch + temp,
  data = GasolineYield,
  family = "kw",
  control = gkw_control(
    method = "nlminb",
    maxit = 1000,
    silent = FALSE, # Show optimization progress
    trace = 1 # Print iteration details
  )
)

# Example 2.3: Fast fitting without standard errors
# Useful for model exploration or large datasets
fit_fast <- gkwreg(
  yield ~ batch + temp,
  data = GasolineYield,
  family = "kw",
  control = gkw_control(hessian = FALSE)
)
# Note: Cannot compute confint() without hessian
coef(fit_fast) # Point estimates still available

# Example 2.4: Custom convergence tolerances
fit_tight <- gkwreg(
  yield ~ batch + temp,
  data = GasolineYield,
  family = "kw",
  control = gkw_control(
    reltol = 1e-10, # Tighter convergence
    maxit = 2000 # More iterations allowed
  )
)

# SECTION 3: Advanced Formula Specifications

# Example 3.1: Different predictors for different parameters
# alpha depends on batch, beta depends on temp
fit_diff <- gkwreg(
  yield ~ batch | temp,
  data = GasolineYield,
  family = "kw"
)
summary(fit_diff)

# Example 3.2: Intercept-only for one parameter
# alpha varies with predictors, beta is constant
fit_partial <- gkwreg(
  yield ~ batch + temp | 1,
  data = GasolineYield,
  family = "kw"
)

# Example 3.3: Complex model with interactions
fit_interact <- gkwreg(
  yield ~ batch * temp | temp + I(temp^2),
  data = GasolineYield,
  family = "kw"
)

# SECTION 4: Working with Different Families

# Example 4.1: Fit multiple families and compare
families <- c("beta", "kw", "ekw", "bkw", "gkw")
fits <- lapply(families, function(fam) {
  gkwreg(yield ~ batch + temp, data = GasolineYield, family = fam)
})
names(fits) <- families

# Compare via information criteria
comparison <- data.frame(
  Family = families,
  LogLik = sapply(fits, logLik),
  AIC = sapply(fits, AIC),
  BIC = sapply(fits, BIC),
  npar = sapply(fits, function(x) x$npar)
)
print(comparison)

# Example 4.2: Formal nested model testing
fit_kw <- gkwreg(yield ~ batch + temp, GasolineYield, family = "kw")
fit_ekw <- gkwreg(yield ~ batch + temp, GasolineYield, family = "ekw")
fit_gkw <- gkwreg(yield ~ batch + temp, GasolineYield, family = "gkw")
anova(fit_kw, fit_ekw, fit_gkw)

# SECTION 5: Link Functions and Scales

# Example 5.1: Custom link functions
fit_links <- gkwreg(
  yield ~ batch + temp,
  data = GasolineYield,
  family = "kw",
  link = list(alpha = "sqrt", beta = "log")
)

# Example 5.2: Custom link scales
# Smaller scale = steeper response curve
fit_scale <- gkwreg(
  yield ~ batch + temp,
  data = GasolineYield,
  family = "kw",
  link_scale = list(alpha = 5, beta = 15)
)

# Example 5.3: Uniform link for all parameters
fit_uniform <- gkwreg(
  yield ~ batch + temp,
  data = GasolineYield,
  family = "kw",
  link = "log" # Single string applied to all
)

# SECTION 6: Prediction and Inference

# Fit model for prediction examples
fit <- gkwreg(yield ~ batch + temp, GasolineYield, family = "kw")

# Example 6.1: Confidence intervals at different levels
confint(fit, level = 0.95) # 95% CI
confint(fit, level = 0.90) # 90% CI
confint(fit, level = 0.99) # 99% CI

# SECTION 7: Diagnostic Plots and Model Checking

fit <- gkwreg(yield ~ batch + temp, GasolineYield, family = "kw")

# Example 7.1: All diagnostic plots (default)
par(mfrow = c(3, 2))
plot(fit, ask = FALSE)

# Example 7.2: Select specific plots
par(mfrow = c(3, 1))
plot(fit, which = c(2, 4, 5)) # Cook's distance, Residuals, Half-normal

# Example 7.3: Using ggplot2 for modern graphics
plot(fit, use_ggplot = TRUE, arrange_plots = TRUE)

# Example 7.4: Customized half-normal plot
par(mfrow = c(1, 1))
plot(fit,
  which = 5,
  type = "quantile",
  nsim = 200, # More simulations for smoother envelope
  level = 0.95
) # 95% confidence envelope

# Example 7.5: Extract diagnostic data programmatically
diagnostics <- plot(fit, save_diagnostics = TRUE)
head(diagnostics$data) # Residuals, Cook's distance, etc.

# SECTION 8: Real Data Example - Food Expenditure

# Load and prepare data
data(FoodExpenditure, package = "betareg")
food_data <- FoodExpenditure
food_data$prop <- food_data$food / food_data$income

# Example 8.1: Basic model
fit_food <- gkwreg(
  prop ~ persons | income,
  data = food_data,
  family = "kw"
)
summary(fit_food)

# Example 8.2: Compare with Beta regression
fit_food_beta <- gkwreg(
  prop ~ persons | income,
  data = food_data,
  family = "beta"
)

# Which fits better?
AIC(fit_food, fit_food_beta)

# Example 8.3: Model diagnostics
par(mfrow = c(3, 1))
plot(fit_food, which = c(2, 5, 6))

# Example 8.4: Interpretation via effects
# How does proportion spent on food change with income?
income_seq <- seq(min(food_data$income), max(food_data$income), length = 50)
pred_data <- data.frame(
  persons = median(food_data$persons),
  income = income_seq
)
pred_food <- predict(fit_food, newdata = pred_data, type = "response")

par(mfrow = c(1, 1))
plot(food_data$income, food_data$prop,
  xlab = "Income", ylab = "Proportion Spent on Food",
  main = "Food Expenditure Pattern"
)
lines(income_seq, pred_food, col = "red", lwd = 2)

# SECTION 9: Simulation Studies

# Example 9.1: Simple Kumaraswamy simulation
set.seed(123)
n <- 500
x1 <- runif(n, -2, 2)
x2 <- rnorm(n)

# True model: log(alpha) = 0.8 + 0.3*x1, log(beta) = 1.2 - 0.2*x2
eta_alpha <- 0.8 + 0.3 * x1
eta_beta <- 1.2 - 0.2 * x2
alpha_true <- exp(eta_alpha)
beta_true <- exp(eta_beta)

# Generate response
y <- rkw(n, alpha = alpha_true, beta = beta_true)
sim_data <- data.frame(y = y, x1 = x1, x2 = x2)

# Fit and check parameter recovery
fit_sim <- gkwreg(y ~ x1 | x2, data = sim_data, family = "kw")

# Compare estimated vs true coefficients
cbind(
  True = c(0.8, 0.3, 1.2, -0.2),
  Estimated = coef(fit_sim),
  SE = fit_sim$se
)

# Example 9.2: Complex simulation with all five parameters
set.seed(2203)
n <- 2000
x <- runif(n, -1, 1)

# True parameters
alpha <- exp(0.5 + 0.3 * x)
beta <- exp(1.0 - 0.2 * x)
gamma <- exp(0.7 + 0.4 * x)
delta <- plogis(0.0 + 0.5 * x) # logit scale
lambda <- exp(-0.3 + 0.2 * x)

# Generate from GKw
y <- rgkw(n,
  alpha = alpha, beta = beta, gamma = gamma,
  delta = delta, lambda = lambda
)
sim_data2 <- data.frame(y = y, x = x)

# Fit GKw model
fit_gkw <- gkwreg(
  y ~ x | x | x | x | x,
  data = sim_data2,
  family = "gkw",
  control = gkw_control(method = "L-BFGS-B", maxit = 2000)
)
summary(fit_gkw)

# SECTION 10: Handling Convergence Issues

# Example 10.1: Try different optimizers
methods <- c("nlminb", "BFGS", "Nelder-Mead", "CG")
fits_methods <- lapply(methods, function(m) {
  tryCatch(
    gkwreg(yield ~ batch + temp, GasolineYield,
      family = "kw",
      control = gkw_control(method = m, silent = TRUE)
    ),
    error = function(e) NULL
  )
})
names(fits_methods) <- methods

# Check which converged
converged <- sapply(fits_methods, function(f) {
  if (is.null(f)) {
    return(FALSE)
  }
  f$convergence
})
print(converged)

# Example 10.2: Verbose mode for debugging
fit_debug <- gkwreg(
  yield ~ batch + temp,
  data = GasolineYield,
  family = "kw",
  control = gkw_control(
    method = "BFGS",
    silent = TRUE,
    trace = 0, # 2, Maximum verbosity
    maxit = 1000
  )
)

# SECTION 11: Memory and Performance Optimization

# Example 11.1: Minimal object for large datasets
fit_minimal <- gkwreg(
  yield ~ batch + temp,
  data = GasolineYield,
  family = "kw",
  model = FALSE, # Don't store model frame
  x = FALSE, # Don't store design matrices
  y = FALSE, # Don't store response
  control = gkw_control(hessian = FALSE) # Skip Hessian
)

# Much smaller object
object.size(fit_minimal)

# Trade-off: Limited post-fitting capabilities
# Can still use: coef(), logLik(), AIC(), BIC()
# Cannot use: predict(), some diagnostics

# Example 11.2: Fast exploratory analysis
# Fit many models quickly without standard errors
formulas <- list(
  yield ~ batch,
  yield ~ temp,
  yield ~ batch + temp,
  yield ~ batch * temp
)

fast_fits <- lapply(formulas, function(f) {
  gkwreg(f, GasolineYield,
    family = "kw",
    control = gkw_control(hessian = FALSE),
    model = FALSE, x = FALSE, y = FALSE
  )
})

# Compare models via AIC
sapply(fast_fits, AIC)

# Refit best model with full inference
best_formula <- formulas[[which.min(sapply(fast_fits, AIC))]]
fit_final <- gkwreg(best_formula, GasolineYield, family = "kw")
summary(fit_final)

# SECTION 12: Model Selection and Comparison

# Example 12.1: Nested model testing
fit1 <- gkwreg(yield ~ 1, GasolineYield, family = "kw")
fit2 <- gkwreg(yield ~ batch, GasolineYield, family = "kw")
fit3 <- gkwreg(yield ~ batch + temp, GasolineYield, family = "kw")

# Likelihood ratio tests
anova(fit1, fit2, fit3)

# Example 12.2: Information criteria table
models <- list(
  "Intercept only" = fit1,
  "Batch effect" = fit2,
  "Batch + Temp" = fit3
)

ic_table <- data.frame(
  Model = names(models),
  df = sapply(models, function(m) m$npar),
  LogLik = sapply(models, logLik),
  AIC = sapply(models, AIC),
  BIC = sapply(models, BIC),
  Delta_AIC = sapply(models, AIC) - min(sapply(models, AIC))
)
print(ic_table)

# Example 12.3: Cross-validation for predictive performance
# 5-fold cross-validation
set.seed(2203)
n <- nrow(GasolineYield)
folds <- sample(rep(1:5, length.out = n))

cv_rmse <- sapply(1:5, function(fold) {
  train <- GasolineYield[folds != fold, ]
  test <- GasolineYield[folds == fold, ]

  fit_train <- gkwreg(yield ~ batch + temp, train,
    family = "kw"
  )
  pred_test <- predict(fit_train, newdata = test, type = "response")

  sqrt(mean((test$yield - pred_test)^2))
})

cat("Cross-validated RMSE:", mean(cv_rmse), "\n")
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("gkwreg", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("logLik.gkwreg")
### * logLik.gkwreg

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: logLik.gkwreg
### Title: Extract Log-Likelihood from Generalized Kumaraswamy Regression
###   Models
### Aliases: logLik.gkwreg

### ** Examples

## No test: 
# Load example data
data(GasolineYield)

# Fit a Kumaraswamy regression model
fit <- gkwreg(yield ~ batch + temp, data = GasolineYield, family = "kw")

# Extract log-likelihood
ll <- logLik(fit)
print(ll)

# Access attributes
cat("Log-likelihood:", as.numeric(ll), "\n")
cat("Parameters:", attr(ll, "df"), "\n")
cat("Observations:", attr(ll, "nobs"), "\n")
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("logLik.gkwreg", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("lrtest")
### * lrtest

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: lrtest
### Title: Likelihood Ratio Test for Nested GKw Models
### Aliases: lrtest

### ** Examples

## No test: 
data(GasolineYield)

# Fit nested models
fit_restricted <- gkwreg(yield ~ temp, data = GasolineYield, family = "kw")
fit_full <- gkwreg(yield ~ batch + temp, data = GasolineYield, family = "kw")

# Likelihood ratio test
lrtest(fit_restricted, fit_full)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("lrtest", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("model.frame.gkwreg")
### * model.frame.gkwreg

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: model.frame.gkwreg
### Title: Extract Model Frame from GKw Regression Model
### Aliases: model.frame.gkwreg

### ** Examples

## No test: 
data(GasolineYield)
fit <- gkwreg(yield ~ batch + temp, data = GasolineYield, family = "kw")
head(model.frame(fit))
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("model.frame.gkwreg", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("model.matrix.gkwreg")
### * model.matrix.gkwreg

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: model.matrix.gkwreg
### Title: Extract Model Matrix from GKw Regression Model
### Aliases: model.matrix.gkwreg

### ** Examples

## No test: 
data(GasolineYield)
fit <- gkwreg(yield ~ batch + temp, data = GasolineYield, family = "kw")
head(model.matrix(fit))
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("model.matrix.gkwreg", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("nobs.gkwreg")
### * nobs.gkwreg

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: nobs.gkwreg
### Title: Number of Observations for GKw Regression Models
### Aliases: nobs.gkwreg

### ** Examples

## No test: 
data(GasolineYield)
fit <- gkwreg(yield ~ batch + temp, data = GasolineYield, family = "kw")
nobs(fit)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("nobs.gkwreg", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot.gkwreg")
### * plot.gkwreg

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot.gkwreg
### Title: Diagnostic Plots for Generalized Kumaraswamy Regression Models
### Aliases: plot.gkwreg
### Keywords: diagnostics hplot methods plot regression

### ** Examples

## No test: 
# EXAMPLE 1: Basic Usage with Default Settings

# Simulate data
library(gkwdist)

set.seed(123)
n <- 200
x1 <- runif(n, -2, 2)
x2 <- rnorm(n)

# True model parameters
alpha_true <- exp(0.7 + 0.3 * x1)
beta_true <- exp(1.2 - 0.2 * x2)

# Generate response
y <- rkw(n, alpha = alpha_true, beta = beta_true)
df <- data.frame(y = y, x1 = x1, x2 = x2)

# Fit model
model <- gkwreg(y ~ x1 | x2, data = df, family = "kw")

# Generate all diagnostic plots with defaults
par(mfrow = c(3, 2))
plot(model, ask = FALSE)

# EXAMPLE 2: Selective Plots with Custom Residual Type

# Focus on key diagnostic plots only
par(mfrow = c(3, 1))
plot(model,
  which = c(2, 4, 5), # Cook's distance, Resid vs LinPred, Half-normal
  type = "pearson"
) # Use Pearson residuals

# Check for influential points (plot 2) and non-linearity (plot 4)
par(mfrow = c(2, 1))
plot(model,
  which = c(2, 4),
  type = "deviance"
)

# EXAMPLE 3: Caption Customization - New Named List Interface

# Customize only specific plot titles (RECOMMENDED NEW WAY)
par(mfrow = c(3, 1))
plot(model,
  which = c(1, 4, 6),
  caption = list(
    "1" = "Time Pattern Check",
    "4" = "Linearity Assessment",
    "6" = "Predictive Accuracy"
  )
)

# Customize subtitle and main title
par(mfrow = c(2, 1))
plot(model,
  which = c(1, 5),
  main = "Model Diagnostics",
  sub.caption = "Kumaraswamy Regression - Training Data",
  caption = list("5" = "Normality Check with 95% Envelope")
)

# Suppress subtitle entirely
par(mfrow = c(3, 2))
plot(model, sub.caption = "")

# EXAMPLE 4: Backward Compatible Caption (Vector Interface)

# OLD WAY - still fully supported
par(mfrow = c(3, 2))
plot(model,
  which = 1:6,
  caption = c(
    "Residual Pattern Analysis",
    "Influence Diagnostics",
    "Leverage Assessment",
    "Linearity Check",
    "Distributional Fit",
    "Prediction Quality"
  )
)

# EXAMPLE 5: ggplot2 Graphics with Theming

# Modern publication-quality plots
plot(model,
  use_ggplot = TRUE,
  arrange_plots = TRUE
)

# With custom theme
plot(model,
  use_ggplot = TRUE,
  theme_fn = ggplot2::theme_bw,
  arrange_plots = TRUE
)

# With classic theme and custom colors (via ...)
plot(model,
  use_ggplot = TRUE,
  theme_fn = ggplot2::theme_classic,
  arrange_plots = TRUE
)

# EXAMPLE 6: Arranged Multi-Panel ggplot2 Display

# Requires gridExtra or ggpubr package
plot(model,
  which = 1:4,
  use_ggplot = TRUE,
  arrange_plots = TRUE, # Arrange in grid
  theme_fn = ggplot2::theme_minimal
)

# Focus plots in 2x2 grid
plot(model,
  which = c(2, 3, 4, 6),
  use_ggplot = TRUE,
  arrange_plots = TRUE,
  caption = list(
    "2" = "Influential Cases",
    "3" = "High Leverage Points"
  )
)

# EXAMPLE 7: Half-Normal Plot Customization

# Higher precision envelope (more simulations)
par(mfrow = c(1, 2))
plot(model,
  which = 5,
  nsim = 500, # More accurate envelope
  level = 0.95
) # 95% confidence level

# Quick envelope for large datasets
plot(model,
  which = 5,
  nsim = 500, # Faster computation
  level = 0.90
)

# EXAMPLE 8: Different Residual Types Comparison

# Compare different residual types
par(mfrow = c(2, 2))
plot(model, which = 4, type = "quantile", main = "Quantile")
plot(model, which = 4, type = "pearson", main = "Pearson")
plot(model, which = 4, type = "deviance", main = "Deviance")
par(mfrow = c(1, 1))

# Quantile residuals for half-normal plot (recommended)
plot(model, which = 5, type = "quantile")

# EXAMPLE 9: Family Comparison Diagnostics

# Compare diagnostics under different distributional assumptions
# Helps assess if alternative family would fit better
par(mfrow = c(2, 2))
plot(model,
  which = c(5, 6),
  family = "kw", # Original family
  main = "Kumaraswamy"
)

plot(model,
  which = c(5, 6),
  family = "beta", # Alternative family
  main = "Beta"
)
par(mfrow = c(1, 1))

# EXAMPLE 10: Large Dataset - Performance Optimization

# Simulate large dataset
set.seed(456)
n_large <- 50000
x1_large <- runif(n_large, -2, 2)
x2_large <- rnorm(n_large)
alpha_large <- exp(0.5 + 0.2 * x1_large)
beta_large <- exp(1.0 - 0.1 * x2_large)
y_large <- rkw(n_large, alpha = alpha_large, beta = beta_large)
df_large <- data.frame(y = y_large, x1 = x1_large, x2 = x2_large)

model_large <- gkwreg(y ~ x1 | x2, data = df_large, family = "kw")

# Optimized plotting for large dataset
par(mfrow = c(2, 2), mar = c(3, 3, 2, 2))
plot(model_large,
  which = c(1, 2, 4, 6), # Skip computationally intensive plot 5
  sample_size = 2000, # Use random sample of 2000 observations
  ask = FALSE
) # Don't prompt

# If half-normal plot needed, reduce simulations
par(mfrow = c(1, 1))
plot(model_large,
  which = 5,
  sample_size = 1000, # Smaller sample
  nsim = 50
) # Fewer simulations

# EXAMPLE 11: Saving Diagnostic Data for Custom Analysis

# Extract diagnostic measures without plotting
par(mfrow = c(1, 1))
diag_data <- plot(model_large,
  which = 1:6,
  save_diagnostics = TRUE
)

# Examine structure
str(diag_data)

# Access diagnostic measures
head(diag_data$data) # Residuals, Cook's distance, leverage, etc.

# Identify influential observations
influential <- which(diag_data$data$cook_dist > diag_data$model_info$cook_threshold)
cat("Influential observations:", head(influential), "\n")

# High leverage points
high_lev <- which(diag_data$data$leverage > diag_data$model_info$leverage_threshold)
cat("High leverage points:", head(high_lev), "\n")

# Custom diagnostic plot using saved data
plot(diag_data$data$fitted, diag_data$data$resid,
  xlab = "Fitted Values", ylab = "Residuals",
  main = "Custom Diagnostic Plot",
  col = ifelse(diag_data$data$cook_dist >
    diag_data$model_info$cook_threshold, "red", "black"),
  pch = 16
)
abline(h = 0, col = "gray", lty = 2)
legend("topright", legend = "Influential", col = "red", pch = 16)

# EXAMPLE 12: Interactive Plotting Control

# ask = TRUE Force prompting between plots (useful for presentations)
# Disable prompting (batch processing)
par(mfrow = c(3, 2))
plot(model,
  which = 1:6,
  ask = FALSE
) # Never prompts

# EXAMPLE 13: Base R Graphics Customization via ...

# Customize point appearance
par(mfrow = c(2, 2))
plot(model,
  which = c(1, 4, 6),
  pch = 16, # Filled circles
  col = "steelblue", # Blue points
  cex = 0.8
) # Smaller points

# Multiple customizations
plot(model,
  which = 2,
  pch = 21, # Circles with border
  col = "black", # Border color
  bg = "lightblue", # Fill color
  cex = 1.2, # Larger points
  lwd = 2
) # Thicker lines

# EXAMPLE 14: Comparing Models

# Fit competing models
model_kw <- gkwreg(y ~ x1 | x2, data = df, family = "kw")
model_beta <- gkwreg(y ~ x1 | x2, data = df, family = "beta")

# Compare diagnostics side-by-side
par(mfrow = c(2, 2))

# Kumaraswamy model
plot(model_kw, which = 5, main = "Kumaraswamy - Half-Normal")
plot(model_kw, which = 6, main = "Kumaraswamy - Pred vs Obs")

# Beta model
plot(model_beta, which = 5, main = "Beta - Half-Normal")
plot(model_beta, which = 6, main = "Beta - Pred vs Obs")

par(mfrow = c(1, 1))
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot.gkwreg", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("predict.gkwreg")
### * predict.gkwreg

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: predict.gkwreg
### Title: Predictions from a Fitted Generalized Kumaraswamy Regression
###   Model
### Aliases: predict.gkwreg
### Keywords: models predict regression

### ** Examples

## No test: 
# Generate a sample dataset (n = 1000)
library(gkwdist)
set.seed(123)
n <- 1000

# Create predictors
x1 <- runif(n, -2, 2)
x2 <- rnorm(n)
x3 <- factor(rbinom(n, 1, 0.4))

# Simulate Kumaraswamy distributed data
# True parameters with specific relationships to predictors
true_alpha <- exp(0.7 + 0.3 * x1)
true_beta <- exp(1.2 - 0.2 * x2 + 0.4 * (x3 == "1"))

# Generate random responses
y <- rkw(n, alpha = true_alpha, beta = true_beta)

# Ensure responses are strictly in (0, 1)
y <- pmax(pmin(y, 1 - 1e-7), 1e-7)

# Create data frame
df <- data.frame(y = y, x1 = x1, x2 = x2, x3 = x3)

# Split into training and test sets
set.seed(456)
train_idx <- sample(n, 800)
train_data <- df[train_idx, ]
test_data <- df[-train_idx, ]

# ====================================================================
# Example 1: Basic usage - Fit a Kumaraswamy model and make predictions
# ====================================================================

# Fit the model
kw_model <- gkwreg(y ~ x1 | x2 + x3, data = train_data, family = "kw")

# Predict mean response for test data
pred_mean <- predict(kw_model, newdata = test_data, type = "response")

# Calculate prediction error
mse <- mean((test_data$y - pred_mean)^2)
cat("Mean Squared Error:", mse, "\n")

# ====================================================================
# Example 2: Different prediction types
# ====================================================================

# Create a grid of values for visualization
x1_grid <- seq(-2, 2, length.out = 100)
grid_data <- data.frame(x1 = x1_grid, x2 = 0, x3 = 0)

# Predict different quantities
pred_mean <- predict(kw_model, newdata = grid_data, type = "response")
pred_var <- predict(kw_model, newdata = grid_data, type = "variance")
pred_params <- predict(kw_model, newdata = grid_data, type = "parameter")
pred_alpha <- predict(kw_model, newdata = grid_data, type = "alpha")
pred_beta <- predict(kw_model, newdata = grid_data, type = "beta")

# Plot predicted mean and parameters against x1
plot(x1_grid, pred_mean,
  type = "l", col = "blue",
  xlab = "x1", ylab = "Predicted Mean", main = "Mean Response vs x1"
)
plot(x1_grid, pred_var,
  type = "l", col = "red",
  xlab = "x1", ylab = "Predicted Variance", main = "Response Variance vs x1"
)
plot(x1_grid, pred_alpha,
  type = "l", col = "purple",
  xlab = "x1", ylab = "Alpha", main = "Alpha Parameter vs x1"
)
plot(x1_grid, pred_beta,
  type = "l", col = "green",
  xlab = "x1", ylab = "Beta", main = "Beta Parameter vs x1"
)

# ====================================================================
# Example 3: Computing densities, CDFs, and quantiles
# ====================================================================

# Select a single observation
obs_data <- test_data[1, ]

# Create a sequence of y values for plotting
y_seq <- seq(0.01, 0.99, length.out = 100)

# Compute density at each y value
dens_values <- predict(kw_model,
  newdata = obs_data,
  type = "density", at = y_seq, elementwise = FALSE
)

# Compute CDF at each y value
cdf_values <- predict(kw_model,
  newdata = obs_data,
  type = "probability", at = y_seq, elementwise = FALSE
)

# Compute quantiles for a sequence of probabilities
prob_seq <- seq(0.1, 0.9, by = 0.1)
quant_values <- predict(kw_model,
  newdata = obs_data,
  type = "quantile", at = prob_seq, elementwise = FALSE
)

# Plot density and CDF
plot(y_seq, dens_values,
  type = "l", col = "blue",
  xlab = "y", ylab = "Density", main = "Predicted PDF"
)
plot(y_seq, cdf_values,
  type = "l", col = "red",
  xlab = "y", ylab = "Cumulative Probability", main = "Predicted CDF"
)

# ====================================================================
# Example 4: Prediction under different distributional assumptions
# ====================================================================

# Fit models with different families
beta_model <- gkwreg(y ~ x1 | x2 + x3, data = train_data, family = "beta")
gkw_model <- gkwreg(y ~ x1 | x2 + x3 | 1 | 1 | x3, data = train_data, family = "gkw")

# Predict means using different families
pred_kw <- predict(kw_model, newdata = test_data, type = "response")
pred_beta <- predict(beta_model, newdata = test_data, type = "response")
pred_gkw <- predict(gkw_model, newdata = test_data, type = "response")

# Calculate MSE for each family
mse_kw <- mean((test_data$y - pred_kw)^2)
mse_beta <- mean((test_data$y - pred_beta)^2)
mse_gkw <- mean((test_data$y - pred_gkw)^2)

cat("MSE by family:\n")
cat("Kumaraswamy:", mse_kw, "\n")
cat("Beta:", mse_beta, "\n")
cat("GKw:", mse_gkw, "\n")

# Compare predictions from different families visually
plot(test_data$y, pred_kw,
  col = "blue", pch = 16,
  xlab = "Observed", ylab = "Predicted", main = "Predicted vs Observed"
)
points(test_data$y, pred_beta, col = "red", pch = 17)
points(test_data$y, pred_gkw, col = "green", pch = 18)
abline(0, 1, lty = 2)
legend("topleft",
  legend = c("Kumaraswamy", "Beta", "GKw"),
  col = c("blue", "red", "green"), pch = c(16, 17, 18)
)

# ====================================================================
# Example 5: Working with linear predictors and link functions
# ====================================================================

# Extract linear predictors and parameter values
lp <- predict(kw_model, newdata = test_data, type = "link")
params <- predict(kw_model, newdata = test_data, type = "parameter")

# Verify that inverse link transformation works correctly
# For Kumaraswamy model, alpha and beta use log links by default
alpha_from_lp <- exp(lp$alpha)
beta_from_lp <- exp(lp$beta)

# Compare with direct parameter predictions
cat("Manual inverse link vs direct parameter prediction:\n")
cat("Alpha difference:", max(abs(alpha_from_lp - params$alpha)), "\n")
cat("Beta difference:", max(abs(beta_from_lp - params$beta)), "\n")

# ====================================================================
# Example 6: Elementwise calculations
# ====================================================================

# Generate probabilities specific to each observation
probs <- runif(nrow(test_data), 0.1, 0.9)

# Calculate quantiles for each observation at its own probability level
quant_elementwise <- predict(kw_model,
  newdata = test_data,
  type = "quantile", at = probs, elementwise = TRUE
)

# Calculate probabilities at each observation's actual value
prob_at_y <- predict(kw_model,
  newdata = test_data,
  type = "probability", at = test_data$y, elementwise = TRUE
)

# Create Q-Q plot
plot(sort(prob_at_y), seq(0, 1, length.out = length(prob_at_y)),
  xlab = "Empirical Probability", ylab = "Theoretical Probability",
  main = "P-P Plot", type = "l"
)
abline(0, 1, lty = 2, col = "red")

# ====================================================================
# Example 7: Predicting for the original data
# ====================================================================

# Fit a model with original data
full_model <- gkwreg(y ~ x1 + x2 + x3 | x1 + x2 + x3, data = df, family = "kw")

# Get fitted values using predict and compare with model's fitted.values
fitted_from_predict <- predict(full_model, type = "response")
fitted_from_model <- full_model$fitted.values

# Compare results
cat(
  "Max difference between predict() and fitted.values:",
  max(abs(fitted_from_predict - fitted_from_model)), "\n"
)

# ====================================================================
# Example 8: Handling missing data
# ====================================================================

# Create test data with some missing values
test_missing <- test_data
test_missing$x1[1:5] <- NA
test_missing$x2[6:10] <- NA

# Predict with different na.action options
pred_na_pass <- tryCatch(
  predict(kw_model, newdata = test_missing, na.action = na.pass),
  error = function(e) rep(NA, nrow(test_missing))
)
pred_na_omit <- tryCatch(
  predict(kw_model, newdata = test_missing, na.action = na.omit),
  error = function(e) rep(NA, nrow(test_missing))
)

# Show which positions have NAs
cat("Rows with missing predictors:", which(is.na(pred_na_pass)), "\n")
cat("Length after na.omit:", length(pred_na_omit), "\n")
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("predict.gkwreg", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("print.gkwreg")
### * print.gkwreg

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: print.gkwreg
### Title: Print Method for Generalized Kumaraswamy Regression Models
### Aliases: print.gkwreg

### ** Examples

## No test: 
data(GasolineYield)
fit <- gkwreg(yield ~ batch + temp, data = GasolineYield, family = "kw")
print(fit)

# With more digits
print(fit, digits = 5)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("print.gkwreg", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("residuals.gkwreg")
### * residuals.gkwreg

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: residuals.gkwreg
### Title: Extract Residuals from a Generalized Kumaraswamy Regression
###   Model
### Aliases: residuals.gkwreg
### Keywords: diagnostics methods regression residuals

### ** Examples

## No test: 
require(gkwreg)
require(gkwdist)

# Example 1: Comprehensive residual analysis for FoodExpenditure
data(FoodExpenditure)
FoodExpenditure$prop <- FoodExpenditure$food / FoodExpenditure$income

fit_kw <- gkwreg(
  prop ~ income + persons | income + persons,
  data = FoodExpenditure,
  family = "kw"
)

# Extract different types of residuals
res_response <- residuals(fit_kw, type = "response")
res_pearson <- residuals(fit_kw, type = "pearson")
res_deviance <- residuals(fit_kw, type = "deviance")
res_quantile <- residuals(fit_kw, type = "quantile")
res_coxsnell <- residuals(fit_kw, type = "cox-snell")

# Summary statistics
residual_summary <- data.frame(
  Type = c("Response", "Pearson", "Deviance", "Quantile", "Cox-Snell"),
  Mean = c(
    mean(res_response), mean(res_pearson),
    mean(res_deviance), mean(res_quantile),
    mean(res_coxsnell)
  ),
  SD = c(
    sd(res_response), sd(res_pearson),
    sd(res_deviance), sd(res_quantile),
    sd(res_coxsnell)
  ),
  Min = c(
    min(res_response), min(res_pearson),
    min(res_deviance), min(res_quantile),
    min(res_coxsnell)
  ),
  Max = c(
    max(res_response), max(res_pearson),
    max(res_deviance), max(res_quantile),
    max(res_coxsnell)
  )
)
print(residual_summary)

# Example 2: Diagnostic plots for model assessment
data(GasolineYield)

fit_ekw <- gkwreg(
  yield ~ batch + temp | temp | batch,
  data = GasolineYield,
  family = "ekw"
)

# Set up plotting grid
par(mfrow = c(2, 3))

# Plot 1: Residuals vs Fitted
fitted_vals <- fitted(fit_ekw)
res_pears <- residuals(fit_ekw, type = "pearson")
plot(fitted_vals, res_pears,
  xlab = "Fitted Values", ylab = "Pearson Residuals",
  main = "Residuals vs Fitted",
  pch = 19, col = rgb(0, 0, 1, 0.5)
)
abline(h = 0, col = "red", lwd = 2, lty = 2)
lines(lowess(fitted_vals, res_pears), col = "blue", lwd = 2)

# Plot 2: Normal QQ-plot (Quantile Residuals)
res_quant <- residuals(fit_ekw, type = "quantile")
qqnorm(res_quant,
  main = "Normal Q-Q Plot (Quantile Residuals)",
  pch = 19, col = rgb(0, 0, 1, 0.5)
)
qqline(res_quant, col = "red", lwd = 2)

# Plot 3: Scale-Location (sqrt of standardized residuals)
plot(fitted_vals, sqrt(abs(res_pears)),
  xlab = "Fitted Values", ylab = expression(sqrt("|Std. Residuals|")),
  main = "Scale-Location",
  pch = 19, col = rgb(0, 0, 1, 0.5)
)
lines(lowess(fitted_vals, sqrt(abs(res_pears))), col = "red", lwd = 2)

# Plot 4: Histogram of Quantile Residuals
hist(res_quant,
  breaks = 15, probability = TRUE,
  xlab = "Quantile Residuals",
  main = "Histogram with Normal Overlay",
  col = "lightblue", border = "white"
)
curve(dnorm(x, mean(res_quant), sd(res_quant)),
  add = TRUE, col = "red", lwd = 2
)

# Plot 5: Cox-Snell Residual Plot
res_cs <- residuals(fit_ekw, type = "cox-snell")
plot(qexp(ppoints(length(res_cs))), sort(res_cs),
  xlab = "Theoretical Exponential Quantiles",
  ylab = "Ordered Cox-Snell Residuals",
  main = "Cox-Snell Residual Plot",
  pch = 19, col = rgb(0, 0, 1, 0.5)
)
abline(0, 1, col = "red", lwd = 2)

# Plot 6: Residuals vs Index
plot(seq_along(res_pears), res_pears,
  xlab = "Observation Index", ylab = "Pearson Residuals",
  main = "Residuals vs Index",
  pch = 19, col = rgb(0, 0, 1, 0.5)
)
abline(h = 0, col = "red", lwd = 2, lty = 2)

par(mfrow = c(1, 1))

# Example 3: Partial residual plots for covariate effects
data(ReadingSkills)

fit_interact <- gkwreg(
  accuracy ~ dyslexia * iq | dyslexia + iq,
  data = ReadingSkills,
  family = "kw"
)

# Partial residuals for IQ effect on alpha parameter
X_alpha <- fit_interact$model_matrices$alpha
iq_col_alpha <- which(colnames(X_alpha) == "iq")

if (length(iq_col_alpha) > 0) {
  res_partial_alpha <- residuals(fit_interact,
    type = "partial",
    parameter = "alpha",
    covariate_idx = iq_col_alpha
  )

  par(mfrow = c(1, 2))

  # Partial residual plot for alpha
  plot(ReadingSkills$iq, res_partial_alpha,
    xlab = "IQ (z-scores)",
    ylab = "Partial Residual (alpha)",
    main = "Effect of IQ on Mean (alpha)",
    pch = 19, col = ReadingSkills$dyslexia
  )
  lines(lowess(ReadingSkills$iq, res_partial_alpha),
    col = "blue", lwd = 2
  )
  legend("topleft",
    legend = c("Control", "Dyslexic"),
    col = c("black", "red"), pch = 19
  )

  # Partial residuals for IQ effect on beta parameter
  X_beta <- fit_interact$model_matrices$beta
  iq_col_beta <- which(colnames(X_beta) == "iq")

  if (length(iq_col_beta) > 0) {
    res_partial_beta <- residuals(fit_interact,
      type = "partial",
      parameter = "beta",
      covariate_idx = iq_col_beta
    )

    plot(ReadingSkills$iq, res_partial_beta,
      xlab = "IQ (z-scores)",
      ylab = "Partial Residual (beta)",
      main = "Effect of IQ on Precision (beta)",
      pch = 19, col = ReadingSkills$dyslexia
    )
    lines(lowess(ReadingSkills$iq, res_partial_beta),
      col = "blue", lwd = 2
    )
  }

  par(mfrow = c(1, 1))
}

# Example 4: Comparing residuals across different families
data(StressAnxiety)

fit_kw_stress <- gkwreg(
  anxiety ~ stress | stress,
  data = StressAnxiety,
  family = "kw"
)

# Quantile residuals under different family assumptions
res_quant_kw <- residuals(fit_kw_stress, type = "quantile", family = "kw")
res_quant_beta <- residuals(fit_kw_stress, type = "quantile", family = "beta")

# Compare normality
par(mfrow = c(1, 2))

qqnorm(res_quant_kw,
  main = "QQ-Plot: Kumaraswamy Residuals",
  pch = 19, col = rgb(0, 0, 1, 0.5)
)
qqline(res_quant_kw, col = "red", lwd = 2)

qqnorm(res_quant_beta,
  main = "QQ-Plot: Beta Residuals",
  pch = 19, col = rgb(0, 0.5, 0, 0.5)
)
qqline(res_quant_beta, col = "red", lwd = 2)

par(mfrow = c(1, 1))

# Formal normality tests
shapiro_kw <- shapiro.test(res_quant_kw)
shapiro_beta <- shapiro.test(res_quant_beta)

cat("\nShapiro-Wilk Test Results:\n")
cat(
  "Kumaraswamy:  W =", round(shapiro_kw$statistic, 4),
  ", p-value =", round(shapiro_kw$p.value, 4), "\n"
)
cat(
  "Beta:         W =", round(shapiro_beta$statistic, 4),
  ", p-value =", round(shapiro_beta$p.value, 4), "\n"
)

# Example 5: Outlier detection using standardized residuals
data(MockJurors)

fit_mc <- gkwreg(
  confidence ~ verdict * conflict | verdict + conflict,
  data = MockJurors,
  family = "mc"
)

res_dev <- residuals(fit_mc, type = "deviance")
res_quant <- residuals(fit_mc, type = "quantile")

# Identify potential outliers (|z| > 2.5)
outlier_idx <- which(abs(res_quant) > 2.5)

if (length(outlier_idx) > 0) {
  cat("\nPotential outliers detected at indices:", outlier_idx, "\n")

  # Display outlier information
  outlier_data <- data.frame(
    Index = outlier_idx,
    Confidence = MockJurors$confidence[outlier_idx],
    Verdict = MockJurors$verdict[outlier_idx],
    Conflict = MockJurors$conflict[outlier_idx],
    Quantile_Residual = round(res_quant[outlier_idx], 3),
    Deviance_Residual = round(res_dev[outlier_idx], 3)
  )
  print(outlier_data)

  # Influence plot
  plot(seq_along(res_quant), res_quant,
    xlab = "Observation Index",
    ylab = "Quantile Residual",
    main = "Outlier Detection: Mock Jurors",
    pch = 19, col = rgb(0, 0, 1, 0.5)
  )
  points(outlier_idx, res_quant[outlier_idx],
    col = "red", pch = 19, cex = 1.5
  )
  abline(
    h = c(-2.5, 0, 2.5), col = c("orange", "black", "orange"),
    lty = c(2, 1, 2), lwd = 2
  )
  legend("topright",
    legend = c("Normal", "Outlier", "±2.5 SD"),
    col = c(rgb(0, 0, 1, 0.5), "red", "orange"),
    pch = c(19, 19, NA),
    lty = c(NA, NA, 2),
    lwd = c(NA, NA, 2)
  )
} else {
  cat("\nNo extreme outliers detected (|z| > 2.5)\n")
}
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("residuals.gkwreg", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("response")
### * response

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: response
### Title: Extract Response Variable from GKw Regression Model
### Aliases: response response.gkwreg

### ** Examples

## No test: 
data(GasolineYield)
fit <- gkwreg(yield ~ batch + temp, data = GasolineYield, family = "kw")
y <- response(fit)
head(y)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("response", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("retinal")
### * retinal

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: retinal
### Title: Intraocular Gas Decay in Retinal Surgery
### Aliases: retinal
### Keywords: datasets

### ** Examples

## No test: 
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

# Interpretation:
# - Beta: Precision varies with both level and time
#   Early measurements may be more variable than late measurements

# Test heteroscedasticity
anova(fit_kw, fit_kw_hetero)

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
summary(fit_mc)

# Model comparison
AIC(fit_kw, fit_kw_hetero, fit_ekw, fit_mc)
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("retinal", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("sdac")
### * sdac

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: sdac
### Title: Autologous Peripheral Blood Stem Cell Transplants Data
### Aliases: sdac
### Keywords: datasets

### ** Examples

## No test: 
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

# Interpretation:
# - Gender effects in both mean and precision
# - Precision may differ between males and females

# Test gender significance
anova(fit_kw, fit_kw_gender)

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
summary(fit_ekw)

# Clinical interpretation:
# - Lambda varies by chemotherapy protocol: Some protocols produce more
#   extreme recovery patterns (very high or very low CD34+ counts)
# - G-CSF-only vs multi-day protocols may differ in tail behavior
# - Important for risk stratification and clinical decision-making

# Test if extreme patterns differ by protocol
anova(fit_kw_gender, fit_ekw)

# Example 4: Interaction between age and protocol
# Protocol effectiveness may vary with patient age
fit_kw_interact <- gkwreg(
  rcd ~ ageadj * chemo |
    age * chemo,
  data = sdac,
  family = "kw"
)
summary(fit_kw_interact)

# Interpretation:
# - Interaction: Does protocol effectiveness decline with age?
# - Critical for personalized treatment selection
## End(No test)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("sdac", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("summary.gkwreg")
### * summary.gkwreg

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: summary.gkwreg
### Title: Summary Method for Generalized Kumaraswamy Regression Models
### Aliases: summary.gkwreg
### Keywords: models regression summary

### ** Examples

## No test: 
set.seed(123)
n <- 100
x1 <- runif(n, -2, 2)
x2 <- rnorm(n)
alpha_coef <- c(0.8, 0.3, -0.2)
beta_coef <- c(1.2, -0.4, 0.1)
eta_alpha <- alpha_coef[1] + alpha_coef[2] * x1 + alpha_coef[3] * x2
eta_beta <- beta_coef[1] + beta_coef[2] * x1 + beta_coef[3] * x2
alpha_true <- exp(eta_alpha)
beta_true <- exp(eta_beta)
# Use stats::rbeta as a placeholder if rkw is unavailable
y <- stats::rbeta(n, shape1 = alpha_true, shape2 = beta_true)
y <- pmax(pmin(y, 1 - 1e-7), 1e-7)
df <- data.frame(y = y, x1 = x1, x2 = x2)

# Fit a Kumaraswamy regression model
kw_reg <- gkwreg(y ~ x1 + x2 | x1 + x2, data = df, family = "kw")

# Generate detailed summary using the summary method
summary_kw <- summary(kw_reg)

# Print the summary object (uses print.summary.gkwreg)
print(summary_kw)

# Extract coefficient table directly from the summary object
coef_table <- coef(summary_kw) # Equivalent to summary_kw$coefficients
print(coef_table)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("summary.gkwreg", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("terms.gkwreg")
### * terms.gkwreg

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: terms.gkwreg
### Title: Extract Terms from GKw Regression Model
### Aliases: terms.gkwreg

### ** Examples

## No test: 
data(GasolineYield)
fit <- gkwreg(yield ~ batch + temp, data = GasolineYield, family = "kw")
terms(fit)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("terms.gkwreg", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("update.gkwreg")
### * update.gkwreg

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: update.gkwreg
### Title: Update and Re-fit a GKw Regression Model
### Aliases: update.gkwreg

### ** Examples

## No test: 
# Load example data
require(gkwreg)

data(GasolineYield)

# EXAMPLE 1: Simple formulas (1 part - alpha only)

m1_0 <- gkwreg(yield ~ 1, data = GasolineYield, family = "kw")
m1_1 <- update(m1_0, . ~ . + temp)
m1_2 <- update(m1_1, . ~ . + batch)
m1_3 <- update(m1_2, . ~ . - temp)

anova(m1_0, m1_1, m1_2)
AIC(m1_0, m1_1, m1_2, m1_3)
BIC(m1_0, m1_1, m1_2, m1_3)

# EXAMPLE 2: Two-part formulas (alpha | beta)

# Start with intercept-only for both
m2_0 <- gkwreg(yield ~ 1 | 1, data = GasolineYield, family = "kw")

# Add temp to alpha
m2_1 <- update(m2_0, . ~ . + temp | .)

# Add batch to beta
m2_2 <- update(m2_1, . ~ . | . + batch)

# Add batch to alpha too
m2_3 <- update(m2_2, . ~ . + batch | .)

anova(m2_0, m2_1, m2_2, m2_3)
AIC(m2_0, m2_1, m2_2, m2_3)

# EXAMPLE 3: Three-part formulas (alpha | beta | gamma)

m3_0 <- gkwreg(yield ~ 1,
  data = GasolineYield,
  family = "gkw",
  control = gkw_control(method = "BFGS", maxit = 2000)
)

m3_1 <- update(m3_0, . ~ . + temp | . | .)
m3_2 <- update(m3_1, . ~ . | . + batch | .)
m3_3 <- update(m3_2, . ~ . | . | . + temp)

anova(m3_0, m3_1, m3_2, m3_3)

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
AIC(fit0, fit1, fit2, fit3, fit4)
BIC(fit0, fit1, fit2, fit3, fit4)

# EXAMPLE 5: Changing other parameters

# Change family
fit_gkw <- update(fit2, family = "gkw")

# Change link function
fit_logit <- update(fit2, link = list(alpha = "logit"))

# View call without fitting
update(fit2, . ~ . | . + temp, evaluate = FALSE)
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("update.gkwreg", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
