# Changelog

## gkwreg 2.1.1

CRAN release: 2025-11-15

## gkwreg 2.1.0

#### Comparative Testing

- **Introduced a dedicated comparative test suite** to validate
  `gkwreg`’s beta family implementation against the reference `betareg`
  package, ensuring numerical accuracy and reliability.

- **Confirmed statistical equivalence** despite different internal
  parameterizations. Tests demonstrate that `gkwreg`’s shape-based
  (`gamma`, `delta+1`) approach produces equivalent statistical models
  to `betareg`’s mean-precision (`mu`, `phi`) approach.

- **Validated key outputs**, showing that log-likelihood, AIC/BIC,
  fitted values, and predictions are virtually identical between the two
  packages when fitting the same beta regression model.

- **Successfully benchmarked `gkwreg` with `family = "beta"`** as a
  robust and reliable alternative for beta regression, yielding the same
  inferential conclusions as the established `betareg` package.

- **Verified consistency across multiple scenarios**, including
  controlled simulations with known parameters and real-world datasets
  (`GasolineYield`, `FoodExpenditure`), ensuring robust performance in
  diverse modeling contexts.

## gkwreg 2.0.0

### Major Changes

#### Package Restructuring

- **Complete package reformulation** following JOSS reviewer feedback to
  reduce complexity and improve maintainability.

- **Distribution functions moved to separate package** `gkwdist`: All
  `d*`, `p*`, `q*`, `r*` density/CDF/quantile/random generation
  functions have been extracted to the companion package `gkwdist` for
  cleaner namespace and reduced dependencies. The `gkwreg` package now
  focuses exclusively on regression modeling.

- **Univariate fitting functions removed**: `gkwfit()`, `gkwgof()`, and
  `gkwfitall()` have been removed to maintain package focus on
  regression. Users needing univariate distribution fitting should use
  the `gkwdist` package directly or standard MLE approaches.

#### Simplified Interface

- **Introduced
  [`gkw_control()`](https://evandeilton.github.io/gkwreg/reference/gkw_control.md)**:
  All technical/optimization parameters (method, start, fixed, hessian,
  maxit, tolerances, etc.) are now consolidated in a dedicated control
  function following the
  [`glm.control()`](https://rdrr.io/r/stats/glm.control.html) design
  pattern. This dramatically simplifies the main
  [`gkwreg()`](https://evandeilton.github.io/gkwreg/reference/gkwreg.md)
  interface.

- **Removed arguments violating separation of concerns** from
  [`gkwreg()`](https://evandeilton.github.io/gkwreg/reference/gkwreg.md):

  - `plot` argument removed (use
    [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method
    instead)
  - `conf.level` argument removed (use
    [`confint()`](https://rdrr.io/r/stats/confint.html) method instead)
  - `profile`, `submodels`, `npoints` arguments removed (focused
    functionality)

- **Streamlined
  [`gkwreg()`](https://evandeilton.github.io/gkwreg/reference/gkwreg.md)
  signature**: Reduced from 15+ arguments to ~12 core arguments, with
  technical options delegated to `control`.

#### Complete S3 Method Implementation

- **Standard methods suite**: Implemented complete S3 methods following
  R conventions:
  - [`print.gkwreg()`](https://evandeilton.github.io/gkwreg/reference/print.gkwreg.md),
    [`summary.gkwreg()`](https://evandeilton.github.io/gkwreg/reference/summary.gkwreg.md),
    [`print.summary.gkwreg()`](https://evandeilton.github.io/gkwreg/reference/print.summary.gkwreg.md)
  - [`coef.gkwreg()`](https://evandeilton.github.io/gkwreg/reference/coef.gkwreg.md),
    [`vcov.gkwreg()`](https://evandeilton.github.io/gkwreg/reference/vcov.gkwreg.md),
    [`fitted.gkwreg()`](https://evandeilton.github.io/gkwreg/reference/fitted.gkwreg.md)
  - [`logLik.gkwreg()`](https://evandeilton.github.io/gkwreg/reference/logLik.gkwreg.md),
    [`AIC.gkwreg()`](https://evandeilton.github.io/gkwreg/reference/AIC.gkwreg.md),
    [`BIC.gkwreg()`](https://evandeilton.github.io/gkwreg/reference/BIC.gkwreg.md),
    [`nobs.gkwreg()`](https://evandeilton.github.io/gkwreg/reference/nobs.gkwreg.md)
  - [`confint.gkwreg()`](https://evandeilton.github.io/gkwreg/reference/confint.gkwreg.md),
    [`residuals.gkwreg()`](https://evandeilton.github.io/gkwreg/reference/residuals.gkwreg.md),
    [`predict.gkwreg()`](https://evandeilton.github.io/gkwreg/reference/predict.gkwreg.md)
  - [`anova.gkwreg()`](https://evandeilton.github.io/gkwreg/reference/anova.gkwreg.md),
    [`print.anova.gkwreg()`](https://evandeilton.github.io/gkwreg/reference/print.anova.gkwreg.md),
    [`lrtest()`](https://evandeilton.github.io/gkwreg/reference/lrtest.md)

#### Enhanced Diagnostics

- **Comprehensive
  [`plot.gkwreg()`](https://evandeilton.github.io/gkwreg/reference/plot.gkwreg.md)
  method** with 6 diagnostic plot types:

  1.  Residuals vs Observation Indices
  2.  Cook’s Distance
  3.  Generalized Leverage vs Fitted Values
  4.  Residuals vs Linear Predictor
  5.  Half-Normal Plot with Simulated Envelope
  6.  Predicted vs Observed Values

- **Dual graphics system support**: Base R graphics (default) or ggplot2
  with automatic grid arrangement via `gridExtra`/`ggpubr`.

- **Advanced customization**: Named-list interface for plot captions
  (partial customization without repeating all titles), theme control,
  sampling for large datasets.

#### Powerful Prediction

- **[`predict.gkwreg()`](https://evandeilton.github.io/gkwreg/reference/predict.gkwreg.md)
  with 9 prediction types**:
  - `"response"`, `"variance"`, `"link"`, `"parameter"`
  - Individual parameters: `"alpha"`, `"beta"`, `"gamma"`, `"delta"`,
    `"lambda"`
  - Distribution functions: `"density"`, `"probability"`, `"quantile"`
- **Element-wise and vectorized modes**: Flexible evaluation via
  `elementwise` argument for CDF/PDF/quantile calculations.

#### Model Comparison Tools

- **Likelihood ratio tests**:
  [`anova.gkwreg()`](https://evandeilton.github.io/gkwreg/reference/anova.gkwreg.md)
  for comparing nested models with automatic ordering and chi-squared
  tests; dedicated
  [`lrtest()`](https://evandeilton.github.io/gkwreg/reference/lrtest.md)
  function for pairwise comparisons.

- **Information criteria**:
  [`AIC.gkwreg()`](https://evandeilton.github.io/gkwreg/reference/AIC.gkwreg.md)
  and
  [`BIC.gkwreg()`](https://evandeilton.github.io/gkwreg/reference/BIC.gkwreg.md)
  with multi-model comparison support returning data frames.

#### Documentation Improvements

- **Extensive Roxygen documentation** for all exported functions with
  detailed examples, mathematical formulas, and usage guidance.

- **Updated README.md** with comprehensive feature overview, quick start
  guide, advanced examples, and ecosystem comparison table.

- **NULL default intelligent behavior**: Several arguments default to
  `NULL` triggering smart auto-configuration (e.g., `sub.caption`,
  `ask`, `theme_fn` in
  [`plot.gkwreg()`](https://evandeilton.github.io/gkwreg/reference/plot.gkwreg.md)).

### Testing Framework

#### Comprehensive Test Suite Added

The package now includes a robust testing framework with **1000+ unit
tests** covering all major functionalities:

##### Core Function Testing

- **[`gkwreg()`](https://evandeilton.github.io/gkwreg/reference/gkwreg.md)**:
  20 tests for model fitting, parameter estimation, formula handling,
  all distribution families, link functions, and convergence
- **[`predict.gkwreg()`](https://evandeilton.github.io/gkwreg/reference/predict.gkwreg.md)**:
  10 tests for predictions, including response means, densities, CDFs,
  quantiles, and parameter extraction
- **[`residuals.gkwreg()`](https://evandeilton.github.io/gkwreg/reference/residuals.gkwreg.md)**:
  10 tests for all residual types (response, Pearson, deviance,
  quantile, standardized, working, partial)
- **[`fitted.gkwreg()`](https://evandeilton.github.io/gkwreg/reference/fitted.gkwreg.md)**:
  10 tests for fitted value extraction and validation

##### S3 Methods Testing

- **[`anova.gkwreg()`](https://evandeilton.github.io/gkwreg/reference/anova.gkwreg.md)**:
  45 tests for model comparisons, likelihood ratio tests, and nested
  model hierarchies
- **Print methods**: Tests for
  [`print.gkwreg()`](https://evandeilton.github.io/gkwreg/reference/print.gkwreg.md)
  and
  [`print.summary.gkwreg()`](https://evandeilton.github.io/gkwreg/reference/print.summary.gkwreg.md)
- **Accessor methods**: Tests for
  [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`nobs()`](https://rdrr.io/r/stats/nobs.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html)
- **Summary method**: Tests for
  [`summary.gkwreg()`](https://evandeilton.github.io/gkwreg/reference/summary.gkwreg.md)
  including coefficient tables, confidence intervals, and fit statistics

##### Test Coverage Includes

- ll 7 distribution families (GKw, BKw, KKw, EKw, MC, Kw, Beta)
- Different link functions and scales
- Edge cases and boundary conditions
- Missing data handling (NA)
- Subset and weight specifications
- Large dataset performance
- Error handling and input validation
- Statistical correctness verification
- Numerical accuracy checks

##### Testing Framework

- Built with `testthat` package
- Uses simulated data from `gkwdist` package
- Tests with real datasets (GasolineYield, FoodExpenditure)
- Reproducible with fixed random seeds

### Minor Improvements

- **Link scaling support**: Added `link_scale` argument to
  [`gkwreg()`](https://evandeilton.github.io/gkwreg/reference/gkwreg.md)
  for controlling transformation intensity.

- **Performance optimizations**: Intelligent caching, sampling support
  for diagnostics on large datasets, optional Hessian computation.

------------------------------------------------------------------------

**Breaking Changes**: Version 2.0.0 introduces breaking changes. Code
using `gkwfit()`, `gkwgof()`, `gkwfitall()`, or distribution functions
([`dgkw()`](https://evandeilton.github.io/gkwdist/reference/dgkw.html),
etc.) must be updated to use the `gkwdist` package or the new
[`gkwreg()`](https://evandeilton.github.io/gkwreg/reference/gkwreg.md)
interface with
[`gkw_control()`](https://evandeilton.github.io/gkwreg/reference/gkw_control.md).
