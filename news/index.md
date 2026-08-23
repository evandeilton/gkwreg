# Changelog

## gkwreg 2.1.18

Resubmission addressing the CRAN pre-test feedback on 2.1.17: the
vignette took about 11 minutes to rebuild, which is more than CRAN can
afford to run regularly.

### Vignette build time

- The Monte Carlo study is now precomputed. The study itself is
  unchanged – still 3 scenarios x 200 replications x 4 models – but the
  summaries are stored in `inst/extdata/vignette-simulations.rds` and
  read by the vignette instead of being recomputed on every build. The
  generating script is `data-raw/vignette-simulations.R`.
- The illustrative model fits stay live, so the vignette still exercises
  the package; they use a single family, which cuts the run-time TMB
  compilation from three models to one.
- Dropped `cache = TRUE` from the vignette: with the study precomputed
  it buys nothing, and it was what produced the stray cache directory
  removed in 2.1.17.

Local vignette build time drops from about 3m30s to 35s.

### Fixed

- Corrected a measurement artefact in the simulation study. Timings were
  taken around each fitting call, so the one-off TMB compilation (~30s
  per family) fell inside the first replicate and was then averaged over
  all 200, inflating every reported `gkwreg` timing by roughly 0.15s.
  Scenario 1 consequently reported Kumaraswamy as slower than `betareg`
  while the text claimed it was faster. The models are now compiled
  before the timed loops.
- Figures quoted in the vignette prose (speed-ups, convergence rates,
  AIC) are now computed inline from the result tables rather than
  written out by hand, so the narrative cannot drift away from the
  results it describes.

------------------------------------------------------------------------

## gkwreg 2.1.17

Maintenance release. No changes to the statistical methods or to the
exported API.

Version 2.1.16 was used for the JOSS review archive on Zenodo and was
never released on CRAN; 2.1.17 is the CRAN update that follows 2.1.14.

### Publication

- The methodology and software are now published in the *Journal of Open
  Source Software*: Lopes and Bonat (2026),
  <https://doi.org/10.21105/joss.08991>.
- Added `inst/CITATION`, so that `citation("gkwreg")` returns the
  peer-reviewed reference, and added that reference to the `Description`
  field.

### Dependencies

- `utils` and `grDevices` are now declared in `Imports`. Both were
  already used via `::`
  ([`utils::modifyList`](https://rdrr.io/r/utils/modifyList.html),
  [`utils::globalVariables`](https://rdrr.io/r/utils/globalVariables.html),
  [`grDevices::dev.hold`](https://rdrr.io/r/grDevices/dev.flush.html),
  [`grDevices::devAskNewPage`](https://rdrr.io/r/grDevices/devAskNewPage.html))
  without being declared.
- The suggested package **betareg** is now used conditionally in the
  vignette and in the comparative test file, as required by the CRAN
  policy on packages listed in `Suggests`.

### Fixed

- Restored genuine UTF-8 characters in the documentation and examples.
  Accented author names and mathematical symbols had regressed to
  literal `<U+XXXX>` escape sequences, affecting the `LossAversion`,
  `ReadingSkills`, `gkwreg`, `anova.gkwreg` and `residuals.gkwreg` help
  pages, the README and the vignette.
- Normalised `DESCRIPTION`, `LICENSE`, `NAMESPACE` and all package
  sources to LF line endings.
- Restored `LICENSE` to the two-line DCF stub that
  `License: MIT + file LICENSE` requires. It had been replaced by the
  full MIT text, which `R CMD check` reports as “License stub is invalid
  DCF”. The full text remains in `LICENSE.md` for GitHub.
- Fixed a typo in
  [`utils::globalVariables()`](https://rdrr.io/r/utils/globalVariables.html),
  where `"dkw dmc"` was a single string instead of two separate entries.

### Packaging

- The knitr cache directory of the vignette is no longer under version
  control and no longer reaches the source tarball; the vignette is
  always rebuilt from scratch.
- Tightened `.Rbuildignore` and `.gitignore` so that build artefacts,
  session files and cache directories cannot leak into the tarball.

### Documentation

- README: added the JOSS badge, removed a stale BibTeX block that
  advertised an outdated version and omitted the second author, and
  fixed the author footer, which was being rendered as a broken table.
- Reworked the `pkgdown` site (light theme, KaTeX math rendering).

------------------------------------------------------------------------

## gkwreg 2.1.14

CRAN release: 2026-01-09

### Fixed

- **clang-san runtime error (integer overflow).** Fixed a
  `static_cast<int>` overflow in the cache-key generation used by the
  TMB models. A `safe_int_cast()` helper now prevents undefined
  behaviour when distribution parameters reach extreme values during
  optimisation. Affects `gkwreg.cpp`, `bkwreg.cpp`, `kkwreg.cpp`,
  `ekwreg.cpp`, `mcreg.cpp` and `kwreg.cpp`.

------------------------------------------------------------------------

## gkwreg 2.1.13

### CRAN Resubmission

Addresses all remaining issues for CRAN acceptance after archival on
2025-11-30.

#### Fixed

- Added `inst/WORDLIST` with ‘Kumaraswamy’ to resolve spelling NOTE
- Added `skip_on_cran()` to all test files to reduce check time from
  18min to ~9min
- Added `cran-comments.md` documenting changes since archival (excluded
  from build via .Rbuildignore)

#### Confirmed

- Cache policy now fully compliant: uses only
  [`tempdir()`](https://rdrr.io/r/base/tempfile.html) for
  session-specific TMB DLL cache
- All `~/.cache/gkwreg` usage completely removed (fixed in v2.1.11)
- Check time now under 10 minutes

------------------------------------------------------------------------

## gkwreg 2.1.12

### CRAN Compliance Fixes

This release addresses all issues that led to package archival on
2025-11-30.

#### Fixed

- **CRITICAL**: Removed `RcppArmadillo` and `RcppEigen` from `Imports`
  field in DESCRIPTION (they remain in `LinkingTo` only, as they are
  used solely for C++ compilation)
- Added single quotes around technical term ‘Kumaraswamy’ in DESCRIPTION
- Removed unnecessary `@import RcppArmadillo` roxygen directive
- Regenerated NAMESPACE to reflect correct imports

#### Note

The cache policy violation (`~/.cache/gkwreg`) reported in version 2.1.6
was already fixed in version 2.1.11. The package now uses ONLY
[`tempdir()`](https://rdrr.io/r/base/tempfile.html) for session-specific
temporary cache, fully complying with CRAN policies.

## gkwreg 2.1.4

#### Minors and documentation

- Fixed the link to the LICENSE.md file.
- Corrected the bold-faced sentences in paper.md.
- Added contribution guidelines.
- Fix README.md equations and other small text issues

Other changes:

- Added new subsections: Distributional Regression Framework and Model
  Diagnostics in paper.md. This makes the paper more comprehensive.
- Removed mentions of the distribution family, as those implementations
  are now in the new package “gkwdisst”.

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
