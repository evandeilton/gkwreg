# Control Parameters for Generalized Kumaraswamy Regression

Auxiliary function for controlling
[`gkwreg()`](https://evandeilton.github.io/gkwreg/reference/gkwreg.md)
fitting process. This function consolidates all technical/advanced
fitting options in one place, keeping the main
[`gkwreg()`](https://evandeilton.github.io/gkwreg/reference/gkwreg.md)
interface clean and user-friendly. Follows the same design pattern as
[`glm.control`](https://rdrr.io/r/stats/glm.control.html),
[`betareg.control`](https://rdrr.io/pkg/betareg/man/betareg.control.html),
and similar control functions in R.

## Usage

``` r
gkw_control(
  method = c("nlminb", "BFGS", "Nelder-Mead", "CG", "SANN", "L-BFGS-B"),
  start = NULL,
  fixed = NULL,
  hessian = TRUE,
  maxit = 500,
  reltol = sqrt(.Machine$double.eps),
  abstol = 0,
  trace = 0,
  silent = TRUE,
  eval.max = 500,
  iter.max = 300,
  step.min = 1e-08,
  step.max = 1,
  x.tol = 1.5e-08,
  rel.tol = sqrt(.Machine$double.eps),
  alpha = 1,
  beta = 0.5,
  gamma = 2,
  warn.1d.NelderMead = TRUE,
  type = 1,
  temp = 10,
  tmax = 10,
  lmm = 5,
  factr = 1e+07,
  pgtol = 0,
  REPORT = NULL,
  fnscale = 1,
  parscale = NULL,
  ndeps = NULL,
  ...
)

# S3 method for class 'gkw_control'
print(x, ...)
```

## Arguments

- method:

  Character string specifying the optimization algorithm. Options:
  `"nlminb"` (default), `"BFGS"`, `"Nelder-Mead"`, `"CG"`, `"SANN"`,
  `"L-BFGS-B"`. If `"nlminb"`, uses
  [`nlminb`](https://rdrr.io/r/stats/nlminb.html); otherwise uses
  [`optim`](https://rdrr.io/r/stats/optim.html) with the specified
  method.

- start:

  Optional named list of starting values for regression coefficients.
  Names should match parameter names (alpha, beta, gamma, delta,
  lambda). If `NULL` (default), starting values are determined
  automatically.

- fixed:

  Optional named list of parameters to hold fixed at specific values
  during estimation. Currently experimental. Default `NULL`.

- hessian:

  Logical. If `TRUE` (default), compute the Hessian matrix via
  [`sdreport`](https://rdrr.io/pkg/TMB/man/sdreport.html) to obtain
  standard errors and variance-covariance matrix. Set to `FALSE` for
  faster fitting when standard errors are not needed.

- maxit:

  Integer. Maximum number of iterations for the optimizer. Default 500
  for derivative-based methods, 10000 for SANN. Increase for difficult
  optimization problems.

- reltol:

  Numeric. Relative convergence tolerance for the optimizer. Default
  `sqrt(.Machine$double.eps)` approx. 1.5e-8. Smaller values require
  tighter convergence but may increase computation time. Used by
  Nelder-Mead, BFGS, and CG methods.

- abstol:

  Numeric. Absolute convergence tolerance. Default 0. Used by some
  optimization methods as an additional stopping criterion.

- trace:

  Integer. Controls verbosity of the optimizer.

  - 0: Silent (default)

  - 1: Print iteration progress

  - 2+: Print detailed diagnostic information (up to 6 for L-BFGS-B)

  Ignored if `silent = TRUE`.

- silent:

  Logical. If `TRUE` (default), suppress all progress messages from TMB
  compilation and optimization. Set to `FALSE` for debugging or to
  monitor long-running fits.

- eval.max:

  Integer. Maximum number of function evaluations (nlminb only).
  Default 500. Increase for difficult optimization problems.

- iter.max:

  Integer. Maximum number of iterations (nlminb only). Default 300.
  Usually less than `eval.max`.

- step.min:

  Numeric. Minimum step length (nlminb only). Default 1e-8. Controls how
  small steps can become before stopping.

- step.max:

  Numeric. Maximum step length (nlminb only). Default 1. Useful for
  preventing overshooting in difficult optimization problems.

- x.tol:

  Numeric. Tolerance for parameter convergence (nlminb only). Default
  1.5e-8. Optimizer stops if parameter changes are smaller than this.

- rel.tol:

  Numeric. Relative tolerance for function value (nlminb only). Default
  `sqrt(.Machine$double.eps)`. Alternative specification of relative
  tolerance.

- alpha:

  Numeric. Reflection factor for Nelder-Mead method. Default 1.0. Only
  used when `method = "Nelder-Mead"`.

- beta:

  Numeric. Contraction factor for Nelder-Mead method. Default 0.5. Only
  used when `method = "Nelder-Mead"`.

- gamma:

  Numeric. Expansion factor for Nelder-Mead method. Default 2.0. Only
  used when `method = "Nelder-Mead"`.

- warn.1d.NelderMead:

  Logical. Whether to warn when Nelder-Mead is used for one-dimensional
  optimization. Default `TRUE`.

- type:

  Integer. Update formula for CG method. Options:

  - 1: Fletcher-Reeves update

  - 2: Polak-Ribiere update

  - 3: Beale-Sorenson update

  Default 1. Only used when `method = "CG"`.

- temp:

  Numeric. Starting temperature for SANN method. Default 10. Only used
  when `method = "SANN"`.

- tmax:

  Integer. Number of function evaluations at each temperature for SANN
  method. Default 10. Only used when `method = "SANN"`.

- lmm:

  Integer. Number of BFGS updates retained in L-BFGS-B method.
  Default 5. Only used when `method = "L-BFGS-B"`.

- factr:

  Numeric. Convergence tolerance factor for L-BFGS-B method. Convergence
  occurs when the reduction in the objective is within this factor of
  the machine tolerance. Default 1e7 (tolerance ~1e-8). Only used when
  `method = "L-BFGS-B"`.

- pgtol:

  Numeric. Tolerance on the projected gradient for L-BFGS-B method.
  Default 0 (check suppressed). Only used when `method = "L-BFGS-B"`.

- REPORT:

  Integer. Frequency of progress reports for BFGS, L-BFGS-B and SANN
  methods when `trace > 0`. Default 10 for BFGS/L-BFGS-B, 100 for SANN.

- fnscale:

  Numeric. Overall scaling to be applied to the function value and
  gradient during optimization. Default 1. If negative, turns the
  problem into a maximization problem.

- parscale:

  Numeric vector. Scaling values for parameters. Optimization is
  performed on par/parscale. Default `rep(1, n_params)`.

- ndeps:

  Numeric vector. Step sizes for finite-difference approximation to the
  gradient. Default 1e-3.

- ...:

  Additional arguments passed to the optimizer. Allows fine-grained
  control without formally adding parameters. Advanced users only.

- x:

  An object of class `"gkw_control"`.

## Value

An object of class `"gkw_control"`, which is a list containing all
control parameters with validated and default-filled values. This object
is passed to
[`gkwreg()`](https://evandeilton.github.io/gkwreg/reference/gkwreg.md)
via the `control` argument.

## Details

This function provides a centralized way to set all technical parameters
for model fitting. It serves several purposes:

- **Clean interface**:
  [`gkwreg()`](https://evandeilton.github.io/gkwreg/reference/gkwreg.md)
  has fewer arguments

- **Organized documentation**: All technical options documented here

- **Input validation**: Parameters validated before fitting

- **Extensibility**: New options can be added without changing
  [`gkwreg()`](https://evandeilton.github.io/gkwreg/reference/gkwreg.md)

- **Backward compatibility**: Old code continues working

**Method-specific parameters:**

Each optimization method accepts different control parameters:

- **Nelder-Mead**: `alpha`, `beta`, `gamma`, `maxit`, `reltol`,
  `abstol`, `trace`, `REPORT`, `warn.1d.NelderMead`

- **BFGS**: `maxit`, `reltol`, `abstol`, `trace`, `REPORT`

- **CG**: `type`, `maxit`, `reltol`, `abstol`, `trace`

- **SANN**: `temp`, `tmax`, `maxit`, `trace`, `REPORT`

- **L-BFGS-B**: `lmm`, `factr`, `pgtol`, `trace`, `REPORT`

**When to use gkw_control():**

Most users never need to adjust these settings. Use `gkw_control()`
when:

- Optimization fails to converge (increase `maxit`, adjust tolerances)

- Debugging fit problems (set `silent = FALSE`, `trace = 1`)

- Comparing optimizers (try `method = "BFGS"` vs `"nlminb"`)

- Fine-tuning performance (disable `hessian` if SEs not needed)

- Using custom starting values (`start = list(...)`)

**Recommended practices:**

- Start with defaults, only adjust if needed

- Increase `maxit` before adjusting tolerances

- Use `trace = 1` to diagnose convergence issues

- Disable `hessian` for speed if only point estimates needed

- Try different `method`s if one fails (BFGS often more robust)

- For L-BFGS-B with bounds, adjust `factr` and `pgtol` if needed

## References

Nocedal, J., & Wright, S. J. (2006). *Numerical Optimization* (2nd ed.).
Springer.

Belisle, C. J. P. (1992). Convergence theorems for a class of simulated
annealing algorithms on R^d. *Journal of Applied Probability*, 29,
885-895.

Byrd, R. H., Lu, P., Nocedal, J. and Zhu, C. (1995). A limited memory
algorithm for bound constrained optimization. *SIAM Journal on
Scientific Computing*, 16, 1190-1208.

## See also

[`gkwreg`](https://evandeilton.github.io/gkwreg/reference/gkwreg.md) for
the main fitting function,
[`nlminb`](https://rdrr.io/r/stats/nlminb.html),
[`optim`](https://rdrr.io/r/stats/optim.html) for optimizer details,
[`betareg.control`](https://rdrr.io/pkg/betareg/man/betareg.control.html)
for similar design pattern.

## Author

Lopes, J. E.

## Examples

``` r
# \donttest{
# Default control (used automatically if not specified)
ctrl <- gkw_control()
print(ctrl)
#> Generalized Kumaraswamy Control Parameters
#> ===========================================
#> 
#> Optimization:
#>   Method:             nlminb 
#>   Max evaluations:    500 
#>   Max iterations:     300 
#>   Relative tolerance: 1.490116e-08 
#>   Absolute tolerance: 0e+00 
#> 
#> Output:
#>   Compute Hessian:    TRUE 
#>   Silent mode:        TRUE 

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
# }
```
