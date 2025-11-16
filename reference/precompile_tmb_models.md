# Precompile All TMB Models

Compiles and caches all TMB models found in inst/tmb/. Useful to run
after package installation to avoid compilation delays during first use.

## Usage

``` r
precompile_tmb_models(pkg_name = "gkwreg", verbose = TRUE)
```

## Arguments

- pkg_name:

  Package name (default: "gkwreg")

- verbose:

  Logical

## Value

Invisibly returns list of compilation results
