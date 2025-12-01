# Check and Compile TMB Model Code

Ensures a TMB model is compiled and loaded. Uses session-specific cache
within tempdir() to avoid recompilation during the same R session.

## Usage

``` r
.check_and_compile_TMB_code(
  dll_name,
  pkg_name = "gkwreg",
  force_recompile = FALSE,
  verbose = FALSE
)
```

## Arguments

- dll_name:

  Character, base name of the C++ file (e.g., "kwreg", "ekwreg")

- pkg_name:

  Character, package name (default: "gkwreg")

- force_recompile:

  Logical, force recompilation (default: FALSE)

- verbose:

  Logical, print diagnostic messages (default: FALSE)

## Value

Invisibly returns a list with DLL information
