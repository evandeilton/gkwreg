# Clean TMB Cache

Removes all cached compiled TMB models. Useful for troubleshooting or
after package updates.

## Usage

``` r
clean_tmb_cache(pkg_name = "gkwreg", verbose = TRUE)
```

## Arguments

- pkg_name:

  Package name (default: "gkwreg")

- verbose:

  Logical

## Value

Invisibly returns number of files deleted
