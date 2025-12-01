# Get Cache Directory for Compiled Models

Creates version-specific cache directory for TMB DLLs within tempdir().
Uses temporary directory to comply with CRAN policies (no user directory
access). Cache is session-specific and will not persist across R
sessions.

## Usage

``` r
.get_cache_dir(pkg_name)
```

## Arguments

- pkg_name:

  Package name

## Value

Path to cache directory within tempdir()
