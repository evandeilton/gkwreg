# Get Cache Directory for Compiled Models

Creates version-specific cache directory for TMB DLLs. Structure:
\<cache_root\>/\<pkg_name\>/tmb_cache/R/\<pkg_version\>/

## Usage

``` r
.get_cache_dir(pkg_name)
```

## Arguments

- pkg_name:

  Package name

## Value

Path to cache directory
