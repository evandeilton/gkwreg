# Check and Compile TMB Model Code with Persistent Cache

This utility function verifies whether the TMB model shared object
(`.so/.dll`) file has already been compiled for a specified DLL name. If
not, it compiles the corresponding C++ file and caches it in a
persistent directory across R sessions.

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

  A character string specifying the base name of the C++ file and the
  resulting DLL. The function assumes the code file is `dll_name.cpp`
  located in the `inst/tmb/` directory of the package.

- pkg_name:

  A character string specifying the package name. Defaults to "gkwreg".

- force_recompile:

  Logical; if `TRUE`, forces recompilation even if a valid compiled file
  exists (default is `FALSE`).

- verbose:

  Logical; if `TRUE`, prints detailed status messages (default is
  `TRUE`).

## Value

Returns (invisibly) a list with information about the compiled model,
including path, normalized path, name, and compilation status. If any
step fails, an error is thrown.

## Details

The function works through the following steps:

1.  Creates a persistent cache directory for storing compiled TMB
    models.

2.  Checks if a compiled file for the specified DLL already exists in
    the cache directory and whether it's up-to-date compared to the
    source code.

3.  If a valid compiled file exists, it loads it directly.

4.  If not, the function locates the corresponding C++ file inside the
    package, compiles it, and stores the result in the cache directory.

5.  Provides diagnostic messages regarding compilation status and
    exported symbols.
