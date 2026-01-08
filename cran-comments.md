## Test environments

* Local: Linux (Zorin OS 18), R 4.5.2
* win-builder: R-devel (pending)
* macOS (GitHub Actions): R-release (pending)

## R CMD check results

0 errors | 0 warnings | 2 notes

Notes:
- Installed size (5.6Mb): expected for package with C++ code and vignettes
- New submission after archival

## Resubmission

This is a resubmission. Package was archived on 2025-11-30 for policy violation 
regarding creation of files in `~/.cache/gkwreg`.

### Changes since last submission (v2.1.13 → v2.1.14):

1. **FIXED**: clang-san runtime error (integer overflow)
   - Fixed `static_cast<int>` overflow in cache key generation for TMB models
   - Added `safe_int_cast()` function to prevent undefined behavior when 
     distribution parameters reach extreme values during optimization
   - Affects: ekwreg.cpp, kwreg.cpp, mcreg.cpp, kkwreg.cpp, bkwreg.cpp, gkwreg.cpp

### Changes from archival (v2.1.12 → v2.1.13):

1. **FIXED**: Cache policy violation
   - Completely removed use of `~/.cache/` and similar user directories
   - Now uses ONLY `tempdir()` for session-specific TMB DLL caching
   - Cache is automatically cleaned on R session exit
   - See `R/gkwreg-compile-tmb.R`, function `.get_cache_dir()` (line 287)

2. **FIXED**: RcppArmadillo and RcppEigen removed from Imports 
   - These packages are now only in LinkingTo (correct for header-only usage)
   
3. **FIXED**: Reduced example and test execution time
   - Total check time now under 10 minutes

4. **NEW**: Added inst/WORDLIST with 'Kumaraswamy' 
   - Technical term for the distribution family

## Downstream dependencies

No strong reverse dependencies.
