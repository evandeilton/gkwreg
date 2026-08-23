## Submission summary

This is a maintenance update of **gkwreg**, currently on CRAN as version 2.1.14
(published 2026-01-09). Version 2.1.17 contains no changes to the statistical
methods and no changes to the exported API.

Note on the version number: 2.1.16 was used for the JOSS review archive on
Zenodo and was never submitted to CRAN, so this update goes from 2.1.14 to
2.1.17.

Changes relevant to CRAN:

* `utils` and `grDevices` are now declared in `Imports`. Both were already used
  via `::` in the package code without being declared.
* The suggested package **betareg** is now used conditionally
  (`requireNamespace()` in the vignette, `skip_if_not_installed()` in the tests),
  as required by the policy on packages listed in `Suggests`.
* A knitr cache directory that had been committed under `vignettes/` no longer
  reaches the tarball, so the vignette is always rebuilt from scratch.
  `.Rbuildignore` was tightened accordingly.
* The accompanying paper has been published in the *Journal of Open Source
  Software*. The reference `<doi:10.21105/joss.08991>` was added to the
  `Description` field, and an `inst/CITATION` file was added so that
  `citation("gkwreg")` returns the peer-reviewed reference.
* Restored genuine UTF-8 characters in several help pages, where accented author
  names and mathematical symbols had regressed to literal `<U+XXXX>` escape
  sequences, and normalised all sources to LF line endings.

## Test environments

* Local: Linux x86_64 (Ubuntu-based), R 4.6.1 -- `R CMD check --as-cran`
* <!-- TODO before submitting: run win-builder (R-devel and R-release) and the
  GitHub Actions matrix (Ubuntu R-release/R-devel/R-oldrel-1, macOS R-release,
  Windows R-release), then list them here. Remove this comment. -->

## R CMD check results

0 errors | 0 warnings | 1 note

The NOTE concerns the installed size, which is expected for a package with
compiled C++/TMB code and one HTML vignette.

## Vignette check time

The single vignette runs a Monte Carlo study and dominates the check time
(roughly 6-12 minutes on the CRAN Linux flavours for the previous version). Its
content is unchanged in this release. Please let us know if a shorter vignette
is preferred, and we will pre-compute the simulation results.

## Downstream dependencies

There are no reverse dependencies on CRAN.
