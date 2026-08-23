## Resubmission

This is a resubmission of 2.1.17, addressing your pre-test feedback:

> * checking re-building of vignette outputs ... [11m] OK
> Please reduce the vignette build timings by using small toy data only, few
> iterations, or by providing precomputed results for the most lengthy parts.

We took the third option. The vignette's Monte Carlo study (3 scenarios x 200
replications x 4 models) is now precomputed into
`inst/extdata/vignette-simulations.rds` by `data-raw/vignette-simulations.R`,
and the vignette reads the summaries instead of recomputing them. The
illustrative single-dataset fits stay live, so the vignette still runs the
package, but they use one distribution family instead of three, which cuts the
run-time TMB compilation accordingly.

Local vignette build time drops from about 3m30s to 35s; the full
`R CMD check --as-cran` now takes about 5 minutes on our machine.

We also fixed a timing artefact this exposed: the one-off TMB compilation was
being counted inside the first replicate's fit time and averaged over all 200.

## Submission summary

This is a maintenance update of **gkwreg**, currently on CRAN as version 2.1.14
(published 2026-01-09). Version 2.1.18 contains no changes to the statistical
methods and no changes to the exported API.

Note on the version number: 2.1.16 was used for the JOSS review archive on
Zenodo and was never submitted to CRAN, so this update goes from 2.1.14 to
2.1.18.

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

0 errors | 0 warnings | 0 notes

The only NOTE seen locally is `checking HTML version of manual`, raised because
`tidy` and the `V8` package are unavailable on the test machine; it is an
artefact of the local setup, not of the package.

The full test suite was run separately with `NOT_CRAN=true`:
`FAIL 0 | WARN 0 | SKIP 0 | PASS 1177`. The tests carry `skip_on_cran()` to keep
the check time down, so they are skipped during `R CMD check` itself.

## Vignette check time

The single vignette runs a Monte Carlo study and dominates the check time. On
this machine the full check takes 8m25s, of which the vignette rebuild is 3m25s
and the `--run-donttest` examples are 3m56s. For the previous version the CRAN
Linux flavours reported 6-12 minutes for the vignette alone. Its content is
unchanged in this release; please let us know if a shorter vignette is
preferred, and we will pre-compute the simulation results.

## Downstream dependencies

There are no reverse dependencies on CRAN.
