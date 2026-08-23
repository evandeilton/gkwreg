## Submission summary

This is a maintenance update of **gkwreg**, currently on CRAN as version 2.1.14
(published 2026-01-09). The previous version passes on all 13 CRAN check
flavours with no ERRORs, WARNINGs or NOTEs.

Version 2.1.16 contains no changes to the statistical methods or to the
exported API. It is a documentation and packaging release:

* The accompanying paper has been published in the *Journal of Open Source
  Software*. The reference `<doi:10.21105/joss.08991>` was added to the
  `Description` field and an `inst/CITATION` file was added, so that
  `citation("gkwreg")` now returns the peer-reviewed reference.
* Restored genuine UTF-8 characters in several help pages, where accented
  author names and mathematical symbols had regressed to literal `<U+XXXX>`
  escape sequences.
* Normalised all package sources to LF line endings.
* Removed a stray knitr cache directory that had been committed under
  `vignettes/`, and tightened `.Rbuildignore` so build artefacts cannot leak
  into the tarball.
* The vignette now uses the suggested package **betareg** conditionally
  (`requireNamespace()`), as required by the CRAN policy on suggested packages.

## Test environments

* Local: Linux x86_64, R 4.6.1
* GitHub Actions: Ubuntu (R-release, R-devel, R-oldrel), macOS (R-release),
  Windows (R-release)
* win-builder: R-devel, R-release

## R CMD check results

0 errors | 0 warnings | 1 note

The remaining NOTE concerns the installed size, which is expected for a package
with compiled C++/TMB code and one vignette.

## Downstream dependencies

There are no reverse dependencies on CRAN.
