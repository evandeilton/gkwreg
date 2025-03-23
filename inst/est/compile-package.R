# usethis::use_package("ggplot2")
# usethis::use_mit_license()
# usethis::use_git()
# usethis::use_rcpp()
# usethis::use_namespace()
# usethis::use_readme_rmd()
# usethis::use_vignette("gkwreg-package", "Generalized Kumaraswamy Regression Models for Bounded Data")

# Defina a flag de roxygen2 antes de documentar
# usethis::use_version()
options(roxygen2.running = TRUE)

## -------------------------------------------------------------------------- ##
# usethis::use_version(which = "dev")
devtools::clean_dll()
styler::style_pkg()
devtools::load_all(recompile = TRUE, quiet = FALSE)
Rcpp::compileAttributes()
roxygen2::roxygenize()
devtools::document()
# devtools::check()

# usethis::use_version()
devtools::install(force = TRUE)


## -------------------------------------------------------------------------- ##
devtools::build()


devtools::build_manual()
# devtools::build_readme()
# devtools::build_vignettes()

## -------------------------------------------------------------------------- ##
# pkgdown::init_site()
pkgdown::build_site()
pkgdown::build_site_github_pages()
