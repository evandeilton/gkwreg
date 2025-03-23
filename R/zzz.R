utils::globalVariables(c(
  "packageVersion", "dev.interactive", "calculateMean", "index", "theoretical", "observed", "lower",
  "upper", "log", "log_scale", "calculateResponseResiduals", "::", ":::", ".gkwreg_env", "get_tmb_info",
  "x", "Theoretical", "Empirical", "value", "loglik"
))


#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
NULL


## usethis namespace: start
#' @importFrom Rcpp sourceCpp evalCpp
#' @import RcppArmadillo
#' @import graphics
## usethis namespace: end
NULL

#' @useDynLib gkwreg, .registration = TRUE
NULL
