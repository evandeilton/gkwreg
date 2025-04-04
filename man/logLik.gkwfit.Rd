% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gkwfit.R
\name{logLik.gkwfit}
\alias{logLik.gkwfit}
\title{Extract Log-Likelihood from a gkwfit Object}
\usage{
\method{logLik}{gkwfit}(object, ...)
}
\arguments{
\item{object}{An object of class \code{"gkwfit"}, typically the result of a call to \code{\link{gkwfit}}.}

\item{...}{Additional arguments (currently ignored).}
}
\value{
An object of class \code{"logLik"}. This is the numeric log-likelihood value
with the following attributes:
\item{df}{The number of estimated parameters in the model (integer).}
\item{nobs}{The number of observations used for fitting the model (integer).}
}
\description{
Extracts the maximized log-likelihood value from a model fitted by \code{\link{gkwfit}}.
It returns an object of class \code{"logLik"}, which includes attributes for the
degrees of freedom (\code{"df"}) and the number of observations (\code{"nobs"}) used in the fit.
}
\details{
This method provides compatibility with standard R functions that operate on
log-likelihood values, such as \code{\link[stats]{AIC}}, \code{\link[stats]{BIC}},
and likelihood ratio tests. It retrieves the log-likelihood stored during the
model fitting process (in \code{object$loglik}) and attaches the required
attributes (\code{object$df} for the number of estimated parameters and
\code{object$nobs} for the number of observations).
}
\examples{
\dontrun{
# Assume 'rkw' and 'gkwfit' functions from your package exist

set.seed(123)
# Generate sample data (use rkw if available, otherwise placeholder)
# kw_data_sample <- rkw(100, alpha = 2.5, beta = 1.5)
kw_data_sample <- runif(100)^(1 / 1.5) # Placeholder if rkw not available
kw_data_sample <- 1 - (1 - kw_data_sample)^(1 / 2.5) # Placeholder
kw_data_sample <- pmax(1e-6, pmin(1 - 1e-6, kw_data_sample)) # Ensure (0,1)

# Fit the model
fit_obj <- gkwfit(data = kw_data_sample, family = "kw", silent = TRUE)

# Check if fit converged and has necessary components
if (fit_obj$convergence == 0 && !is.null(fit_obj$loglik)) {
  # Extract the log-likelihood object
  ll <- logLik(fit_obj)

  # Print the value (uses print.logLik)
  print(ll)

  # Show attributes
  attributes(ll)

  # Use with standard functions (requires the corresponding AIC/BIC methods
  # for gkwfit OR relies on them calling this logLik method)
  AIC(fit_obj) # Should work if AIC.gkwfit calls logLik(fit_obj)
  BIC(fit_obj) # Should work if BIC.gkwfit calls logLik(fit_obj)

  # Extract components directly
  logLik_value <- as.numeric(ll)
  degrees_freedom <- attr(ll, "df")
  num_observations <- attr(ll, "nobs")
  cat("LogLik:", logLik_value, "\n")
  cat("df:", degrees_freedom, "\n")
  cat("nobs:", num_observations, "\n")
} else {
  print("Fit did not converge or is missing log-likelihood.")
}
}

}
\seealso{
\code{\link{gkwfit}}, \code{\link[stats]{AIC}}, \code{\link[stats]{BIC}}, \code{\link[stats]{logLik}}
}
\author{
Lopes, J. E.
}
\keyword{methods}
\keyword{models}
\keyword{utility}
