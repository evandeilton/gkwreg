% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gkwfit.R
\name{BIC.gkwfit}
\alias{BIC.gkwfit}
\title{Calculate Bayesian Information Criterion (BIC) for gkwfit Objects}
\usage{
\method{BIC}{gkwfit}(object, ...)
}
\arguments{
\item{object}{An object of class \code{"gkwfit"}, typically the result of a call
to \code{\link{gkwfit}}.}

\item{...}{Optionally, more fitted model objects of class \code{"gkwfit"}.}
}
\value{
\itemize{
\item If only one \code{object} is provided: A single numeric value, the calculated BIC.
\item If multiple objects are provided: A \code{data.frame} with rows corresponding
to the models and columns for the degrees of freedom (\code{df}) and the
calculated BIC value (named \code{BIC}). The data frame is sorted in
ascending order based on the BIC values. Row names are generated from the
deparsed calls or the names of the arguments passed to BIC.
}
}
\description{
Computes the Bayesian Information Criterion (BIC), sometimes called the
Schwarz criterion (SIC), for one or more fitted model objects of class \code{"gkwfit"}.
}
\details{
This function calculates the BIC based on the formula
\eqn{-2 \times \log Likelihood + \log(n) \times df}, where \eqn{n} is the number
of observations and \eqn{df} represents the number of estimated parameters in the
model (degrees of freedom).

It relies on the \code{\link{logLik.gkwfit}} method to extract the log-likelihood,
the degrees of freedom (\code{df}), and the number of observations (\code{nobs})
for each model. Ensure that \code{logLik.gkwfit} is defined and returns a valid
\code{"logLik"} object with appropriate attributes.

When comparing multiple models fitted to the \strong{same data}, the model with the
lower BIC value is generally preferred, as BIC tends to penalize model complexity
more heavily than AIC for larger sample sizes. The function returns a sorted
data frame to facilitate this comparison when multiple objects are provided. A
warning is issued if models were fitted to different numbers of observations.
}
\examples{
\donttest{

set.seed(2203)
y <- rkw(1000, alpha = 2.5, beta = 1.5)

# Fit different models to the same data
fit1_kw <- gkwfit(y, family = "kw", silent = TRUE)
fit2_bkw <- gkwfit(y, family = "bkw", silent = TRUE)
fit3_gkw <- gkwfit(y, family = "gkw", silent = TRUE)

# Calculate BIC for a single model
bic1 <- BIC(fit1_kw)
print(bic1)

# Compare BIC values for multiple models
bic_comparison <- c(BIC(fit1_kw), BIC(fit2_bkw), BIC(fit3_gkw))
print(bic_comparison)
}

}
\seealso{
\code{\link{gkwfit}}, \code{\link[stats]{BIC}}, \code{\link{logLik.gkwfit}}, \code{\link{AIC.gkwfit}}
}
\author{
Lopes, J. E. (with refinements)
}
\keyword{methods}
\keyword{models}
\keyword{stats}
