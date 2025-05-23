% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gkwreg.R
\name{.prepare_tmb_params}
\alias{.prepare_tmb_params}
\title{Prepare TMB Parameters for GKw Regression}
\usage{
.prepare_tmb_params(model_data, family, param_names, fixed, param_positions)
}
\arguments{
\item{model_data}{List of model data.}

\item{family}{Family name.}

\item{param_names}{Names of parameters.}

\item{fixed}{List of fixed parameters.}

\item{param_positions}{Parameter position mapping for the family.}
}
\value{
A list with TMB parameters.
}
\description{
Prepare TMB Parameters for GKw Regression
}
\keyword{internal}
