% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gkwreg.R
\name{.calculate_half_normal_data}
\alias{.calculate_half_normal_data}
\title{Calculate half-normal plot data with envelope}
\usage{
.calculate_half_normal_data(
  resid_vec,
  idx,
  nsim,
  level,
  param_mat,
  param_vectors,
  type,
  family
)
}
\arguments{
\item{resid_vec}{Vector of residuals}

\item{idx}{Vector of observation indices}

\item{nsim}{Number of simulations for envelope}

\item{level}{Confidence level for envelope}

\item{param_mat}{Matrix of calculated parameters}

\item{param_vectors}{List of parameter vectors}

\item{type}{Character string specifying residual type}

\item{family}{Character string specifying distribution family}
}
\value{
A data frame with half-normal plot data
}
\description{
Calculate half-normal plot data with envelope
}
\keyword{internal}
