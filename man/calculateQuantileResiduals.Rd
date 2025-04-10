% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{calculateQuantileResiduals}
\alias{calculateQuantileResiduals}
\title{Calculate Quantile Residuals}
\usage{
calculateQuantileResiduals(y, params, family = "gkw")
}
\arguments{
\item{y}{NumericVector of observations.}

\item{params}{NumericMatrix with parameters (columns: alpha, beta, gamma, delta, lambda).}

\item{family}{String specifying the distribution family (default: "gkw").}
}
\value{
NumericVector of quantile residuals.
}
\description{
Computes quantile residuals by transforming the cumulative distribution function (CDF) values to the standard normal quantiles.
}
