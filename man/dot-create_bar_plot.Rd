% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gkwgof.R
\name{.create_bar_plot}
\alias{.create_bar_plot}
\title{Create Bar Plot for Model Comparison}
\usage{
.create_bar_plot(data, title, metrics, higher_better = FALSE)
}
\arguments{
\item{data}{Data frame containing the comparison data}

\item{title}{Plot title}

\item{metrics}{Vector of metric names to include in the plot}

\item{higher_better}{Logical; if TRUE, higher values indicate better fit}
}
\value{
A ggplot2 object
}
\description{
Create Bar Plot for Model Comparison
}
\keyword{internal}
