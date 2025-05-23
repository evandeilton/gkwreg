% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{llgkw}
\alias{llgkw}
\title{Negative Log-Likelihood for the Generalized Kumaraswamy Distribution}
\usage{
llgkw(par, data)
}
\arguments{
\item{par}{A numeric vector of length 5 containing the distribution parameters
in the order: \code{alpha} (\eqn{\alpha > 0}), \code{beta} (\eqn{\beta > 0}),
\code{gamma} (\eqn{\gamma > 0}), \code{delta} (\eqn{\delta \ge 0}),
\code{lambda} (\eqn{\lambda > 0}).}

\item{data}{A numeric vector of observations. All values must be strictly
between 0 and 1 (exclusive).}
}
\value{
Returns a single \code{double} value representing the negative
log-likelihood (\eqn{-\ell(\theta|\mathbf{x})}). Returns a large positive
value (e.g., \code{Inf}) if any parameter values in \code{par} are invalid
according to their constraints, or if any value in \code{data} is not in
the interval (0, 1).
}
\description{
Computes the negative log-likelihood function for the five-parameter
Generalized Kumaraswamy (GKw) distribution given a vector of observations.
This function is designed for use in optimization routines (e.g., maximum
likelihood estimation).
}
\details{
The probability density function (PDF) of the GKw distribution is given in
\code{\link{dgkw}}. The log-likelihood function \eqn{\ell(\theta)} for a sample
\eqn{\mathbf{x} = (x_1, \dots, x_n)} is:
\deqn{
\ell(\theta | \mathbf{x}) = n\ln(\lambda\alpha\beta) - n\ln B(\gamma,\delta+1) +
  \sum_{i=1}^{n} [(\alpha-1)\ln(x_i) + (\beta-1)\ln(v_i) + (\gamma\lambda-1)\ln(w_i) + \delta\ln(z_i)]
}
where \eqn{\theta = (\alpha, \beta, \gamma, \delta, \lambda)}, \eqn{B(a,b)}
is the Beta function (\code{\link[base]{beta}}), and:
\itemize{
\item \eqn{v_i = 1 - x_i^{\alpha}}
\item \eqn{w_i = 1 - v_i^{\beta} = 1 - (1-x_i^{\alpha})^{\beta}}
\item \eqn{z_i = 1 - w_i^{\lambda} = 1 - [1-(1-x_i^{\alpha})^{\beta}]^{\lambda}}
}
This function computes \eqn{-\ell(\theta|\mathbf{x})}.

Numerical stability is prioritized using:
\itemize{
\item \code{\link[base]{lbeta}} function for the log-Beta term.
\item Log-transformations of intermediate terms (\eqn{v_i, w_i, z_i}) and
use of \code{\link[base]{log1p}} where appropriate to handle values
close to 0 or 1 accurately.
\item Checks for invalid parameters and data.
}
}
\examples{
\donttest{
# Generate sample data from a known GKw distribution
set.seed(123)
true_par <- c(alpha = 2, beta = 3, gamma = 1.0, delta = 0.5, lambda = 0.5)
sample_data <- rgkw(100, alpha = true_par[1], beta = true_par[2],
                   gamma = true_par[3], delta = true_par[4], lambda = true_par[5])
hist(sample_data, breaks = 20, main = "Sample GKw Data")

# --- Maximum Likelihood Estimation using optim ---
# Initial parameter guess (can be crucial)
start_par <- c(1.5, 2.5, 1.2, 0.3, 0.6)

# Perform optimization (minimizing negative log-likelihood)
# Ensure data is passed correctly to llgkw
mle_result <- stats::optim(par = start_par,
                           fn = llgkw,
                           method = "BFGS", # Method supporting bounds might be safer
                           hessian = TRUE,
                           data = sample_data)

# Check convergence and results
if (mle_result$convergence == 0) {
  print("Optimization converged successfully.")
  mle_par <- mle_result$par
  print("Estimated parameters:")
  print(mle_par)
  print("True parameters:")
  print(true_par)

  # Standard errors from Hessian (optional)
  # fisher_info <- solve(mle_result$hessian) # Need positive definite Hessian
  # std_errors <- sqrt(diag(fisher_info))
  # print("Approximate Standard Errors:")
  # print(std_errors)

} else {
  warning("Optimization did not converge!")
  print(mle_result$message)
}

# --- Compare numerical and analytical derivatives (if available) ---
# Requires the 'numDeriv' package and analytical functions 'grgkw', 'hsgkw'
if (requireNamespace("numDeriv", quietly = TRUE) &&
    exists("grgkw") && exists("hsgkw") && mle_result$convergence == 0) {

  cat("\nComparing Derivatives at MLE estimates:\n")

  # Numerical derivatives
  num_grad <- numDeriv::grad(func = llgkw, x = mle_par, data = sample_data)
  num_hess <- numDeriv::hessian(func = llgkw, x = mle_par, data = sample_data)

  # Analytical derivatives (assuming they exist)
  # Note: grgkw/hsgkw might compute derivatives of log-likelihood,
  # while llgkw is negative log-likelihood. Adjust signs if needed.
  # Assuming grgkw/hsgkw compute derivatives of NEGATIVE log-likelihood here:
  ana_grad <- grgkw(par = mle_par, data = sample_data)
  ana_hess <- hsgkw(par = mle_par, data = sample_data)

  # Check differences (should be small if analytical functions are correct)
  cat("Difference between numerical and analytical gradient:\n")
  print(summary(abs(num_grad - ana_grad)))

  cat("Difference between numerical and analytical Hessian:\n")
  print(summary(abs(num_hess - ana_hess)))

} else {
   cat("\nSkipping derivative comparison.\n")
   cat("Requires 'numDeriv' package and functions 'grgkw', 'hsgkw'.\n")
}

}

}
\references{
Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized
distributions. \emph{Journal of Statistical Computation and Simulation}

Kumaraswamy, P. (1980). A generalized probability density function for
double-bounded random processes. \emph{Journal of Hydrology}, \emph{46}(1-2), 79-88.
}
\seealso{
\code{\link{dgkw}}, \code{\link{pgkw}}, \code{\link{qgkw}}, \code{\link{rgkw}},
\code{\link{grgkw}}, \code{\link{hsgkw}} (gradient and Hessian functions, if available),
\code{\link[stats]{optim}}, \code{\link[base]{lbeta}}, \code{\link[base]{log1p}}
}
\author{
Lopes, J. E.
}
\keyword{distribution}
\keyword{likelihood}
\keyword{optimize}
