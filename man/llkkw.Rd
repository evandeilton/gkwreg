% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{llkkw}
\alias{llkkw}
\title{Negative Log-Likelihood for the kkw Distribution}
\usage{
llkkw(par, data)
}
\arguments{
\item{par}{A numeric vector of length 4 containing the distribution parameters
in the order: \code{alpha} (\eqn{\alpha > 0}), \code{beta} (\eqn{\beta > 0}),
\code{delta} (\eqn{\delta \ge 0}), \code{lambda} (\eqn{\lambda > 0}).}

\item{data}{A numeric vector of observations. All values must be strictly
between 0 and 1 (exclusive).}
}
\value{
Returns a single \code{double} value representing the negative
log-likelihood (\eqn{-\ell(\theta|\mathbf{x})}). Returns \code{Inf}
if any parameter values in \code{par} are invalid according to their
constraints, or if any value in \code{data} is not in the interval (0, 1).
}
\description{
Computes the negative log-likelihood function for the Kumaraswamy-Kumaraswamy
(kkw) distribution with parameters \code{alpha} (\eqn{\alpha}), \code{beta}
(\eqn{\beta}), \code{delta} (\eqn{\delta}), and \code{lambda} (\eqn{\lambda}),
given a vector of observations. This distribution is a special case of the
Generalized Kumaraswamy (GKw) distribution where \eqn{\gamma = 1}.
}
\details{
The kkw distribution is the GKw distribution (\code{\link{dgkw}}) with \eqn{\gamma=1}.
Its probability density function (PDF) is:
\deqn{
f(x | \theta) = (\delta + 1) \lambda \alpha \beta x^{\alpha - 1} (1 - x^\alpha)^{\beta - 1} \bigl[1 - (1 - x^\alpha)^\beta\bigr]^{\lambda - 1} \bigl\{1 - \bigl[1 - (1 - x^\alpha)^\beta\bigr]^\lambda\bigr\}^{\delta}
}
for \eqn{0 < x < 1} and \eqn{\theta = (\alpha, \beta, \delta, \lambda)}.
The log-likelihood function \eqn{\ell(\theta | \mathbf{x})} for a sample
\eqn{\mathbf{x} = (x_1, \dots, x_n)} is \eqn{\sum_{i=1}^n \ln f(x_i | \theta)}:
\deqn{
\ell(\theta | \mathbf{x}) = n[\ln(\delta+1) + \ln(\lambda) + \ln(\alpha) + \ln(\beta)]
+ \sum_{i=1}^{n} [(\alpha-1)\ln(x_i) + (\beta-1)\ln(v_i) + (\lambda-1)\ln(w_i) + \delta\ln(z_i)]
}
where:
\itemize{
\item \eqn{v_i = 1 - x_i^{\alpha}}
\item \eqn{w_i = 1 - v_i^{\beta} = 1 - (1-x_i^{\alpha})^{\beta}}
\item \eqn{z_i = 1 - w_i^{\lambda} = 1 - [1-(1-x_i^{\alpha})^{\beta}]^{\lambda}}
}
This function computes and returns the \emph{negative} log-likelihood, \eqn{-\ell(\theta|\mathbf{x})},
suitable for minimization using optimization routines like \code{\link[stats]{optim}}.
Numerical stability is maintained similarly to \code{\link{llgkw}}.
}
\examples{
\donttest{
# Assuming existence of rkkw, grkkw, hskkw functions for kkw distribution

# Generate sample data from a known kkw distribution
set.seed(123)
true_par_kkw <- c(alpha = 2, beta = 3, delta = 1.5, lambda = 0.5)
# Use rkkw if it exists, otherwise use rgkw with gamma=1
if (exists("rkkw")) {
  sample_data_kkw <- rkkw(100, alpha = true_par_kkw[1], beta = true_par_kkw[2],
                         delta = true_par_kkw[3], lambda = true_par_kkw[4])
} else {
  sample_data_kkw <- rgkw(100, alpha = true_par_kkw[1], beta = true_par_kkw[2],
                         gamma = 1, delta = true_par_kkw[3], lambda = true_par_kkw[4])
}
hist(sample_data_kkw, breaks = 20, main = "kkw(2, 3, 1.5, 0.5) Sample")

# --- Maximum Likelihood Estimation using optim ---
# Initial parameter guess
start_par_kkw <- c(1.5, 2.5, 1.0, 0.6)

# Perform optimization (minimizing negative log-likelihood)
mle_result_kkw <- stats::optim(par = start_par_kkw,
                               fn = llkkw, # Use the kkw neg-log-likelihood
                               method = "BFGS",
                               hessian = TRUE,
                               data = sample_data_kkw)

# Check convergence and results
if (mle_result_kkw$convergence == 0) {
  print("Optimization converged successfully.")
  mle_par_kkw <- mle_result_kkw$par
  print("Estimated kkw parameters:")
  print(mle_par_kkw)
  print("True kkw parameters:")
  print(true_par_kkw)
} else {
  warning("Optimization did not converge!")
  print(mle_result_kkw$message)
}

# --- Compare numerical and analytical derivatives (if available) ---
# Requires 'numDeriv' package and analytical functions 'grkkw', 'hskkw'
if (mle_result_kkw$convergence == 0 &&
    requireNamespace("numDeriv", quietly = TRUE) &&
    exists("grkkw") && exists("hskkw")) {

  cat("\nComparing Derivatives at kkw MLE estimates:\n")

  # Numerical derivatives of llkkw
  num_grad_kkw <- numDeriv::grad(func = llkkw, x = mle_par_kkw, data = sample_data_kkw)
  num_hess_kkw <- numDeriv::hessian(func = llkkw, x = mle_par_kkw, data = sample_data_kkw)

  # Analytical derivatives (assuming they return derivatives of negative LL)
  ana_grad_kkw <- grkkw(par = mle_par_kkw, data = sample_data_kkw)
  ana_hess_kkw <- hskkw(par = mle_par_kkw, data = sample_data_kkw)

  # Check differences
  cat("Max absolute difference between gradients:\n")
  print(max(abs(num_grad_kkw - ana_grad_kkw)))
  cat("Max absolute difference between Hessians:\n")
  print(max(abs(num_hess_kkw - ana_hess_kkw)))

} else {
   cat("\nSkipping derivative comparison for kkw.\n")
   cat("Requires convergence, 'numDeriv' package and functions 'grkkw', 'hskkw'.\n")
}

}

}
\references{
Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized
distributions. \emph{Journal of Statistical Computation and Simulation},

Kumaraswamy, P. (1980). A generalized probability density function for
double-bounded random processes. \emph{Journal of Hydrology}, \emph{46}(1-2), 79-88.
}
\seealso{
\code{\link{llgkw}} (parent distribution negative log-likelihood),
\code{\link{dkkw}}, \code{\link{pkkw}}, \code{\link{qkkw}}, \code{\link{rkkw}},
\code{\link{grkkw}} (gradient, if available),
\code{\link{hskkw}} (Hessian, if available),
\code{\link[stats]{optim}}
}
\author{
Lopes, J. E.
}
\keyword{distribution}
\keyword{likelihood}
\keyword{optimize}
