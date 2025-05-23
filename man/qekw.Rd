% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{qekw}
\alias{qekw}
\title{Quantile Function of the Exponentiated Kumaraswamy (EKw) Distribution}
\usage{
qekw(p, alpha, beta, lambda, lower_tail = TRUE, log_p = FALSE)
}
\arguments{
\item{p}{Vector of probabilities (values between 0 and 1).}

\item{alpha}{Shape parameter \code{alpha} > 0. Can be a scalar or a vector.
Default: 1.0.}

\item{beta}{Shape parameter \code{beta} > 0. Can be a scalar or a vector.
Default: 1.0.}

\item{lambda}{Shape parameter \code{lambda} > 0 (exponent parameter).
Can be a scalar or a vector. Default: 1.0.}

\item{lower_tail}{Logical; if \code{TRUE} (default), probabilities are \eqn{p = P(X \le q)},
otherwise, probabilities are \eqn{p = P(X > q)}.}

\item{log_p}{Logical; if \code{TRUE}, probabilities \code{p} are given as
\eqn{\log(p)}. Default: \code{FALSE}.}
}
\value{
A vector of quantiles corresponding to the given probabilities \code{p}.
The length of the result is determined by the recycling rule applied to
the arguments (\code{p}, \code{alpha}, \code{beta}, \code{lambda}).
Returns:
\itemize{
\item \code{0} for \code{p = 0} (or \code{p = -Inf} if \code{log_p = TRUE},
when \code{lower_tail = TRUE}).
\item \code{1} for \code{p = 1} (or \code{p = 0} if \code{log_p = TRUE},
when \code{lower_tail = TRUE}).
\item \code{NaN} for \code{p < 0} or \code{p > 1} (or corresponding log scale).
\item \code{NaN} for invalid parameters (e.g., \code{alpha <= 0},
\code{beta <= 0}, \code{lambda <= 0}).
}
Boundary return values are adjusted accordingly for \code{lower_tail = FALSE}.
}
\description{
Computes the quantile function (inverse CDF) for the Exponentiated
Kumaraswamy (EKw) distribution with parameters \code{alpha} (\eqn{\alpha}),
\code{beta} (\eqn{\beta}), and \code{lambda} (\eqn{\lambda}).
It finds the value \code{q} such that \eqn{P(X \le q) = p}. This distribution
is a special case of the Generalized Kumaraswamy (GKw) distribution where
\eqn{\gamma = 1} and \eqn{\delta = 0}.
}
\details{
The quantile function \eqn{Q(p)} is the inverse of the CDF \eqn{F(q)}. The CDF
for the EKw (\eqn{\gamma=1, \delta=0}) distribution is \eqn{F(q) = [1 - (1 - q^\alpha)^\beta ]^\lambda}
(see \code{\link{pekw}}). Inverting this equation for \eqn{q} yields the
quantile function:
\deqn{
Q(p) = \left\{ 1 - \left[ 1 - p^{1/\lambda} \right]^{1/\beta} \right\}^{1/\alpha}
}
The function uses this closed-form expression and correctly handles the
\code{lower_tail} and \code{log_p} arguments by transforming \code{p}
appropriately before applying the formula. This is equivalent to the general
GKw quantile function (\code{\link{qgkw}}) evaluated with \eqn{\gamma=1, \delta=0}.
}
\examples{
\donttest{
# Example values
p_vals <- c(0.1, 0.5, 0.9)
alpha_par <- 2.0
beta_par <- 3.0
lambda_par <- 1.5

# Calculate quantiles
quantiles <- qekw(p_vals, alpha_par, beta_par, lambda_par)
print(quantiles)

# Calculate quantiles for upper tail probabilities P(X > q) = p
quantiles_upper <- qekw(p_vals, alpha_par, beta_par, lambda_par,
                        lower_tail = FALSE)
print(quantiles_upper)
# Check: qekw(p, ..., lt=F) == qekw(1-p, ..., lt=T)
print(qekw(1 - p_vals, alpha_par, beta_par, lambda_par))

# Calculate quantiles from log probabilities
log_p_vals <- log(p_vals)
quantiles_logp <- qekw(log_p_vals, alpha_par, beta_par, lambda_par,
                       log_p = TRUE)
print(quantiles_logp)
# Check: should match original quantiles
print(quantiles)

# Compare with qgkw setting gamma = 1, delta = 0
quantiles_gkw <- qgkw(p_vals, alpha = alpha_par, beta = beta_par,
                     gamma = 1.0, delta = 0.0, lambda = lambda_par)
print(paste("Max difference:", max(abs(quantiles - quantiles_gkw)))) # Should be near zero

# Verify inverse relationship with pekw
p_check <- 0.75
q_calc <- qekw(p_check, alpha_par, beta_par, lambda_par)
p_recalc <- pekw(q_calc, alpha_par, beta_par, lambda_par)
print(paste("Original p:", p_check, " Recalculated p:", p_recalc))
# abs(p_check - p_recalc) < 1e-9 # Should be TRUE

# Boundary conditions
print(qekw(c(0, 1), alpha_par, beta_par, lambda_par)) # Should be 0, 1
print(qekw(c(-Inf, 0), alpha_par, beta_par, lambda_par, log_p = TRUE)) # Should be 0, 1
}

}
\references{
Nadarajah, S., Cordeiro, G. M., & Ortega, E. M. (2012). The exponentiated
Kumaraswamy distribution. \emph{Journal of the Franklin Institute}, \emph{349}(3),

Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized
distributions. \emph{Journal of Statistical Computation and Simulation},

Kumaraswamy, P. (1980). A generalized probability density function for
double-bounded random processes. \emph{Journal of Hydrology}, \emph{46}(1-2), 79-88.
}
\seealso{
\code{\link{qgkw}} (parent distribution quantile function),
\code{\link{dekw}}, \code{\link{pekw}}, \code{\link{rekw}} (other EKw functions),
\code{\link[stats]{qunif}}
}
\author{
Lopes, J. E.
}
\keyword{distribution}
\keyword{quantile}
