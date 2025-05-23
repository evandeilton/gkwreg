% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{qbeta_}
\alias{qbeta_}
\title{Quantile Function of the Beta Distribution (gamma, delta+1 Parameterization)}
\usage{
qbeta_(p, gamma, delta, lower_tail = TRUE, log_p = FALSE)
}
\arguments{
\item{p}{Vector of probabilities (values between 0 and 1).}

\item{gamma}{First shape parameter (\code{shape1}), \eqn{\gamma > 0}. Can be a
scalar or a vector. Default: 1.0.}

\item{delta}{Second shape parameter is \code{delta + 1} (\code{shape2}), requires
\eqn{\delta \ge 0} so that \code{shape2 >= 1}. Can be a scalar or a vector.
Default: 0.0 (leading to \code{shape2 = 1}).}

\item{lower_tail}{Logical; if \code{TRUE} (default), probabilities are \eqn{p = P(X \le q)},
otherwise, probabilities are \eqn{p = P(X > q)}.}

\item{log_p}{Logical; if \code{TRUE}, probabilities \code{p} are given as
\eqn{\log(p)}. Default: \code{FALSE}.}
}
\value{
A vector of quantiles corresponding to the given probabilities \code{p}.
The length of the result is determined by the recycling rule applied to
the arguments (\code{p}, \code{gamma}, \code{delta}).
Returns:
\itemize{
\item \code{0} for \code{p = 0} (or \code{p = -Inf} if \code{log_p = TRUE},
when \code{lower_tail = TRUE}).
\item \code{1} for \code{p = 1} (or \code{p = 0} if \code{log_p = TRUE},
when \code{lower_tail = TRUE}).
\item \code{NaN} for \code{p < 0} or \code{p > 1} (or corresponding log scale).
\item \code{NaN} for invalid parameters (e.g., \code{gamma <= 0},
\code{delta < 0}).
}
Boundary return values are adjusted accordingly for \code{lower_tail = FALSE}.
}
\description{
Computes the quantile function (inverse CDF) for the standard Beta
distribution, using a parameterization common in generalized distribution
families. It finds the value \code{q} such that \eqn{P(X \le q) = p}. The
distribution is parameterized by \code{gamma} (\eqn{\gamma}) and \code{delta}
(\eqn{\delta}), corresponding to the standard Beta distribution with shape
parameters \code{shape1 = gamma} and \code{shape2 = delta + 1}.
}
\details{
This function computes the quantiles of a Beta distribution with parameters
\code{shape1 = gamma} and \code{shape2 = delta + 1}. It is equivalent to
calling \code{stats::qbeta(p, shape1 = gamma, shape2 = delta + 1,
lower.tail = lower_tail, log.p = log_p)}.

This distribution arises as a special case of the five-parameter
Generalized Kumaraswamy (GKw) distribution (\code{\link{qgkw}}) obtained
by setting \eqn{\alpha = 1}, \eqn{\beta = 1}, and \eqn{\lambda = 1}.
It is therefore also equivalent to the McDonald (Mc)/Beta Power distribution
(\code{\link{qmc}}) with \eqn{\lambda = 1}.

The function likely calls R's underlying \code{qbeta} function but ensures
consistent parameter recycling and handling within the C++ environment,
matching the style of other functions in the related families. Boundary
conditions (p=0, p=1) are handled explicitly.
}
\examples{
\donttest{
# Example values
p_vals <- c(0.1, 0.5, 0.9)
gamma_par <- 2.0 # Corresponds to shape1
delta_par <- 3.0 # Corresponds to shape2 - 1
shape1 <- gamma_par
shape2 <- delta_par + 1

# Calculate quantiles using qbeta_
quantiles <- qbeta_(p_vals, gamma_par, delta_par)
print(quantiles)

# Compare with stats::qbeta
quantiles_stats <- stats::qbeta(p_vals, shape1 = shape1, shape2 = shape2)
print(paste("Max difference vs stats::qbeta:", max(abs(quantiles - quantiles_stats))))

# Compare with qgkw setting alpha=1, beta=1, lambda=1
quantiles_gkw <- qgkw(p_vals, alpha = 1.0, beta = 1.0, gamma = gamma_par,
                      delta = delta_par, lambda = 1.0)
print(paste("Max difference vs qgkw:", max(abs(quantiles - quantiles_gkw))))

# Compare with qmc setting lambda=1
quantiles_mc <- qmc(p_vals, gamma = gamma_par, delta = delta_par, lambda = 1.0)
print(paste("Max difference vs qmc:", max(abs(quantiles - quantiles_mc))))

# Calculate quantiles for upper tail
quantiles_upper <- qbeta_(p_vals, gamma_par, delta_par, lower_tail = FALSE)
print(quantiles_upper)
print(stats::qbeta(p_vals, shape1, shape2, lower.tail = FALSE))

# Calculate quantiles from log probabilities
log_p_vals <- log(p_vals)
quantiles_logp <- qbeta_(log_p_vals, gamma_par, delta_par, log_p = TRUE)
print(quantiles_logp)
print(stats::qbeta(log_p_vals, shape1, shape2, log.p = TRUE))

# Verify inverse relationship with pbeta_
p_check <- 0.75
q_calc <- qbeta_(p_check, gamma_par, delta_par)
p_recalc <- pbeta_(q_calc, gamma_par, delta_par)
print(paste("Original p:", p_check, " Recalculated p:", p_recalc))
# abs(p_check - p_recalc) < 1e-9 # Should be TRUE

# Boundary conditions
print(qbeta_(c(0, 1), gamma_par, delta_par)) # Should be 0, 1
print(qbeta_(c(-Inf, 0), gamma_par, delta_par, log_p = TRUE)) # Should be 0, 1

}

}
\references{
Johnson, N. L., Kotz, S., & Balakrishnan, N. (1995). \emph{Continuous Univariate
Distributions, Volume 2} (2nd ed.). Wiley.

Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized
distributions. \emph{Journal of Statistical Computation and Simulation},
}
\seealso{
\code{\link[stats]{qbeta}} (standard R implementation),
\code{\link{qgkw}} (parent distribution quantile function),
\code{\link{qmc}} (McDonald/Beta Power quantile function),
\code{dbeta_}, \code{pbeta_}, \code{rbeta_} (other functions for this parameterization, if they exist).
}
\author{
Lopes, J. E.
}
\keyword{beta}
\keyword{distribution}
\keyword{quantile}
