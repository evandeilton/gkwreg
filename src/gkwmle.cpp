// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <limits>

using namespace Rcpp;
using namespace arma;

/*
* ===========================================================================
* NUMERIC STABILITY AUXILIARY FUNCTIONS
* ===========================================================================
* These functions ensure accurate numerical calculations even in extreme
* situations, near distribution boundaries, or with very small/large values.
*/

// Constants for numeric stability and precision
static const double EPSILON      = std::numeric_limits<double>::epsilon();
static const double DBL_MIN_SAFE = std::numeric_limits<double>::min() * 10.0;
static const double LOG_DBL_MIN  = std::log(DBL_MIN_SAFE);
static const double LOG_DBL_MAX  = std::log(std::numeric_limits<double>::max() / 10.0);
static const double LN2          = 0.6931471805599453;
static const double SQRT_EPSILON = 1.4901161193847656e-08; // sqrt(EPSILON)

/**
* log1mexp(u) calculates log(1 - exp(u)) with enhanced numerical stability
*
* This function is crucial for accurate calculations when u is negative and
* close to zero, where direct computation would suffer catastrophic cancellation.
* Uses different approximation methods depending on the range of u.
*
* @param u A negative value (log(x) where x < 1)
* @return log(1 - exp(u)), or NaN if u > 0
*/
inline double log1mexp(double u) {
// Input validation - u must be non-positive
if (u > 0.0) {
  return R_NaN;  // log(1 - exp(positive)) would yield log of negative number
}

// For u in (-ln(2), 0], use log(-expm1(u)) for better accuracy
if (u > -LN2) {
  return std::log(-std::expm1(u));
}

// For u <= -ln(2), use log1p(-exp(u)) for better accuracy
return std::log1p(-std::exp(u));
}

/**
* log1pexp(x) calculates log(1 + exp(x)) with protection against overflow
*
* This function handles various regimes of x with appropriate approximations
* to maintain numerical stability.
*
* @param x Input value
* @return log(1 + exp(x)) calculated with numerical stability
*/
inline double log1pexp(double x) {
if (x > 700.0)    return x;                      // For very large x, log(1+exp(x)) ≈ x
if (x > 18.0)     return x + std::log1p(std::exp(-x)); // For large x
if (x > -37.0)    return std::log1p(std::exp(x));      // For moderate x
return std::exp(x);                              // For very negative x, where exp(x) ≈ 0
}

/**
* safe_log(x) computes log(x) with protection against invalid inputs
*
* @param x Input value
* @return log(x) or appropriate limiting value for x <= 0 or very small x
*/
inline double safe_log(double x) {
if (x <= 0.0)        return R_NegInf;   // Log of non-positive number
if (x < DBL_MIN_SAFE) return LOG_DBL_MIN; // Prevent underflow
return std::log(x);
}

/**
* safe_exp(x) computes exp(x) with protection against overflow/underflow
*
* @param x Input value
* @return exp(x) or appropriate limiting value for extreme x
*/
inline double safe_exp(double x) {
if (x > LOG_DBL_MAX) return R_PosInf;  // Prevent overflow
if (x < LOG_DBL_MIN) return 0.0;       // Prevent underflow
return std::exp(x);
}

/**
* safe_pow(x, y) computes x^y with robust error handling
*
* Handles special cases like negative base with non-integer exponent,
* and uses logarithmic transformation for numerical stability with positive base.
*
* @param x Base value
* @param y Exponent value
* @return x^y calculated with numerical stability
*/
inline double safe_pow(double x, double y) {
// Handle special cases
if (x == 0.0) {
  if (y > 0.0)  return 0.0;
  if (y == 0.0) return 1.0;   // 0^0 convention
  return R_PosInf;            // 0^negative is undefined/infinity
}
if (x == 1.0 || y == 0.0) return 1.0;
if (y == 1.0) return x;

// Check for negative base with non-integer exponent (undefined in real domain)
if (x < 0.0 && std::abs(y - std::round(y)) > SQRT_EPSILON) {
  return R_NaN;
}

// For positive base, compute via logarithm for better numerical stability
if (x > 0.0) {
  double lx = std::log(x);
  double log_result = y * lx;
  return safe_exp(log_result);
}

// For negative base with integer exponent, use standard pow
return std::pow(x, y);
}

/**
* Vector version of log1mexp for element-wise operations on arma::vec
*
* @param u Vector of input values
* @return Vector of log(1 - exp(u)) values
*/
inline arma::vec vec_log1mexp(const arma::vec& u) {
arma::vec result(u.n_elem);
for (size_t i = 0; i < u.n_elem; ++i) {
  result[i] = log1mexp(u[i]);
}
return result;
}

/**
* Vector version of log1pexp for element-wise operations on arma::vec
*
* @param x Vector of input values
* @return Vector of log(1 + exp(x)) values
*/
inline arma::vec vec_log1pexp(const arma::vec& x) {
arma::vec result(x.n_elem);
for (size_t i = 0; i < x.n_elem; ++i) {
  result[i] = log1pexp(x[i]);
}
return result;
}

/**
* Vector version of safe_log for element-wise operations on arma::vec
*
* @param x Vector of input values
* @return Vector of safe_log(x) values
*/
inline arma::vec vec_safe_log(const arma::vec& x) {
arma::vec result(x.n_elem);
for (size_t i = 0; i < x.n_elem; ++i) {
  result[i] = safe_log(x[i]);
}
return result;
}

/**
* Vector version of safe_exp for element-wise operations on arma::vec
*
* @param x Vector of input values
* @return Vector of safe_exp(x) values
*/
inline arma::vec vec_safe_exp(const arma::vec& x) {
arma::vec result(x.n_elem);
for (size_t i = 0; i < x.n_elem; ++i) {
  result[i] = safe_exp(x[i]);
}
return result;
}

/**
* Vector version of safe_pow for element-wise operations
*
* @param x Vector of base values
* @param y Single exponent value
* @return Vector of x[i]^y values
*/
inline arma::vec vec_safe_pow(const arma::vec& x, double y) {
arma::vec result(x.n_elem);
for (size_t i = 0; i < x.n_elem; ++i) {
  result[i] = safe_pow(x[i], y);
}
return result;
}

/**
* Vector version of safe_pow with vector exponents
*
* @param x Vector of base values
* @param y Vector of exponent values (must match size of x)
* @return Vector of x[i]^y[i] values
*/
inline arma::vec vec_safe_pow(const arma::vec& x, const arma::vec& y) {
if (x.n_elem != y.n_elem) {
  Rcpp::stop("Vectors must have same length in vec_safe_pow");
}

arma::vec result(x.n_elem);
for (size_t i = 0; i < x.n_elem; ++i) {
  result[i] = safe_pow(x[i], y[i]);
}
return result;
}

/**
* Checks if GKw parameters are in the valid domain
*
* Verifies that all parameters satisfy the constraints:
* alpha > 0, beta > 0, gamma > 0, delta >= 0, lambda > 0
*
* With strict=true, also enforces reasonable bounds to avoid numerical issues.
*
* @param alpha Shape parameter
* @param beta Shape parameter
* @param gamma Shape parameter
* @param delta Shape parameter
* @param lambda Shape parameter
* @param strict Whether to enforce additional bounds for numerical stability
* @return true if parameters are valid, false otherwise
*/
inline bool check_pars(double alpha,
                     double beta,
                     double gamma,
                     double delta,
                     double lambda,
                     bool strict = false) {
// Basic parameter constraints
if (alpha <= 0.0 || beta <= 0.0 || gamma <= 0.0 || delta < 0.0 || lambda <= 0.0) {
  return false;
}

// Optional stricter constraints to avoid numerical issues
if (strict) {
  const double MIN_PARAM = 1e-5;
  const double MAX_PARAM = 1e5;

  if (alpha < MIN_PARAM || beta < MIN_PARAM || gamma < MIN_PARAM || lambda < MIN_PARAM) {
    return false;
  }
  if (alpha > MAX_PARAM || beta > MAX_PARAM || gamma > MAX_PARAM ||
      delta > MAX_PARAM || lambda > MAX_PARAM) {
    return false;
  }
}
return true;
}

/**
* Vector version of parameter checker for GKw distribution
*
* Checks all combinations of parameter values for validity.
*
* @param alpha Vector of alpha values
* @param beta Vector of beta values
* @param gamma Vector of gamma values
* @param delta Vector of delta values
* @param lambda Vector of lambda values
* @param strict Whether to enforce additional bounds for numerical stability
* @return arma::uvec of boolean values indicating parameter validity
*/
inline arma::uvec check_pars_vec(const arma::vec& alpha,
                               const arma::vec& beta,
                               const arma::vec& gamma,
                               const arma::vec& delta,
                               const arma::vec& lambda,
                               bool strict = false) {
// Find maximum length for broadcasting
size_t n = std::max({alpha.n_elem, beta.n_elem, gamma.n_elem,
                    delta.n_elem, lambda.n_elem});

arma::uvec valid(n, arma::fill::ones);

for (size_t i = 0; i < n; ++i) {
  // Get parameter values with proper cycling/broadcasting
  double a = alpha[i % alpha.n_elem];
  double b = beta[i % beta.n_elem];
  double g = gamma[i % gamma.n_elem];
  double d = delta[i % delta.n_elem];
  double l = lambda[i % lambda.n_elem];

  valid[i] = check_pars(a, b, g, d, l, strict);
}

return valid;
}

/*
* ===========================================================================
* PRIMARY FUNCTIONS FOR GENERALIZED KUMARASWAMY DISTRIBUTION
* ===========================================================================
*/

//' @title Density Function for Generalized Kumaraswamy Distribution
//'
//' @description
//' Calculates the probability density function (PDF) of the Generalized Kumaraswamy
//' distribution GKw(α, β, γ, δ, λ) for values in the open interval (0,1). Returns
//' either f(x; θ) or log(f(x; θ)) depending on the log_prob parameter.
//'
//' @param x Vector of quantiles where the density will be evaluated.
//' @param alpha Shape parameter α > 0 (scalar or vector). Controls the left tail behavior. Default: 1.0.
//' @param beta Shape parameter β > 0 (scalar or vector). Controls the right tail behavior. Default: 1.0.
//' @param gamma Shape parameter γ > 0 (scalar or vector). Affects the central shape. Default: 1.0.
//' @param delta Shape parameter δ ≥ 0 (scalar or vector). Introduces additional flexibility. Default: 0.0.
//' @param lambda Shape parameter λ > 0 (scalar or vector). Controls overall dispersion. Default: 1.0.
//' @param log_prob Logical; if TRUE, probabilities are returned as log(f(x)). Default: FALSE.
//'
//' @return Vector of density values corresponding to each element in x. If x contains
//'         values outside (0,1), those positions will return 0 (or -Inf if log_prob=TRUE).
//'
//' @details
//' The probability density function of the GKw distribution as derived by Carrasco et al. (2010) is:
//'
//' \deqn{
//' f(x; \alpha, \beta, \gamma, \delta, \lambda) =
//'   \frac{\lambda \alpha \beta x^{\alpha-1}(1-x^{\alpha})^{\beta-1}}
//'        {B(\gamma, \delta+1)}
//'   [1-(1-x^{\alpha})^{\beta}]^{\gamma\lambda-1}
//'   [1-[1-(1-x^{\alpha})^{\beta}]^{\lambda}]^{\delta}
//' }
//'
//' for x in (0,1), where B(γ, δ+1) is the beta function.
//'
//' This distribution arises from the construction:
//'
//' \deqn{f(x; \theta) = g_2(G_1(x; \alpha, \beta); \gamma, \delta, \lambda) \cdot g_1(x; \alpha, \beta)}
//'
//' where G_1 is the CDF of the Kumaraswamy distribution, g_1 is its PDF, and g_2 is the
//' generalized beta density of first kind.
//'
//' The function handles several edge cases:
//' \itemize{
//'   \item Returns 0 (or -Inf if log_prob=TRUE) for x outside (0,1)
//'   \item Implements numerical stability precautions for x very close to 0 or 1
//'   \item Returns 0 (or -Inf) for parameter combinations that would cause computational issues
//' }
//'
//' The GKw distribution includes several special cases:
//' \itemize{
//'   \item Kumaraswamy (Kw): γ=1, δ=0, λ=1
//'   \item McDonald distribution: α=1, β=1
//'   \item Beta distribution: α=1, β=1, λ=1
//'   \item Beta-Kumaraswamy (BKw): λ=1
//'   \item Kumaraswamy-Kumaraswamy: γ=1
//'   \item Exponentiated Kumaraswamy (EKw): γ=1, δ=0
//'   \item Beta Power (BP): α=1, β=1
//' }
//'
//' @references
//' Kumaraswamy, P. (1980). A generalized probability density function for double-bounded random processes.
//' Journal of Hydrology, 46(1-2), 79-88.
//'
//' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized distributions.
//' Journal of Statistical Computation and Simulation, 81(7), 883-898.
//'
//' @examples
//' \dontrun{
//' # Simple density evaluation at a point
//' dgkw(0.5, 2, 3, 1, 0, 1)
//'
//' # Plot the PDF for various parameter sets
//' x <- seq(0.01, 0.99, by = 0.01)
//'
//' # Standard Kumaraswamy (γ=1, δ=0, λ=1)
//' plot(x, dgkw(x, 2, 3, 1, 0, 1), type = "l",
//'      main = "GKw Densities", ylab = "f(x)", col = "blue")
//'
//' # Beta equivalent (α=1, β=1, λ=1)
//' lines(x, dgkw(x, 1, 1, 2, 3, 1), col = "red")
//'
//' # Exponentiated Kumaraswamy (γ=1, δ=0)
//' lines(x, dgkw(x, 2, 3, 1, 0, 2), col = "green")
//'
//' # Vectorized parameter example
//' alphas <- c(0.5, 1.5, 3.0)
//' # Returns 3 density values for the same x
//' dgkw(0.5, alphas, 2, 1, 0, 1)
//' }
//'
//' @export
// [[Rcpp::export]]
arma::vec dgkw(
   const arma::vec& x,
   const Rcpp::NumericVector& alpha = Rcpp::NumericVector::create(1.0),
   const Rcpp::NumericVector& beta = Rcpp::NumericVector::create(1.0),
   const Rcpp::NumericVector& gamma = Rcpp::NumericVector::create(1.0),
   const Rcpp::NumericVector& delta = Rcpp::NumericVector::create(0.0),
   const Rcpp::NumericVector& lambda = Rcpp::NumericVector::create(1.0),
   bool log_prob = false
) {
 // Convert NumericVector to arma::vec
 arma::vec alpha_vec(alpha.begin(), alpha.size());
 arma::vec beta_vec(beta.begin(), beta.size());
 arma::vec gamma_vec(gamma.begin(), gamma.size());
 arma::vec delta_vec(delta.begin(), delta.size());
 arma::vec lambda_vec(lambda.begin(), lambda.size());

 // Find the maximum length for broadcasting
 size_t n = std::max({x.n_elem, alpha_vec.n_elem, beta_vec.n_elem,
                     gamma_vec.n_elem, delta_vec.n_elem, lambda_vec.n_elem});

 // Initialize result vector
 arma::vec result(n);
 if (log_prob) {
   result.fill(R_NegInf);
 } else {
   result.fill(0.0);
 }

 // Process each element
 for (size_t i = 0; i < n; ++i) {
   // Get parameter values with broadcasting/recycling
   double a = alpha_vec[i % alpha_vec.n_elem];
   double b = beta_vec[i % beta_vec.n_elem];
   double g = gamma_vec[i % gamma_vec.n_elem];
   double d = delta_vec[i % delta_vec.n_elem];
   double l = lambda_vec[i % lambda_vec.n_elem];
   double xi = x[i % x.n_elem];

   // Validate parameters
   if (!check_pars(a, b, g, d, l)) {
     Rcpp::warning("dgkw: invalid parameters at index %d (alpha,beta,gamma>0, delta>=0, lambda>0)", i+1);
     continue;
   }

   // Check if x is within (0,1)
   if (xi <= 0.0 || xi >= 1.0 || !R_finite(xi)) {
     continue;
   }

   // Numerical stability: avoid calculations very close to 0 or 1
   double x_near_zero = std::pow(SQRT_EPSILON, 1.0/a);
   double x_near_one = 1.0 - std::pow(SQRT_EPSILON, 1.0/a);

   if (xi < x_near_zero || xi > x_near_one) {
     continue;
   }

   // Precalculate common terms used in the PDF formula
   double log_beta_val = R::lgammafn(g) + R::lgammafn(d + 1.0) - R::lgammafn(g + d + 1.0);
   double log_const = std::log(l) + std::log(a) + std::log(b) - log_beta_val;
   double gamma_lambda = g * l;

   // Calculate x^α
   double log_xi = std::log(xi);
   double log_x_alpha = a * log_xi;
   double x_alpha = std::exp(log_x_alpha);

   // Check if x^α < 1 for numerical stability
   if (x_alpha >= 1.0 - SQRT_EPSILON) {
     continue;
   }

   // Calculate (1 - x^α)
   double log_one_minus_x_alpha = log1mexp(log_x_alpha);
   if (!R_finite(log_one_minus_x_alpha)) {
     continue;
   }

   // Calculate (1 - x^α)^β
   double log_one_minus_x_alpha_beta = b * log_one_minus_x_alpha;

   // Calculate 1 - (1 - x^α)^β
   double log_term1 = log1mexp(log_one_minus_x_alpha_beta);
   if (!R_finite(log_term1)) {
     continue;
   }

   // Calculate [1-(1-x^α)^β]^λ
   double log_term1_lambda = l * log_term1;

   // Calculate 1 - [1-(1-x^α)^β]^λ
   double log_term2 = log1mexp(log_term1_lambda);
   if (!R_finite(log_term2)) {
     continue;
   }

   // Assemble the full log-density expression
   double logdens = log_const +
     (a - 1.0) * log_xi +
     (b - 1.0) * log_one_minus_x_alpha +
     (gamma_lambda - 1.0) * log_term1 +
     d * log_term2;

   // Check for invalid result
   if (!R_finite(logdens)) {
     continue;
   }

   // Return log-density or density as requested
   result(i) = log_prob ? logdens : safe_exp(logdens);
 }

 return result;
}


//' @title Cumulative Distribution Function for Generalized Kumaraswamy Distribution
//'
//' @description
//' Calculates the cumulative distribution function (CDF) of the Generalized Kumaraswamy
//' distribution GKw(α, β, γ, δ, λ) for values in the interval (0,1). Returns F(q) = P(X ≤ q)
//' by default, with options for upper tail probability P(X > q) and log transformations.
//'
//' @param q Vector of quantiles where the CDF will be evaluated.
//' @param alpha Shape parameter α > 0 (scalar or vector). Controls the left tail behavior. Default: 1.0.
//' @param beta Shape parameter β > 0 (scalar or vector). Controls the right tail behavior. Default: 1.0.
//' @param gamma Shape parameter γ > 0 (scalar or vector). Affects the central shape. Default: 1.0.
//' @param delta Shape parameter δ ≥ 0 (scalar or vector). Introduces additional flexibility. Default: 0.0.
//' @param lambda Shape parameter λ > 0 (scalar or vector). Controls overall dispersion. Default: 1.0.
//' @param lower_tail Logical; if TRUE (default), probabilities are P(X ≤ q), otherwise P(X > q).
//' @param log_p Logical; if TRUE, probabilities are returned as log(p). Default: FALSE.
//'
//' @return Vector of probabilities corresponding to each element in q.
//'
//' @details
//' The cumulative distribution function for the GKw distribution as derived in
//' Carrasco et al. (2010) is:
//'
//' \deqn{
//' F(q; \alpha, \beta, \gamma, \delta, \lambda) =
//'   I_{[1-(1-q^{\alpha})^{\beta}]^{\lambda}}(\gamma, \delta+1)
//' }
//'
//' where I_x(a,b) is the regularized incomplete beta function defined as:
//'
//' \deqn{
//' I_x(a, b) = \frac{B_x(a, b)}{B(a, b)}
//' }
//'
//' With B_x(a, b) being the incomplete beta function:
//'
//' \deqn{
//' B_x(a, b) = \int_0^x t^{a-1}(1-t)^{b-1} dt
//' }
//'
//' The GKw distribution can be derived from the following construction:
//' Let \eqn{G_1(x; \alpha, \beta)} be the two-parameter Kumaraswamy CDF and \eqn{g_2(x; \gamma, \delta, \lambda)}
//' be the generalized beta density of first kind. Then:
//'
//' \deqn{
//' F(x; \theta) = \int_0^{G_1(x; \alpha, \beta)} g_2(t; \gamma, \delta, \lambda) dt
//' }
//'
//' The function implements specialized numerical techniques for stability from the paper,
//' carefully handling boundary cases and extreme parameter values to ensure proper behavior
//' throughout the parameter space.
//'
//' @references
//' Kumaraswamy, P. (1980). A generalized probability density function for double-bounded random processes.
//' Journal of Hydrology, 46(1-2), 79-88.
//'
//' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized distributions.
//' Journal of Statistical Computation and Simulation, 81(7), 883-898.
//'
//' @examples
//' \dontrun{
//' # Simple CDF evaluation
//' pgkw(0.5, 2, 3, 1, 0, 1)
//'
//' # Use of vectorized parameters
//' alphas <- c(0.5, 1.0, 2.0)
//' betas <- c(1.0, 2.0, 3.0)
//' pgkw(0.5, alphas, betas)
//'
//' # Comparing special cases
//' x <- seq(0.01, 0.99, by = 0.01)
//' # Standard Kumaraswamy
//' pkw <- pgkw(x, 2, 3, 1, 0, 1)
//' # Beta distribution
//' pbeta_equiv <- pgkw(x, 1, 1, 2, 3, 1)
//' }
//'
//' @export
// [[Rcpp::export]]
arma::vec pgkw(
   const arma::vec& q,
   const Rcpp::NumericVector& alpha = Rcpp::NumericVector::create(1.0),
   const Rcpp::NumericVector& beta = Rcpp::NumericVector::create(1.0),
   const Rcpp::NumericVector& gamma = Rcpp::NumericVector::create(1.0),
   const Rcpp::NumericVector& delta = Rcpp::NumericVector::create(0.0),
   const Rcpp::NumericVector& lambda = Rcpp::NumericVector::create(1.0),
   bool lower_tail = true,
   bool log_p = false
) {
 // Convert NumericVector to arma::vec
 arma::vec alpha_vec(alpha.begin(), alpha.size());
 arma::vec beta_vec(beta.begin(), beta.size());
 arma::vec gamma_vec(gamma.begin(), gamma.size());
 arma::vec delta_vec(delta.begin(), delta.size());
 arma::vec lambda_vec(lambda.begin(), lambda.size());

 // Find maximum length for broadcasting
 size_t n = std::max({q.n_elem, alpha_vec.n_elem, beta_vec.n_elem,
                     gamma_vec.n_elem, delta_vec.n_elem, lambda_vec.n_elem});

 // Initialize result vector
 arma::vec result(n);

 // Process each element
 for (size_t i = 0; i < n; ++i) {
   // Get parameter values with broadcasting/recycling
   double a = alpha_vec[i % alpha_vec.n_elem];
   double b = beta_vec[i % beta_vec.n_elem];
   double g = gamma_vec[i % gamma_vec.n_elem];
   double d = delta_vec[i % delta_vec.n_elem];
   double l = lambda_vec[i % lambda_vec.n_elem];
   double qi = q[i % q.n_elem];

   // Check parameter validity
   if (!check_pars(a, b, g, d, l)) {
     result(i) = NA_REAL;
     Rcpp::warning("pgkw: invalid parameters at index %d (alpha,beta,gamma>0, delta>=0, lambda>0)", i+1);
     continue;
   }

   // Check domain boundaries
   if (!R_finite(qi) || qi <= 0.0) {
     result(i) = lower_tail ? (log_p ? R_NegInf : 0.0) : (log_p ? 0.0 : 1.0);
     continue;
   }

   if (qi >= 1.0) {
     result(i) = lower_tail ? (log_p ? 0.0 : 1.0) : (log_p ? R_NegInf : 0.0);
     continue;
   }

   // Compute CDF using stable numerical methods
   double log_qi = std::log(qi);
   double log_qi_alpha = a * log_qi;
   double qi_alpha = safe_exp(log_qi_alpha);

   // Calculate (1 - q^alpha) with numerical stability
   double one_minus_qi_alpha;
   if (qi_alpha < 0.5) {
     one_minus_qi_alpha = 1.0 - qi_alpha;
   } else {
     // If close to 1, use expm1 for better precision
     one_minus_qi_alpha = -std::expm1(log_qi_alpha);
   }

   // Boundary checks
   if (one_minus_qi_alpha <= 0.0) {
     result(i) = lower_tail ? (log_p ? 0.0 : 1.0) : (log_p ? R_NegInf : 0.0);
     continue;
   }

   if (one_minus_qi_alpha >= 1.0) {
     result(i) = lower_tail ? (log_p ? R_NegInf : 0.0) : (log_p ? 0.0 : 1.0);
     continue;
   }

   // Calculate (1 - q^alpha)^beta
   double log_oma = std::log(one_minus_qi_alpha);
   double log_oma_beta = b * log_oma;
   double oma_beta = safe_exp(log_oma_beta);

   // Calculate 1 - (1 - q^alpha)^beta
   double term = 1.0 - oma_beta;

   // Boundary checks
   if (term <= 0.0) {
     result(i) = lower_tail ? (log_p ? R_NegInf : 0.0) : (log_p ? 0.0 : 1.0);
     continue;
   }

   if (term >= 1.0) {
     result(i) = lower_tail ? (log_p ? 0.0 : 1.0) : (log_p ? R_NegInf : 0.0);
     continue;
   }

   // Calculate [1 - (1 - q^alpha)^beta]^lambda
   double log_term = std::log(term);
   double log_y = l * log_term;
   double y = safe_exp(log_y);

   // Boundary checks
   if (y <= 0.0) {
     result(i) = lower_tail ? (log_p ? R_NegInf : 0.0) : (log_p ? 0.0 : 1.0);
     continue;
   }

   if (y >= 1.0) {
     result(i) = lower_tail ? (log_p ? 0.0 : 1.0) : (log_p ? R_NegInf : 0.0);
     continue;
   }

   // Use pbeta for the final calculation
   double prob = R::pbeta(y, g, d + 1.0, /*lower_tail=*/true, /*log_p=*/false);

   // Adjust for upper tail if requested
   if (!lower_tail) {
     prob = 1.0 - prob;
   }

   // Convert to log scale if requested
   if (log_p) {
     if (prob <= 0.0) {
       prob = R_NegInf;
     } else if (prob >= 1.0) {
       prob = 0.0;
     } else {
       prob = std::log(prob);
     }
   }

   result(i) = prob;
 }

 return result;
}
//' @title Quantile Function for Generalized Kumaraswamy Distribution
//'
//' @description
//' Calculates the quantile function (inverse CDF) of the Generalized Kumaraswamy
//' distribution GKw(α, β, γ, δ, λ). For a given probability p, returns the value x
//' such that P(X ≤ x) = p, where X follows a GKw distribution with the specified parameters.
//'
//' @param p Vector of probabilities for which to compute quantiles. Values should be in (0,1),
//'        or in (-Inf,0] if log_p=TRUE.
//' @param alpha Shape parameter α > 0 (scalar or vector). Controls the left tail behavior. Default: 1.0.
//' @param beta Shape parameter β > 0 (scalar or vector). Controls the right tail behavior. Default: 1.0.
//' @param gamma Shape parameter γ > 0 (scalar or vector). Affects the central shape. Default: 1.0.
//' @param delta Shape parameter δ ≥ 0 (scalar or vector). Introduces additional flexibility. Default: 0.0.
//' @param lambda Shape parameter λ > 0 (scalar or vector). Controls overall dispersion. Default: 1.0.
//' @param lower_tail Logical; if TRUE (default), probabilities are P(X ≤ x), otherwise P(X > x).
//' @param log_p Logical; if TRUE, probabilities are provided as log(p). Default: FALSE.
//'
//' @return Vector of quantiles corresponding to each probability in p. Returns:
//'   \itemize{
//'     \item 0 for probabilities ≤ 0 (or when log(p) = -Inf if log_p=TRUE)
//'     \item 1 for probabilities ≥ 1 (or when log(p) = 0 if log_p=TRUE)
//'     \item NA for invalid inputs (including p > 1 when log_p=TRUE, which would be log(p) > 0)
//'     \item NA for invalid parameters (any of α, β, γ ≤ 0, δ < 0, or λ ≤ 0)
//'   }
//'
//' @details
//' According to Carrasco et al. (2010), the quantile function for the GKw distribution is
//' derived by inverting the CDF \deqn{F(x) = I_{(1-(1-x^α)^β)^λ}(γ, δ+1)}:
//'
//' \deqn{
//' Q(p) = \{1 - (1 - (I^{-1}_{p}(\gamma, \delta+1))^{1/\lambda})^{1/\beta}\}^{1/\alpha}
//' }
//'
//' where \eqn{I^{-1}_{p}(a, b)} is the inverse regularized incomplete beta function.
//'
//' The implementation follows these steps to compute the quantile:
//' \enumerate{
//'   \item Find y = qbeta(p, γ, δ+1)
//'   \item Compute \deqn{v = y^(1/λ)}
//'   \item Compute \deqn{tmp = 1 - v}
//'   \item Compute \deqn{tmp2 = tmp^(1/β)}
//'   \item Compute \deqn{q = (1 - tmp2)^(1/α)}
//' }
//'
//' Special care is taken at each step to handle numerical issues:
//' \itemize{
//'   \item Boundary cases (p ≤ 0 or p ≥ 1) are handled directly
//'   \item Intermediate results are checked for validity at each transformation step
//'   \item Special cases where parameters equal 1 are optimized for speed and accuracy
//' }
//'
//' This implementation supports vectorized parameters, allowing any combination of
//' vector and scalar inputs for all distribution parameters, with proper broadcasting
//' for vectors of different lengths.
//'
//' @references
//' Kumaraswamy, P. (1980). A generalized probability density function for double-bounded random processes.
//' Journal of Hydrology, 46(1-2), 79-88.
//'
//' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized distributions.
//' Journal of Statistical Computation and Simulation, 81(7), 883-898.
//'
//' @examples
//' \dontrun{
//' # Basic quantile calculation
//' qgkw(0.5, 2, 3, 1, 0, 1)  # Median of GKw(2,3,1,0,1)
//'
//' # Computing multiple quantiles
//' probs <- c(0.1, 0.25, 0.5, 0.75, 0.9)
//' qgkw(probs, 2, 3, 1, 0, 1)  # Various percentiles
//'
//' # Upper tail quantiles
//' qgkw(0.1, 2, 3, 1, 0, 1, lower_tail = FALSE)  # 90th percentile
//'
//' # Log probabilities
//' qgkw(log(0.5), 2, 3, 1, 0, 1, log_p = TRUE)
//'
//' # Vectorized parameters
//' alphas <- c(0.5, 1.0, 2.0)
//' betas <- c(1.0, 2.0, 3.0)
//' qgkw(0.5, alphas, betas)  # Will return 3 quantile values
//'
//' # Verify inverse relationship with pgkw
//' p <- 0.75
//' x <- qgkw(p, 2, 3, 1, 0, 1)
//' pgkw(x, 2, 3, 1, 0, 1)  # Should be approximately 0.75
//' }
//'
//' @export
// [[Rcpp::export]]
arma::vec qgkw(
   const arma::vec& p,
   const Rcpp::NumericVector& alpha = Rcpp::NumericVector::create(1.0),
   const Rcpp::NumericVector& beta = Rcpp::NumericVector::create(1.0),
   const Rcpp::NumericVector& gamma = Rcpp::NumericVector::create(1.0),
   const Rcpp::NumericVector& delta = Rcpp::NumericVector::create(0.0),
   const Rcpp::NumericVector& lambda = Rcpp::NumericVector::create(1.0),
   bool lower_tail = true,
   bool log_p = false
) {
 // Convert NumericVector to arma::vec
 arma::vec alpha_vec(alpha.begin(), alpha.size());
 arma::vec beta_vec(beta.begin(), beta.size());
 arma::vec gamma_vec(gamma.begin(), gamma.size());
 arma::vec delta_vec(delta.begin(), delta.size());
 arma::vec lambda_vec(lambda.begin(), lambda.size());

 // Find maximum length for broadcasting
 size_t n = std::max({p.n_elem, alpha_vec.n_elem, beta_vec.n_elem,
                     gamma_vec.n_elem, delta_vec.n_elem, lambda_vec.n_elem});

 // Initialize result vector
 arma::vec result(n);

 // Process each element
 for (size_t i = 0; i < n; ++i) {
   // Get parameter values with broadcasting/recycling
   double a = alpha_vec[i % alpha_vec.n_elem];
   double b = beta_vec[i % beta_vec.n_elem];
   double g = gamma_vec[i % gamma_vec.n_elem];
   double d = delta_vec[i % delta_vec.n_elem];
   double l = lambda_vec[i % lambda_vec.n_elem];
   double pp = p[i % p.n_elem];

   // Validate parameters
   if (!check_pars(a, b, g, d, l)) {
     result(i) = NA_REAL;
     Rcpp::warning("qgkw: invalid parameters at index %d (alpha,beta,gamma>0, delta>=0, lambda>0)", i+1);
     continue;
   }

   // Process log_p and lower_tail
   if (log_p) {
     if (pp > 0.0) {
       // log(p) > 0 implies p > 1, which is invalid
       result(i) = NA_REAL;
       continue;
     }
     pp = std::exp(pp);  // Convert from log-scale
   }

   if (!lower_tail) {
     pp = 1.0 - pp;  // Convert from upper tail to lower tail
   }

   // Check probability bounds
   if (!R_finite(pp) || pp < 0.0) {
     result(i) = 0.0;  // For p ≤ 0, quantile is 0
     continue;
   }

   if (pp > 1.0) {
     result(i) = 1.0;  // For p > 1, quantile is 1
     continue;
   }

   if (pp <= 0.0) {
     result(i) = 0.0;  // For p = 0, quantile is 0
     continue;
   }

   if (pp >= 1.0) {
     result(i) = 1.0;  // For p = 1, quantile is 1
     continue;
   }

   // Step 1: Find y = qbeta(p, γ, δ+1)
   double y = R::qbeta(pp, g, d + 1.0, /*lower_tail=*/true, /*log_p=*/false);

   // Check for boundary conditions
   if (y <= 0.0) {
     result(i) = 0.0;
     continue;
   }

   if (y >= 1.0) {
     result(i) = 1.0;
     continue;
   }

   // Step 2: Compute v = y^(1/λ)
   double v;
   if (l == 1.0) {
     v = y;  // Optimization for λ=1
   } else {
     v = safe_pow(y, 1.0/l);
   }

   // Step 3: Compute tmp = 1 - v
   double tmp = 1.0 - v;

   // Check for boundary conditions
   if (tmp <= 0.0) {
     result(i) = 1.0;
     continue;
   }

   if (tmp >= 1.0) {
     result(i) = 0.0;
     continue;
   }

   // Step 4: Compute tmp2 = tmp^(1/β)
   double tmp2;
   if (b == 1.0) {
     tmp2 = tmp;  // Optimization for β=1
   } else {
     tmp2 = safe_pow(tmp, 1.0/b);
   }

   // Check for boundary conditions
   if (tmp2 <= 0.0) {
     result(i) = 1.0;
     continue;
   }

   if (tmp2 >= 1.0) {
     result(i) = 0.0;
     continue;
   }

   // Step 5: Compute q = (1 - tmp2)^(1/α)
   double one_minus_tmp2 = 1.0 - tmp2;
   double qq;

   if (one_minus_tmp2 <= 0.0) {
     qq = 0.0;
   } else if (one_minus_tmp2 >= 1.0) {
     qq = 1.0;
   } else if (a == 1.0) {
     qq = one_minus_tmp2;  // Optimization for α=1
   } else {
     qq = safe_pow(one_minus_tmp2, 1.0/a);
   }

   // Final boundary check to ensure result is in (0,1)
   if (qq < 0.0) {
     qq = 0.0;
   } else if (qq > 1.0) {
     qq = 1.0;
   }

   result(i) = qq;
 }

 return result;
}


//' @title Random Number Generation for Generalized Kumaraswamy Distribution
//'
//' @description
//' Generates random deviates from the Generalized Kumaraswamy distribution GKw(α, β, γ, δ, λ).
//' The implementation uses an efficient transformation method based on the Beta distribution
//' as derived in Carrasco et al. (2010).
//'
//' @param n Number of random values to generate. Must be a positive integer.
//' @param alpha Shape parameter α > 0 (scalar or vector). Controls the left tail behavior. Default: 1.0.
//' @param beta Shape parameter β > 0 (scalar or vector). Controls the right tail behavior. Default: 1.0.
//' @param gamma Shape parameter γ > 0 (scalar or vector). Affects the central shape. Default: 1.0.
//' @param delta Shape parameter δ ≥ 0 (scalar or vector). Introduces additional flexibility. Default: 0.0.
//' @param lambda Shape parameter λ > 0 (scalar or vector). Controls overall dispersion. Default: 1.0.
//'
//' @return Vector of length n containing random values from the GKw distribution.
//'   If any parameters are invalid, the function produces an error.
//'   Returns NA for parameter combinations that are invalid.
//'
//' @details
//' According to Carrasco et al. (2010), if V ~ Beta(γ, δ+1), then:
//'
//' \deqn{
//' X = \{1 - [1 - V^{1/\lambda}]^{1/\beta}\}^{1/\alpha}
//' }
//'
//' follows a GKw(α, β, γ, δ, λ) distribution.
//'
//' The random generation algorithm implements this transformation method:
//' \enumerate{
//'   \item Generate V ~ Beta(γ, δ+1)
//'   \item Compute v = V^(1/λ)
//'   \item Compute tmp = 1 - v
//'   \item Compute tmp2 = tmp^(1/β)
//'   \item Compute x = (1 - tmp2)^(1/α)
//' }
//'
//' This implementation includes several optimizations:
//' \itemize{
//'   \item Special cases for α=1, β=1, or λ=1 are handled directly to improve efficiency
//'   \item Boundary cases are checked at each step to maintain numerical stability
//'   \item Safe power transformations prevent numerical issues with extreme values
//'   \item Full support for vectorized parameters with appropriate broadcasting
//' }
//'
//' @references
//' Kumaraswamy, P. (1980). A generalized probability density function for double-bounded random processes.
//' Journal of Hydrology, 46(1-2), 79-88.
//'
//' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized distributions.
//' Journal of Statistical Computation and Simulation, 81(7), 883-898.
//'
//' @examples
//' \dontrun{
//' # Generate 1000 random values from a GKw(2,3,1,0,1)
//' x <- rgkw(1000, 2, 3, 1, 0, 1)
//'
//' # Histogram of generated values
//' hist(x, breaks = 30, probability = TRUE,
//'      main = "Histogram of GKw(2,3,1,0,1) Random Values")
//'
//' # Add the theoretical density curve
//' curve(dgkw(x, 2, 3, 1, 0, 1), add = TRUE, col = "red", lwd = 2)
//'
//' # Comparing empirical and theoretical quantiles
//' prob_points <- seq(0.1, 0.9, by = 0.1)
//' theo_quantiles <- qgkw(prob_points, 2, 3, 1, 0, 1)
//' emp_quantiles <- quantile(x, prob_points)
//' plot(theo_quantiles, emp_quantiles,
//'      main = "Q-Q Plot for GKw(2,3,1,0,1)",
//'      xlab = "Theoretical Quantiles",
//'      ylab = "Empirical Quantiles")
//' abline(0, 1, col = "blue")
//'
//' # Using vectorized parameters
//' alphas <- c(0.5, 1.0, 2.0)
//' samples <- rgkw(300, alphas, 2, 1, 0, 1)  # 100 samples from each configuration
//' hist(samples, breaks = 30, probability = TRUE,
//'      main = "Mixed GKw Samples with Different Alpha Values")
//' }
//'
//' @export
// [[Rcpp::export]]
arma::vec rgkw(
   int n,
   const Rcpp::NumericVector& alpha = Rcpp::NumericVector::create(1.0),
   const Rcpp::NumericVector& beta = Rcpp::NumericVector::create(1.0),
   const Rcpp::NumericVector& gamma = Rcpp::NumericVector::create(1.0),
   const Rcpp::NumericVector& delta = Rcpp::NumericVector::create(0.0),
   const Rcpp::NumericVector& lambda = Rcpp::NumericVector::create(1.0)
) {
 // Convert NumericVector to arma::vec
 arma::vec alpha_vec(alpha.begin(), alpha.size());
 arma::vec beta_vec(beta.begin(), beta.size());
 arma::vec gamma_vec(gamma.begin(), gamma.size());
 arma::vec delta_vec(delta.begin(), delta.size());
 arma::vec lambda_vec(lambda.begin(), lambda.size());

 // Input validation for n
 if (n <= 0) {
   Rcpp::stop("rgkw: n must be a positive integer");
   return arma::vec();
 }

 // Count of parameter combinations for vectorization
 size_t k = std::max({alpha_vec.n_elem, beta_vec.n_elem, gamma_vec.n_elem,
                     delta_vec.n_elem, lambda_vec.n_elem});

 // Initialize result vector
 arma::vec result(n);

 // Process each element
 for (int i = 0; i < n; ++i) {
   // Index for parameter combination (cycling through k combinations)
   size_t idx = i % k;

   // Get parameter values with broadcasting/recycling
   double a = alpha_vec[idx % alpha_vec.n_elem];
   double b = beta_vec[idx % beta_vec.n_elem];
   double g = gamma_vec[idx % gamma_vec.n_elem];
   double d = delta_vec[idx % delta_vec.n_elem];
   double l = lambda_vec[idx % lambda_vec.n_elem];

   // Validate parameters
   if (!check_pars(a, b, g, d, l)) {
     result(i) = NA_REAL;
     Rcpp::warning("rgkw: invalid parameters at index %d (alpha,beta,gamma>0, delta>=0, lambda>0)", idx+1);
     continue;
   }

   // Generate Beta(γ, δ+1) random value
   double vi = R::rbeta(g, d + 1.0);

   // Check for boundary conditions
   if (vi <= 0.0) {
     result(i) = 0.0;
     continue;
   }

   if (vi >= 1.0) {
     result(i) = 1.0;
     continue;
   }

   // Compute v = V^(1/λ)
   double vl;
   if (l == 1.0) {
     vl = vi;  // Optimization for λ=1
   } else {
     vl = safe_pow(vi, 1.0/l);
   }

   // Compute tmp = 1 - v
   double tmp = 1.0 - vl;

   // Check for boundary conditions
   if (tmp <= 0.0) {
     result(i) = 1.0;
     continue;
   }

   if (tmp >= 1.0) {
     result(i) = 0.0;
     continue;
   }

   // Compute tmp2 = tmp^(1/β)
   double tmp2;
   if (b == 1.0) {
     tmp2 = tmp;  // Optimization for β=1
   } else {
     tmp2 = safe_pow(tmp, 1.0/b);
   }

   // Check for boundary conditions
   if (tmp2 <= 0.0) {
     result(i) = 1.0;
     continue;
   }

   if (tmp2 >= 1.0) {
     result(i) = 0.0;
     continue;
   }

   // Compute x = (1 - tmp2)^(1/α)
   double one_minus_tmp2 = 1.0 - tmp2;
   double xx;

   if (one_minus_tmp2 <= 0.0) {
     xx = 0.0;
   } else if (one_minus_tmp2 >= 1.0) {
     xx = 1.0;
   } else if (a == 1.0) {
     xx = one_minus_tmp2;  // Optimization for α=1
   } else {
     xx = safe_pow(one_minus_tmp2, 1.0/a);
   }

   // Final boundary check to ensure result is in (0,1)
   if (xx < 0.0) {
     xx = 0.0;
   } else if (xx > 1.0) {
     xx = 1.0;
   }

   result(i) = xx;
 }

 return result;
}





//' @title Log-Likelihood Function for Generalized Kumaraswamy Distribution
//'
//' @description
//' Calculates the negative log-likelihood for the Generalized Kumaraswamy distribution
//' GKw(α, β, γ, δ, λ) given a vector of observations. This function is primarily used
//' in maximum likelihood estimation and model fitting procedures.
//'
//' @param par NumericVector of length 5 containing parameters (α, β, γ, δ, λ) in that order.
//'        All parameters must be positive.
//' @param data NumericVector of observations, where all values must be in the open interval (0,1).
//'
//' @return The negative log-likelihood value as a double, or -Inf if any parameters or data values
//'         are invalid (parameters ≤ 0 or data outside (0,1)).
//'
//' @details
//' The probability density function of the GKw distribution is:
//'
//' \deqn{
//' f(x; \alpha, \beta, \gamma, \delta, \lambda) =
//'   \frac{\lambda \alpha \beta x^{\alpha-1}(1-x^{\alpha})^{\beta-1}}
//'        {B(\gamma, \delta+1)}
//'   [1-(1-x^{\alpha})^{\beta}]^{\gamma\lambda-1}
//'   [1-[1-(1-x^{\alpha})^{\beta}]^{\lambda}]^{\delta}
//' }
//'
//' The corresponding log-likelihood function is:
//'
//' \deqn{
//' \ell(\theta) = n\ln(\lambda\alpha\beta) - n\ln B(\gamma,\delta+1) +
//'   \sum_{i=1}^{n} [(\alpha-1)\ln(x_i) + (\beta-1)\ln(v_i) + (\gamma\lambda-1)\ln(w_i) + \delta\ln(z_i)]
//' }
//'
//' where:
//' \itemize{
//'   \item \deqn{v_i = 1 - x_i^{\alpha}}
//'   \item \deqn{w_i = 1 - v_i^{\beta} = 1 - (1-x_i^{\alpha})^{\beta}}
//'   \item \deqn{z_i = 1 - w_i^{\lambda} = 1 - (1-(1-x_i^{\alpha})^{\beta})^{\lambda}}
//' }
//'
//' The implementation uses several techniques for numerical stability:
//' \itemize{
//'   \item R's lbeta function for calculating log(B(γ, δ+1))
//'   \item log1p for calculating log(1-x) when x is close to zero
//'   \item Careful handling of intermediate values to avoid underflow/overflow
//' }
//'
//' Note that this function returns the negative log-likelihood (rather than the log-likelihood),
//' making it suitable for minimization in optimization procedures.
//'
//' @examples
//' \dontrun{
//' # Generate sample data from a GKw distribution
//' set.seed(123)
//' x <- rgkw(100, 2, 3, 1.0, 0.5, 0.5)
//' hist(x, breaks = 20, main = "GKw(2, 3, 1.0, 0.5, 0.5) Sample")
//'
//' # Use in optimization with Hessian-based methods
//' result <- optim(c(0.5, 0.5, 0.5, 0.5, 0.5), llgkw, method = "BFGS",
//'                 hessian = TRUE, data = x)
//'
//' # Compare numerical and analytical derivatives
//' num_grad <- numDeriv::grad(llgkw, x = result$par, data = x)
//' num_hess <- numDeriv::hessian(llgkw, x = result$par, data = x)
//'
//' ana_grad <- grgkw(result$par, data = x)
//' ana_hess <- hsgkw(result$par, data = x)
//'
//' # Check differences (should be very small)
//' round(num_grad - ana_grad, 4)
//' round(num_hess - ana_hess, 4)
//'
//' }
//'
//' @seealso
//' \code{\link[gkwreg]{dgkw}} for the GKw density function,
//' \code{\link[gkwreg]{pgkw}} for the GKw cumulative distribution function,
//' \code{\link[gkwreg]{hsgkw}} for the analytic Hessian of the GKw log-likelihood,
//'
//' @references
//' Kumaraswamy, P. (1980). A generalized probability density function for double-bounded random processes.
//' Journal of Hydrology, 46(1-2), 79-88.
//'
//' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized distributions.
//' Journal of Statistical Computation and Simulation, 81(7), 883-898.
//'
//' @export
// [[Rcpp::export]]
double llgkw(const Rcpp::NumericVector& par, const Rcpp::NumericVector& data) {
 // Parameter extraction
 double alpha = par[0];   // Shape parameter α > 0
 double beta = par[1];    // Shape parameter β > 0
 double gamma = par[2];   // Shape parameter γ > 0
 double delta = par[3];   // Shape parameter δ > 0
 double lambda = par[4];  // Shape parameter λ > 0

 // Parameter validation - all parameters must be positive
 if (alpha <= 0 || beta <= 0 || gamma <= 0 || delta <= 0 || lambda <= 0) {
   return R_NegInf;  // Return negative infinity for invalid parameters
 }

 // Convert data to arma::vec for more efficient operations
 // Use aliasing (false) to avoid copying the data
 // arma::vec x(data.begin(), data.size(), false);
 arma::vec x = Rcpp::as<arma::vec>(data);

 // Data validation - all values must be in the range (0,1)
 if (arma::any(x <= 0) || arma::any(x >= 1)) {
   return R_NegInf;  // Return negative infinity for invalid data
 }

 int n = x.n_elem;  // Sample size

 // Calculate log of Beta function for constant term
 // Use R's lbeta function for better numerical stability
 double log_beta_term = R::lbeta(gamma, delta + 1);

 // Calculate the constant term: n*log(λαβ/B(γ,δ+1))
 double constant_term = n * (std::log(lambda) + std::log(alpha) + std::log(beta) - log_beta_term);

 // Calculate log(x) and sum (α-1)*log(x) terms
 arma::vec log_x = arma::log(x);
 double term1 = arma::sum((alpha - 1.0) * log_x);

 // Calculate v = 1-x^α and sum (β-1)*log(v) terms
 arma::vec x_alpha = arma::pow(x, alpha);
 arma::vec v = 1.0 - x_alpha;
 arma::vec log_v = arma::log(v);
 double term2 = arma::sum((beta - 1.0) * log_v);

 // Calculate w = 1-v^β = 1-(1-x^α)^β and sum (γλ-1)*log(w) terms
 arma::vec v_beta = arma::pow(v, beta);
 arma::vec w = 1.0 - v_beta;

 // Handle numerical stability for log(w) when w is close to zero
 arma::vec log_w(n);
 for (int i = 0; i < n; i++) {
   if (w(i) < 1e-10) {
     // Use log1p for numerical stability: log(w) = log(1-v^β) = log1p(-v^β)
     log_w(i) = std::log1p(-v_beta(i));
   } else {
     log_w(i) = std::log(w(i));
   }
 }
 double term3 = arma::sum((gamma * lambda - 1.0) * log_w);

 // Calculate z = 1-w^λ = 1-[1-(1-x^α)^β]^λ and sum δ*log(z) terms
 arma::vec w_lambda = arma::pow(w, lambda);
 arma::vec z = 1.0 - w_lambda;

 // Handle numerical stability for log(z) when z is close to zero
 arma::vec log_z(n);
 for (int i = 0; i < n; i++) {
   if (z(i) < 1e-10) {
     // Use log1p for numerical stability: log(z) = log(1-w^λ) = log1p(-w^λ)
     log_z(i) = std::log1p(-w_lambda(i));
   } else {
     log_z(i) = std::log(z(i));
   }
 }
 double term4 = arma::sum(delta * log_z);

 // Return final minus-log-likelihood: constant term + sum of all individual terms
 return -(constant_term + term1 + term2 + term3 + term4);
}



//' @title Gradient Function for Generalized Kumaraswamy Log-Likelihood
//'
//' @description
//' Calculates the gradient vector (partial derivatives) of the negative log-likelihood
//' function for the Generalized Kumaraswamy (GKw) distribution. This function provides
//' the exact gradient needed for efficient optimization in maximum likelihood estimation.
//'
//' @param par NumericVector of length 5 containing parameters (α, β, γ, δ, λ) in that order.
//'        All parameters must be positive.
//' @param data NumericVector of observations, where all values must be in the open interval (0,1).
//'
//' @return NumericVector of length 5 containing the gradient components (partial derivatives)
//'         of the negative log-likelihood with respect to each parameter (α, β, γ, δ, λ).
//'         Returns a vector of NaN values if any parameters or data values are invalid.
//'
//' @details
//' The gradient vector contains the following partial derivatives of the negative log-likelihood:
//'
//' \deqn{
//' \frac{\partial \ell}{\partial \alpha} = \frac{n}{\alpha} + \sum_{i=1}^{n}\log(x_i) -
//' \sum_{i=1}^{n}\left[x_i^{\alpha} \log(x_i) \left(\frac{\beta-1}{v_i} -
//' \frac{(\gamma\lambda-1) \beta v_i^{\beta-1}}{w_i} +
//' \frac{\delta \lambda \beta v_i^{\beta-1} w_i^{\lambda-1}}{z_i}\right)\right]
//' }
//'
//' \deqn{
//' \frac{\partial \ell}{\partial \beta} = \frac{n}{\beta} + \sum_{i=1}^{n}\log(v_i) -
//' \sum_{i=1}^{n}\left[v_i^{\beta} \log(v_i) \left(\frac{\gamma\lambda-1}{w_i} -
//' \frac{\delta \lambda w_i^{\lambda-1}}{z_i}\right)\right]
//' }
//'
//' \deqn{
//' \frac{\partial \ell}{\partial \gamma} = -n[\psi(\gamma) - \psi(\gamma+\delta+1)] +
//' \lambda\sum_{i=1}^{n}\log(w_i)
//' }
//'
//' \deqn{
//' \frac{\partial \ell}{\partial \delta} = -n[\psi(\delta+1) - \psi(\gamma+\delta+1)] +
//' \sum_{i=1}^{n}\log(z_i)
//' }
//'
//' \deqn{
//' \frac{\partial \ell}{\partial \lambda} = \frac{n}{\lambda} +
//' \gamma\sum_{i=1}^{n}\log(w_i) - \delta\sum_{i=1}^{n}\frac{w_i^{\lambda}\log(w_i)}{z_i}
//' }
//'
//' where:
//' \itemize{
//'   \item \deqn{v_i = 1 - x_i^{\alpha}}
//'   \item \deqn{w_i = 1 - v_i^{\beta} = 1 - (1-x_i^{\alpha})^{\beta}}
//'   \item \deqn{z_i = 1 - w_i^{\lambda} = 1 - (1-(1-x_i^{\alpha})^{\beta})^{\lambda}}
//'   \item \deqn{\psi} is the digamma function (derivative of the log-gamma function)
//' }
//'
//' The implementation includes several numerical safeguards:
//' \itemize{
//'   \item Parameter and data validation with appropriate error handling
//'   \item Clamping of intermediate values to avoid numerical underflow/overflow
//'   \item Efficient vector operations using Armadillo C++ library
//' }
//'
//' The returned gradient is negated to align with minimization of negative log-likelihood
//' in optimization routines.
//'
//' @examples
//' \dontrun{
//' # Generate sample data from a GKw distribution
//' set.seed(123)
//' x <- rgkw(100, 2, 3, 1.0, 0.5, 0.5)
//' hist(x, breaks = 20, main = "GKw(2, 3, 1.0, 0.5, 0.5) Sample")
//'
//' # Use in optimization with Hessian-based methods
//' result <- optim(c(0.5, 0.5, 0.5, 0.5, 0.5), llgkw, method = "BFGS",
//'                 hessian = TRUE, data = x)
//'
//' # Compare numerical and analytical derivatives
//' num_grad <- numDeriv::grad(llgkw, x = result$par, data = x)
//' num_hess <- numDeriv::hessian(llgkw, x = result$par, data = x)
//'
//' ana_grad <- grgkw(result$par, data = x)
//' ana_hess <- hsgkw(result$par, data = x)
//'
//' # Check differences (should be very small)
//' round(num_grad - ana_grad, 4)
//' round(num_hess - ana_hess, 4)
//'
//' }
//'
//' @seealso
//' \code{\link[gkwreg]{llgkw}} for the negative log-likelihood function,
//' \code{\link[gkwreg]{hsgkw}} for the Hessian matrix of the GKw log-likelihood,
//' \code{\link[gkwreg]{dgkw}} for the GKw density function,
//'
//' @references
//' Kumaraswamy, P. (1980). A generalized probability density function for double-bounded random processes.
//' Journal of Hydrology, 46(1-2), 79-88.
//'
//' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized distributions.
//' Journal of Statistical Computation and Simulation, 81(7), 883-898.
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector grgkw(const Rcpp::NumericVector& par, const Rcpp::NumericVector& data) {
 // Parameter extraction
 double alpha = par[0];   // Shape parameter α > 0
 double beta = par[1];    // Shape parameter β > 0
 double gamma = par[2];   // Shape parameter γ > 0
 double delta = par[3];   // Shape parameter δ > 0
 double lambda = par[4];  // Shape parameter λ > 0

 // Parameter validation
 if (alpha <= 0 || beta <= 0 || gamma <= 0 || delta <= 0 || lambda <= 0) {
   Rcpp::NumericVector grad(5, R_NaN);
   return grad;
 }

 // Data conversion and validation
 arma::vec x = Rcpp::as<arma::vec>(data);

 if (arma::any(x <= 0) || arma::any(x >= 1)) {
   Rcpp::NumericVector grad(5, R_NaN);
   return grad;
 }

 int n = x.n_elem;  // Sample size

 // Initialize gradient vector
 Rcpp::NumericVector grad(5, 0.0);

 // Small constant to avoid numerical issues
 double eps = std::numeric_limits<double>::epsilon() * 100;

 // Compute transformations and intermediate values
 arma::vec log_x = arma::log(x);                // log(x_i)
 arma::vec x_alpha = arma::pow(x, alpha);       // x_i^α
 arma::vec x_alpha_log_x = x_alpha % log_x;     // x_i^α * log(x_i)

 // v_i = 1 - x_i^α
 arma::vec v = 1.0 - x_alpha;
 v = arma::clamp(v, eps, 1.0 - eps);            // Prevent numerical issues

 arma::vec log_v = arma::log(v);                // log(v_i)
 arma::vec v_beta_m1 = arma::pow(v, beta - 1.0); // v_i^(β-1)
 arma::vec v_beta = arma::pow(v, beta);          // v_i^β
 arma::vec v_beta_log_v = v_beta % log_v;        // v_i^β * log(v_i)

 // w_i = 1 - v_i^β = 1 - (1-x_i^α)^β
 arma::vec w = 1.0 - v_beta;
 w = arma::clamp(w, eps, 1.0 - eps);            // Prevent numerical issues

 arma::vec log_w = arma::log(w);                // log(w_i)
 arma::vec w_lambda_m1 = arma::pow(w, lambda - 1.0); // w_i^(λ-1)
 arma::vec w_lambda = arma::pow(w, lambda);          // w_i^λ
 arma::vec w_lambda_log_w = w_lambda % log_w;        // w_i^λ * log(w_i)

 // z_i = 1 - w_i^λ = 1 - [1-(1-x_i^α)^β]^λ
 arma::vec z = 1.0 - w_lambda;
 z = arma::clamp(z, eps, 1.0 - eps);            // Prevent numerical issues

 arma::vec log_z = arma::log(z);                // log(z_i)

 // Calculate partial derivatives for each parameter (for log-likelihood)

 // ∂ℓ/∂α = n/α + Σᵢlog(xᵢ) - Σᵢ[xᵢ^α * log(xᵢ) * ((β-1)/vᵢ - (γλ-1) * β * vᵢ^(β-1) / wᵢ + δ * λ * β * vᵢ^(β-1) * wᵢ^(λ-1) / zᵢ)]
 double d_alpha = n / alpha + arma::sum(log_x);

 // Calculate the complex term in the α gradient
 arma::vec alpha_term2 = (beta - 1.0) / v;                // (β-1)/v_i
 arma::vec alpha_term3 = (gamma * lambda - 1.0) * beta * v_beta_m1 / w;  // (γλ-1) * β * v_i^(β-1) / w_i
 arma::vec alpha_term4 = delta * lambda * beta * v_beta_m1 % w_lambda_m1 / z;  // δ * λ * β * v_i^(β-1) * w_i^(λ-1) / z_i

 d_alpha -= arma::sum(x_alpha_log_x % (alpha_term2 - alpha_term3 + alpha_term4));

 // ∂ℓ/∂β = n/β + Σᵢlog(vᵢ) - Σᵢ[vᵢ^β * log(vᵢ) * ((γλ-1) / wᵢ - δ * λ * wᵢ^(λ-1) / zᵢ)]
 double d_beta = n / beta + arma::sum(log_v);

 // Calculate the complex term in the β gradient
 arma::vec beta_term2 = (gamma * lambda - 1.0) / w;       // (γλ-1) / w_i
 arma::vec beta_term3 = delta * lambda * w_lambda_m1 / z; // δ * λ * w_i^(λ-1) / z_i

 d_beta -= arma::sum(v_beta_log_v % (beta_term2 - beta_term3));

 // ∂ℓ/∂γ = -n[ψ(γ) - ψ(γ+δ+1)] + λΣᵢlog(wᵢ)
 double d_gamma = -n * (R::digamma(gamma) - R::digamma(gamma + delta + 1)) + lambda * arma::sum(log_w);

 // ∂ℓ/∂δ = -n[ψ(δ+1) - ψ(γ+δ+1)] + Σᵢlog(zᵢ)
 double d_delta = -n * (R::digamma(delta + 1) - R::digamma(gamma + delta + 1)) + arma::sum(log_z);

 // ∂ℓ/∂λ = n/λ + γΣᵢlog(wᵢ) - δΣᵢ[(wᵢ^λ*log(wᵢ))/zᵢ]
 double d_lambda = n / lambda + gamma * arma::sum(log_w) - delta * arma::sum(w_lambda_log_w / z);

 // Since we're optimizing negative log-likelihood, negate all derivatives
 grad[0] = -d_alpha;
 grad[1] = -d_beta;
 grad[2] = -d_gamma;
 grad[3] = -d_delta;
 grad[4] = -d_lambda;

 return grad;
}

//' @title Analytic Hessian Matrix for Generalized Kumaraswamy Distribution
//'
//' @description
//' Computes the analytic Hessian matrix of the log-likelihood function for
//' the Generalized Kumaraswamy (GKw) distribution. This function provides
//' exact second derivatives needed for optimization and inference.
//'
//' @param par Numeric vector of length 5 containing the parameters
//'        (α, β, γ, δ, λ) in that order. All parameters must be positive.
//' @param data Numeric vector of observations, where all values must be
//'        in the open interval (0,1).
//'
//' @return A 5×5 numeric matrix representing the Hessian of the negative
//'         log-likelihood function. If parameters or data are invalid
//'         (parameters ≤ 0 or data outside (0,1)), returns a matrix of
//'         NaN values.
//'
//' @details
//' The log-likelihood for the generalized Kumaraswamy distribution is:
//'
//' \deqn{
//' \ell(\theta) = n \ln(\lambda) + n \ln(\alpha) + n \ln(\beta) - n \ln B(\gamma, \delta+1)
//' + (\alpha-1) \sum \ln(x_i)
//' + (\beta-1) \sum \ln(1 - x_i^\alpha)
//' + (\gamma\lambda - 1) \sum \ln\{1 - (1 - x_i^\alpha)^\beta\}
//' + \delta \sum \ln\{1 - \{1 - (1 - x_i^\alpha)^\beta\}^\lambda\}
//' }
//'
//' where B refers to the Beta function.
//'
//' The implementation computes all second derivatives analytically for each term.
//' For computational efficiency, the following transformations are used:
//' \itemize{
//'   \item \deqn{A = x^α} and derivatives
//'   \item \deqn{v = 1 - A}
//'   \item \deqn{w = 1 - v^β}
//'   \item \deqn{z = 1 - w^λ}
//' }
//'
//' The returned Hessian matrix has the following structure:
//' \itemize{
//'   \item Rows/columns 1-5 correspond to α, β, γ, δ, λ respectively
//'   \item The matrix is symmetric (as expected for a Hessian)
//'   \item The matrix represents second derivatives of the negative log-likelihood
//' }
//'
//' This function is implemented in C++ for computational efficiency.
//'
//' @examples
//' \dontrun{
//' # Generate sample data from a GKw distribution
//' set.seed(123)
//' x <- rgkw(100, 2, 3, 1.0, 0.5, 0.5)
//' hist(x, breaks = 20, main = "GKw(2, 3, 1.0, 0.5, 0.5) Sample")
//'
//' # Use in optimization with Hessian-based methods
//' result <- optim(c(0.5, 0.5, 0.5, 0.5, 0.5), llgkw, method = "BFGS",
//'                 hessian = TRUE, data = x)
//'
//' # Compare numerical and analytical derivatives
//' num_grad <- numDeriv::grad(llgkw, x = result$par, data = x)
//' num_hess <- numDeriv::hessian(llgkw, x = result$par, data = x)
//'
//' ana_grad <- grgkw(result$par, data = x)
//' ana_hess <- hsgkw(result$par, data = x)
//'
//' # Check differences (should be very small)
//' round(num_grad - ana_grad, 4)
//' round(num_hess - ana_hess, 4)
//'
//' }
//'
//' @seealso
//' \code{\link[gkwreg]{dgkw}} for the GKw density function,
//' \code{\link[gkwreg]{gkwreg}} for fitting GKw regression models,
//' \code{\link[gkwreg]{pgkw}} for the GKw cumulative distribution function
//'
//' @references
//' Kumaraswamy, P. (1980). A generalized probability density function for double-bounded random processes.
//' Journal of Hydrology, 46(1-2), 79-88.
//'
//' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized distributions.
//' Journal of Statistical Computation and Simulation, 81(7), 883-898.
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix hsgkw(const Rcpp::NumericVector& par, const Rcpp::NumericVector& data) {
 // Parameter extraction
 double alpha  = par[0];   // θ[0] = α
 double beta   = par[1];   // θ[1] = β
 double gamma  = par[2];   // θ[2] = γ
 double delta  = par[3];   // θ[3] = δ
 double lambda = par[4];   // θ[4] = λ

 // Simple parameter validation (all > 0)
 if(alpha <= 0 || beta <= 0 || gamma <= 0 || delta <= 0 || lambda <= 0) {
   Rcpp::NumericMatrix nanH(5,5);
   nanH.fill(R_NaN);
   return nanH;
 }

 // Data conversion and basic validation
 arma::vec x = Rcpp::as<arma::vec>(data);
 if(arma::any(x <= 0) || arma::any(x >= 1)) {
   Rcpp::NumericMatrix nanH(5,5);
   nanH.fill(R_NaN);
   return nanH;
 }

 int n = x.n_elem;  // sample size

 // Initialize Hessian matrix H (of ℓ(θ)) as 5x5
 arma::mat H(5,5, arma::fill::zeros);

 // --- CONSTANT TERMS (do not depend on x) ---
 // L1: n ln(λ)  => d²/dλ² = -n/λ²
 H(4,4) += - n/(lambda*lambda);
 // L2: n ln(α)  => d²/dα² = -n/α²
 H(0,0) += - n/(alpha*alpha);
 // L3: n ln(β)  => d²/dβ² = -n/β²
 H(1,1) += - n/(beta*beta);
 // L4: - n ln[B(γ, δ+1)]
 //   d²/dγ² = -n [ψ₁(γ) - ψ₁(γ+δ+1)]  where ψ₁ is the trigamma function
 H(2,2) += - n * ( R::trigamma(gamma) - R::trigamma(gamma+delta+1) );
 //   d²/dδ² = -n [ψ₁(δ+1) - ψ₁(γ+δ+1)]
 H(3,3) += - n * ( R::trigamma(delta+1) - R::trigamma(gamma+delta+1) );
 //   Mixed derivative (γ,δ): = n ψ₁(γ+δ+1)
 H(2,3) += n * R::trigamma(gamma+delta+1);
 H(3,2) = H(2,3);
 // L5: (α-1) Σ ln(x_i)  --> contributes only to first derivatives

 // Accumulators for mixed derivatives with λ
 double acc_gamma_lambda = 0.0;  // Sum of ln(w)
 double acc_delta_lambda = 0.0;  // Sum of dz_dlambda / z
 double acc_alpha_lambda = 0.0;  // For α,λ contributions
 double acc_beta_lambda = 0.0;   // For β,λ contributions

 // --- TERMS THAT INVOLVE THE OBSERVATIONS ---
 // Loop over each observation to accumulate contributions from:
 // L6: (β-1) Σ ln(v), where v = 1 - x^α
 // L7: (γλ-1) Σ ln(w), where w = 1 - v^β
 // L8: δ Σ ln(z), where z = 1 - w^λ
 for (int i = 0; i < n; i++) {
   double xi    = x(i);
   double ln_xi = std::log(xi);

   // -- Compute A = x^α and its derivatives --
   double A = std::pow(xi, alpha);                  // A = x^α
   double dA_dalpha = A * ln_xi;                    // dA/dα = x^α ln(x)
   double d2A_dalpha2 = A * ln_xi * ln_xi;          // d²A/dα² = x^α (ln(x))²

   // -- v = 1 - A and its derivatives --
   double v = 1.0 - A;                              // v = 1 - x^α
   double ln_v = std::log(v);                       // ln(v)
   double dv_dalpha = -dA_dalpha;                   // dv/dα = -dA/dα = -x^α ln(x)
   double d2v_dalpha2 = -d2A_dalpha2;               // d²v/dα² = -d²A/dα² = -x^α (ln(x))²

   // --- L6: (β-1) ln(v) ---
   // Second derivative w.r.t. α: (β-1)*[(d²v/dα²*v - (dv/dα)²)/v²]
   double d2L6_dalpha2 = (beta - 1.0) * ((d2v_dalpha2 * v - dv_dalpha * dv_dalpha) / (v*v));
   // Mixed derivative: d²L6/(dα dβ) = d/dβ[(β-1)*(dv_dalpha/v)] = (dv_dalpha/v)
   double d2L6_dalpha_dbeta = dv_dalpha / v;

   // --- L7: (γλ - 1) ln(w), where w = 1 - v^β ---
   double v_beta = std::pow(v, beta);               // v^β
   double w = 1.0 - v_beta;                         // w = 1 - v^β
   double ln_w = std::log(w);                       // ln(w)
   // Derivative of w w.r.t. v: dw/dv = -β * v^(β-1)
   double dw_dv = -beta * std::pow(v, beta - 1.0);
   // Chain rule: dw/dα = dw/dv * dv/dα
   double dw_dalpha = dw_dv * dv_dalpha;
   // Second derivative w.r.t. α for L7:
   // d²/dα² ln(w) = [d²w/dα² * w - (dw/dα)²] / w²
   // Computing d²w/dα²:
   //   dw/dα = -β * v^(β-1)*dv_dalpha,
   //   d²w/dα² = -β * [(β-1)*v^(β-2)*(dv_dalpha)² + v^(β-1)*d²v_dalpha²]
   double d2w_dalpha2 = -beta * ((beta - 1.0) * std::pow(v, beta-2.0) * (dv_dalpha * dv_dalpha)
                                   + std::pow(v, beta-1.0) * d2v_dalpha2);
   double d2L7_dalpha2 = (gamma * lambda - 1.0) * ((d2w_dalpha2 * w - (dw_dalpha * dw_dalpha)) / (w*w));
   // Derivative w.r.t. β: d/dβ ln(w). Note: d/dβ(v^β) = v^β ln(v) => d/dβ w = -v^β ln(v)
   double dw_dbeta = -v_beta * ln_v;
   // Second derivative w.r.t. β for L7:
   // d²/dβ² ln(w) = [d²w/dβ² * w - (dw/dβ)²]/w², where d²w/dβ² = -v^β (ln(v))²
   double d2w_dbeta2 = -v_beta * (ln_v * ln_v);
   double d2L7_dbeta2 = (gamma * lambda - 1.0) * ((d2w_dbeta2 * w - (dw_dbeta * dw_dbeta))/(w*w));
   // Mixed derivative L7 (α,β): d²/(dα dβ) ln(w) =
   //   = d/dβ[(dw_dalpha)/w] = (d/dβ dw_dalpha)/w - (dw_dalpha*dw_dbeta)/(w*w)
   // Approximate d/dβ dw_dalpha:
   double d_dw_dalpha_dbeta = -std::pow(v, beta-1.0) * (1.0 + beta * ln_v) * dv_dalpha;
   double d2L7_dalpha_dbeta = (gamma * lambda - 1.0) * ((d_dw_dalpha_dbeta / w) - (dw_dalpha * dw_dbeta)/(w*w));

   // --- L8: δ ln(z), where z = 1 - w^λ ---
   double w_lambda_val = std::pow(w, lambda);       // w^λ
   double z = 1.0 - w_lambda_val;                   // z = 1 - w^λ
   // Derivative w.r.t. α: dz/dα = -λ * w^(λ-1) * dw/dα
   double dz_dalpha = -lambda * std::pow(w, lambda-1.0) * dw_dalpha;
   // Second derivative w.r.t. α for L8:
   // d²z/dα² = -λ * [(λ-1)*w^(λ-2)*(dw/dα)² + w^(λ-1)*d²w/dα²]
   double d2z_dalpha2 = -lambda * ((lambda - 1.0) * std::pow(w, lambda-2.0) * (dw_dalpha*dw_dalpha)
                                     + std::pow(w, lambda-1.0) * d2w_dalpha2);
   double d2L8_dalpha2 = delta * ((d2z_dalpha2 * z - dz_dalpha*dz_dalpha)/(z*z));

   // Derivative w.r.t. β: dz/dβ = -λ * w^(λ-1) * dw/dβ
   double dz_dbeta = -lambda * std::pow(w, lambda-1.0) * dw_dbeta;
   // Second derivative w.r.t. β for L8:
   // d²z/dβ² = -λ * [(λ-1)*w^(λ-2)*(dw/dβ)² + w^(λ-1)*d²w/dβ²]
   double d2z_dbeta2 = -lambda * ((lambda - 1.0) * std::pow(w, lambda-2.0) * (dw_dbeta*dw_dbeta)
                                    + std::pow(w, lambda-1.0) * d2w_dbeta2);
   double d2L8_dbeta2 = delta * ((d2z_dbeta2 * z - dz_dbeta*dz_dbeta)/(z*z));

   // Mixed derivative L8 (α,β): d²/(dα dβ) ln(z)
   // = (d/dβ dz_dalpha)/z - (dz_dalpha*dz_dbeta)/(z*z)
   // Approximate d/dβ dz_dalpha = -λ * [(λ-1)*w^(λ-2)*(dw_dβ*dw_dα) + w^(λ-1)*(d/dβ dw_dalpha)]
   double d_dw_dalpha_dbeta_2 = -lambda * ((lambda - 1.0) * std::pow(w, lambda-2.0) * dw_dbeta * dw_dalpha
                                             + std::pow(w, lambda-1.0) * d_dw_dalpha_dbeta);
   double d2L8_dalpha_dbeta = delta * ((d_dw_dalpha_dbeta_2 / z) - (dz_dalpha*dz_dbeta)/(z*z));

   // Derivatives of L8 with respect to λ:
   // d/dλ ln(z) = (1/z)*dz/dλ, with dz/dλ = -w^λ ln(w)
   double dz_dlambda = -w_lambda_val * ln_w;
   // d²/dλ² ln(z) = [d²z/dλ² * z - (dz_dlambda)²]/z² (assuming w constant in λ)
   double d2z_dlambda2 = -w_lambda_val * (ln_w * ln_w);
   double d2L8_dlambda2 = delta * ((d2z_dlambda2 * z - dz_dlambda*dz_dlambda)/(z*z));

   // Mixed derivative L8 (α,λ): d²/(dα dλ) ln(z) = (d/dα dz_dλ)/z - (dz_dλ*dz_dalpha)/(z*z)
   // Correct formula: sum of two terms, not multiplication
   double d_dalpha_dz_dlambda = -std::pow(w, lambda-1.0) * dw_dalpha -
     lambda * ln_w * std::pow(w, lambda-1.0) * dw_dalpha;
   double d2L8_dalpha_dlambda = delta * ((d_dalpha_dz_dlambda / z) - (dz_dlambda*dz_dalpha)/(z*z));

   // Mixed derivative L8 (β,λ): d²/(dβ dλ) ln(z) = (d/dβ dz_dλ)/z - (dz_dlambda*dz_dbeta)/(z*z)
   // Correct formula: sum of two terms, not multiplication
   double d_dbeta_dz_dlambda = -std::pow(w, lambda-1.0) * dw_dbeta -
     lambda * ln_w * std::pow(w, lambda-1.0) * dw_dbeta;
   double d2L8_dbeta_dlambda = delta * ((d_dbeta_dz_dlambda / z) - (dz_dlambda*dz_dbeta)/(z*z));

   // --- ACCUMULATING CONTRIBUTIONS TO THE HESSIAN MATRIX ---
   // Index: 0 = α, 1 = β, 2 = γ, 3 = δ, 4 = λ

   // H(α,α): sum of L2, L6, L7, and L8 (constants already added)
   H(0,0) += d2L6_dalpha2 + d2L7_dalpha2 + d2L8_dalpha2;

   // H(α,β): mixed from L6, L7, and L8
   H(0,1) += d2L6_dalpha_dbeta + d2L7_dalpha_dbeta + d2L8_dalpha_dbeta;
   H(1,0) = H(0,1);

   // H(β,β): contributions from L3, L7, and L8
   H(1,1) += d2L7_dbeta2 + d2L8_dbeta2;

   // H(λ,λ): contains L1 and L8 (L1 already added)
   H(4,4) += d2L8_dlambda2;

   // H(γ,α): from L7 - derivative of ln(w) in α multiplied by λ factor of (γλ-1)
   H(2,0) += lambda * (dw_dalpha / w);
   H(0,2) = H(2,0);

   // H(γ,β): from L7 - derivative of ln(w) in β multiplied by λ
   H(2,1) += lambda * (dw_dbeta / w);
   H(1,2) = H(2,1);

   // H(δ,α): L8 - mixed derivative: d/dα ln(z)
   H(3,0) += dz_dalpha / z;
   H(0,3) = H(3,0);

   // H(δ,β): L8 - d/dβ ln(z)
   H(3,1) += dz_dbeta / z;
   H(1,3) = H(3,1);

   // Accumulating terms for mixed derivatives with λ
   // (α,λ): Term from L7 (γ contribution) + term from L8 (δ contribution)
   double term1_alpha_lambda = gamma * (dw_dalpha / w);
   double term2_alpha_lambda = d2L8_dalpha_dlambda;
   acc_alpha_lambda += term1_alpha_lambda + term2_alpha_lambda;

   // (β,λ): Term from L7 (γ contribution) + term from L8 (δ contribution)
   double term1_beta_lambda = gamma * (dw_dbeta / w);
   double term2_beta_lambda = d2L8_dbeta_dlambda;
   acc_beta_lambda += term1_beta_lambda + term2_beta_lambda;

   // (γ,λ): Contribution from L7 (γλ-1)*ln(w)
   acc_gamma_lambda += ln_w;

   // (δ,λ): Contribution from L8 δ*ln(z)
   acc_delta_lambda += dz_dlambda / z;
 } // end of loop

 // Applying mixed derivatives with λ
 // Note: All signs are positive for log-likelihood (not negative log-likelihood)

 // H(α,λ): Positive sign for log-likelihood
 H(0,4) = acc_alpha_lambda;
 H(4,0) = H(0,4);

 // H(β,λ): Positive sign for log-likelihood
 H(1,4) = acc_beta_lambda;
 H(4,1) = H(1,4);

 // H(γ,λ): Positive sign for log-likelihood
 H(2,4) = acc_gamma_lambda;
 H(4,2) = H(2,4);

 // H(δ,λ): Positive sign for log-likelihood
 H(3,4) = acc_delta_lambda;
 H(4,3) = H(3,4);

 // Returns the analytic Hessian matrix of the log-likelihood
 return Rcpp::wrap(-H);
}


// //' @title Analytic Hessian Matrix for Generalized Kumaraswamy Distribution
// //'
// //' @description
// //' Computes the analytic Hessian matrix of the log-likelihood function for
// //' the Generalized Kumaraswamy (GKw) distribution. This function provides
// //' exact second derivatives needed for optimization and inference.
// //'
// //' @param par Numeric vector of length 5 containing the parameters
// //'        (α, β, γ, δ, λ) in that order. All parameters must be positive.
// //' @param data Numeric vector of observations, where all values must be
// //'        in the open interval (0,1).
// //'
// //' @return A 5×5 numeric matrix representing the Hessian of the negative
// //'         log-likelihood function. If parameters or data are invalid
// //'         (parameters ≤ 0 or data outside (0,1)), returns a matrix of
// //'         NaN values.
// //'
// //' @details
// //' The log-likelihood for the generalized Kumaraswamy distribution is:
// //'
// //' \deqn{
// //' \ell(\theta) = n \ln(\lambda) + n \ln(\alpha) + n \ln(\beta) - n \ln B(\gamma, \delta+1)
// //' + (\alpha-1) \sum \ln(x_i)
// //' + (\beta-1) \sum \ln(1 - x_i^\alpha)
// //' + (\gamma\lambda - 1) \sum \ln\{1 - (1 - x_i^\alpha)^\beta\}
// //' + \delta \sum \ln\{1 - \{1 - (1 - x_i^\alpha)^\beta\}^\lambda\}
// //' }
// //'
// //' where B refers to the Beta function.
// //'
// //' The implementation computes all second derivatives analytically for each term.
// //' For computational efficiency, the following transformations are used:
// //' \itemize{
// //'   \item \deqn{A = x^α} and derivatives
// //'   \item \deqn{v = 1 - A}
// //'   \item \deqn{w = 1 - v^β}
// //'   \item \deqn{z = 1 - w^λ}
// //' }
// //'
// //' The returned Hessian matrix has the following structure:
// //' \itemize{
// //'   \item Rows/columns 1-5 correspond to α, β, γ, δ, λ respectively
// //'   \item The matrix is symmetric (as expected for a Hessian)
// //'   \item The matrix represents second derivatives of the negative log-likelihood
// //' }
// //'
// //' This function is implemented in C++ for computational efficiency.
// //'
// //' @examples
// //' \dontrun{
// //' # Generate sample data from a GKw distribution
// //' set.seed(123)
// //' x <- rgkw(100, 2, 3, 1.0, 0.5, 0.5)
// //' hist(x, breaks = 20, main = "GKw(2, 3, 1.0, 0.5, 0.5) Sample")
// //'
// //' # Use in optimization with Hessian-based methods
// //' result <- optim(c(0.5, 0.5, 0.5, 0.5, 0.5), llgkw, method = "BFGS",
// //'                 hessian = TRUE, data = x)
// //'
// //' # Compare numerical and analytical derivatives
// //' num_grad <- numDeriv::grad(llgkw, x = result$par, data = x)
// //' num_hess <- numDeriv::hessian(llgkw, x = result$par, data = x)
// //'
// //' ana_grad <- grgkw(result$par, data = x)
// //' ana_hess <- hsgkw(result$par, data = x)
// //'
// //' # Check differences (should be very small)
// //' round(num_grad - ana_grad, 4)
// //' round(num_hess - ana_hess, 4)
// //'
// //' }
// //'
// //' @seealso
// //' \code{\link[gkwreg]{dgkw}} for the GKw density function,
// //' \code{\link[gkwreg]{gkwreg}} for fitting GKw regression models,
// //' \code{\link[gkwreg]{pgkw}} for the GKw cumulative distribution function
// //'
// //' @references
// //' Kumaraswamy, P. (1980). A generalized probability density function for double-bounded random processes.
// //' Journal of Hydrology, 46(1-2), 79-88.
// //'
// //' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized distributions.
// //' Journal of Statistical Computation and Simulation, 81(7), 883-898.
// //'
// //' @export
// // [[Rcpp::export]]
// Rcpp::NumericMatrix hsgkw(const Rcpp::NumericVector& par, const Rcpp::NumericVector& data) {
//  // Parameter extraction
//  double alpha  = par[0];   // θ[0] = α
//  double beta   = par[1];   // θ[1] = β
//  double gamma  = par[2];   // θ[2] = γ
//  double delta  = par[3];   // θ[3] = δ
//  double lambda = par[4];   // θ[4] = λ
//
//  // Simple parameter validation (all > 0)
//  if(alpha <= 0 || beta <= 0 || gamma <= 0 || delta <= 0 || lambda <= 0) {
//    Rcpp::NumericMatrix nanH(5,5);
//    nanH.fill(R_NaN);
//    return nanH;
//  }
//
//  // Data conversion and basic validation
//  arma::vec x = Rcpp::as<arma::vec>(data);
//  if(arma::any(x <= 0) || arma::any(x >= 1)) {
//    Rcpp::NumericMatrix nanH(5,5);
//    nanH.fill(R_NaN);
//    return nanH;
//  }
//
//  int n = x.n_elem;  // sample size
//
//  // Initialize Hessian matrix H (of ℓ(θ)) as 5x5
//  arma::mat H(5,5, arma::fill::zeros);
//
//  // --- CONSTANT TERMS (do not depend on x) ---
//  // L1: n ln(λ)  => d²/dλ² = -n/λ²
//  H(4,4) += - n/(lambda*lambda);
//  // L2: n ln(α)  => d²/dα² = -n/α²
//  H(0,0) += - n/(alpha*alpha);
//  // L3: n ln(β)  => d²/dβ² = -n/β²
//  H(1,1) += - n/(beta*beta);
//  // L4: - n ln[B(γ, δ+1)]
//  //   d²/dγ² = -n [ψ₁(γ) - ψ₁(γ+δ+1)]  where ψ₁ is the trigamma function
//  H(2,2) += - n * ( R::trigamma(gamma) - R::trigamma(gamma+delta+1) );
//  //   d²/dδ² = -n [ψ₁(δ+1) - ψ₁(γ+δ+1)]
//  H(3,3) += - n * ( R::trigamma(delta+1) - R::trigamma(gamma+delta+1) );
//  //   Mixed derivative (γ,δ): = n ψ₁(γ+δ+1)
//  H(2,3) += n * R::trigamma(gamma+delta+1);
//  H(3,2) = H(2,3);
//  // L5: (α-1) Σ ln(x_i)  --> contributes only to first derivatives
//
//  // Accumulators for mixed derivatives with λ
//  double acc_gamma_lambda = 0.0;  // Sum of ln(w)
//  double acc_delta_lambda = 0.0;  // Sum of dz_dlambda / z
//  double acc_alpha_lambda = 0.0;  // For α,λ contributions
//  double acc_beta_lambda = 0.0;   // For β,λ contributions
//
//  // --- TERMS THAT INVOLVE THE OBSERVATIONS ---
//  // Loop over each observation to accumulate contributions from:
//  // L6: (β-1) Σ ln(v), where v = 1 - x^α
//  // L7: (γλ-1) Σ ln(w), where w = 1 - v^β
//  // L8: δ Σ ln(z), where z = 1 - w^λ
//  for (int i = 0; i < n; i++) {
//    double xi    = x(i);
//    double ln_xi = std::log(xi);
//
//    // -- Compute A = x^α and its derivatives --
//    double A = std::pow(xi, alpha);                  // A = x^α
//    double dA_dalpha = A * ln_xi;                    // dA/dα = x^α ln(x)
//    double d2A_dalpha2 = A * ln_xi * ln_xi;          // d²A/dα² = x^α (ln(x))²
//
//    // -- v = 1 - A and its derivatives --
//    double v = 1.0 - A;                              // v = 1 - x^α
//    double ln_v = std::log(v);                       // ln(v)
//    double dv_dalpha = -dA_dalpha;                   // dv/dα = -dA/dα = -x^α ln(x)
//    double d2v_dalpha2 = -d2A_dalpha2;               // d²v/dα² = -d²A/dα² = -x^α (ln(x))²
//
//    // --- L6: (β-1) ln(v) ---
//    // First derivative w.r.t. α: (β-1) * (1/v)*dv_dalpha
//    double dL6_dalpha = (beta - 1.0) * (dv_dalpha / v);
//    // Second derivative w.r.t. α: (β-1)*[(d²v/dα²*v - (dv/dα)²)/v²]
//    double d2L6_dalpha2 = (beta - 1.0) * ((d2v_dalpha2 * v - dv_dalpha * dv_dalpha) / (v*v));
//    // First derivative w.r.t. β: dL6/dβ = ln(v)
//    double dL6_dbeta = ln_v;
//    // Mixed derivative: d²L6/(dα dβ) = d/dβ[(β-1)*(dv_dalpha/v)] = (dv_dalpha/v)
//    double d2L6_dalpha_dbeta = dv_dalpha / v;
//
//    // --- L7: (γλ - 1) ln(w), where w = 1 - v^β ---
//    double v_beta = std::pow(v, beta);               // v^β
//    double w = 1.0 - v_beta;                         // w = 1 - v^β
//    double ln_w = std::log(w);                       // ln(w)
//    // Derivative of w w.r.t. v: dw/dv = -β * v^(β-1)
//    double dw_dv = -beta * std::pow(v, beta - 1.0);
//    // Chain rule: dw/dα = dw/dv * dv/dα
//    double dw_dalpha = dw_dv * dv_dalpha;
//    // First derivative w.r.t. α: d/dα ln(w) = (1/w)*dw_dalpha
//    double dL7_dalpha = (gamma * lambda - 1.0) * (dw_dalpha / w);
//    // Second derivative w.r.t. α for L7:
//    // d²/dα² ln(w) = [d²w/dα² * w - (dw/dα)²] / w²
//    // Computing d²w/dα²:
//    //   dw/dα = -β * v^(β-1)*dv_dalpha,
//    //   d²w/dα² = -β * [(β-1)*v^(β-2)*(dv_dalpha)² + v^(β-1)*d²v_dalpha²]
//    double d2w_dalpha2 = -beta * ((beta - 1.0) * std::pow(v, beta-2.0) * (dv_dalpha * dv_dalpha)
//                                    + std::pow(v, beta-1.0) * d2v_dalpha2);
//    double d2L7_dalpha2 = (gamma * lambda - 1.0) * ((d2w_dalpha2 * w - (dw_dalpha * dw_dalpha)) / (w*w));
//    // Derivative w.r.t. β: d/dβ ln(w). Note: d/dβ(v^β) = v^β ln(v) => d/dβ w = -v^β ln(v)
//    double dw_dbeta = -v_beta * ln_v;
//    double dL7_dbeta = (gamma * lambda - 1.0) * (dw_dbeta / w);
//    // Second derivative w.r.t. β for L7:
//    // d²/dβ² ln(w) = [d²w/dβ² * w - (dw/dβ)²]/w², where d²w/dβ² = -v^β (ln(v))²
//    double d2w_dbeta2 = -v_beta * (ln_v * ln_v);
//    double d2L7_dbeta2 = (gamma * lambda - 1.0) * ((d2w_dbeta2 * w - (dw_dbeta * dw_dbeta))/(w*w));
//    // Mixed derivative L7 (α,β): d²/(dα dβ) ln(w) =
//    //   = d/dβ[(dw_dalpha)/w] = (d/dβ dw_dalpha)/w - (dw_dalpha*dw_dbeta)/(w*w)
//    // Approximate d/dβ dw_dalpha:
//    double d_dw_dalpha_dbeta = -std::pow(v, beta-1.0) * (1.0 + beta * ln_v) * dv_dalpha;
//    double d2L7_dalpha_dbeta = (gamma * lambda - 1.0) * ((d_dw_dalpha_dbeta / w) - (dw_dalpha * dw_dbeta)/(w*w));
//
//    // --- L8: δ ln(z), where z = 1 - w^λ ---
//    double w_lambda_val = std::pow(w, lambda);       // w^λ
//    double z = 1.0 - w_lambda_val;                   // z = 1 - w^λ
//    double ln_z = std::log(z);                       // ln(z)
//    // Derivative w.r.t. α: dz/dα = -λ * w^(λ-1) * dw/dα
//    double dz_dalpha = -lambda * std::pow(w, lambda-1.0) * dw_dalpha;
//    double dL8_dalpha = delta * (dz_dalpha / z);
//    // Second derivative w.r.t. α for L8:
//    // d²z/dα² = -λ * [(λ-1)*w^(λ-2)*(dw/dα)² + w^(λ-1)*d²w/dα²]
//    double d2w_dalpha2_dummy = d2w_dalpha2; // already calculated for L7
//    double d2z_dalpha2 = -lambda * ((lambda - 1.0) * std::pow(w, lambda-2.0) * (dw_dalpha*dw_dalpha)
//                                      + std::pow(w, lambda-1.0) * d2w_dalpha2_dummy);
//    double d2L8_dalpha2 = delta * ((d2z_dalpha2 * z - dz_dalpha*dz_dalpha)/(z*z));
//
//    // Derivative w.r.t. β: dz/dβ = -λ * w^(λ-1) * dw/dβ
//    double dz_dbeta = -lambda * std::pow(w, lambda-1.0) * dw_dbeta;
//    double dL8_dbeta = delta * (dz_dbeta / z);
//    // Second derivative w.r.t. β for L8:
//    // d²z/dβ² = -λ * [(λ-1)*w^(λ-2)*(dw/dβ)² + w^(λ-1)*d²w/dβ²]
//    double d2z_dbeta2 = -lambda * ((lambda - 1.0) * std::pow(w, lambda-2.0) * (dw_dbeta*dw_dbeta)
//                                     + std::pow(w, lambda-1.0) * d2w_dbeta2);
//    double d2L8_dbeta2 = delta * ((d2z_dbeta2 * z - dz_dbeta*dz_dbeta)/(z*z));
//
//    // Mixed derivative L8 (α,β): d²/(dα dβ) ln(z)
//    // = (d/dβ dz_dalpha)/z - (dz_dalpha*dz_dbeta)/(z*z)
//    // Approximate d/dβ dz_dalpha = -λ * [(λ-1)*w^(λ-2)*(dw_dβ*dw_dα) + w^(λ-1)*(d/dβ dw_dalpha)]
//    double d_dw_dalpha_dbeta_2 = -lambda * ((lambda - 1.0) * std::pow(w, lambda-2.0) * dw_dbeta * dw_dalpha
//                                              + std::pow(w, lambda-1.0) * d_dw_dalpha_dbeta);
//    double d2L8_dalpha_dbeta = delta * ((d_dw_dalpha_dbeta_2 / z) - (dz_dalpha*dz_dbeta)/(z*z));
//
//    // Derivatives of L8 with respect to λ:
//    // d/dλ ln(z) = (1/z)*dz/dλ, with dz/dλ = -w^λ ln(w)
//    double dz_dlambda = -w_lambda_val * ln_w;
//    double dL8_dlambda = delta * (dz_dlambda / z);
//    // d²/dλ² ln(z) = [d²z/dλ² * z - (dz_dlambda)²]/z² (assuming w constant in λ)
//    double d2z_dlambda2 = -w_lambda_val * (ln_w * ln_w);
//    double d2L8_dlambda2 = delta * ((d2z_dlambda2 * z - dz_dlambda*dz_dlambda)/(z*z));
//
//    // Mixed derivative L8 (α,λ): d²/(dα dλ) ln(z) = (d/dα dz_dλ)/z - (dz_dλ*dz_dalpha)/(z*z)
//    // Correct formula: sum of two terms, not multiplication
//    double d_dalpha_dz_dlambda = -std::pow(w, lambda-1.0) * dw_dalpha -
//      lambda * ln_w * std::pow(w, lambda-1.0) * dw_dalpha;
//    double d2L8_dalpha_dlambda = delta * ((d_dalpha_dz_dlambda / z) - (dz_dlambda*dz_dalpha)/(z*z));
//
//    // Mixed derivative L8 (β,λ): d²/(dβ dλ) ln(z) = (d/dβ dz_dλ)/z - (dz_dlambda*dz_dbeta)/(z*z)
//    // Correct formula: sum of two terms, not multiplication
//    double d_dbeta_dz_dlambda = -std::pow(w, lambda-1.0) * dw_dbeta -
//      lambda * ln_w * std::pow(w, lambda-1.0) * dw_dbeta;
//    double d2L8_dbeta_dlambda = delta * ((d_dbeta_dz_dlambda / z) - (dz_dlambda*dz_dbeta)/(z*z));
//
//    // --- ACCUMULATING CONTRIBUTIONS TO THE HESSIAN MATRIX ---
//    // Index: 0 = α, 1 = β, 2 = γ, 3 = δ, 4 = λ
//
//    // H(α,α): sum of L2, L6, L7, and L8 (constants already added)
//    H(0,0) += d2L6_dalpha2 + d2L7_dalpha2 + d2L8_dalpha2;
//
//    // H(α,β): mixed from L6, L7, and L8
//    H(0,1) += d2L6_dalpha_dbeta + d2L7_dalpha_dbeta + d2L8_dalpha_dbeta;
//    H(1,0) = H(0,1);
//
//    // H(β,β): contributions from L3, L7, and L8
//    H(1,1) += d2L7_dbeta2 + d2L8_dbeta2;
//
//    // H(λ,λ): contains L1 and L8 (L1 already added)
//    H(4,4) += d2L8_dlambda2;
//
//    // H(γ,α): from L7 - derivative of ln(w) in α multiplied by λ factor of (γλ-1)
//    H(2,0) += lambda * (dw_dalpha / w);
//    H(0,2) = H(2,0);
//
//    // H(γ,β): from L7 - derivative of ln(w) in β multiplied by λ
//    H(2,1) += lambda * (dw_dbeta / w);
//    H(1,2) = H(2,1);
//
//    // H(δ,α): L8 - mixed derivative: d/dα ln(z)
//    H(3,0) += dz_dalpha / z;
//    H(0,3) = H(3,0);
//
//    // H(δ,β): L8 - d/dβ ln(z)
//    H(3,1) += dz_dbeta / z;
//    H(1,3) = H(3,1);
//
//    // Accumulating terms for mixed derivatives with λ
//    // (α,λ): Term from L7 (γ contribution) + term from L8 (δ contribution)
//    double term1_alpha_lambda = gamma * (dw_dalpha / w);
//    double term2_alpha_lambda = d2L8_dalpha_dlambda;
//    acc_alpha_lambda += term1_alpha_lambda + term2_alpha_lambda;
//
//    // (β,λ): Term from L7 (γ contribution) + term from L8 (δ contribution)
//    double term1_beta_lambda = gamma * (dw_dbeta / w);
//    double term2_beta_lambda = d2L8_dbeta_dlambda;
//    acc_beta_lambda += term1_beta_lambda + term2_beta_lambda;
//
//    // (γ,λ): Contribution from L7 (γλ-1)*ln(w)
//    acc_gamma_lambda += ln_w;
//
//    // (δ,λ): Contribution from L8 δ*ln(z)
//    acc_delta_lambda += dz_dlambda / z;
//  } // end of loop
//
//  // Applying mixed derivatives with λ
//  // Note: All signs are positive for log-likelihood (not negative log-likelihood)
//
//  // H(α,λ): Positive sign for log-likelihood
//  H(0,4) = acc_alpha_lambda;
//  H(4,0) = H(0,4);
//
//  // H(β,λ): Positive sign for log-likelihood
//  H(1,4) = acc_beta_lambda;
//  H(4,1) = H(1,4);
//
//  // H(γ,λ): Positive sign for log-likelihood
//  H(2,4) = acc_gamma_lambda;
//  H(4,2) = H(2,4);
//
//  // H(δ,λ): Positive sign for log-likelihood
//  H(3,4) = acc_delta_lambda;
//  H(4,3) = H(3,4);
//
//  // Returns the analytic Hessian matrix of the log-likelihood
//  return Rcpp::wrap(-H);
// }






/*
----------------------------------------------------------------------------
REUSE OF NUMERIC STABILITY FUNCTIONS AND CHECKS
----------------------------------------------------------------------------
NOTE: We assume the following inline functions are already available, exactly as in gkwdist.cpp:
- log1mexp(double)
- log1pexp(double)
- safe_log(double)
- safe_exp(double)
- safe_pow(double,double)
etc.

The KwKw distribution here sets γ=1 in the GKw(α,β,γ=1,δ,λ). So we only keep α>0, β>0,
δ≥0, λ>0, ignoring γ. We'll define a small parameter checker:

check_kkw_pars(alpha,beta,delta,lambda) => boolean
*/

// -----------------------------------------------------------------------------
// Parameter Checker for KwKw Distribution
// KwKw(α, β, 1, δ, λ) => alpha>0, beta>0, delta≥0, λ>0
// -----------------------------------------------------------------------------
inline bool check_kkw_pars(double alpha,
                         double beta,
                         double delta,
                         double lambda,
                         bool strict = false) {
if (alpha <= 0.0 || beta <= 0.0 || delta < 0.0 || lambda <= 0.0) {
  return false;
}
if (strict) {
  const double MINP = 1e-5;
  const double MAXP = 1e5;
  if (alpha < MINP || beta < MINP || lambda < MINP) {
    return false;
  }
  if (alpha > MAXP || beta > MAXP || delta > MAXP || lambda > MAXP) {
    return false;
  }
}
return true;
}

/*
----------------------------------------------------------------------------
KUMARASWAMY-KUMARASWAMY DISTRIBUTION  KwKw(α, β, 1, δ, λ)
----------------------------------------------------------------------------
PDF:
f(x) = λ α β (δ+1) x^(α-1) (1 - x^α)^(β-1)
[1 - (1 - x^α)^β]^(λ - 1)
{1 - [1 - (1 - x^α)^β]^λ}^δ ,  0 < x < 1.

CDF:
F(x) = 1 - { 1 - [1 - (1 - x^α)^β]^λ }^(δ+1).

QUANTILE (inverse CDF):
Solve F(x)=p => x = ...
We get
y = [1 - (1 - x^α)^β]^λ,
F(x)=1 - (1-y)^(δ+1),
=> (1-y)^(δ+1) = 1-p
=> y = 1 - (1-p)^(1/(δ+1))
=> (1 - x^α)^β = 1 - y
=> x^α = 1 - (1-y)^(1/β)
=> x = [1 - (1-y)^(1/β)]^(1/α).
with y = 1 - (1-p)^(1/(δ+1)) all raised to 1/λ in the general GKw, but here it's directly [1 - (1-p)^(1/(δ+1))] since γ=1. Actually we must be consistent with the formula from the article. Let's confirm carefully:

The table in the user's message says:
F(x) = 1 - [1 - (1 - x^α)^β]^λ)^(δ+1),
=> (1 - [1-(1-x^α)^β]^λ)^(δ+1) = 1 - p
=> 1 - [1-(1-x^α)^β]^λ = 1 - p^(1/(δ+1)) is not correct. We must do it carefully:

F(x)=1 - [1 - y]^ (δ+1) with y=[1-(1 - x^α)^β]^λ.
=> [1 - y]^(δ+1) = 1 - p => 1-y = (1 - p)^(1/(δ+1)) => y=1 - (1-p)^(1/(δ+1)).
Then y^(1/λ) if we had a general GKw, but here "y" itself is already [1-(1-x^α)^β]^λ. So to invert that we do y^(1/λ). Indeed, so that part is needed because the exponent λ is still free. So let's define:

y = [1 - (1 - x^α)^β]^λ
=> y^(1/λ) = 1 - (1 - x^α)^β
=> (1 - x^α)^β = 1 - y^(1/λ)

Then (1 - x^α)= [1 - y^(1/λ)]^(1/β).
=> x^α=1 - [1 - y^(1/λ)]^(1/β).
=> x=[1 - [1 - y^(1/λ)]^(1/β)]^(1/α).

So the quantile formula is indeed:

QKwKw(p)= [ 1 - [ 1 - ( 1 - (1-p)^(1/(δ+1)) )^(1/λ }^(1/β ) ]^(1/α).

We'll code it carefully.

RNG:
In the user's table, the recommended approach is:
V ~ Uniform(0,1)
U = 1 - (1 - V)^(1/(δ+1))    ( that is the portion for the (δ+1) exponent )
X= [1 - [1 - U^(1/λ}]^(1/β)]^(1/α)

LOG-LIKELIHOOD:
log f(x) = log(λ) + log(α) + log(β) + log(δ+1)
+ (α-1)*log(x)
+ (β-1)*log(1 - x^α)
+ (λ-1)*log(1 - (1 - x^α)^β)
+ δ* log(1 - [1-(1 - x^α)^β]^λ).
Then sum across data, multiply n to the constants, etc.
*/


// -----------------------------------------------------------------------------
// 1) dkkw: PDF of KwKw
// -----------------------------------------------------------------------------

//' @title Density of the KwKw Distribution
//'
//' @description
//' Computes the PDF of the KwKw(\eqn{\alpha, \beta, \delta, \lambda}) distribution, which
//' is the GKw distribution restricted to \eqn{\gamma=1}.
//'
//' @param x Vector of quantiles in \eqn{(0,1)}.
//' @param alpha \eqn{\alpha > 0}.
//' @param beta \eqn{\beta > 0}.
//' @param delta \eqn{\delta \ge 0}.
//' @param lambda \eqn{\lambda > 0}.
//' @param log_prob Logical; if \code{TRUE}, returns log-PDF, else PDF.
//'
//' @return Vector of densities (or log-densities) of the same length as the broadcast of inputs.
//'
//' @details
//' The PDF is
//' \deqn{
//'   f(x) = \lambda \,\alpha \,\beta \,(\delta + 1)\; x^{\alpha - 1}\, (1 - x^\alpha)^{\beta - 1}\,
//'          \bigl[1 - (1 - x^\alpha)^\beta\bigr]^{\lambda - 1}\,
//'          \bigl\{1 - \bigl[1 - (1 - x^\alpha)^\beta\bigr]^\lambda\bigr\}^{\delta},
//'   \quad 0 < x < 1.
//' }
//'
//' @export
// [[Rcpp::export]]
arma::vec dkkw(
   const arma::vec& x,
   const Rcpp::NumericVector& alpha,
   const Rcpp::NumericVector& beta,
   const Rcpp::NumericVector& delta,
   const Rcpp::NumericVector& lambda,
   bool log_prob = false
) {
 arma::vec a_vec(alpha.begin(), alpha.size());
 arma::vec b_vec(beta.begin(), beta.size());
 arma::vec d_vec(delta.begin(), delta.size());
 arma::vec l_vec(lambda.begin(), lambda.size());

 // broadcast
 size_t N = std::max({x.n_elem,
                     a_vec.n_elem,
                     b_vec.n_elem,
                     d_vec.n_elem,
                     l_vec.n_elem});

 arma::vec out(N);
 out.fill(log_prob ? R_NegInf : 0.0);

 for (size_t i = 0; i < N; ++i) {
   double a = a_vec[i % a_vec.n_elem];
   double b = b_vec[i % b_vec.n_elem];
   double dd = d_vec[i % d_vec.n_elem];
   double ll = l_vec[i % l_vec.n_elem];
   double xx = x[i % x.n_elem];

   if (!check_kkw_pars(a, b, dd, ll)) {
     // invalid => 0 or -Inf
     continue;
   }
   // domain
   if (xx <= 0.0 || xx >= 1.0 || !R_finite(xx)) {
     continue;
   }

   // Precompute logs for PDF
   // log(f(x)) = log(λ) + log(α) + log(β) + log(δ+1)
   //            + (α-1)*log(x)
   //            + (β-1)*log(1 - x^α)
   //            + (λ-1)*log(1 - (1 - x^α)^β)
   //            + δ* log(1 - [1 - (1 - x^α)^β]^λ)
   double logCst = std::log(ll) + std::log(a) + std::log(b) + std::log(dd + 1.0);

   // x^alpha
   double lx = std::log(xx);
   double log_xalpha = a * lx;  // log(x^alpha)= alpha*log(x)
   // 1 - x^alpha
   double log_1_minus_xalpha = log1mexp(log_xalpha); // stable
   if (!R_finite(log_1_minus_xalpha)) {
     continue;
   }

   // (β-1)* log(1 - x^alpha)
   double term1 = (b - 1.0) * log_1_minus_xalpha;

   // let A = (1 - x^alpha)^β => logA = b * log_1_minus_xalpha
   double logA = b * log_1_minus_xalpha;
   double log_1_minusA = log1mexp(logA); // stable => log(1 - A)
   if (!R_finite(log_1_minusA)) {
     continue;
   }
   // (λ-1)* log(1 - A)
   double term2 = (ll - 1.0) * log_1_minusA;

   // let B = [1 - (1 - x^alpha)^β]^λ => logB = λ*log(1 - A)
   double logB = ll * log_1_minusA;
   double log_1_minus_B = log1mexp(logB);
   if (!R_finite(log_1_minus_B)) {
     continue;
   }
   // δ * log( 1 - B )
   double term3 = dd * log_1_minus_B;

   double log_pdf = logCst
   + (a - 1.0)*lx
   + term1
   + term2
   + term3;

   if (log_prob) {
     out(i) = log_pdf;
   } else {
     out(i) = std::exp(log_pdf);
   }
 }

 return out;
}


// -----------------------------------------------------------------------------
// 2) pkkw: CDF of KwKw
// -----------------------------------------------------------------------------

//' @title CDF of the KwKw Distribution
//'
//' @description
//' Computes the CDF for KwKw(\eqn{\alpha, \beta, \delta, \lambda}), i.e. P(X≤x).
//'
//' @param q Vector of quantiles in \eqn{(0,1)}.
//' @param alpha \eqn{\alpha > 0}.
//' @param beta \eqn{\beta > 0}.
//' @param delta \eqn{\delta \ge 0}.
//' @param lambda \eqn{\lambda > 0}.
//' @param lower_tail Logical; if TRUE (default), returns F(q)=P(X≤q), else 1-F(q).
//' @param log_p Logical; if TRUE, returns log of the probability.
//'
//' @details
//' The CDF is
//' \deqn{
//'   F(x) = 1 - \bigl\{1 - \bigl[1 - (1 - x^\alpha)^\beta\bigr]^\lambda\bigr\}^{\delta + 1}.
//' }
//'
//' @export
// [[Rcpp::export]]
arma::vec pkkw(
   const arma::vec& q,
   const Rcpp::NumericVector& alpha,
   const Rcpp::NumericVector& beta,
   const Rcpp::NumericVector& delta,
   const Rcpp::NumericVector& lambda,
   bool lower_tail = true,
   bool log_p = false
) {
 arma::vec a_vec(alpha.begin(), alpha.size());
 arma::vec b_vec(beta.begin(), beta.size());
 arma::vec d_vec(delta.begin(), delta.size());
 arma::vec l_vec(lambda.begin(), lambda.size());

 size_t N = std::max({q.n_elem,
                     a_vec.n_elem,
                     b_vec.n_elem,
                     d_vec.n_elem,
                     l_vec.n_elem});

 arma::vec out(N);

 for (size_t i = 0; i < N; ++i) {
   double a = a_vec[i % a_vec.n_elem];
   double b = b_vec[i % b_vec.n_elem];
   double dd = d_vec[i % d_vec.n_elem];
   double ll = l_vec[i % l_vec.n_elem];
   double xx = q[i % q.n_elem];

   if (!check_kkw_pars(a, b, dd, ll)) {
     out(i) = NA_REAL;
     continue;
   }

   // boundaries
   if (!R_finite(xx) || xx <= 0.0) {
     // F(0) = 0
     double val0 = (lower_tail ? 0.0 : 1.0);
     out(i) = log_p ? std::log(val0) : val0;
     continue;
   }
   if (xx >= 1.0) {
     // F(1)=1
     double val1 = (lower_tail ? 1.0 : 0.0);
     out(i) = log_p ? std::log(val1) : val1;
     continue;
   }

   // x^alpha
   double lx = std::log(xx);
   double log_xalpha = a * lx;
   double xalpha = std::exp(log_xalpha);

   double one_minus_xalpha = 1.0 - xalpha;
   if (one_minus_xalpha <= 0.0) {
     // near 1 => F ~ 1
     double val1 = (lower_tail ? 1.0 : 0.0);
     out(i) = log_p ? std::log(val1) : val1;
     continue;
   }
   // (1 - x^alpha)^beta => ...
   double vbeta = std::pow(one_minus_xalpha, b);
   double y = 1.0 - vbeta;  // [1-(1 - x^alpha)^β]
   if (y <= 0.0) {
     // => F=0
     double val0 = (lower_tail ? 0.0 : 1.0);
     out(i) = log_p ? std::log(val0) : val0;
     continue;
   }
   if (y >= 1.0) {
     // => F=1
     double val1 = (lower_tail ? 1.0 : 0.0);
     out(i) = log_p ? std::log(val1) : val1;
     continue;
   }

   double ylambda = std::pow(y, ll);   // [1-(1-x^alpha)^β]^λ
   if (ylambda <= 0.0) {
     // => F=0
     double val0 = (lower_tail ? 0.0 : 1.0);
     out(i) = log_p ? std::log(val0) : val0;
     continue;
   }
   if (ylambda >= 1.0) {
     // => F=1
     double val1 = (lower_tail ? 1.0 : 0.0);
     out(i) = log_p ? std::log(val1) : val1;
     continue;
   }

   double outer = 1.0 - ylambda; // 1 - ...
   double cdfval = 1.0 - std::pow(outer, dd+1.0);

   if (!lower_tail) {
     cdfval = 1.0 - cdfval;
   }
   if (log_p) {
     cdfval = std::log(cdfval);
   }
   out(i) = cdfval;
 }

 return out;
}


// -----------------------------------------------------------------------------
// 3) qkkw: Quantile of KwKw
// -----------------------------------------------------------------------------

//' @title Quantile Function of the KwKw Distribution
//'
//' @description
//' Computes the quantiles for KwKw(\eqn{\alpha, \beta, \delta, \lambda}).
//'
//' @param p Vector of probabilities in (0,1).
//' @param alpha \eqn{\alpha > 0}.
//' @param beta \eqn{\beta > 0}.
//' @param delta \eqn{\delta \ge 0}.
//' @param lambda \eqn{\lambda > 0}.
//' @param lower_tail Logical; if TRUE, p is F(x). If FALSE, p is 1-F(x).
//' @param log_p Logical; if TRUE, p is given as log(p).
//'
//' @details
//' Invert the CDF
//' \eqn{
//'   F(x)= 1 - [1 - (1 - x^α)^β]^λ)^(δ+1).
//' }
//' The resulting formula is
//' \deqn{
//'   Q(p)= \Bigl\{1 - \Bigl[1 - \Bigl(1 - p\Bigr)^{\frac{1}{\delta+1}}\Bigr]^{\frac{1}{\lambda}}\Bigr\}^{\frac{1}{\beta}}\Bigr\}^{\frac{1}{\alpha}}.
//' }
//'
//' @export
// [[Rcpp::export]]
arma::vec qkkw(
   const arma::vec& p,
   const Rcpp::NumericVector& alpha,
   const Rcpp::NumericVector& beta,
   const Rcpp::NumericVector& delta,
   const Rcpp::NumericVector& lambda,
   bool lower_tail = true,
   bool log_p = false
) {
 arma::vec a_vec(alpha.begin(), alpha.size());
 arma::vec b_vec(beta.begin(), beta.size());
 arma::vec d_vec(delta.begin(), delta.size());
 arma::vec l_vec(lambda.begin(), lambda.size());

 size_t N = std::max({p.n_elem,
                     a_vec.n_elem,
                     b_vec.n_elem,
                     d_vec.n_elem,
                     l_vec.n_elem});

 arma::vec out(N);

 for (size_t i = 0; i < N; ++i) {
   double a = a_vec[i % a_vec.n_elem];
   double b = b_vec[i % b_vec.n_elem];
   double dd = d_vec[i % d_vec.n_elem];
   double ll = l_vec[i % l_vec.n_elem];
   double pp = p[i % p.n_elem];

   if (!check_kkw_pars(a, b, dd, ll)) {
     out(i) = NA_REAL;
     continue;
   }

   // Convert p if log_p
   if (log_p) {
     if (pp > 0.0) {
       // log(p)>0 => p>1 => invalid
       out(i) = NA_REAL;
       continue;
     }
     pp = std::exp(pp);
   }
   // if upper tail
   if (!lower_tail) {
     pp = 1.0 - pp;
   }

   // boundaries
   if (pp <= 0.0) {
     out(i) = 0.0;
     continue;
   }
   if (pp >= 1.0) {
     out(i) = 1.0;
     continue;
   }

   // formula:
   // F(x)=p => 1 - [1 - (1 - x^α)^β]^λ]^(δ+1) = p
   // => [1 - (1 - x^α)^β]^λ = 1 - (1-p)^(1/(δ+1))
   double tmp1 = 1.0 - std::pow(1.0 - pp, 1.0/(dd+1.0));
   if (tmp1 < 0.0)  tmp1=0.0;
   if (tmp1>1.0)    tmp1=1.0;  // safety

   // let T= tmp1^(1/λ)
   double T;
   if (ll==1.0) {
     T=tmp1;
   } else {
     T=std::pow(tmp1, 1.0/ll);
   }
   // => (1 - x^α)^β = 1 - T
   double M=1.0 - T;
   if (M<0.0)  M=0.0;
   if (M>1.0)  M=1.0;
   // => 1 - x^α= M^(1/β)
   double Mpow= std::pow(M, 1.0/b);
   double xalpha=1.0 - Mpow;
   if (xalpha<0.0)  xalpha=0.0;
   if (xalpha>1.0)  xalpha=1.0;

   // x= xalpha^(1/α) => actually => x= [1 - M^(1/β)]^(1/α)
   double xx;
   if (a==1.0) {
     xx=xalpha;
   } else {
     xx=std::pow(xalpha, 1.0/a);
   }

   if (xx<0.0)  xx=0.0;
   if (xx>1.0)  xx=1.0;

   out(i)= xx;
 }

 return out;
}


// -----------------------------------------------------------------------------
// 4) rkkw: RNG for KwKw
// -----------------------------------------------------------------------------

//' @title Random Generation from KwKw Distribution
//'
//' @description
//' Generates \code{n} samples from the KwKw(\eqn{\alpha, \beta, \delta, \lambda}) distribution.
//'
//' @param n Integer number of samples.
//' @param alpha \eqn{\alpha > 0}.
//' @param beta \eqn{\beta > 0}.
//' @param delta \eqn{\delta \ge 0}.
//' @param lambda \eqn{\lambda > 0}.
//'
//' @details
//' The table suggests: if \eqn{V ~ Unif(0,1)}, define
//' \eqn{U = 1 - (1 - V)^{1/(\delta+1)}}, then
//' \eqn{X = \{1 - [1 - U^{1/\lambda}]^{1/\beta}\}^{1/\alpha}}.
//'
//' @export
// [[Rcpp::export]]
arma::vec rkkw(
   int n,
   const Rcpp::NumericVector& alpha,
   const Rcpp::NumericVector& beta,
   const Rcpp::NumericVector& delta,
   const Rcpp::NumericVector& lambda
) {
 if (n<=0) {
   Rcpp::stop("rkkw: n must be positive");
 }

 arma::vec a_vec(alpha.begin(), alpha.size());
 arma::vec b_vec(beta.begin(), beta.size());
 arma::vec d_vec(delta.begin(), delta.size());
 arma::vec l_vec(lambda.begin(), lambda.size());

 size_t k= std::max({a_vec.n_elem, b_vec.n_elem, d_vec.n_elem, l_vec.n_elem});
 arma::vec out(n);

 for (int i=0; i<n; i++) {
   size_t idx= i % k;
   double a= a_vec[idx % a_vec.n_elem];
   double b= b_vec[idx % b_vec.n_elem];
   double dd= d_vec[idx % d_vec.n_elem];
   double ll= l_vec[idx % l_vec.n_elem];

   if (!check_kkw_pars(a,b,dd,ll)) {
     out(i)= NA_REAL;
     Rcpp::warning("rkkw: invalid parameters at index %d", i+1);
     continue;
   }

   double V= R::runif(0.0,1.0);
   // U=1 - (1 - V)^(1/(δ+1))
   double U = 1.0 - std::pow(1.0 - V, 1.0/(dd+1.0));
   if (U<0.0)  U=0.0;
   if (U>1.0)  U=1.0;

   // x = {1 - [1 - U^(1/λ}]^(1/β)}^(1/α)
   double u_pow;
   if (ll==1.0) {
     u_pow=U;
   } else {
     u_pow=std::pow(U, 1.0/ll);
   }
   double bracket= 1.0- u_pow;
   if (bracket<0.0) bracket=0.0;
   if (bracket>1.0) bracket=1.0;
   double bracket2= std::pow(bracket, 1.0/b);
   double xalpha= 1.0 - bracket2;
   if (xalpha<0.0) xalpha=0.0;
   if (xalpha>1.0) xalpha=1.0;

   double xx;
   if (a==1.0) {
     xx=xalpha;
   } else {
     xx= std::pow(xalpha, 1.0/a);
     if (!R_finite(xx) || xx<0.0) xx=0.0;
     if (xx>1.0) xx=1.0;
   }
   out(i)=xx;
 }

 return out;
}


//' @title Negative Log-Likelihood of KwKw Distribution
//'
//' @description
//' Computes the negative log-likelihood for \eqn{\mathrm{KwKw}(\alpha,\beta,\delta,\lambda)}
//' given a data vector in (0,1).
//'
//' @param par NumericVector of length 4, \eqn{( \alpha,\beta,\delta,\lambda )}.
//' @param data NumericVector of observations, each in \eqn{(0,1)}.
//'
//' @details
//' The PDF is
//' \deqn{
//'   f(x) = \lambda\,\alpha\,\beta\,(\delta+1)\, x^{\alpha-1}\,(1-x^\alpha)^{\beta-1}\,
//'          [1 - (1 - x^\alpha)^\beta]^{\lambda-1}\,\{1 - [1 - (1-x^\alpha)^\beta]^\lambda\}^\delta.
//' }
//' The log-likelihood is the sum of \eqn{\log(f(x_i))}. This function returns the negative
//' log-likelihood (i.e. \eqn{-\ell(\theta)}).
//'
//' @return A single numeric value. \code{Inf} if invalid data or parameters.
//'
//' @examples
//' \dontrun{
//' # Generate sample data from a KKw distribution
//' set.seed(123)
//' x <- rkkw(100, 2, 3, 1.5, 0.5)
//' hist(x, breaks = 20, main = "KKw(2, 3, 1.5, 0.5) Sample")
//'
//' # Use in optimization with Hessian-based methods
//' result <- optim(c(0.5, 0.5, 0.5, 0.5), llkkw, method = "BFGS",
//'                 hessian = TRUE, data = x)
//'
//' # Compare numerical and analytical derivatives
//' num_grad <- numDeriv::grad(llkkw, x = result$par, data = x)
//' num_hess <- numDeriv::hessian(llkkw, x = result$par, data = x)
//'
//' ana_grad <- grkkw(result$par, data = x)
//' ana_hess <- hskkw(result$par, data = x)
//'
//' # Check differences (should be very small)
//' round(num_grad - ana_grad, 4)
//' round(num_hess - ana_hess, 4)
//'
//' }
//'
//' @export
// [[Rcpp::export]]
double llkkw(const Rcpp::NumericVector& par,
            const Rcpp::NumericVector& data) {
 if (par.size()<4) {
   return R_PosInf; // invalid
 }
 double a= par[0];
 double b= par[1];
 double dd= par[2];
 double ll= par[3];

 if (!check_kkw_pars(a,b,dd,ll)) {
   return R_PosInf; // param invalid => +Inf for negative log-likelihood
 }

 arma::vec x = Rcpp::as<arma::vec>(data);
 if (x.n_elem<1) {
   return R_PosInf;
 }
 // check data in (0,1)
 if (arma::any(x <= 0.0) || arma::any(x >= 1.0)) {
   return R_PosInf;
 }

 int n= x.n_elem;
 // sum of log-likelihood

 // constant part: n * [ log(λ)+ log(α)+ log(β)+ log(δ+1) ]
 double cst = n*( std::log(ll) + std::log(a) + std::log(b) + std::log(dd+1.0) );

 // sum parts
 // sum((α-1)*log(x))
 arma::vec lx = arma::log(x);
 double term1= (a-1.0)*arma::sum(lx);

 // sum((β-1)* log(1 - x^α))
 arma::vec xalpha = arma::pow(x,a);
 arma::vec log1mxalpha = arma::log(1.0 - xalpha);
 double term2= (b-1.0)*arma::sum(log1mxalpha);

 // sum((λ-1)* log(1 - (1 - x^α)^β))
 arma::vec vbeta = arma::pow(1.0 - xalpha, b);
 arma::vec one_minus_vbeta = 1.0 - vbeta;
 // safe log
 arma::vec log_omv = arma::log(one_minus_vbeta);
 double term3= (ll-1.0)*arma::sum(log_omv);

 // sum( δ* log(1 - [1-(1-x^α)^β]^λ ) )
 // let w= one_minus_vbeta => [1-w]^λ => so 1- w^λ is the factor
 arma::vec wlambda= arma::pow(one_minus_vbeta, ll);
 arma::vec one_minus_wlambda= 1.0 - wlambda;
 arma::vec log_1mw= arma::log(one_minus_wlambda);
 double term4= dd* arma::sum(log_1mw);

 double loglike = cst + term1 + term2 + term3 + term4;
 // negative
 return -loglike;
}


//' @title Gradient Function for Kumaraswamy-Kumaraswamy Log-Likelihood
//'
//' @description
//' Calculates the gradient vector (partial derivatives) of the negative log-likelihood
//' function for the Kumaraswamy-Kumaraswamy (KwKw) distribution. This function provides
//' the exact gradient needed for efficient optimization in maximum likelihood estimation.
//' The KwKw is a submodel of GKw with γ = 1 fixed.
//'
//' @param par NumericVector of length 4 containing parameters (α, β, δ, λ) in that order.
//'        All parameters must be positive.
//' @param data NumericVector of observations, where all values must be in the open interval (0,1).
//'
//' @return NumericVector of length 4 containing the gradient components (partial derivatives)
//'         of the negative log-likelihood with respect to each parameter (α, β, δ, λ).
//'         Returns a vector of NaN values if any parameters or data values are invalid.
//'
//' @details
//' The gradient vector contains the following partial derivatives of the negative log-likelihood:
//'
//' \deqn{
//' \frac{\partial \ell}{\partial \alpha} = \frac{n}{\alpha} + \sum_{i=1}^{n}\log(x_i) -
//' (\beta-1)\sum_{i=1}^{n}\frac{x_i^{\alpha}\log(x_i)}{1-x_i^{\alpha}} +
//' (\lambda-1)\sum_{i=1}^{n}\frac{\beta(1-x_i^{\alpha})^{\beta-1}x_i^{\alpha}\log(x_i)}{1-(1-x_i^{\alpha})^{\beta}} -
//' \delta\sum_{i=1}^{n}\frac{\lambda[1-(1-x_i^{\alpha})^{\beta}]^{\lambda-1}\beta(1-x_i^{\alpha})^{\beta-1}x_i^{\alpha}\log(x_i)}{1-[1-(1-x_i^{\alpha})^{\beta}]^{\lambda}}
//' }
//'
//' \deqn{
//' \frac{\partial \ell}{\partial \beta} = \frac{n}{\beta} + \sum_{i=1}^{n}\log(1-x_i^{\alpha}) -
//' (\lambda-1)\sum_{i=1}^{n}\frac{(1-x_i^{\alpha})^{\beta}\log(1-x_i^{\alpha})}{1-(1-x_i^{\alpha})^{\beta}} +
//' \delta\sum_{i=1}^{n}\frac{\lambda[1-(1-x_i^{\alpha})^{\beta}]^{\lambda-1}(1-x_i^{\alpha})^{\beta}\log(1-x_i^{\alpha})}{1-[1-(1-x_i^{\alpha})^{\beta}]^{\lambda}}
//' }
//'
//' \deqn{
//' \frac{\partial \ell}{\partial \delta} = \frac{n}{\delta+1} +
//' \sum_{i=1}^{n}\log(1-[1-(1-x_i^{\alpha})^{\beta}]^{\lambda})
//' }
//'
//' \deqn{
//' \frac{\partial \ell}{\partial \lambda} = \frac{n}{\lambda} +
//' \sum_{i=1}^{n}\log(1-(1-x_i^{\alpha})^{\beta}) -
//' \delta\sum_{i=1}^{n}\frac{[1-(1-x_i^{\alpha})^{\beta}]^{\lambda}\log(1-(1-x_i^{\alpha})^{\beta})}{1-[1-(1-x_i^{\alpha})^{\beta}]^{\lambda}}
//' }
//'
//' where:
//' \itemize{
//'   \item \deqn{v_i = 1 - x_i^{\alpha}}
//'   \item \deqn{w_i = 1 - v_i^{\beta} = 1 - (1-x_i^{\alpha})^{\beta}}
//'   \item \deqn{z_i = 1 - w_i^{\lambda} = 1 - (1-(1-x_i^{\alpha})^{\beta})^{\lambda}}
//' }
//'
//' The implementation includes several numerical safeguards:
//' \itemize{
//'   \item Parameter and data validation with appropriate error handling
//'   \item Clamping of intermediate values to avoid numerical underflow/overflow
//'   \item Efficient vector operations using Armadillo C++ library
//' }
//'
//' The returned gradient is negated to align with minimization of negative log-likelihood
//' in optimization routines.
//'
//' @examples
//' \dontrun{
//' # Generate sample data from a KKw distribution
//' set.seed(123)
//' x <- rkkw(100, 2, 3, 1.5, 0.5)
//' hist(x, breaks = 20, main = "KKw(2, 3, 1.5, 0.5) Sample")
//'
//' # Use in optimization with Hessian-based methods
//' result <- optim(c(0.5, 0.5, 0.5, 0.5), llkkw, method = "BFGS",
//'                 hessian = TRUE, data = x)
//'
//' # Compare numerical and analytical derivatives
//' num_grad <- numDeriv::grad(llkkw, x = result$par, data = x)
//' num_hess <- numDeriv::hessian(llkkw, x = result$par, data = x)
//'
//' ana_grad <- grkkw(result$par, data = x)
//' ana_hess <- hskkw(result$par, data = x)
//'
//' # Check differences (should be very small)
//' round(num_grad - ana_grad, 4)
//' round(num_hess - ana_hess, 4)
//' }
//'
//' @seealso
//' \code{\link[gkwreg]{llkkw}} for the negative log-likelihood function,
//' \code{\link[gkwreg]{hskkw}} for the Hessian matrix of the KwKw log-likelihood,
//' \code{\link[gkwreg]{dkkw}} for the KwKw density function,
//'
//' @references
//' Kumaraswamy, P. (1980). A generalized probability density function for double-bounded random processes.
//' Journal of Hydrology, 46(1-2), 79-88.
//'
//' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized distributions.
//' Journal of Statistical Computation and Simulation, 81(7), 883-898.
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector grkkw(const Rcpp::NumericVector& par, const Rcpp::NumericVector& data) {
 // Parameter extraction
 double alpha = par[0];   // Shape parameter α > 0
 double beta = par[1];    // Shape parameter β > 0
 double delta = par[2];   // Shape parameter δ > 0
 double lambda = par[3];  // Shape parameter λ > 0

 // Parameter validation
 if (alpha <= 0 || beta <= 0 || delta <= 0 || lambda <= 0) {
   Rcpp::NumericVector grad(4, R_NaN);
   return grad;
 }

 // Data conversion and validation
 arma::vec x = Rcpp::as<arma::vec>(data);

 if (arma::any(x <= 0) || arma::any(x >= 1)) {
   Rcpp::NumericVector grad(4, R_NaN);
   return grad;
 }

 int n = x.n_elem;  // Sample size

 // Initialize gradient vector
 Rcpp::NumericVector grad(4, 0.0);

 // Small constant to avoid numerical issues
 double eps = std::numeric_limits<double>::epsilon() * 100;

 // Compute transformations and intermediate values
 arma::vec log_x = arma::log(x);                // log(x_i)
 arma::vec x_alpha = arma::pow(x, alpha);       // x_i^α
 arma::vec x_alpha_log_x = x_alpha % log_x;     // x_i^α * log(x_i)

 // v_i = 1 - x_i^α
 arma::vec v = 1.0 - x_alpha;
 v = arma::clamp(v, eps, 1.0 - eps);            // Prevent numerical issues

 arma::vec log_v = arma::log(v);                // log(v_i)
 arma::vec v_beta_m1 = arma::pow(v, beta - 1.0); // v_i^(β-1)
 arma::vec v_beta = arma::pow(v, beta);          // v_i^β
 arma::vec v_beta_log_v = v_beta % log_v;        // v_i^β * log(v_i)

 // w_i = 1 - v_i^β = 1 - (1-x_i^α)^β
 arma::vec w = 1.0 - v_beta;
 w = arma::clamp(w, eps, 1.0 - eps);            // Prevent numerical issues

 arma::vec log_w = arma::log(w);                // log(w_i)
 arma::vec w_lambda_m1 = arma::pow(w, lambda - 1.0); // w_i^(λ-1)
 arma::vec w_lambda = arma::pow(w, lambda);          // w_i^λ
 arma::vec w_lambda_log_w = w_lambda % log_w;        // w_i^λ * log(w_i)

 // z_i = 1 - w_i^λ = 1 - [1-(1-x_i^α)^β]^λ
 arma::vec z = 1.0 - w_lambda;
 z = arma::clamp(z, eps, 1.0 - eps);            // Prevent numerical issues

 arma::vec log_z = arma::log(z);                // log(z_i)

 // Calculate partial derivatives for each parameter (for log-likelihood)

 // ∂ℓ/∂α = n/α + Σᵢlog(xᵢ) - (β-1)Σᵢ[xᵢ^α*log(xᵢ)/vᵢ] + (λ-1)Σᵢ[β*vᵢ^(β-1)*xᵢ^α*log(xᵢ)/wᵢ] - δΣᵢ[λ*wᵢ^(λ-1)*β*vᵢ^(β-1)*xᵢ^α*log(xᵢ)/zᵢ]
 double d_alpha = n / alpha + arma::sum(log_x);

 // Calculate terms separately with correct signs
 arma::vec alpha_term1 = (beta - 1.0) / v;                       // (β-1)/v_i
 arma::vec alpha_term2 = (lambda - 1.0) * beta * v_beta_m1 / w;  // (λ-1)*β*v_i^(β-1)/w_i
 arma::vec alpha_term3 = delta * lambda * beta * v_beta_m1 % w_lambda_m1 / z;  // δ*λ*β*v_i^(β-1)*w_i^(λ-1)/z_i

 // Combine terms with correct signs
 d_alpha -= arma::sum(x_alpha_log_x % alpha_term1);  // Subtract term with (β-1)
 d_alpha += arma::sum(x_alpha_log_x % alpha_term2);  // Add term with (λ-1)
 d_alpha -= arma::sum(x_alpha_log_x % alpha_term3);  // Subtract term with δ

 // ∂ℓ/∂β = n/β + Σᵢlog(vᵢ) - (λ-1)Σᵢ[vᵢ^β*log(vᵢ)/wᵢ] + δΣᵢ[λ*wᵢ^(λ-1)*vᵢ^β*log(vᵢ)/zᵢ]
 double d_beta = n / beta + arma::sum(log_v);

 // Calculate terms separately with correct signs
 arma::vec beta_term1 = (lambda - 1.0) / w;       // (λ-1)/w_i
 arma::vec beta_term2 = delta * lambda * w_lambda_m1 / z; // δ*λ*w_i^(λ-1)/z_i

 // CORRECTED: These signs were wrong in the previous implementation
 d_beta -= arma::sum(v_beta_log_v % beta_term1);  // Subtract term with (λ-1)
 d_beta += arma::sum(v_beta_log_v % beta_term2);  // Add term with δ

 // ∂ℓ/∂δ = n/(δ+1) + Σᵢlog(zᵢ)
 double d_delta = n / (delta + 1.0) + arma::sum(log_z);

 // ∂ℓ/∂λ = n/λ + Σᵢlog(wᵢ) - δΣᵢ[(wᵢ^λ*log(wᵢ))/zᵢ]
 double d_lambda = n / lambda + arma::sum(log_w) - delta * arma::sum(w_lambda_log_w / z);

 // Since we're optimizing negative log-likelihood, negate all derivatives
 grad[0] = -d_alpha;
 grad[1] = -d_beta;
 grad[2] = -d_delta;
 grad[3] = -d_lambda;

 return grad;
}



//' @title Analytic Hessian Matrix for Kumaraswamy-Kumaraswamy Distribution
//'
//' @description
//' Computes the analytic Hessian matrix of the log-likelihood function for
//' the Kumaraswamy-Kumaraswamy (KKw) distribution. This function provides
//' exact second derivatives needed for optimization and inference.
//'
//' @param par Numeric vector of length 4 containing the parameters
//'        (α, β, δ, λ) in that order. All parameters must be positive.
//' @param data Numeric vector of observations, where all values must be
//'        in the open interval (0,1).
//'
//' @return A 4×4 numeric matrix representing the Hessian of the negative
//'         log-likelihood function. If parameters or data are invalid
//'         (parameters ≤ 0 or data outside (0,1)), returns a matrix of
//'         NaN values.
//'
//' @details
//' The log-likelihood for the Kumaraswamy-Kumaraswamy distribution is:
//'
//' \deqn{
//' \ell(\theta) = n \ln(\lambda) + n \ln(\alpha) + n \ln(\beta) + n \ln(\delta+1)
//' + (\alpha-1) \sum \ln(x_i)
//' + (\beta-1) \sum \ln(1 - x_i^\alpha)
//' + (\lambda-1) \sum \ln\{1 - (1 - x_i^\alpha)^\beta\}
//' + \delta \sum \ln\{1 - \{1 - (1 - x_i^\alpha)^\beta\}^\lambda\}
//' }
//'
//' where γ is fixed at 1 for this distribution.
//'
//' The implementation computes all second derivatives analytically for each term.
//' For computational efficiency, the following transformations are used:
//' \itemize{
//'   \item \deqn{A = x^α} and derivatives
//'   \item \deqn{v = 1 - A}
//'   \item \deqn{w = 1 - v^β}
//'   \item \deqn{z = 1 - w^λ}
//' }
//'
//' The returned Hessian matrix has the following structure:
//' \itemize{
//'   \item Rows/columns 1-4 correspond to α, β, δ, λ respectively
//'   \item The matrix is symmetric (as expected for a Hessian)
//'   \item The matrix represents second derivatives of the negative log-likelihood
//' }
//'
//' This function is implemented in C++ for computational efficiency.
//'
//' @examples
//' \dontrun{
//' # Generate sample data from a KKw distribution
//' set.seed(123)
//' x <- rkkw(100, 2, 3, 1.5, 0.5)
//' hist(x, breaks = 20, main = "KKw(2, 3, 1.5, 0.5) Sample")
//'
//' # Use in optimization with Hessian-based methods
//' result <- optim(c(0.5, 0.5, 0.5, 0.5), llkkw, method = "BFGS",
//'                 hessian = TRUE, data = x)
//'
//' # Compare numerical and analytical derivatives
//' num_grad <- numDeriv::grad(llkkw, x = result$par, data = x)
//' num_hess <- numDeriv::hessian(llkkw, x = result$par, data = x)
//'
//' ana_grad <- grkkw(result$par, data = x)
//' ana_hess <- hskkw(result$par, data = x)
//'
//' # Check differences (should be very small)
//' round(num_grad - ana_grad, 4)
//' round(num_hess - ana_hess, 4)
//' }
//'
//'
//' @references
//' Kumaraswamy, P. (1980). A generalized probability density function for double-bounded random processes.
//' Journal of Hydrology, 46(1-2), 79-88.
//'
//' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized distributions.
//' Journal of Statistical Computation and Simulation, 81(7), 883-898.
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix hskkw(const Rcpp::NumericVector& par, const Rcpp::NumericVector& data) {
 // Parameter extraction
 double alpha  = par[0];   // θ[0] = α
 double beta   = par[1];   // θ[1] = β
 double delta  = par[2];   // θ[2] = δ
 double lambda = par[3];   // θ[3] = λ

 // Simple parameter validation (all > 0)
 if(alpha <= 0 || beta <= 0 || delta <= 0 || lambda <= 0) {
   Rcpp::NumericMatrix nanH(4,4);
   nanH.fill(R_NaN);
   return nanH;
 }

 // Data conversion and basic validation
 arma::vec x = Rcpp::as<arma::vec>(data);
 if(arma::any(x <= 0) || arma::any(x >= 1)) {
   Rcpp::NumericMatrix nanH(4,4);
   nanH.fill(R_NaN);
   return nanH;
 }

 int n = x.n_elem;  // sample size

 // Initialize Hessian matrix H (of ℓ(θ)) as 4x4
 arma::mat H(4,4, arma::fill::zeros);

 // --- CONSTANT TERMS (do not depend on x) ---
 // L1: n ln(λ)  => d²/dλ² = -n/λ²
 H(3,3) += -n/(lambda*lambda);
 // L2: n ln(α)  => d²/dα² = -n/α²
 H(0,0) += -n/(alpha*alpha);
 // L3: n ln(β)  => d²/dβ² = -n/β²
 H(1,1) += -n/(beta*beta);
 // L4: n ln(δ+1) => d²/dδ² = -n/(δ+1)²
 H(2,2) += -n/std::pow(delta+1.0, 2.0);

 // Accumulators for mixed derivatives with λ
 double acc_delta_lambda = 0.0;  // Sum of dz_dlambda / z
 double acc_alpha_lambda = 0.0;  // For α,λ contributions
 double acc_beta_lambda = 0.0;   // For β,λ contributions

 // --- TERMS THAT INVOLVE THE OBSERVATIONS ---
 // Loop over each observation to accumulate contributions from:
 // L5: (α-1) Σ ln(x_i)  --> contributes only to first derivatives
 // L6: (β-1) Σ ln(v), where v = 1 - x^α
 // L7: (λ-1) Σ ln(w), where w = 1 - v^β
 // L8: δ Σ ln(z), where z = 1 - w^λ
 for (int i = 0; i < n; i++) {
   double xi    = x(i);
   double ln_xi = std::log(xi);

   // -- Compute A = x^α and its derivatives --
   double A = std::pow(xi, alpha);                  // A = x^α
   double dA_dalpha = A * ln_xi;                    // dA/dα = x^α ln(x)
   double d2A_dalpha2 = A * ln_xi * ln_xi;          // d²A/dα² = x^α (ln(x))²

   // -- v = 1 - A and its derivatives --
   double v = 1.0 - A;                              // v = 1 - x^α
   double ln_v = std::log(v);                       // ln(v)
   double dv_dalpha = -dA_dalpha;                   // dv/dα = -dA/dα = -x^α ln(x)
   double d2v_dalpha2 = -d2A_dalpha2;               // d²v/dα² = -d²A/dα² = -x^α (ln(x))²

   // --- L6: (β-1) ln(v) ---
   // Second derivative w.r.t. α: (β-1)*[(d²v/dα²*v - (dv/dα)²)/v²]
   double d2L6_dalpha2 = (beta - 1.0) * ((d2v_dalpha2 * v - dv_dalpha * dv_dalpha) / (v*v));
   // Mixed derivative: d²L6/(dα dβ) = d/dβ[(β-1)*(dv_dalpha/v)] = (dv_dalpha/v)
   double d2L6_dalpha_dbeta = dv_dalpha / v;

   // --- L7: (λ - 1) ln(w), where w = 1 - v^β ---
   double v_beta = std::pow(v, beta);               // v^β
   double w = 1.0 - v_beta;                         // w = 1 - v^β
   double ln_w = std::log(w);                       // ln(w)
   // Derivative of w w.r.t. v: dw/dv = -β * v^(β-1)
   double dw_dv = -beta * std::pow(v, beta - 1.0);
   // Chain rule: dw/dα = dw/dv * dv/dα
   double dw_dalpha = dw_dv * dv_dalpha;
   // Second derivative w.r.t. α for L7:
   // d²/dα² ln(w) = [d²w/dα² * w - (dw/dα)²] / w²
   // Computing d²w/dα²:
   //   dw/dα = -β * v^(β-1)*dv_dalpha,
   //   d²w/dα² = -β * [(β-1)*v^(β-2)*(dv_dalpha)² + v^(β-1)*d²v_dalpha²]
   double d2w_dalpha2 = -beta * ((beta - 1.0) * std::pow(v, beta-2.0) * (dv_dalpha * dv_dalpha)
                                   + std::pow(v, beta-1.0) * d2v_dalpha2);
   double d2L7_dalpha2 = (lambda - 1.0) * ((d2w_dalpha2 * w - (dw_dalpha * dw_dalpha)) / (w*w));
   // Derivative w.r.t. β: d/dβ ln(w). Note: d/dβ(v^β) = v^β ln(v) => d/dβ w = -v^β ln(v)
   double dw_dbeta = -v_beta * ln_v;
   // Second derivative w.r.t. β for L7:
   // d²/dβ² ln(w) = [d²w/dβ² * w - (dw/dβ)²]/w², where d²w/dβ² = -v^β (ln(v))²
   double d2w_dbeta2 = -v_beta * (ln_v * ln_v);
   double d2L7_dbeta2 = (lambda - 1.0) * ((d2w_dbeta2 * w - (dw_dbeta * dw_dbeta))/(w*w));
   // Mixed derivative L7 (α,β): d²/(dα dβ) ln(w) =
   //   = d/dβ[(dw_dalpha)/w] = (d/dβ dw_dalpha)/w - (dw_dalpha*dw_dbeta)/(w*w)
   double d_dw_dalpha_dbeta = -std::pow(v, beta-1.0) * (1.0 + beta * ln_v) * dv_dalpha;
   double d2L7_dalpha_dbeta = (lambda - 1.0) * ((d_dw_dalpha_dbeta / w) - (dw_dalpha * dw_dbeta)/(w*w));

   // --- L8: δ ln(z), where z = 1 - w^λ ---
   double w_lambda_val = std::pow(w, lambda);       // w^λ
   double z = 1.0 - w_lambda_val;                   // z = 1 - w^λ
   // Derivative w.r.t. α: dz/dα = -λ * w^(λ-1) * dw/dα
   double dz_dalpha = -lambda * std::pow(w, lambda-1.0) * dw_dalpha;
   // Second derivative w.r.t. α for L8:
   // d²z/dα² = -λ * [(λ-1)*w^(λ-2)*(dw/dα)² + w^(λ-1)*d²w/dα²]
   double d2z_dalpha2 = -lambda * ((lambda - 1.0) * std::pow(w, lambda-2.0) * (dw_dalpha*dw_dalpha)
                                     + std::pow(w, lambda-1.0) * d2w_dalpha2);
   double d2L8_dalpha2 = delta * ((d2z_dalpha2 * z - dz_dalpha*dz_dalpha)/(z*z));

   // Derivative w.r.t. β: dz/dβ = -λ * w^(λ-1) * dw/dβ
   double dz_dbeta = -lambda * std::pow(w, lambda-1.0) * dw_dbeta;
   // Second derivative w.r.t. β for L8:
   // d²z/dβ² = -λ * [(λ-1)*w^(λ-2)*(dw/dβ)² + w^(λ-1)*d²w/dβ²]
   double d2z_dbeta2 = -lambda * ((lambda - 1.0) * std::pow(w, lambda-2.0) * (dw_dbeta*dw_dbeta)
                                    + std::pow(w, lambda-1.0) * d2w_dbeta2);
   double d2L8_dbeta2 = delta * ((d2z_dbeta2 * z - dz_dbeta*dz_dbeta)/(z*z));

   // Mixed derivative L8 (α,β): d²/(dα dβ) ln(z)
   // = (d/dβ dz_dalpha)/z - (dz_dalpha*dz_dbeta)/(z*z)
   double d_dw_dalpha_dbeta_2 = -lambda * ((lambda - 1.0) * std::pow(w, lambda-2.0) * dw_dbeta * dw_dalpha
                                             + std::pow(w, lambda-1.0) * d_dw_dalpha_dbeta);
   double d2L8_dalpha_dbeta = delta * ((d_dw_dalpha_dbeta_2 / z) - (dz_dalpha*dz_dbeta)/(z*z));

   // Derivatives of L8 with respect to λ:
   // d/dλ ln(z) = (1/z)*dz/dλ, with dz/dλ = -w^λ ln(w)
   double dz_dlambda = -w_lambda_val * ln_w;
   // d²/dλ² ln(z) = [d²z/dλ² * z - (dz_dlambda)²]/z² (assuming w constant in λ)
   double d2z_dlambda2 = -w_lambda_val * (ln_w * ln_w);
   double d2L8_dlambda2 = delta * ((d2z_dlambda2 * z - dz_dlambda*dz_dlambda)/(z*z));

   // Mixed derivative L8 (α,λ): d²/(dα dλ) ln(z) = (d/dα dz_dλ)/z - (dz_dλ*dz_dalpha)/(z*z)
   double d_dalpha_dz_dlambda = -lambda * std::pow(w, lambda-1.0) * dw_dalpha * ln_w -
     w_lambda_val * (dw_dalpha / w);
   double d2L8_dalpha_dlambda = delta * ((d_dalpha_dz_dlambda / z) - (dz_dlambda*dz_dalpha)/(z*z));

   // Mixed derivative L8 (β,λ): d²/(dβ dλ) ln(z) = (d/dβ dz_dλ)/z - (dz_dlambda*dz_dbeta)/(z*z)
   double d_dbeta_dz_dlambda = -lambda * std::pow(w, lambda-1.0) * dw_dbeta * ln_w -
     w_lambda_val * (dw_dbeta / w);
   double d2L8_dbeta_dlambda = delta * ((d_dbeta_dz_dlambda / z) - (dz_dlambda*dz_dbeta)/(z*z));

   // Mixed derivative L8 (δ,λ): d²/(dδ dλ) ln(z) = d/dλ [δ * (dz_dlambda / z)] = dz_dlambda / z
   double d2L8_ddelta_dlambda = dz_dlambda / z;

   // --- ACCUMULATING CONTRIBUTIONS TO THE HESSIAN MATRIX ---
   // Index: 0 = α, 1 = β, 2 = δ, 3 = λ

   // H(α,α): sum of L2, L6, L7, and L8 (constants already added)
   H(0,0) += d2L6_dalpha2 + d2L7_dalpha2 + d2L8_dalpha2;

   // H(α,β): mixed from L6, L7, and L8
   H(0,1) += d2L6_dalpha_dbeta + d2L7_dalpha_dbeta + d2L8_dalpha_dbeta;
   H(1,0) = H(0,1);

   // H(β,β): contributions from L3, L7, and L8
   H(1,1) += d2L7_dbeta2 + d2L8_dbeta2;

   // H(λ,λ): contributions from L1 and L8 (L1 already added)
   H(3,3) += d2L8_dlambda2;

   // H(α,δ): From L8 - d/dα ln(z)
   H(0,2) += dz_dalpha / z;
   H(2,0) = H(0,2);

   // H(β,δ): From L8 - d/dβ ln(z)
   H(1,2) += dz_dbeta / z;
   H(2,1) = H(1,2);

   // Accumulating terms for mixed derivatives with λ
   // The term from L7 is Σ ln(w) derived w.r.t. α and Σ ln(w) derived w.r.t. β
   acc_alpha_lambda += (dw_dalpha / w) + d2L8_dalpha_dlambda;
   acc_beta_lambda += (dw_dbeta / w) + d2L8_dbeta_dlambda;
   acc_delta_lambda += d2L8_ddelta_dlambda;

 } // end of loop

 // Applying mixed derivatives with λ
 // H(α,λ):
 H(0,3) = acc_alpha_lambda;
 H(3,0) = H(0,3);

 // H(β,λ):
 H(1,3) = acc_beta_lambda;
 H(3,1) = H(1,3);

 // H(δ,λ):
 H(2,3) = acc_delta_lambda;
 H(3,2) = H(2,3);

 // Returns the analytic Hessian matrix of the negative log-likelihood
 return Rcpp::wrap(-H);
}

// //' @title Analytic Hessian Matrix for Kumaraswamy-Kumaraswamy Distribution
// //'
// //' @description
// //' Computes the analytic Hessian matrix of the log-likelihood function for
// //' the Kumaraswamy-Kumaraswamy (KKw) distribution. This function provides
// //' exact second derivatives needed for optimization and inference.
// //'
// //' @param par Numeric vector of length 4 containing the parameters
// //'        (α, β, δ, λ) in that order. All parameters must be positive.
// //' @param data Numeric vector of observations, where all values must be
// //'        in the open interval (0,1).
// //'
// //' @return A 4×4 numeric matrix representing the Hessian of the negative
// //'         log-likelihood function. If parameters or data are invalid
// //'         (parameters ≤ 0 or data outside (0,1)), returns a matrix of
// //'         NaN values.
// //'
// //' @details
// //' The log-likelihood for the Kumaraswamy-Kumaraswamy distribution is:
// //'
// //' \deqn{
// //' \ell(\theta) = n \ln(\lambda) + n \ln(\alpha) + n \ln(\beta) + n \ln(\delta+1)
// //' + (\alpha-1) \sum \ln(x_i)
// //' + (\beta-1) \sum \ln(1 - x_i^\alpha)
// //' + (\lambda-1) \sum \ln\{1 - (1 - x_i^\alpha)^\beta\}
// //' + \delta \sum \ln\{1 - \{1 - (1 - x_i^\alpha)^\beta\}^\lambda\}
// //' }
// //'
// //' where γ is fixed at 1 for this distribution.
// //'
// //' The implementation computes all second derivatives analytically for each term.
// //' For computational efficiency, the following transformations are used:
// //' \itemize{
// //'   \item \deqn{A = x^α} and derivatives
// //'   \item \deqn{v = 1 - A}
// //'   \item \deqn{w = 1 - v^β}
// //'   \item \deqn{z = 1 - w^λ}
// //' }
// //'
// //' The returned Hessian matrix has the following structure:
// //' \itemize{
// //'   \item Rows/columns 1-4 correspond to α, β, δ, λ respectively
// //'   \item The matrix is symmetric (as expected for a Hessian)
// //'   \item The matrix represents second derivatives of the negative log-likelihood
// //' }
// //'
// //' This function is implemented in C++ for computational efficiency.
// //'
// //' @examples
// //' \dontrun{
// //' # Generate sample data from a KKw distribution
// //' set.seed(123)
// //' x <- rkkw(100, 2, 3, 1.5, 0.5)
// //' hist(x, breaks = 20, main = "KKw(2, 3, 1.5, 0.5) Sample")
// //'
// //' # Use in optimization with Hessian-based methods
// //' result <- optim(c(0.5, 0.5, 0.5, 0.5), llkkw, method = "BFGS",
// //'                 hessian = TRUE, data = x)
// //'
// //' # Compare numerical and analytical derivatives
// //' num_grad <- numDeriv::grad(llkkw, x = result$par, data = x)
// //' num_hess <- numDeriv::hessian(llkkw, x = result$par, data = x)
// //'
// //' ana_grad <- grkkw(result$par, data = x)
// //' ana_hess <- hskkw(result$par, data = x)
// //'
// //' # Check differences (should be very small)
// //' round(num_grad - ana_grad, 4)
// //' round(num_hess - ana_hess, 4)
// //' }
// //'
// //'
// //' @references
// //' Kumaraswamy, P. (1980). A generalized probability density function for double-bounded random processes.
// //' Journal of Hydrology, 46(1-2), 79-88.
// //'
// //' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized distributions.
// //' Journal of Statistical Computation and Simulation, 81(7), 883-898.
// //'
// //' @export
// // [[Rcpp::export]]
// Rcpp::NumericMatrix hskkw(const Rcpp::NumericVector& par, const Rcpp::NumericVector& data) {
//  // Parameter extraction
//  double alpha  = par[0];   // θ[0] = α
//  double beta   = par[1];   // θ[1] = β
//  double delta  = par[2];   // θ[2] = δ
//  double lambda = par[3];   // θ[3] = λ
//
//  // Simple parameter validation (all > 0)
//  if(alpha <= 0 || beta <= 0 || delta <= 0 || lambda <= 0) {
//    Rcpp::NumericMatrix nanH(4,4);
//    nanH.fill(R_NaN);
//    return nanH;
//  }
//
//  // Data conversion and basic validation
//  arma::vec x = Rcpp::as<arma::vec>(data);
//  if(arma::any(x <= 0) || arma::any(x >= 1)) {
//    Rcpp::NumericMatrix nanH(4,4);
//    nanH.fill(R_NaN);
//    return nanH;
//  }
//
//  int n = x.n_elem;  // sample size
//
//  // Initialize Hessian matrix H (of ℓ(θ)) as 4x4
//  arma::mat H(4,4, arma::fill::zeros);
//
//  // --- CONSTANT TERMS (do not depend on x) ---
//  // L1: n ln(λ)  => d²/dλ² = -n/λ²
//  H(3,3) += -n/(lambda*lambda);
//  // L2: n ln(α)  => d²/dα² = -n/α²
//  H(0,0) += -n/(alpha*alpha);
//  // L3: n ln(β)  => d²/dβ² = -n/β²
//  H(1,1) += -n/(beta*beta);
//  // L4: n ln(δ+1) => d²/dδ² = -n/(δ+1)²
//  H(2,2) += -n/std::pow(delta+1.0, 2.0);
//
//  // Accumulators for mixed derivatives with λ
//  double acc_delta_lambda = 0.0;  // Sum of dz_dlambda / z
//  double acc_alpha_lambda = 0.0;  // For α,λ contributions
//  double acc_beta_lambda = 0.0;   // For β,λ contributions
//
//  // --- TERMS THAT INVOLVE THE OBSERVATIONS ---
//  // Loop over each observation to accumulate contributions from:
//  // L5: (α-1) Σ ln(x_i)  --> contributes only to first derivatives
//  // L6: (β-1) Σ ln(v), where v = 1 - x^α
//  // L7: (λ-1) Σ ln(w), where w = 1 - v^β
//  // L8: δ Σ ln(z), where z = 1 - w^λ
//  for (int i = 0; i < n; i++) {
//    double xi    = x(i);
//    double ln_xi = std::log(xi);
//
//    // -- Compute A = x^α and its derivatives --
//    double A = std::pow(xi, alpha);                  // A = x^α
//    double dA_dalpha = A * ln_xi;                    // dA/dα = x^α ln(x)
//    double d2A_dalpha2 = A * ln_xi * ln_xi;          // d²A/dα² = x^α (ln(x))²
//
//    // -- v = 1 - A and its derivatives --
//    double v = 1.0 - A;                              // v = 1 - x^α
//    double ln_v = std::log(v);                       // ln(v)
//    double dv_dalpha = -dA_dalpha;                   // dv/dα = -dA/dα = -x^α ln(x)
//    double d2v_dalpha2 = -d2A_dalpha2;               // d²v/dα² = -d²A/dα² = -x^α (ln(x))²
//
//    // --- L6: (β-1) ln(v) ---
//    // First derivative w.r.t. α: (β-1) * (1/v)*dv_dalpha
//    double dL6_dalpha = (beta - 1.0) * (dv_dalpha / v);
//    // Second derivative w.r.t. α: (β-1)*[(d²v/dα²*v - (dv/dα)²)/v²]
//    double d2L6_dalpha2 = (beta - 1.0) * ((d2v_dalpha2 * v - dv_dalpha * dv_dalpha) / (v*v));
//    // First derivative w.r.t. β: dL6/dβ = ln(v)
//    double dL6_dbeta = ln_v;
//    // Mixed derivative: d²L6/(dα dβ) = d/dβ[(β-1)*(dv_dalpha/v)] = (dv_dalpha/v)
//    double d2L6_dalpha_dbeta = dv_dalpha / v;
//
//    // --- L7: (λ - 1) ln(w), where w = 1 - v^β ---
//    double v_beta = std::pow(v, beta);               // v^β
//    double w = 1.0 - v_beta;                         // w = 1 - v^β
//    double ln_w = std::log(w);                       // ln(w)
//    // Derivative of w w.r.t. v: dw/dv = -β * v^(β-1)
//    double dw_dv = -beta * std::pow(v, beta - 1.0);
//    // Chain rule: dw/dα = dw/dv * dv/dα
//    double dw_dalpha = dw_dv * dv_dalpha;
//    // First derivative w.r.t. α: d/dα ln(w) = (1/w)*dw_dalpha
//    double dL7_dalpha = (lambda - 1.0) * (dw_dalpha / w);
//    // Second derivative w.r.t. α for L7:
//    // d²/dα² ln(w) = [d²w/dα² * w - (dw/dα)²] / w²
//    // Computing d²w/dα²:
//    //   dw/dα = -β * v^(β-1)*dv_dalpha,
//    //   d²w/dα² = -β * [(β-1)*v^(β-2)*(dv_dalpha)² + v^(β-1)*d²v_dalpha²]
//    double d2w_dalpha2 = -beta * ((beta - 1.0) * std::pow(v, beta-2.0) * (dv_dalpha * dv_dalpha)
//                                    + std::pow(v, beta-1.0) * d2v_dalpha2);
//    double d2L7_dalpha2 = (lambda - 1.0) * ((d2w_dalpha2 * w - (dw_dalpha * dw_dalpha)) / (w*w));
//    // Derivative w.r.t. β: d/dβ ln(w). Note: d/dβ(v^β) = v^β ln(v) => d/dβ w = -v^β ln(v)
//    double dw_dbeta = -v_beta * ln_v;
//    double dL7_dbeta = (lambda - 1.0) * (dw_dbeta / w);
//    // Second derivative w.r.t. β for L7:
//    // d²/dβ² ln(w) = [d²w/dβ² * w - (dw/dβ)²]/w², where d²w/dβ² = -v^β (ln(v))²
//    double d2w_dbeta2 = -v_beta * (ln_v * ln_v);
//    double d2L7_dbeta2 = (lambda - 1.0) * ((d2w_dbeta2 * w - (dw_dbeta * dw_dbeta))/(w*w));
//    // Mixed derivative L7 (α,β): d²/(dα dβ) ln(w) =
//    //   = d/dβ[(dw_dalpha)/w] = (d/dβ dw_dalpha)/w - (dw_dalpha*dw_dbeta)/(w*w)
//    // Approximate d/dβ dw_dalpha:
//    double d_dw_dalpha_dbeta = -std::pow(v, beta-1.0) * (1.0 + beta * ln_v) * dv_dalpha;
//    double d2L7_dalpha_dbeta = (lambda - 1.0) * ((d_dw_dalpha_dbeta / w) - (dw_dalpha * dw_dbeta)/(w*w));
//
//    // --- L8: δ ln(z), where z = 1 - w^λ ---
//    double w_lambda_val = std::pow(w, lambda);       // w^λ
//    double z = 1.0 - w_lambda_val;                   // z = 1 - w^λ
//    double ln_z = std::log(z);                       // ln(z)
//    // Derivative w.r.t. α: dz/dα = -λ * w^(λ-1) * dw/dα
//    double dz_dalpha = -lambda * std::pow(w, lambda-1.0) * dw_dalpha;
//    double dL8_dalpha = delta * (dz_dalpha / z);
//    // Second derivative w.r.t. α for L8:
//    // d²z/dα² = -λ * [(λ-1)*w^(λ-2)*(dw/dα)² + w^(λ-1)*d²w/dα²]
//    double d2w_dalpha2_dummy = d2w_dalpha2; // already calculated for L7
//    double d2z_dalpha2 = -lambda * ((lambda - 1.0) * std::pow(w, lambda-2.0) * (dw_dalpha*dw_dalpha)
//                                      + std::pow(w, lambda-1.0) * d2w_dalpha2_dummy);
//    double d2L8_dalpha2 = delta * ((d2z_dalpha2 * z - dz_dalpha*dz_dalpha)/(z*z));
//
//    // Derivative w.r.t. β: dz/dβ = -λ * w^(λ-1) * dw/dβ
//    double dz_dbeta = -lambda * std::pow(w, lambda-1.0) * dw_dbeta;
//    double dL8_dbeta = delta * (dz_dbeta / z);
//    // Second derivative w.r.t. β for L8:
//    // d²z/dβ² = -λ * [(λ-1)*w^(λ-2)*(dw/dβ)² + w^(λ-1)*d²w/dβ²]
//    double d2z_dbeta2 = -lambda * ((lambda - 1.0) * std::pow(w, lambda-2.0) * (dw_dbeta*dw_dbeta)
//                                     + std::pow(w, lambda-1.0) * d2w_dbeta2);
//    double d2L8_dbeta2 = delta * ((d2z_dbeta2 * z - dz_dbeta*dz_dbeta)/(z*z));
//
//    // Mixed derivative L8 (α,β): d²/(dα dβ) ln(z)
//    // = (d/dβ dz_dalpha)/z - (dz_dalpha*dz_dbeta)/(z*z)
//    // Approximate d/dβ dz_dalpha = -λ * [(λ-1)*w^(λ-2)*(dw_dβ*dw_dα) + w^(λ-1)*(d/dβ dw_dalpha)]
//    double d_dw_dalpha_dbeta_2 = -lambda * ((lambda - 1.0) * std::pow(w, lambda-2.0) * dw_dbeta * dw_dalpha
//                                              + std::pow(w, lambda-1.0) * d_dw_dalpha_dbeta);
//    double d2L8_dalpha_dbeta = delta * ((d_dw_dalpha_dbeta_2 / z) - (dz_dalpha*dz_dbeta)/(z*z));
//
//    // Derivatives of L8 with respect to λ:
//    // d/dλ ln(z) = (1/z)*dz/dλ, with dz/dλ = -w^λ ln(w)
//    double dz_dlambda = -w_lambda_val * ln_w;
//    double dL8_dlambda = delta * (dz_dlambda / z);
//    // d²/dλ² ln(z) = [d²z/dλ² * z - (dz_dlambda)²]/z² (assuming w constant in λ)
//    double d2z_dlambda2 = -w_lambda_val * (ln_w * ln_w);
//    double d2L8_dlambda2 = delta * ((d2z_dlambda2 * z - dz_dlambda*dz_dlambda)/(z*z));
//
//    // Mixed derivative L8 (α,λ): d²/(dα dλ) ln(z) = (d/dα dz_dλ)/z - (dz_dλ*dz_dalpha)/(z*z)
//    // The correct formula for d/dα[dz_dλ] = d/dα[-w^λ*ln(w)]
//    // = -λ*w^(λ-1)*dw/dα*ln(w) - w^λ*(1/w)*dw/dα
//    double d_dalpha_dz_dlambda = -lambda * std::pow(w, lambda-1.0) * dw_dalpha * ln_w -
//      w_lambda_val * (dw_dalpha / w);
//    double d2L8_dalpha_dlambda = delta * ((d_dalpha_dz_dlambda / z) - (dz_dlambda*dz_dalpha)/(z*z));
//
//    // Mixed derivative L8 (β,λ): d²/(dβ dλ) ln(z) = (d/dβ dz_dλ)/z - (dz_dlambda*dz_dbeta)/(z*z)
//    // The correct formula for d/dβ[dz_dλ] = d/dβ[-w^λ*ln(w)]
//    // = -λ*w^(λ-1)*dw/dβ*ln(w) - w^λ*(1/w)*dw/dβ
//    double d_dbeta_dz_dlambda = -lambda * std::pow(w, lambda-1.0) * dw_dbeta * ln_w -
//      w_lambda_val * (dw_dbeta / w);
//    double d2L8_dbeta_dlambda = delta * ((d_dbeta_dz_dlambda / z) - (dz_dlambda*dz_dbeta)/(z*z));
//
//    // Mixed derivative L8 (δ,λ): d²/(dδ dλ) ln(z) = d/dλ [δ * (dz_dlambda / z)] = dz_dlambda / z
//    double d2L8_ddelta_dlambda = dz_dlambda / z;
//
//    // --- ACCUMULATING CONTRIBUTIONS TO THE HESSIAN MATRIX ---
//    // Index: 0 = α, 1 = β, 2 = δ, 3 = λ
//
//    // H(α,α): sum of L2, L6, L7, and L8 (constants already added)
//    H(0,0) += d2L6_dalpha2 + d2L7_dalpha2 + d2L8_dalpha2;
//
//    // H(α,β): mixed from L6, L7, and L8
//    H(0,1) += d2L6_dalpha_dbeta + d2L7_dalpha_dbeta + d2L8_dalpha_dbeta;
//    H(1,0) = H(0,1);
//
//    // H(β,β): contributions from L3, L7, and L8
//    H(1,1) += d2L7_dbeta2 + d2L8_dbeta2;
//
//    // H(λ,λ): contributions from L1 and L8 (L1 already added)
//    H(3,3) += d2L8_dlambda2;
//
//    // H(α,δ): From L8 - d/dα ln(z)
//    H(0,2) += dz_dalpha / z;
//    H(2,0) = H(0,2);
//
//    // H(β,δ): From L8 - d/dβ ln(z)
//    H(1,2) += dz_dbeta / z;
//    H(2,1) = H(1,2);
//
//    // Accumulating terms for mixed derivatives with λ
//    // The term from L7 is Σ ln(w) derived w.r.t. α and Σ ln(w) derived w.r.t. β
//    acc_alpha_lambda += (dw_dalpha / w) + d2L8_dalpha_dlambda;
//    acc_beta_lambda += (dw_dbeta / w) + d2L8_dbeta_dlambda;
//    acc_delta_lambda += d2L8_ddelta_dlambda;
//
//  } // end of loop
//
//  // Applying mixed derivatives with λ
//  // H(α,λ):
//  H(0,3) = acc_alpha_lambda;
//  H(3,0) = H(0,3);
//
//  // H(β,λ):
//  H(1,3) = acc_beta_lambda;
//  H(3,1) = H(1,3);
//
//  // H(δ,λ):
//  H(2,3) = acc_delta_lambda;
//  H(3,2) = H(2,3);
//
//  // Returns the analytic Hessian matrix of the negative log-likelihood
//  return Rcpp::wrap(-H);
// }





/*
----------------------------------------------------------------------------
REUSE OF THE NUMERIC STABILITY FUNCTIONS AND CHECKS
----------------------------------------------------------------------------
NOTE: We assume the following inline functions are already implemented
and available in the compilation environment, as requested:
- log1mexp(double)
- log1pexp(double)
- safe_log(double)
- safe_exp(double)
- safe_pow(double, double)
- check_pars(...) or any analogous checker
- etc.

For clarity, below we define a small parameter-check function specifically
for the Beta-Kumaraswamy (BKw) distribution, which needs alpha>0, beta>0,
gamma>0, delta>=0. We do NOT have lambda in BKw because it's fixed at 1.
*/

// -----------------------------------------------------------------------------
// Parameter Checker for Beta-Kumaraswamy (BKw) distribution
// -----------------------------------------------------------------------------
inline bool check_bkw_pars(double alpha,
                         double beta,
                         double gamma,
                         double delta,
                         bool strict = false) {
if (alpha <= 0.0 || beta <= 0.0 || gamma <= 0.0 || delta < 0.0) {
  return false;
}

if (strict) {
  // Optional stricter numeric bounds
  const double MIN_PARAM = 1e-5;
  const double MAX_PARAM = 1e5;
  if (alpha < MIN_PARAM || beta < MIN_PARAM || gamma < MIN_PARAM || delta > MAX_PARAM) {
    return false;
  }
  if (alpha > MAX_PARAM || beta > MAX_PARAM || gamma > MAX_PARAM) {
    return false;
  }
}
return true;
}

/*
----------------------------------------------------------------------------
BETA-KUMARASWAMY (BKw) DISTRIBUTION
----------------------------------------------------------------------------
PDF:
f(x; α, β, γ, δ) = (α β / B(γ, δ+1)) x^(α-1) (1 - x^α)^( β(δ+1) - 1 )
[ 1 - (1 - x^α)^β ]^(γ - 1)

CDF:
F(x; α, β, γ, δ) = I_{ [1 - (1 - x^α)^β ] } ( γ, δ + 1 )

QUANTILE:
Q(p; α, β, γ, δ) = { 1 - [1 - ( I^{-1}_{p}(γ, δ+1) ) ]^(1/β) }^(1/α)
(But see transformations step-by-step in code for numeric stability.)

RNG:
If V ~ Beta(γ, δ+1) then
X = { 1 - [1 - V ]^(1/β) }^(1/α)

LOG-LIKELIHOOD:
ℓ(θ) = n log(α β) - n log B(γ, δ+1)
+ Σ { (α-1) log(x_i) + [β(δ+1)-1] log(1 - x_i^α) + (γ - 1) log( 1 - (1 - x_i^α)^β ) }

This file defines:
- dbkw()  : density
- pbkw()  : cumulative distribution
- qbkw()  : quantile
- rbkw()  : random generation
- llbkw() : negative log-likelihood
*/


// -----------------------------------------------------------------------------
// 1) dbkw: PDF of Beta-Kumaraswamy
// -----------------------------------------------------------------------------

//' @title Density Function for Beta-Kumaraswamy Distribution
//'
//' @description
//' Computes the PDF of the Beta-Kumaraswamy distribution BKw(α, β, γ, δ) for x in (0, 1).
//' Optionally returns the log of the PDF.
//'
//' @param x Vector of quantiles in (0,1).
//' @param alpha Shape parameter α > 0.
//' @param beta Shape parameter β > 0.
//' @param gamma Shape parameter γ > 0.
//' @param delta Shape parameter δ ≥ 0.
//' @param log_prob Logical; if TRUE, returns log(f(x)).
//'
//' @return A vector of length max(length(x), lengths of parameters) with computed densities
//'         or log-densities. Values of x outside (0,1) yield 0 or -Inf (log-scale).
//'
//' @details
//' The PDF is given by
//' \deqn{
//'   f(x) = \frac{\alpha \beta}{B(\gamma, \delta+1)} \; x^{\alpha - 1} \; \bigl(1 - x^\alpha\bigr)^{\beta(\delta+1) - 1}\;
//'          \bigl[1 - \bigl(1 - x^\alpha\bigr)^\beta\bigr]^{\gamma - 1}, \quad 0<x<1.
//' }
//'
//' Note that BKw is a special case of the Generalized Kumaraswamy (GKw) with \eqn{\lambda=1}.
//'
//' @export
// [[Rcpp::export]]
arma::vec dbkw(
   const arma::vec& x,
   const Rcpp::NumericVector& alpha,
   const Rcpp::NumericVector& beta,
   const Rcpp::NumericVector& gamma,
   const Rcpp::NumericVector& delta,
   bool log_prob = false
) {
 // Convert to arma::vec
 arma::vec alpha_vec(alpha.begin(), alpha.size());
 arma::vec beta_vec(beta.begin(), beta.size());
 arma::vec gamma_vec(gamma.begin(), gamma.size());
 arma::vec delta_vec(delta.begin(), delta.size());

 // Broadcast length
 size_t n = std::max({x.n_elem,
                     alpha_vec.n_elem,
                     beta_vec.n_elem,
                     gamma_vec.n_elem,
                     delta_vec.n_elem});

 // Result
 arma::vec result(n);
 result.fill(log_prob ? R_NegInf : 0.0);

 for (size_t i = 0; i < n; ++i) {
   double a = alpha_vec[i % alpha_vec.n_elem];
   double b = beta_vec[i % beta_vec.n_elem];
   double g = gamma_vec[i % gamma_vec.n_elem];
   double d = delta_vec[i % delta_vec.n_elem];
   double xx = x[i % x.n_elem];

   // Check parameter validity
   if (!check_bkw_pars(a, b, g, d)) {
     // Invalid params => density = 0 or -Inf
     continue;
   }

   // Outside (0,1) => density=0 or log_density=-Inf
   if (xx <= 0.0 || xx >= 1.0 || !R_finite(xx)) {
     continue;
   }

   // PDF formula
   // f(x) = (alpha*beta / B(gamma, delta+1)) *
   //        x^(alpha-1) * (1 - x^alpha)^(beta*(delta+1) - 1) *
   //        [1 - (1 - x^alpha)^beta]^(gamma - 1)

   // Precompute log_B = lbeta(g, d+1)
   double logB = R::lbeta(g, d + 1.0);
   double log_const = std::log(a) + std::log(b) - logB;

   double lx = std::log(xx);
   double xalpha = a * lx;                    // log(x^alpha) = a * log(x)
   double log_1_minus_xalpha = log1mexp(xalpha);

   // (beta*(delta+1) - 1) * log(1 - x^alpha)
   double exponent1 = b * (d + 1.0) - 1.0;
   double term1 = exponent1 * log_1_minus_xalpha;

   // [1 - (1 - x^alpha)^beta]^(gamma - 1)
   // log(1 - (1 - x^alpha)^beta) = log1mexp( b * log(1 - x^alpha) )
   double log_1_minus_xalpha_beta = b * log_1_minus_xalpha;
   double log_bracket = log1mexp(log_1_minus_xalpha_beta);
   double exponent2 = g - 1.0;
   double term2 = exponent2 * log_bracket;

   // Full log pdf
   double log_pdf = log_const +
     (a - 1.0) * lx +
     term1 +
     term2;

   if (log_prob) {
     result(i) = log_pdf;
   } else {
     // exp safely
     result(i) = std::exp(log_pdf);
   }
 }

 return result;
}


// -----------------------------------------------------------------------------
// 2) pbkw: CDF of Beta-Kumaraswamy
// -----------------------------------------------------------------------------

//' @title CDF of Beta-Kumaraswamy Distribution
//'
//' @description
//' Computes the cumulative distribution function for BKw(α, β, γ, δ).
//'
//' @param q Vector of quantiles in (0,1).
//' @param alpha Shape parameter α > 0.
//' @param beta Shape parameter β > 0.
//' @param gamma Shape parameter γ > 0.
//' @param delta Shape parameter δ ≥ 0.
//' @param lower_tail Logical; if TRUE (default), returns F(q)=P(X≤q), else 1-F(q).
//' @param log_p Logical; if TRUE, returns log-probabilities.
//'
//' @details
//' CDF: \eqn{F(x) = I_{1 - (1 - x^α)^β}(\gamma, \delta+1)}, where \eqn{I_z(a,b)} is the
//' regularized incomplete beta function.
//'
//' @return A vector of probabilities, matching the length of the broadcast inputs.
//'
//' @export
// [[Rcpp::export]]
arma::vec pbkw(
   const arma::vec& q,
   const Rcpp::NumericVector& alpha,
   const Rcpp::NumericVector& beta,
   const Rcpp::NumericVector& gamma,
   const Rcpp::NumericVector& delta,
   bool lower_tail = true,
   bool log_p = false
) {
 // Convert
 arma::vec alpha_vec(alpha.begin(), alpha.size());
 arma::vec beta_vec(beta.begin(), beta.size());
 arma::vec gamma_vec(gamma.begin(), gamma.size());
 arma::vec delta_vec(delta.begin(), delta.size());

 // Broadcast
 size_t n = std::max({q.n_elem,
                     alpha_vec.n_elem,
                     beta_vec.n_elem,
                     gamma_vec.n_elem,
                     delta_vec.n_elem});

 arma::vec res(n);

 for (size_t i = 0; i < n; ++i) {
   double a = alpha_vec[i % alpha_vec.n_elem];
   double b = beta_vec[i % beta_vec.n_elem];
   double g = gamma_vec[i % gamma_vec.n_elem];
   double d = delta_vec[i % delta_vec.n_elem];
   double xx = q[i % q.n_elem];

   if (!check_bkw_pars(a, b, g, d)) {
     res(i) = NA_REAL;
     continue;
   }

   if (!R_finite(xx) || xx <= 0.0) {
     // x=0 => F=0
     double prob0 = lower_tail ? 0.0 : 1.0;
     res(i) = log_p ? std::log(prob0) : prob0;
     continue;
   }

   if (xx >= 1.0) {
     // x=1 => F=1
     double prob1 = lower_tail ? 1.0 : 0.0;
     res(i) = log_p ? std::log(prob1) : prob1;
     continue;
   }

   // We want z = 1 - (1 - x^alpha)^beta
   double lx = std::log(xx);
   double xalpha = std::exp(a * lx);
   double one_minus_xalpha = 1.0 - xalpha;

   if (one_minus_xalpha <= 0.0) {
     // F(x) ~ 1 if x^alpha>=1
     double prob1 = lower_tail ? 1.0 : 0.0;
     res(i) = log_p ? std::log(prob1) : prob1;
     continue;
   }

   double temp = 1.0 - std::pow(one_minus_xalpha, b);
   if (temp <= 0.0) {
     double prob0 = lower_tail ? 0.0 : 1.0;
     res(i) = log_p ? std::log(prob0) : prob0;
     continue;
   }

   if (temp >= 1.0) {
     double prob1 = lower_tail ? 1.0 : 0.0;
     res(i) = log_p ? std::log(prob1) : prob1;
     continue;
   }

   // Then F(x) = pbeta(temp, gamma, delta+1, TRUE, FALSE)
   double val = R::pbeta(temp, g, d+1.0, true, false); // F
   if (!lower_tail) {
     val = 1.0 - val;
   }
   if (log_p) {
     val = std::log(val);
   }
   res(i) = val;
 }

 return res;
}


// -----------------------------------------------------------------------------
// 3) qbkw: QUANTILE of Beta-Kumaraswamy
// -----------------------------------------------------------------------------

//' @title Quantile Function of Beta-Kumaraswamy Distribution
//'
//' @description
//' Computes the quantile function for BKw(α, β, γ, δ). For p in (0,1), returns x such that
//' P(X ≤ x) = p under the BKw distribution.
//'
//' @param p Vector of probabilities in (0,1).
//' @param alpha Shape parameter α > 0.
//' @param beta Shape parameter β > 0.
//' @param gamma Shape parameter γ > 0.
//' @param delta Shape parameter δ ≥ 0.
//' @param lower_tail Logical; if TRUE (default), p is F(x)=P(X≤x), otherwise p=1-F(x).
//' @param log_p Logical; if TRUE, p is given as log(p).
//'
//' @details
//' Inversion approach using \eqn{Q(p) = \{1 - [1 - qbeta(p, γ, δ+1)]^(1/β)\}^{1/\alpha}}.
//'
//' @return A vector of quantiles of the same length as the broadcast of p and parameters.
//'
//' @export
// [[Rcpp::export]]
arma::vec qbkw(
   const arma::vec& p,
   const Rcpp::NumericVector& alpha,
   const Rcpp::NumericVector& beta,
   const Rcpp::NumericVector& gamma,
   const Rcpp::NumericVector& delta,
   bool lower_tail = true,
   bool log_p = false
) {
 arma::vec alpha_vec(alpha.begin(), alpha.size());
 arma::vec beta_vec(beta.begin(), beta.size());
 arma::vec gamma_vec(gamma.begin(), gamma.size());
 arma::vec delta_vec(delta.begin(), delta.size());

 size_t n = std::max({p.n_elem,
                     alpha_vec.n_elem,
                     beta_vec.n_elem,
                     gamma_vec.n_elem,
                     delta_vec.n_elem});

 arma::vec res(n);

 for (size_t i = 0; i < n; ++i) {
   double a = alpha_vec[i % alpha_vec.n_elem];
   double b = beta_vec[i % beta_vec.n_elem];
   double g = gamma_vec[i % gamma_vec.n_elem];
   double d = delta_vec[i % delta_vec.n_elem];
   double pp = p[i % p.n_elem];

   if (!check_bkw_pars(a, b, g, d)) {
     res(i) = NA_REAL;
     continue;
   }

   // Convert from log_p if needed
   if (log_p) {
     if (pp > 0.0) {
       // log(p) > 0 => p>1 => invalid
       res(i) = NA_REAL;
       continue;
     }
     pp = std::exp(pp);
   }
   // Convert if upper tail
   if (!lower_tail) {
     pp = 1.0 - pp;
   }

   // Check boundaries
   if (pp <= 0.0) {
     res(i) = 0.0;
     continue;
   } else if (pp >= 1.0) {
     res(i) = 1.0;
     continue;
   }

   // We do: y = qbeta(pp, gamma, delta+1)
   double y = R::qbeta(pp, g, d+1.0, true, false);
   if (y <= 0.0) {
     res(i) = 0.0;
     continue;
   } else if (y >= 1.0) {
     res(i) = 1.0;
     continue;
   }

   // Then x = {1 - [1 - y]^(1/b)}^(1/a)
   double part = 1.0 - y;
   if (part <= 0.0) {
     res(i) = 1.0;
     continue;
   } else if (part >= 1.0) {
     res(i) = 0.0;
     continue;
   }

   double inner = std::pow(part, 1.0/b);
   double xval = 1.0 - inner;
   if (xval < 0.0)  xval = 0.0;
   if (xval > 1.0)  xval = 1.0;

   if (a == 1.0) {
     // small optimization
     res(i) = xval;
   } else {
     double qv = std::pow(xval, 1.0/a);
     if (qv < 0.0)      qv = 0.0;
     else if (qv > 1.0) qv = 1.0;
     res(i) = qv;
   }
 }

 return res;
}


// -----------------------------------------------------------------------------
// 4) rbkw: RNG for Beta-Kumaraswamy
// -----------------------------------------------------------------------------

//' @title Random Generation from Beta-Kumaraswamy Distribution
//'
//' @description
//' Generates n samples from BKw(α, β, γ, δ) using the transformation:
//'   \eqn{V ~ Beta(γ, δ+1) => X = {1 - (1 - V)^(1/β)}^(1/α)}.
//'
//' @param n Integer number of observations to generate.
//' @param alpha Shape parameter α > 0.
//' @param beta Shape parameter β > 0.
//' @param gamma Shape parameter γ > 0.
//' @param delta Shape parameter δ ≥ 0.
//'
//' @return A numeric vector of length n with random draws from the BKw distribution.
//'
//' @export
// [[Rcpp::export]]
arma::vec rbkw(
   int n,
   const Rcpp::NumericVector& alpha,
   const Rcpp::NumericVector& beta,
   const Rcpp::NumericVector& gamma,
   const Rcpp::NumericVector& delta
) {
 if (n <= 0) {
   Rcpp::stop("rbkw: n must be positive");
 }

 arma::vec alpha_vec(alpha.begin(), alpha.size());
 arma::vec beta_vec(beta.begin(), beta.size());
 arma::vec gamma_vec(gamma.begin(), gamma.size());
 arma::vec delta_vec(delta.begin(), delta.size());

 size_t k = std::max({alpha_vec.n_elem,
                     beta_vec.n_elem,
                     gamma_vec.n_elem,
                     delta_vec.n_elem});

 arma::vec out(n);

 for (int i = 0; i < n; ++i) {
   size_t idx = i % k;
   double a = alpha_vec[idx % alpha_vec.n_elem];
   double b = beta_vec[idx % beta_vec.n_elem];
   double g = gamma_vec[idx % gamma_vec.n_elem];
   double d = delta_vec[idx % delta_vec.n_elem];

   if (!check_bkw_pars(a, b, g, d)) {
     out(i) = NA_REAL;
     Rcpp::warning("rbkw: invalid parameters at index %d", i+1);
     continue;
   }

   // V ~ Beta(g, d+1)
   double V = R::rbeta(g, d + 1.0);
   // X = {1 - (1 - V)^(1/b)}^(1/a)
   double one_minus_V = 1.0 - V;
   if (one_minus_V <= 0.0) {
     out(i) = 1.0;
     continue;
   }
   if (one_minus_V >= 1.0) {
     out(i) = 0.0;
     continue;
   }

   double temp = std::pow(one_minus_V, 1.0/b);
   double xval = 1.0 - temp;
   if (xval < 0.0)  xval = 0.0;
   if (xval > 1.0)  xval = 1.0;

   if (a == 1.0) {
     out(i) = xval;
   } else {
     double rv = std::pow(xval, 1.0/a);
     if (rv < 0.0) rv = 0.0;
     if (rv > 1.0) rv = 1.0;
     out(i) = rv;
   }
 }

 return out;
}


// -----------------------------------------------------------------------------
// 5) llbkw: Negative Log-Likelihood for Beta-Kumaraswamy
// -----------------------------------------------------------------------------

//' @title Negative Log-Likelihood for Beta-Kumaraswamy Distribution
//'
//' @description
//' Computes the negative log-likelihood of the BKw(α, β, γ, δ) model given data in (0,1).
//'
//' @param par NumericVector of length 4, \eqn{(α, β, γ, δ)}.
//' @param data NumericVector of observations, each in (0,1).
//'
//' @details
//' The PDF is
//' \deqn{
//'   f(x)=\frac{\alpha \beta}{B(\gamma, \delta+1)} x^{\alpha-1} (1-x^\alpha)^{\beta(\delta+1)-1} [1-(1-x^\alpha)^\beta]^{\gamma-1}.
//' }
//' The log-likelihood \eqn{\ell(\theta)} is the sum of log(f(x_i)) for i=1..n. We return
//' the negative log-likelihood for convenience in optimization.
//'
//' @return A single numeric value (the negative log-likelihood). Returns \code{-Inf} if
//'         parameters are invalid or if any data point is outside (0,1).
//'
//' @export
// [[Rcpp::export]]
double llbkw(const Rcpp::NumericVector& par,
            const Rcpp::NumericVector& data) {
 if (par.size() < 4) {
   return R_NegInf;
 }
 double a = par[0];  // alpha>0
 double b = par[1];  // beta>0
 double g = par[2];  // gamma>0
 double d = par[3];  // delta>=0

 if (!check_bkw_pars(a, b, g, d)) {
   return R_NegInf;
 }

 // Convert data
 arma::vec x = Rcpp::as<arma::vec>(data);
 // Check data
 if (arma::any(x <= 0.0) || arma::any(x >= 1.0)) {
   return R_NegInf;
 }

 int n = x.n_elem;
 double logB = R::lbeta(g, d + 1.0);

 // sum( log(f(x)) ) = n*log(alpha*beta) - n*logB + ...
 double ll_const = n * (std::log(a) + std::log(b) - logB);

 // (alpha - 1)*sum(log(x))
 arma::vec lx = arma::log(x);
 double sum1 = (a - 1.0) * arma::sum(lx);

 // exponent1 = beta*(delta+1) - 1
 double exp1 = b*(d+1.0) - 1.0;
 // sum( exp1 * log(1 - x^alpha) )
 arma::vec xalpha = arma::pow(x, a);
 arma::vec log_1_xalpha = arma::log(1.0 - xalpha);
 double sum2 = exp1 * arma::sum(log_1_xalpha);

 // exponent2 = (gamma - 1)
 // sum( exponent2* log( 1 - (1 - x^alpha)^beta ) )
 // (1 - x^alpha)^beta => vbeta
 arma::vec vbeta = arma::pow(1.0 - xalpha, b);
 arma::vec one_minus_vbeta = 1.0 - vbeta;
 // safe log(1.0 - vbeta)
 arma::vec log_omv(vbeta.n_elem);
 for (int i = 0; i < (int)vbeta.n_elem; i++) {
   double val = one_minus_vbeta[i];
   if (val <= 0.0) {
     return R_NegInf;  // invalid => -Inf
   }
   log_omv[i] = std::log(val);
 }
 double sum3 = (g - 1.0) * arma::sum(log_omv);

 double ll = ll_const + sum1 + sum2 + sum3;

 // negative log-likelihood
 return -ll;
}


//' @title Gradient Function for Beta-Kumaraswamy Log-Likelihood
//'
//' @description
//' Calculates the gradient vector (partial derivatives) of the negative log-likelihood
//' function for the Beta-Kumaraswamy (BKw) distribution. This function provides
//' the exact gradient needed for efficient optimization in maximum likelihood estimation.
//' The BKw is a submodel of GKw with λ = 1 fixed.
//'
//' @param par NumericVector of length 4 containing parameters (α, β, γ, δ) in that order.
//'        All parameters must be positive.
//' @param data NumericVector of observations, where all values must be in the open interval (0,1).
//'
//' @return NumericVector of length 4 containing the gradient components (partial derivatives)
//'         of the negative log-likelihood with respect to each parameter (α, β, γ, δ).
//'         Returns a vector of NaN values if any parameters or data values are invalid.
//'
//' @details
//' The gradient vector contains the following partial derivatives of the negative log-likelihood:
//'
//' \deqn{
//' \frac{\partial \ell}{\partial \alpha} = \frac{n}{\alpha} + \sum_{i=1}^{n}\log(x_i) -
//' \sum_{i=1}^{n}\left[x_i^{\alpha} \log(x_i) \left(\frac{\beta(\delta+1)-1}{v_i} -
//' \frac{(\gamma-1) \beta v_i^{\beta-1}}{w_i}\right)\right]
//' }
//'
//' \deqn{
//' \frac{\partial \ell}{\partial \beta} = \frac{n}{\beta} + (\delta+1)\sum_{i=1}^{n}\log(v_i) -
//' \sum_{i=1}^{n}\left[v_i^{\beta} \log(v_i) \frac{\gamma-1}{w_i}\right]
//' }
//'
//' \deqn{
//' \frac{\partial \ell}{\partial \gamma} = -n[\psi(\gamma) - \psi(\gamma+\delta+1)] +
//' \sum_{i=1}^{n}\log(w_i)
//' }
//'
//' \deqn{
//' \frac{\partial \ell}{\partial \delta} = -n[\psi(\delta+1) - \psi(\gamma+\delta+1)] +
//' \beta\sum_{i=1}^{n}\log(v_i)
//' }
//'
//' where:
//' \itemize{
//'   \item \deqn{v_i = 1 - x_i^{\alpha}}
//'   \item \deqn{w_i = 1 - v_i^{\beta} = 1 - (1-x_i^{\alpha})^{\beta}}
//'   \item \deqn{\psi} is the digamma function (derivative of the log-gamma function)
//' }
//'
//' The implementation includes several numerical safeguards:
//' \itemize{
//'   \item Parameter and data validation with appropriate error handling
//'   \item Clamping of intermediate values to avoid numerical underflow/overflow
//'   \item Efficient vector operations using Armadillo C++ library
//' }
//'
//' The returned gradient is negated to align with minimization of negative log-likelihood
//' in optimization routines.
//'
//' @examples
//' \dontrun{
//' # Generate sample data from a BKw distribution
//' set.seed(123)
//' x <- rbkw(100, 2, 3, 1, 0.5)
//' hist(x, breaks = 20, main = "BKw(2, 3, 1, 0.5) Sample")
//'
//' # Use in optimization with Hessian-based methods
//' result <- optim(c(0.5, 0.5, 0.5, 0.5), llbkw, method = "BFGS",
//'                 hessian = TRUE, data = x)
//'
//' # Compare numerical and analytical derivatives
//' num_grad <- numDeriv::grad(llbkw, x = result$par, data = x)
//' num_hess <- numDeriv::hessian(llbkw, x = result$par, data = x)
//'
//' ana_grad <- grbkw(result$par, data = x)
//' ana_hess <- hsbkw(result$par, data = x)
//'
//' # Check differences (should be very small)
//' round(num_grad - ana_grad, 4)
//' round(num_hess - ana_hess, 4)
//'
//' }
//'
//' @references
//' Kumaraswamy, P. (1980). A generalized probability density function for double-bounded random processes.
//' Journal of Hydrology, 46(1-2), 79-88.
//'
//' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized distributions.
//' Journal of Statistical Computation and Simulation, 81(7), 883-898.
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector grbkw(const Rcpp::NumericVector& par, const Rcpp::NumericVector& data) {
 // Parameter extraction
 double alpha = par[0];   // Shape parameter α > 0
 double beta = par[1];    // Shape parameter β > 0
 double gamma = par[2];   // Shape parameter γ > 0
 double delta = par[3];   // Shape parameter δ > 0

 // Parameter validation
 if (alpha <= 0 || beta <= 0 || gamma <= 0 || delta <= 0) {
   Rcpp::NumericVector grad(4, R_NaN);
   return grad;
 }

 // Data conversion and validation
 arma::vec x = Rcpp::as<arma::vec>(data);

 if (arma::any(x <= 0) || arma::any(x >= 1)) {
   Rcpp::NumericVector grad(4, R_NaN);
   return grad;
 }

 int n = x.n_elem;  // Sample size

 // Initialize gradient vector
 Rcpp::NumericVector grad(4, 0.0);

 // Small constant to avoid numerical issues
 double eps = std::numeric_limits<double>::epsilon() * 100;

 // Compute transformations and intermediate values
 arma::vec log_x = arma::log(x);                // log(x_i)
 arma::vec x_alpha = arma::pow(x, alpha);       // x_i^α
 arma::vec x_alpha_log_x = x_alpha % log_x;     // x_i^α * log(x_i)

 // v_i = 1 - x_i^α
 arma::vec v = 1.0 - x_alpha;
 v = arma::clamp(v, eps, 1.0 - eps);            // Prevent numerical issues

 arma::vec log_v = arma::log(v);                // log(v_i)
 arma::vec v_beta_m1 = arma::pow(v, beta - 1.0); // v_i^(β-1)
 arma::vec v_beta = arma::pow(v, beta);          // v_i^β
 arma::vec v_beta_log_v = v_beta % log_v;        // v_i^β * log(v_i)

 // w_i = 1 - v_i^β = 1 - (1-x_i^α)^β
 arma::vec w = 1.0 - v_beta;
 w = arma::clamp(w, eps, 1.0 - eps);            // Prevent numerical issues

 arma::vec log_w = arma::log(w);                // log(w_i)

 // Calculate partial derivatives for each parameter (for log-likelihood)

 // ∂ℓ/∂α = n/α + Σᵢlog(xᵢ) - Σᵢ[xᵢ^α * log(xᵢ) * ((β(δ+1)-1)/vᵢ - (γ-1) * β * vᵢ^(β-1) / wᵢ)]
 double d_alpha = n / alpha + arma::sum(log_x);

 // Calculate the terms in the α gradient
 arma::vec alpha_term1 = (beta * (delta + 1) - 1.0) / v;           // (β(δ+1)-1)/v_i
 arma::vec alpha_term2 = (gamma - 1.0) * beta * v_beta_m1 / w;     // (γ-1) * β * v_i^(β-1) / w_i

 d_alpha -= arma::sum(x_alpha_log_x % (alpha_term1 - alpha_term2));

 // ∂ℓ/∂β = n/β + (δ+1)Σᵢlog(vᵢ) - Σᵢ[vᵢ^β * log(vᵢ) * ((γ-1) / wᵢ)]
 double d_beta = n / beta + (delta + 1) * arma::sum(log_v);

 // Calculate the term in the β gradient
 arma::vec beta_term = (gamma - 1.0) / w;                          // (γ-1) / w_i

 d_beta -= arma::sum(v_beta_log_v % beta_term);

 // ∂ℓ/∂γ = -n[ψ(γ) - ψ(γ+δ+1)] + Σᵢlog(wᵢ)
 double d_gamma = -n * (R::digamma(gamma) - R::digamma(gamma + delta + 1)) + arma::sum(log_w);

 // ∂ℓ/∂δ = -n[ψ(δ+1) - ψ(γ+δ+1)] + βΣᵢlog(vᵢ)
 double d_delta = -n * (R::digamma(delta + 1) - R::digamma(gamma + delta + 1)) + beta * arma::sum(log_v);

 // Since we're optimizing negative log-likelihood, negate all derivatives
 grad[0] = -d_alpha;
 grad[1] = -d_beta;
 grad[2] = -d_gamma;
 grad[3] = -d_delta;

 return grad;
}



//' @title Analytic Hessian Matrix for Beta-Kumaraswamy Distribution
//'
//' @description
//' Computes the analytic Hessian matrix of the log-likelihood function for
//' the Beta-Kumaraswamy (BKw) distribution. This function provides
//' exact second derivatives needed for optimization and inference.
//'
//' @param par Numeric vector of length 4 containing the parameters
//'        (α, β, γ, δ) in that order. All parameters must be positive.
//' @param data Numeric vector of observations, where all values must be
//'        in the open interval (0,1).
//'
//' @return A 4×4 numeric matrix representing the Hessian of the negative
//'         log-likelihood function. If parameters or data are invalid
//'         (parameters ≤ 0 or data outside (0,1)), returns a matrix of
//'         NaN values.
//'
//' @details
//' The log-likelihood for the Beta-Kumaraswamy distribution is:
//'
//' \deqn{
//' \ell(\theta) = n \ln(\alpha) + n \ln(\beta) - n \ln B(\gamma, \delta+1)
//' + (\alpha-1) \sum \ln(x_i)
//' + (\beta(\delta+1)-1) \sum \ln(1 - x_i^\alpha)
//' + (\gamma-1) \sum \ln\{1 - (1 - x_i^\alpha)^\beta\}
//' }
//'
//' where λ is fixed at 1 for this distribution.
//'
//' The implementation computes all second derivatives analytically for each term.
//' For computational efficiency, the following transformations are used:
//' \itemize{
//'   \item \deqn{A = x^α} and derivatives
//'   \item \deqn{v = 1 - A}
//'   \item \deqn{w = 1 - v^β}
//' }
//'
//' The returned Hessian matrix has the following structure:
//' \itemize{
//'   \item Rows/columns 1-4 correspond to α, β, γ, δ respectively
//'   \item The matrix is symmetric (as expected for a Hessian)
//'   \item The matrix represents second derivatives of the negative log-likelihood
//' }
//'
//' This function is implemented in C++ for computational efficiency.
//'
//' @examples
//' \dontrun{
//' # Generate sample data from a BKw distribution
//' set.seed(123)
//' x <- rbkw(100, 2, 3, 1, 0.5)
//' hist(x, breaks = 20, main = "BKw(2, 3, 1, 0.5) Sample")
//'
//' # Use in optimization with Hessian-based methods
//' result <- optim(c(0.5, 0.5, 0.5, 0.5), llbkw, method = "BFGS",
//'                 hessian = TRUE, data = x)
//'
//' # Compare numerical and analytical derivatives
//' num_grad <- numDeriv::grad(llbkw, x = result$par, data = x)
//' num_hess <- numDeriv::hessian(llbkw, x = result$par, data = x)
//'
//' ana_grad <- grbkw(result$par, data = x)
//' ana_hess <- hsbkw(result$par, data = x)
//'
//' # Check differences (should be very small)
//' round(num_grad - ana_grad, 4)
//' round(num_hess - ana_hess, 4)
//'
//' }
//'
//' @references
//' Kumaraswamy, P. (1980). A generalized probability density function for double-bounded random processes.
//' Journal of Hydrology, 46(1-2), 79-88.
//'
//' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized distributions.
//' Journal of Statistical Computation and Simulation, 81(7), 883-898.
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix hsbkw(const Rcpp::NumericVector& par, const Rcpp::NumericVector& data) {
 // Parameter extraction
 double alpha  = par[0];   // θ[0] = α
 double beta   = par[1];   // θ[1] = β
 double gamma  = par[2];   // θ[2] = γ
 double delta  = par[3];   // θ[3] = δ

 // Simple parameter validation (all > 0)
 if(alpha <= 0 || beta <= 0 || gamma <= 0 || delta <= 0) {
   Rcpp::NumericMatrix nanH(4,4);
   nanH.fill(R_NaN);
   return nanH;
 }

 // Data conversion and basic validation
 arma::vec x = Rcpp::as<arma::vec>(data);
 if(arma::any(x <= 0) || arma::any(x >= 1)) {
   Rcpp::NumericMatrix nanH(4,4);
   nanH.fill(R_NaN);
   return nanH;
 }

 int n = x.n_elem;  // sample size

 // Initialize Hessian matrix H (of ℓ(θ)) as 4x4
 arma::mat H(4,4, arma::fill::zeros);

 // --- CONSTANT TERMS (do not depend on x) ---
 // L1: n ln(α)  => d²/dα² = -n/α²
 H(0,0) += -n/(alpha*alpha);
 // L2: n ln(β)  => d²/dβ² = -n/β²
 H(1,1) += -n/(beta*beta);
 // L3: -n ln[B(γ, δ+1)]
 //   d²/dγ² = -n [ψ₁(γ) - ψ₁(γ+δ+1)]  where ψ₁ is the trigamma function
 H(2,2) += -n * (R::trigamma(gamma) - R::trigamma(gamma+delta+1));
 //   d²/dδ² = -n [ψ₁(δ+1) - ψ₁(γ+δ+1)]
 H(3,3) += -n * (R::trigamma(delta+1) - R::trigamma(gamma+delta+1));
 //   Mixed derivative (γ,δ): = n ψ₁(γ+δ+1)
 H(2,3) += n * R::trigamma(gamma+delta+1);
 H(3,2) = H(2,3);

 // --- TERMS THAT INVOLVE THE OBSERVATIONS ---
 // Loop over each observation to accumulate contributions from:
 // L4: (α-1) Σ ln(x_i)  --> contributes only to first derivatives
 // L5: (β(δ+1)-1) Σ ln(v), where v = 1 - x^α
 // L6: (γ-1) Σ ln(w), where w = 1 - v^β
 for (int i = 0; i < n; i++) {
   double xi    = x(i);
   double ln_xi = std::log(xi);

   // -- Compute A = x^α and its derivatives --
   double A = std::pow(xi, alpha);                  // A = x^α
   double dA_dalpha = A * ln_xi;                    // dA/dα = x^α ln(x)
   double d2A_dalpha2 = A * ln_xi * ln_xi;          // d²A/dα² = x^α (ln(x))²

   // -- v = 1 - A and its derivatives --
   double v = 1.0 - A;                              // v = 1 - x^α
   double ln_v = std::log(v);                       // ln(v)
   double dv_dalpha = -dA_dalpha;                   // dv/dα = -dA/dα = -x^α ln(x)
   double d2v_dalpha2 = -d2A_dalpha2;               // d²v/dα² = -d²A/dα² = -x^α (ln(x))²

   // --- L5: (β(δ+1)-1) ln(v) ---
   double beta_delta_factor = beta * (delta + 1.0) - 1.0;
   // Second derivative w.r.t. α: (β(δ+1)-1)*[(d²v/dα²*v - (dv/dα)²)/v²]
   double d2L5_dalpha2 = beta_delta_factor * ((d2v_dalpha2 * v - dv_dalpha * dv_dalpha) / (v*v));
   // Mixed derivative: d²L5/(dα dβ) = d/dβ[(β(δ+1)-1)*(dv_dalpha/v)] = (δ+1)*(dv_dalpha/v)
   double d2L5_dalpha_dbeta = (delta + 1.0) * (dv_dalpha / v);
   // Mixed derivative: d²L5/(dα dδ) = d/dδ[(β(δ+1)-1)*(dv_dalpha/v)] = β*(dv_dalpha/v)
   double d2L5_dalpha_ddelta = beta * (dv_dalpha / v);
   // Mixed derivative: d²L5/(dβ dδ) = d/dδ[(δ+1)*ln(v)] = ln(v)
   double d2L5_dbeta_ddelta = ln_v;

   // --- L6: (γ - 1) ln(w), where w = 1 - v^β ---
   double v_beta = std::pow(v, beta);              // v^β
   double w = 1.0 - v_beta;                        // w = 1 - v^β

   // Derivative of w w.r.t. v: dw/dv = -β * v^(β-1)
   double dw_dv = -beta * std::pow(v, beta - 1.0);
   // Chain rule: dw/dα = dw/dv * dv/dα
   double dw_dalpha = dw_dv * dv_dalpha;
   // Second derivative w.r.t. α for L6:
   // d²/dα² ln(w) = [d²w/dα² * w - (dw/dα)²] / w²
   // Computing d²w/dα²:
   //   dw/dα = -β * v^(β-1)*dv_dalpha,
   //   d²w/dα² = -β * [(β-1)*v^(β-2)*(dv_dalpha)² + v^(β-1)*d²v_dalpha²]
   double d2w_dalpha2 = -beta * ((beta - 1.0) * std::pow(v, beta-2.0) * (dv_dalpha * dv_dalpha)
                                   + std::pow(v, beta-1.0) * d2v_dalpha2);
   double d2L6_dalpha2 = (gamma - 1.0) * ((d2w_dalpha2 * w - (dw_dalpha * dw_dalpha)) / (w*w));

   // Derivative w.r.t. β: d/dβ ln(w). Note: d/dβ(v^β) = v^β ln(v) => d/dβ w = -v^β ln(v)
   double dw_dbeta = -v_beta * ln_v;
   // Second derivative w.r.t. β for L6:
   // d²/dβ² ln(w) = [d²w/dβ² * w - (dw/dβ)²]/w², where d²w/dβ² = -v^β (ln(v))²
   double d2w_dbeta2 = -v_beta * (ln_v * ln_v);
   double d2L6_dbeta2 = (gamma - 1.0) * ((d2w_dbeta2 * w - (dw_dbeta * dw_dbeta))/(w*w));

   // Mixed derivative L6 (α,β): d²/(dα dβ) ln(w) =
   //   = d/dβ[(dw_dalpha)/w] = (d/dβ dw_dalpha)/w - (dw_dalpha*dw_dbeta)/(w*w)
   double d_dw_dalpha_dbeta = -std::pow(v, beta-1.0) * (1.0 + beta * ln_v) * dv_dalpha;
   double d2L6_dalpha_dbeta = (gamma - 1.0) * ((d_dw_dalpha_dbeta / w) - (dw_dalpha * dw_dbeta)/(w*w));

   // Mixed derivatives with γ
   double d2L6_dalpha_dgamma = dw_dalpha / w;
   double d2L6_dbeta_dgamma = dw_dbeta / w;

   // --- ACCUMULATING CONTRIBUTIONS TO THE HESSIAN MATRIX ---
   // Index: 0 = α, 1 = β, 2 = γ, 3 = δ

   // H(α,α): sum of L1, L5, and L6 (constants already added)
   H(0,0) += d2L5_dalpha2 + d2L6_dalpha2;

   // H(β,β): contributions from L2, L5, and L6
   H(1,1) += d2L6_dbeta2;

   // H(α,β): mixed from L5 and L6
   H(0,1) += d2L5_dalpha_dbeta + d2L6_dalpha_dbeta;
   H(1,0) = H(0,1);

   // H(α,γ): mixed from L6
   H(0,2) += d2L6_dalpha_dgamma;
   H(2,0) = H(0,2);

   // H(α,δ): mixed from L5
   H(0,3) += d2L5_dalpha_ddelta;
   H(3,0) = H(0,3);

   // H(β,γ): mixed from L6
   H(1,2) += d2L6_dbeta_dgamma;
   H(2,1) = H(1,2);

   // H(β,δ): mixed from L5
   H(1,3) += d2L5_dbeta_ddelta;
   H(3,1) = H(1,3);

 } // end of loop

 // Returns the analytic Hessian matrix of the negative log-likelihood
 return Rcpp::wrap(-H);
}


// //' @title Analytic Hessian Matrix for Beta-Kumaraswamy Distribution
// //'
// //' @description
// //' Computes the analytic Hessian matrix of the log-likelihood function for
// //' the Beta-Kumaraswamy (BKw) distribution. This function provides
// //' exact second derivatives needed for optimization and inference.
// //'
// //' @param par Numeric vector of length 4 containing the parameters
// //'        (α, β, γ, δ) in that order. All parameters must be positive.
// //' @param data Numeric vector of observations, where all values must be
// //'        in the open interval (0,1).
// //'
// //' @return A 4×4 numeric matrix representing the Hessian of the negative
// //'         log-likelihood function. If parameters or data are invalid
// //'         (parameters ≤ 0 or data outside (0,1)), returns a matrix of
// //'         NaN values.
// //'
// //' @details
// //' The log-likelihood for the Beta-Kumaraswamy distribution is:
// //'
// //' \deqn{
// //' \ell(\theta) = n \ln(\alpha) + n \ln(\beta) - n \ln B(\gamma, \delta+1)
// //' + (\alpha-1) \sum \ln(x_i)
// //' + (\beta(\delta+1)-1) \sum \ln(1 - x_i^\alpha)
// //' + (\gamma-1) \sum \ln\{1 - (1 - x_i^\alpha)^\beta\}
// //' }
// //'
// //' where λ is fixed at 1 for this distribution.
// //'
// //' The implementation computes all second derivatives analytically for each term.
// //' For computational efficiency, the following transformations are used:
// //' \itemize{
// //'   \item \deqn{A = x^α} and derivatives
// //'   \item \deqn{v = 1 - A}
// //'   \item \deqn{w = 1 - v^β}
// //' }
// //'
// //' The returned Hessian matrix has the following structure:
// //' \itemize{
// //'   \item Rows/columns 1-4 correspond to α, β, γ, δ respectively
// //'   \item The matrix is symmetric (as expected for a Hessian)
// //'   \item The matrix represents second derivatives of the negative log-likelihood
// //' }
// //'
// //' This function is implemented in C++ for computational efficiency.
// //'
// //' @examples
// //' \dontrun{
// //' # Generate sample data from a BKw distribution
// //' set.seed(123)
// //' x <- rbkw(100, 2, 3, 1, 0.5)
// //' hist(x, breaks = 20, main = "BKw(2, 3, 1, 0.5) Sample")
// //'
// //' # Use in optimization with Hessian-based methods
// //' result <- optim(c(0.5, 0.5, 0.5, 0.5), llbkw, method = "BFGS",
// //'                 hessian = TRUE, data = x)
// //'
// //' # Compare numerical and analytical derivatives
// //' num_grad <- numDeriv::grad(llbkw, x = result$par, data = x)
// //' num_hess <- numDeriv::hessian(llbkw, x = result$par, data = x)
// //'
// //' ana_grad <- grbkw(result$par, data = x)
// //' ana_hess <- hsbkw(result$par, data = x)
// //'
// //' # Check differences (should be very small)
// //' round(num_grad - ana_grad, 4)
// //' round(num_hess - ana_hess, 4)
// //'
// //' }
// //'
// //' @references
// //' Kumaraswamy, P. (1980). A generalized probability density function for double-bounded random processes.
// //' Journal of Hydrology, 46(1-2), 79-88.
// //'
// //' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized distributions.
// //' Journal of Statistical Computation and Simulation, 81(7), 883-898.
// //'
// //' @export
// // [[Rcpp::export]]
// Rcpp::NumericMatrix hsbkw(const Rcpp::NumericVector& par, const Rcpp::NumericVector& data) {
//  // Parameter extraction
//  double alpha  = par[0];   // θ[0] = α
//  double beta   = par[1];   // θ[1] = β
//  double gamma  = par[2];   // θ[2] = γ
//  double delta  = par[3];   // θ[3] = δ
//
//  // Fixed parameter for BKw distribution
//  const double lambda = 1.0;  // λ=1 for BKw
//
//  // Simple parameter validation (all > 0)
//  if(alpha <= 0 || beta <= 0 || gamma <= 0 || delta <= 0) {
//    Rcpp::NumericMatrix nanH(4,4);
//    nanH.fill(R_NaN);
//    return nanH;
//  }
//
//  // Data conversion and basic validation
//  arma::vec x = Rcpp::as<arma::vec>(data);
//  if(arma::any(x <= 0) || arma::any(x >= 1)) {
//    Rcpp::NumericMatrix nanH(4,4);
//    nanH.fill(R_NaN);
//    return nanH;
//  }
//
//  int n = x.n_elem;  // sample size
//
//  // Initialize Hessian matrix H (of ℓ(θ)) as 4x4
//  arma::mat H(4,4, arma::fill::zeros);
//
//  // --- CONSTANT TERMS (do not depend on x) ---
//  // L1: n ln(α)  => d²/dα² = -n/α²
//  H(0,0) += -n/(alpha*alpha);
//  // L2: n ln(β)  => d²/dβ² = -n/β²
//  H(1,1) += -n/(beta*beta);
//  // L3: -n ln[B(γ, δ+1)]
//  //   d²/dγ² = -n [ψ₁(γ) - ψ₁(γ+δ+1)]  where ψ₁ is the trigamma function
//  H(2,2) += -n * (R::trigamma(gamma) - R::trigamma(gamma+delta+1));
//  //   d²/dδ² = -n [ψ₁(δ+1) - ψ₁(γ+δ+1)]
//  H(3,3) += -n * (R::trigamma(delta+1) - R::trigamma(gamma+delta+1));
//  //   Mixed derivative (γ,δ): = n ψ₁(γ+δ+1)
//  H(2,3) += n * R::trigamma(gamma+delta+1);
//  H(3,2) = H(2,3);
//
//  // --- TERMS THAT INVOLVE THE OBSERVATIONS ---
//  // Loop over each observation to accumulate contributions from:
//  // L4: (α-1) Σ ln(x_i)  --> contributes only to first derivatives
//  // L5: (β(δ+1)-1) Σ ln(v), where v = 1 - x^α
//  // L6: (γ-1) Σ ln(w), where w = 1 - v^β
//  for (int i = 0; i < n; i++) {
//    double xi    = x(i);
//    double ln_xi = std::log(xi);
//
//    // -- Compute A = x^α and its derivatives --
//    double A = std::pow(xi, alpha);                  // A = x^α
//    double dA_dalpha = A * ln_xi;                    // dA/dα = x^α ln(x)
//    double d2A_dalpha2 = A * ln_xi * ln_xi;          // d²A/dα² = x^α (ln(x))²
//
//    // -- v = 1 - A and its derivatives --
//    double v = 1.0 - A;                              // v = 1 - x^α
//    double ln_v = std::log(v);                       // ln(v)
//    double dv_dalpha = -dA_dalpha;                   // dv/dα = -dA/dα = -x^α ln(x)
//    double d2v_dalpha2 = -d2A_dalpha2;               // d²v/dα² = -d²A/dα² = -x^α (ln(x))²
//
//    // --- L5: (β(δ+1)-1) ln(v) ---
//    double beta_delta_factor = beta * (delta + 1.0) - 1.0;
//    // First derivative w.r.t. α: (β(δ+1)-1) * (1/v)*dv_dalpha
//    double dL5_dalpha = beta_delta_factor * (dv_dalpha / v);
//    // Second derivative w.r.t. α: (β(δ+1)-1)*[(d²v/dα²*v - (dv/dα)²)/v²]
//    double d2L5_dalpha2 = beta_delta_factor * ((d2v_dalpha2 * v - dv_dalpha * dv_dalpha) / (v*v));
//    // First derivative w.r.t. β: dL5/dβ = (δ+1) * ln(v)
//    double dL5_dbeta = (delta + 1.0) * ln_v;
//    // First derivative w.r.t. δ: dL5/dδ = β * ln(v)
//    double dL5_ddelta = beta * ln_v;
//    // Mixed derivative: d²L5/(dα dβ) = d/dβ[(β(δ+1)-1)*(dv_dalpha/v)] = (δ+1)*(dv_dalpha/v)
//    double d2L5_dalpha_dbeta = (delta + 1.0) * (dv_dalpha / v);
//    // Mixed derivative: d²L5/(dα dδ) = d/dδ[(β(δ+1)-1)*(dv_dalpha/v)] = β*(dv_dalpha/v)
//    double d2L5_dalpha_ddelta = beta * (dv_dalpha / v);
//    // Mixed derivative: d²L5/(dβ dδ) = d/dδ[(δ+1)*ln(v)] = ln(v)
//    double d2L5_dbeta_ddelta = ln_v;
//
//    // --- L6: (γ - 1) ln(w), where w = 1 - v^β ---
//    double v_beta = std::pow(v, beta);              // v^β
//    double w = 1.0 - v_beta;                        // w = 1 - v^β
//    double ln_w = std::log(w);                      // ln(w)
//    // Derivative of w w.r.t. v: dw/dv = -β * v^(β-1)
//    double dw_dv = -beta * std::pow(v, beta - 1.0);
//    // Chain rule: dw/dα = dw/dv * dv/dα
//    double dw_dalpha = dw_dv * dv_dalpha;
//    // First derivative w.r.t. α: d/dα ln(w) = (1/w)*dw_dalpha
//    double dL6_dalpha = (gamma - 1.0) * (dw_dalpha / w);
//    // Second derivative w.r.t. α for L6:
//    // d²/dα² ln(w) = [d²w/dα² * w - (dw/dα)²] / w²
//    // Computing d²w/dα²:
//    //   dw/dα = -β * v^(β-1)*dv_dalpha,
//    //   d²w/dα² = -β * [(β-1)*v^(β-2)*(dv_dalpha)² + v^(β-1)*d²v_dalpha²]
//    double d2w_dalpha2 = -beta * ((beta - 1.0) * std::pow(v, beta-2.0) * (dv_dalpha * dv_dalpha)
//                                    + std::pow(v, beta-1.0) * d2v_dalpha2);
//    double d2L6_dalpha2 = (gamma - 1.0) * ((d2w_dalpha2 * w - (dw_dalpha * dw_dalpha)) / (w*w));
//    // Derivative w.r.t. β: d/dβ ln(w). Note: d/dβ(v^β) = v^β ln(v) => d/dβ w = -v^β ln(v)
//    double dw_dbeta = -v_beta * ln_v;
//    double dL6_dbeta = (gamma - 1.0) * (dw_dbeta / w);
//    // Second derivative w.r.t. β for L6:
//    // d²/dβ² ln(w) = [d²w/dβ² * w - (dw/dβ)²]/w², where d²w/dβ² = -v^β (ln(v))²
//    double d2w_dbeta2 = -v_beta * (ln_v * ln_v);
//    double d2L6_dbeta2 = (gamma - 1.0) * ((d2w_dbeta2 * w - (dw_dbeta * dw_dbeta))/(w*w));
//    // Mixed derivative L6 (α,β): d²/(dα dβ) ln(w) =
//    //   = d/dβ[(dw_dalpha)/w] = (d/dβ dw_dalpha)/w - (dw_dalpha*dw_dbeta)/(w*w)
//    // Approximate d/dβ dw_dalpha:
//    double d_dw_dalpha_dbeta = -std::pow(v, beta-1.0) * (1.0 + beta * ln_v) * dv_dalpha;
//    double d2L6_dalpha_dbeta = (gamma - 1.0) * ((d_dw_dalpha_dbeta / w) - (dw_dalpha * dw_dbeta)/(w*w));
//    // First derivative w.r.t. γ: dL6/dγ = ln(w)
//    double dL6_dgamma = ln_w;
//    // Mixed derivatives with γ
//    double d2L6_dalpha_dgamma = dw_dalpha / w;
//    double d2L6_dbeta_dgamma = dw_dbeta / w;
//
//    // --- ACCUMULATING CONTRIBUTIONS TO THE HESSIAN MATRIX ---
//    // Index: 0 = α, 1 = β, 2 = γ, 3 = δ
//
//    // H(α,α): sum of L1, L5, and L6 (constants already added)
//    H(0,0) += d2L5_dalpha2 + d2L6_dalpha2;
//
//    // H(β,β): contributions from L2, L5, and L6
//    H(1,1) += d2L6_dbeta2;
//
//    // H(α,β): mixed from L5 and L6
//    H(0,1) += d2L5_dalpha_dbeta + d2L6_dalpha_dbeta;
//    H(1,0) = H(0,1);
//
//    // H(α,γ): mixed from L6
//    H(0,2) += d2L6_dalpha_dgamma;
//    H(2,0) = H(0,2);
//
//    // H(α,δ): mixed from L5
//    H(0,3) += d2L5_dalpha_ddelta;
//    H(3,0) = H(0,3);
//
//    // H(β,γ): mixed from L6
//    H(1,2) += d2L6_dbeta_dgamma;
//    H(2,1) = H(1,2);
//
//    // H(β,δ): mixed from L5
//    H(1,3) += d2L5_dbeta_ddelta;
//    H(3,1) = H(1,3);
//
//  } // end of loop
//
//  // Returns the analytic Hessian matrix of the negative log-likelihood
//  return Rcpp::wrap(-H);
// }




/*
----------------------------------------------------------------------------
NUMERIC STABILITY FUNCTIONS AND CHECKS
----------------------------------------------------------------------------
NOTE: We assume the following inline functions already exist in the compilation
environment (similarly to gkwdist.cpp), so we do NOT redefine them here:
- log1mexp(double)
- log1pexp(double)
- safe_log(double)
- safe_exp(double)
- safe_pow(double, double)
etc.

Also, we define here a small parameter checker for the "Exponentiated Kumaraswamy"
(EKw) family, interpreted as a sub-model of GKw with gamma=1 and delta=0,
leaving three free parameters: (alpha, beta, lambda).
*/

/*
----------------------------------------------------------------------------
EXPONENTIATED KUMARASWAMY (EKw) DISTRIBUTION
----------------------------------------------------------------------------

We interpret EKw(α, β, λ) as the GKw distribution with gamma=1 and delta=0.

* PDF:
f(x) = λ * α * β * x^(α-1) * (1 - x^α)^(β - 1) *
[1 - (1 - x^α)^β ]^(λ - 1),    for 0 < x < 1.

* CDF:
F(x) = [1 - (1 - x^α)^β ]^λ,         for 0 < x < 1.

* QUANTILE:
If p = F(x) = [1 - (1 - x^α)^β]^λ, then
p^(1/λ) = 1 - (1 - x^α)^β
(1 - x^α)^β = 1 - p^(1/λ)
x^α = 1 - [1 - p^(1/λ)]^(1/β)
x = {1 - [1 - p^(1/λ)]^(1/β)}^(1/α).

* RNG:
We can generate via the quantile method: U ~ Uniform(0,1), X= Q(U).

X = Q(U) = {1 - [1 - U^(1/λ)]^(1/β)}^(1/α).

* LOG-LIKELIHOOD:
The log-density for observation x in (0,1):
log f(x) = log(λ) + log(α) + log(β)
+ (α-1)*log(x)
+ (β-1)*log(1 - x^α)
+ (λ-1)*log(1 - (1 - x^α)^β).

Summation of log-likelihood over all x. We return negative of that for 'llekw'.
*/


// -----------------------------------------------------------------------------
// Parameter checker for EKw distribution
// EKw(α, β, λ):  all must be > 0
// -----------------------------------------------------------------------------
inline bool check_ekw_pars(double alpha, double beta, double lambda, bool strict=false) {
if (alpha <= 0.0 || beta <= 0.0 || lambda <= 0.0) {
  return false;
}
if (strict) {
  const double MINP = 1e-6;
  const double MAXP = 1e6;
  if (alpha < MINP || beta < MINP || lambda < MINP)  return false;
  if (alpha > MAXP || beta > MAXP || lambda > MAXP)  return false;
}
return true;
}

// -----------------------------------------------------------------------------
// 1) dekw: PDF of Exponentiated Kumaraswamy
// -----------------------------------------------------------------------------

//' @title Density for Exponentiated Kumaraswamy Distribution
//'
//' @description
//' Computes the PDF of EKw(α, β, λ), where \eqn{\alpha, β, λ > 0}.
//'
//' @param x Vector of quantiles in (0,1).
//' @param alpha Shape parameter \eqn{\alpha > 0}.
//' @param beta Shape parameter \eqn{\beta > 0}.
//' @param lambda Shape parameter \eqn{\lambda > 0}.
//' @param log_prob Logical; if TRUE, returns log-density, else density.
//'
//' @details
//' The PDF is
//' \deqn{
//'   f(x) = \lambda\,\alpha\,\beta \; x^{\alpha-1} (1 - x^α)^{\beta-1} \;
//'          [1 - (1 - x^α)^\beta ]^{\lambda - 1}, \quad 0<x<1.
//' }
//' This follows from the GKw distribution with \eqn{\gamma=1,\;\delta=0}.
//'
//' @return A vector of the same length as the broadcast of (x, alpha, beta, lambda).
//'
//' @export
// [[Rcpp::export]]
arma::vec dekw(
   const arma::vec& x,
   const Rcpp::NumericVector& alpha,
   const Rcpp::NumericVector& beta,
   const Rcpp::NumericVector& lambda,
   bool log_prob = false
) {
 arma::vec a_vec(alpha.begin(), alpha.size());
 arma::vec b_vec(beta.begin(), beta.size());
 arma::vec l_vec(lambda.begin(), lambda.size());

 size_t N = std::max({ x.n_elem, a_vec.n_elem, b_vec.n_elem, l_vec.n_elem });
 arma::vec out(N);
 out.fill(log_prob ? R_NegInf : 0.0);

 for (size_t i=0; i<N; i++) {
   double a = a_vec[i % a_vec.n_elem];
   double b = b_vec[i % b_vec.n_elem];
   double l = l_vec[i % l_vec.n_elem];
   double xx = x[i % x.n_elem];

   if (!check_ekw_pars(a, b, l)) {
     // invalid => PDF=0 or logPDF=-Inf
     continue;
   }
   // domain check
   if (xx <= 0.0 || xx >= 1.0 || !R_finite(xx)) {
     continue;
   }

   // log f(x) = log(lambda) + log(a) + log(b) + (a-1)*log(x)
   //            + (b-1)*log(1 - x^a)
   //            + (lambda-1)*log(1 - (1 - x^a)^b)
   double ll  = std::log(l);
   double la  = std::log(a);
   double lb  = std::log(b);
   double lx  = std::log(xx);

   double xalpha = a*lx; // log(x^a)
   double log_1_xalpha = log1mexp(xalpha); // log(1 - x^a)
   if (!R_finite(log_1_xalpha)) {
     continue;
   }

   double term2 = (b - 1.0)*log_1_xalpha; // (b-1)* log(1 - x^a)

   // let A= (1 - x^a)^b => logA= b*log_1_xalpha
   double logA = b*log_1_xalpha;
   double log_1_minus_A = log1mexp(logA); // log(1 - A)
   if (!R_finite(log_1_minus_A)) {
     continue;
   }
   double term3 = (l - 1.0)* log_1_minus_A;

   double log_pdf = ll + la + lb
   + (a - 1.0)* lx
   + term2
   + term3;

   if (log_prob) {
     out(i)= log_pdf;
   } else {
     out(i)= std::exp(log_pdf);
   }
 }

 return out;
}

// -----------------------------------------------------------------------------
// 2) pekw: CDF of Exponentiated Kumaraswamy
// -----------------------------------------------------------------------------

//' @title CDF for Exponentiated Kumaraswamy Distribution
//'
//' @description
//' Computes F(x)= P(X ≤ x) for EKw(α, β, λ).
//'
//' @param q Vector of quantiles in (0,1).
//' @param alpha Shape parameter \eqn{\alpha > 0}.
//' @param beta Shape parameter \eqn{\beta > 0}.
//' @param lambda Shape parameter \eqn{\lambda > 0}.
//' @param lower_tail Logical; if TRUE, returns F(q); else 1-F(q).
//' @param log_p Logical; if TRUE, returns log(F(q)) or log(1-F(q)).
//'
//' @details
//' The CDF is
//' \deqn{
//'   F(x)= [1 - (1 - x^α)^β ]^λ, \quad 0<x<1.
//' }
//'
//' @return A vector of probabilities of the same length as the broadcast of inputs.
//'
//' @export
// [[Rcpp::export]]
arma::vec pekw(
   const arma::vec& q,
   const Rcpp::NumericVector& alpha,
   const Rcpp::NumericVector& beta,
   const Rcpp::NumericVector& lambda,
   bool lower_tail = true,
   bool log_p = false
) {
 arma::vec a_vec(alpha.begin(), alpha.size());
 arma::vec b_vec(beta.begin(), beta.size());
 arma::vec l_vec(lambda.begin(), lambda.size());

 size_t N = std::max({ q.n_elem, a_vec.n_elem, b_vec.n_elem, l_vec.n_elem });
 arma::vec out(N);

 for (size_t i=0; i<N; i++) {
   double a = a_vec[i % a_vec.n_elem];
   double b = b_vec[i % b_vec.n_elem];
   double l = l_vec[i % l_vec.n_elem];
   double xx = q[i % q.n_elem];

   if (!check_ekw_pars(a, b, l)) {
     out(i)= NA_REAL;
     continue;
   }

   // boundary
   if (!R_finite(xx) || xx <= 0.0) {
     double val0 = (lower_tail ? 0.0 : 1.0);
     out(i) = (log_p ? std::log(val0) : val0);
     continue;
   }
   if (xx >= 1.0) {
     double val1 = (lower_tail ? 1.0 : 0.0);
     out(i) = (log_p ? std::log(val1) : val1);
     continue;
   }

   // F(x)= [1 - (1 - x^a)^b]^lambda
   double lx = std::log(xx);
   double xalpha = std::exp(a*lx);
   double omx = 1.0 - xalpha;         // (1 - x^α)
   if (omx <= 0.0) {
     // => F=1
     double val1 = (lower_tail ? 1.0 : 0.0);
     out(i) = (log_p ? std::log(val1) : val1);
     continue;
   }
   double t = 1.0 - std::pow(omx, b);
   if (t <= 0.0) {
     // => F=0
     double val0 = (lower_tail ? 0.0 : 1.0);
     out(i) = (log_p ? std::log(val0) : val0);
     continue;
   }
   if (t >= 1.0) {
     // => F=1
     double val1 = (lower_tail ? 1.0 : 0.0);
     out(i) = (log_p ? std::log(val1) : val1);
     continue;
   }
   double val = std::pow(t, l);
   // F(x)=val => if not lower tail => 1-val
   if (!lower_tail) {
     val = 1.0 - val;
   }
   if (log_p) {
     val = std::log(val);
   }
   out(i) = val;
 }

 return out;
}


// -----------------------------------------------------------------------------
// 3) qekw: Quantile of Exponentiated Kumaraswamy
// -----------------------------------------------------------------------------

//' @title Quantile Function for Exponentiated Kumaraswamy
//'
//' @description
//' Computes the quantiles x = Q(p) for EKw(α, β, λ).
//'
//' @param p Vector of probabilities in (0,1).
//' @param alpha Shape parameter \eqn{\alpha > 0}.
//' @param beta Shape parameter \eqn{\beta > 0}.
//' @param lambda Shape parameter \eqn{\lambda > 0}.
//' @param lower_tail Logical; if TRUE, p is F(x). If FALSE, p=1-F(x).
//' @param log_p Logical; if TRUE, p is given as log(p).
//'
//' @details
//' The inverse of \eqn{F(x)= [1 - (1 - x^\alpha)^\beta ]^\lambda} is
//' \deqn{
//'   Q(p)= \Bigl\{1 - \bigl[1 - p^{1/\lambda}\bigr]^{1/\beta}\Bigr\}^{1/\alpha}.
//' }
//' if lower_tail=TRUE. Adjust accordingly if \code{lower_tail=FALSE} or \code{log_p=TRUE}.
//'
//' @export
// [[Rcpp::export]]
arma::vec qekw(
   const arma::vec& p,
   const Rcpp::NumericVector& alpha,
   const Rcpp::NumericVector& beta,
   const Rcpp::NumericVector& lambda,
   bool lower_tail = true,
   bool log_p = false
) {
 arma::vec a_vec(alpha.begin(), alpha.size());
 arma::vec b_vec(beta.begin(), beta.size());
 arma::vec l_vec(lambda.begin(), lambda.size());

 size_t N = std::max({ p.n_elem, a_vec.n_elem, b_vec.n_elem, l_vec.n_elem });
 arma::vec out(N);

 for (size_t i=0; i<N; i++){
   double a = a_vec[i % a_vec.n_elem];
   double b = b_vec[i % b_vec.n_elem];
   double l = l_vec[i % l_vec.n_elem];
   double pp = p[i % p.n_elem];

   if (!check_ekw_pars(a, b, l)) {
     out(i) = NA_REAL;
     continue;
   }

   // handle log_p
   if (log_p) {
     if (pp > 0.0) {
       // log(p)>0 => p>1 => invalid
       out(i) = NA_REAL;
       continue;
     }
     pp = std::exp(pp);
   }
   // handle tail
   if (!lower_tail) {
     pp = 1.0 - pp;
   }

   // boundaries
   if (pp <= 0.0) {
     out(i) = 0.0;
     continue;
   }
   if (pp >= 1.0) {
     out(i) = 1.0;
     continue;
   }

   // Q(p)= {1 - [1 - p^(1/λ)]^(1/β)}^(1/α)
   double step1 = std::pow(pp, 1.0/l);          // p^(1/λ)
   double step2 = 1.0 - step1;                  // 1 - p^(1/λ)
   if (step2 < 0.0) step2 = 0.0;
   double step3 = std::pow(step2, 1.0/b);       // [1 - p^(1/λ)]^(1/β)
   double step4 = 1.0 - step3;                  // 1 - ...
   if (step4 < 0.0) step4 = 0.0;

   double x;
   if (a == 1.0) {
     x = step4;
   } else {
     x = std::pow(step4, 1.0/a);
     if (x < 0.0) x = 0.0;
     if (x > 1.0) x = 1.0;
   }

   out(i) = x;
 }

 return out;
}


// -----------------------------------------------------------------------------
// 4) rekw: RNG for Exponentiated Kumaraswamy
// -----------------------------------------------------------------------------

//' @title Random Generation for Exponentiated Kumaraswamy Distribution
//'
//' @description
//' Generates n samples from EKw(α, β, λ).
//'
//' @param n Integer number of observations.
//' @param alpha Shape parameter \eqn{\alpha > 0}.
//' @param beta Shape parameter \eqn{\beta > 0}.
//' @param lambda Shape parameter \eqn{\lambda > 0}.
//'
//' @details
//' Implementation via the quantile method:  X = Q(U),  U ~ Uniform(0,1).
//' The quantile is
//' \deqn{
//'   Q(u)= \Bigl\{1 - \bigl[1 - u^{1/\lambda}\bigr]^{1/\beta}\Bigr\}^{1/\alpha}.
//' }
//'
//' @export
// [[Rcpp::export]]
arma::vec rekw(
   int n,
   const Rcpp::NumericVector& alpha,
   const Rcpp::NumericVector& beta,
   const Rcpp::NumericVector& lambda
) {
 if (n <= 0) {
   Rcpp::stop("rekw: n must be positive");
 }

 arma::vec a_vec(alpha.begin(), alpha.size());
 arma::vec b_vec(beta.begin(), beta.size());
 arma::vec l_vec(lambda.begin(), lambda.size());

 size_t k = std::max({ a_vec.n_elem, b_vec.n_elem, l_vec.n_elem });
 arma::vec out(n);

 for (int i=0; i<n; i++){
   size_t idx = i % k;
   double a = a_vec[idx % a_vec.n_elem];
   double b = b_vec[idx % b_vec.n_elem];
   double l = l_vec[idx % l_vec.n_elem];

   if (!check_ekw_pars(a, b, l)) {
     out(i) = NA_REAL;
     Rcpp::warning("rekw: invalid parameters at index %d", i+1);
     continue;
   }

   double U = R::runif(0.0, 1.0);
   // X = Q(U)
   double step1 = std::pow(U, 1.0/l);
   double step2 = 1.0 - step1;
   if (step2 < 0.0) step2 = 0.0;
   double step3 = std::pow(step2, 1.0/b);
   double step4 = 1.0 - step3;
   if (step4 < 0.0) step4 = 0.0;

   double x;
   if (a == 1.0) {
     x = step4;
   } else {
     x = std::pow(step4, 1.0/a);
     if (!R_finite(x) || x < 0.0) x = 0.0;
     if (x > 1.0) x = 1.0;
   }

   out(i) = x;
 }

 return out;
}


// -----------------------------------------------------------------------------
// 5) llekw: Negative Log-Likelihood of EKw
// -----------------------------------------------------------------------------

//' @title Negative Log-Likelihood for Exponentiated Kumaraswamy
//'
//' @description
//' Computes the negative log-likelihood of EKw(α, β, λ), for data in (0,1).
//'
//' @param par NumericVector of length 3, (α, β, λ).
//' @param data NumericVector of observations, each in (0,1).
//'
//' @details
//' The PDF is
//' \deqn{
//'  f(x)= \lambda\,\alpha\,\beta \, x^{\alpha-1}\,(1-x^\alpha)^{\beta-1}\,[1-(1-x^\alpha)^\beta]^{\lambda-1}.
//' }
//' The log-likelihood is the sum of log(f(x_i)). We return the negative log-likelihood.
//'
//' @return A double value: the negative log-likelihood. \code{Inf} if invalid.
//'
//' @examples
//' \dontrun{
//' # Generate sample data from an EKw distribution
//' set.seed(123)
//' x <- rekw(100, 2, 3, 0.5)
//' hist(x, breaks = 20, main = "EKw(2, 3, 0.5) Sample")
//'
//' # Use in optimization with Hessian-based methods
//' result <- optim(c(0.5, 0.5, 0.5), llekw, method = "BFGS",
//'                 hessian = TRUE, data = x)
//'
//' # Compare numerical and analytical derivatives
//' num_grad <- numDeriv::grad(llekw, x = result$par, data = x)
//' num_hess <- numDeriv::hessian(llekw, x = result$par, data = x)
//'
//' ana_grad <- grekw(result$par, data = x)
//' ana_hess <- hsekw(result$par, data = x)
//'
//' # Check differences (should be very small)
//' round(num_grad - ana_grad, 4)
//' round(num_hess - ana_hess, 4)
//' }
//'
//' @export
// [[Rcpp::export]]
double llekw(const Rcpp::NumericVector& par,
            const Rcpp::NumericVector& data) {
 if (par.size() < 3) {
   return R_PosInf; // no param
 }
 double a = par[0];
 double b = par[1];
 double l = par[2];

 if (!check_ekw_pars(a, b, l)) {
   return R_PosInf; // invalid => +Inf
 }

 arma::vec x = Rcpp::as<arma::vec>(data);
 if (x.n_elem < 1) {
   return R_PosInf;
 }
 // check data
 if (arma::any(x <= 0.0) || arma::any(x >= 1.0)) {
   return R_PosInf;
 }

 int n = x.n_elem;
 // constant part: n* [ log(l) + log(a) + log(b) ]
 double cst = n*( std::log(l) + std::log(a) + std::log(b) );

 // sum( (a-1)* log(x_i ) )
 arma::vec lx = arma::log(x);
 double sum1 = (a-1.0)* arma::sum(lx);

 // sum( (b-1)* log(1- x^a ) )
 arma::vec xalpha = arma::pow(x, a);
 arma::vec log_1_xalpha = arma::log(1.0 - xalpha);
 double sum2 = (b-1.0)* arma::sum(log_1_xalpha);

 // sum( (l-1)* log(1- (1- x^a)^b ) )
 arma::vec vbeta = arma::pow((1.0 - xalpha), b);
 arma::vec one_minus_vbeta = 1.0 - vbeta;
 arma::vec log_1_mv = arma::log(one_minus_vbeta);
 double sum3 = (l-1.0)* arma::sum(log_1_mv);

 double loglike = cst + sum1 + sum2 + sum3;
 // negative
 return -loglike;
}

// -----------------------------------------------------------------------------
// 6) grekw: Gradient of Negative Log-Likelihood of EKw
// -----------------------------------------------------------------------------

//' @title Gradient Function for Exponentiated Kumaraswamy Log-Likelihood
//'
//' @description
//' Calculates the gradient vector (partial derivatives) of the negative log-likelihood
//' function for the Exponentiated Kumaraswamy (EKw) distribution. This function provides
//' the exact gradient needed for efficient optimization in maximum likelihood estimation.
//'
//' @param par NumericVector of length 3 containing parameters (α, β, λ) in that order.
//'        All parameters must be positive.
//' @param data NumericVector of observations, where all values must be in the open interval (0,1).
//'
//' @return NumericVector of length 3 containing the gradient components (partial derivatives)
//'         of the negative log-likelihood with respect to each parameter (α, β, λ).
//'         Returns a vector of NaN values if any parameters or data values are invalid.
//'
//' @details
//' The gradient vector contains the following partial derivatives of the negative log-likelihood:
//'
//' \deqn{
//' \frac{\partial \ell}{\partial \alpha} = \frac{n}{\alpha} + \sum_{i=1}^{n}\log(x_i) -
//' \sum_{i=1}^{n}\left[x_i^{\alpha} \log(x_i) \left(\frac{\beta-1}{v_i} -
//' \frac{(\lambda-1) \beta v_i^{\beta-1}}{w_i}\right)\right]
//' }
//'
//' \deqn{
//' \frac{\partial \ell}{\partial \beta} = \frac{n}{\beta} + \sum_{i=1}^{n}\log(v_i) -
//' \sum_{i=1}^{n}\left[v_i^{\beta} \log(v_i) \frac{\lambda-1}{w_i}\right]
//' }
//'
//' \deqn{
//' \frac{\partial \ell}{\partial \lambda} = \frac{n}{\lambda} +
//' \sum_{i=1}^{n}\log(w_i)
//' }
//'
//' where:
//' \itemize{
//'   \item \deqn{v_i = 1 - x_i^{\alpha}}
//'   \item \deqn{w_i = 1 - v_i^{\beta} = 1 - (1-x_i^{\alpha})^{\beta}}
//' }
//' @examples
//' \dontrun{
//' # Generate sample data from an EKw distribution
//' set.seed(123)
//' x <- rekw(100, 2, 3, 0.5)
//' hist(x, breaks = 20, main = "EKw(2, 3, 0.5) Sample")
//'
//' # Use in optimization with Hessian-based methods
//' result <- optim(c(0.5, 0.5, 0.5), llekw, method = "BFGS",
//'                 hessian = TRUE, data = x)
//'
//' # Compare numerical and analytical derivatives
//' num_grad <- numDeriv::grad(llekw, x = result$par, data = x)
//' num_hess <- numDeriv::hessian(llekw, x = result$par, data = x)
//'
//' ana_grad <- grekw(result$par, data = x)
//' ana_hess <- hsekw(result$par, data = x)
//'
//' # Check differences (should be very small)
//' round(num_grad - ana_grad, 4)
//' round(num_hess - ana_hess, 4)
//' }
//'
//' @references
//' Kumaraswamy, P. (1980). A generalized probability density function for double-bounded random processes.
//' Journal of Hydrology, 46(1-2), 79-88.
//'
//' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized distributions.
//' Journal of Statistical Computation and Simulation, 81(7), 883-898.
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector grekw(const Rcpp::NumericVector& par, const Rcpp::NumericVector& data) {
 // Parameter extraction
 double alpha = par[0];   // Shape parameter α > 0
 double beta = par[1];    // Shape parameter β > 0
 double lambda = par[2];  // Shape parameter λ > 0

 // Parameter validation
 if (alpha <= 0 || beta <= 0 || lambda <= 0) {
   Rcpp::NumericVector grad(3, R_NaN);
   return grad;
 }

 // Data conversion and validation
 arma::vec x = Rcpp::as<arma::vec>(data);

 if (arma::any(x <= 0) || arma::any(x >= 1)) {
   Rcpp::NumericVector grad(3, R_NaN);
   return grad;
 }

 int n = x.n_elem;  // Sample size

 // Initialize gradient vector
 Rcpp::NumericVector grad(3, 0.0);

 // Small constant to avoid numerical issues
 double eps = std::numeric_limits<double>::epsilon() * 100;

 // Compute transformations and intermediate values
 arma::vec log_x = arma::log(x);                // log(x_i)
 arma::vec x_alpha = arma::pow(x, alpha);       // x_i^α
 arma::vec x_alpha_log_x = x_alpha % log_x;     // x_i^α * log(x_i)

 // v_i = 1 - x_i^α
 arma::vec v = 1.0 - x_alpha;
 v = arma::clamp(v, eps, 1.0 - eps);            // Prevent numerical issues

 arma::vec log_v = arma::log(v);                // log(v_i)
 arma::vec v_beta_m1 = arma::pow(v, beta - 1.0); // v_i^(β-1)
 arma::vec v_beta = arma::pow(v, beta);          // v_i^β
 arma::vec v_beta_log_v = v_beta % log_v;        // v_i^β * log(v_i)

 // w_i = 1 - v_i^β = 1 - (1-x_i^α)^β
 arma::vec w = 1.0 - v_beta;
 w = arma::clamp(w, eps, 1.0 - eps);            // Prevent numerical issues

 arma::vec log_w = arma::log(w);                // log(w_i)

 // Calculate partial derivatives for each parameter (for log-likelihood)

 // ∂ℓ/∂α = n/α + Σᵢlog(xᵢ) - Σᵢ[xᵢ^α * log(xᵢ) * ((β-1)/vᵢ - (λ-1) * β * vᵢ^(β-1) / wᵢ)]
 double d_alpha = n / alpha + arma::sum(log_x);

 // Calculate the complex term in the α gradient
 arma::vec alpha_term1 = (beta - 1.0) / v;                       // (β-1)/v_i
 arma::vec alpha_term2 = (lambda - 1.0) * beta * v_beta_m1 / w;  // (λ-1) * β * v_i^(β-1) / w_i

 d_alpha -= arma::sum(x_alpha_log_x % (alpha_term1 - alpha_term2));

 // ∂ℓ/∂β = n/β + Σᵢlog(vᵢ) - Σᵢ[vᵢ^β * log(vᵢ) * ((λ-1) / wᵢ)]
 double d_beta = n / beta + arma::sum(log_v);

 // Calculate the term in the β gradient
 arma::vec beta_term = (lambda - 1.0) / w;       // (λ-1) / w_i

 d_beta -= arma::sum(v_beta_log_v % beta_term);

 // ∂ℓ/∂λ = n/λ + Σᵢlog(wᵢ)
 double d_lambda = n / lambda + arma::sum(log_w);

 // Since we're optimizing negative log-likelihood, negate all derivatives
 grad[0] = -d_alpha;
 grad[1] = -d_beta;
 grad[2] = -d_lambda;

 return grad;
}


//' @title Analytic Hessian Matrix for Exponentiated Kumaraswamy Distribution
//'
//' @description
//' Computes the analytic Hessian matrix of the log-likelihood function for
//' the Exponentiated Kumaraswamy (EKw) distribution. This function provides
//' exact second derivatives needed for optimization and inference.
//'
//' @param par Numeric vector of length 3 containing the parameters
//'        (α, β, λ) in that order. All parameters must be positive.
//' @param data Numeric vector of observations, where all values must be
//'        in the open interval (0,1).
//'
//' @return A 3×3 numeric matrix representing the Hessian of the negative
//'         log-likelihood function. If parameters or data are invalid
//'         (parameters ≤ 0 or data outside (0,1)), returns a matrix of
//'         NaN values.
//'
//' @details
//' The log-likelihood for the Exponentiated Kumaraswamy distribution is:
//'
//' \deqn{
//' \ell(\theta) = n \ln(\lambda) + n \ln(\alpha) + n \ln(\beta)
//' + (\alpha-1) \sum \ln(x_i)
//' + (\beta-1) \sum \ln(1 - x_i^\alpha)
//' + (\lambda-1) \sum \ln\{1 - (1 - x_i^\alpha)^\beta\}
//' }
//'
//' The implementation computes all second derivatives analytically for each term.
//' For computational efficiency, the following transformations are used:
//' \itemize{
//'   \item \deqn{A = x^α} and derivatives
//'   \item \deqn{v = 1 - A}
//'   \item \deqn{w = 1 - v^β}
//' }
//'
//' The returned Hessian matrix has the following structure:
//' \itemize{
//'   \item Rows/columns 1-3 correspond to α, β, λ respectively
//'   \item The matrix is symmetric (as expected for a Hessian)
//'   \item The matrix represents second derivatives of the negative log-likelihood
//' }
//'
//' This function is implemented in C++ for computational efficiency.
//'
//' @examples
//' \dontrun{
//' # Generate sample data from an EKw distribution
//' set.seed(123)
//' x <- rekw(100, 2, 3, 0.5)
//' hist(x, breaks = 20, main = "EKw(2, 3, 0.5) Sample")
//'
//' # Use in optimization with Hessian-based methods
//' result <- optim(c(0.5, 0.5, 0.5), llekw, method = "BFGS",
//'                 hessian = TRUE, data = x)
//'
//' # Compare numerical and analytical derivatives
//' num_grad <- numDeriv::grad(llekw, x = result$par, data = x)
//' num_hess <- numDeriv::hessian(llekw, x = result$par, data = x)
//'
//' ana_grad <- grekw(result$par, data = x)
//' ana_hess <- hsekw(result$par, data = x)
//'
//' # Check differences (should be very small)
//' round(num_grad - ana_grad, 4)
//' round(num_hess - ana_hess, 4)
//'
//' }
//'
//'
//' @references
//' Kumaraswamy, P. (1980). A generalized probability density function for double-bounded random processes.
//' Journal of Hydrology, 46(1-2), 79-88.
//'
//' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized distributions.
//' Journal of Statistical Computation and Simulation, 81(7), 883-898.
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix hsekw(const Rcpp::NumericVector& par, const Rcpp::NumericVector& data) {
 // Parameter extraction - EKw has only 3 parameters
 double alpha  = par[0];   // θ[0] = α
 double beta   = par[1];   // θ[1] = β
 double lambda = par[2];   // θ[2] = λ

 // Simple parameter validation (all > 0)
 if(alpha <= 0 || beta <= 0 || lambda <= 0) {
   Rcpp::NumericMatrix nanH(3,3);
   nanH.fill(R_NaN);
   return nanH;
 }

 // Data conversion and basic validation
 arma::vec x = Rcpp::as<arma::vec>(data);
 if(arma::any(x <= 0) || arma::any(x >= 1)) {
   Rcpp::NumericMatrix nanH(3,3);
   nanH.fill(R_NaN);
   return nanH;
 }

 int n = x.n_elem;  // sample size

 // Initialize Hessian matrix H (of ℓ(θ)) as 3x3
 arma::mat H(3,3, arma::fill::zeros);

 // --- CONSTANT TERMS (do not depend on x) ---
 // L1: n ln(λ)  => d²/dλ² = -n/λ²
 H(2,2) += - n/(lambda*lambda);
 // L2: n ln(α)  => d²/dα² = -n/α²
 H(0,0) += - n/(alpha*alpha);
 // L3: n ln(β)  => d²/dβ² = -n/β²
 H(1,1) += - n/(beta*beta);

 // --- TERMS THAT INVOLVE THE OBSERVATIONS ---
 // Loop over each observation to accumulate contributions from:
 // L4: (α-1) Σ ln(x_i) --> contributes only to first derivatives
 // L5: (β-1) Σ ln(v), where v = 1 - x^α
 // L6: (λ-1) Σ ln(w), where w = 1 - v^β
 for (int i = 0; i < n; i++) {
   double xi    = x(i);
   double ln_xi = std::log(xi);

   // -- Compute A = x^α and its derivatives --
   double A = std::pow(xi, alpha);                  // A = x^α
   double dA_dalpha = A * ln_xi;                    // dA/dα = x^α ln(x)
   double d2A_dalpha2 = A * ln_xi * ln_xi;          // d²A/dα² = x^α (ln(x))²

   // -- v = 1 - A and its derivatives --
   double v = 1.0 - A;                              // v = 1 - x^α
   double ln_v = std::log(v);                       // ln(v)
   double dv_dalpha = -dA_dalpha;                   // dv/dα = -dA/dα = -x^α ln(x)
   double d2v_dalpha2 = -d2A_dalpha2;               // d²v/dα² = -d²A/dα² = -x^α (ln(x))²

   // --- L5: (β-1) ln(v) ---
   // Second derivative w.r.t. α: (β-1)*[(d²v/dα²*v - (dv/dα)²)/v²]
   double d2L5_dalpha2 = (beta - 1.0) * ((d2v_dalpha2 * v - dv_dalpha * dv_dalpha) / (v*v));
   // Mixed derivative: d²L5/(dα dβ) = d/dβ[(β-1)*(dv_dalpha/v)] = (dv_dalpha/v)
   double d2L5_dalpha_dbeta = dv_dalpha / v;

   // --- L6: (λ - 1) ln(w), where w = 1 - v^β ---
   double v_beta = std::pow(v, beta);               // v^β
   double w = 1.0 - v_beta;                         // w = 1 - v^β

   // Derivative of w w.r.t. v: dw/dv = -β * v^(β-1)
   double dw_dv = -beta * std::pow(v, beta - 1.0);
   // Chain rule: dw/dα = dw/dv * dv/dα
   double dw_dalpha = dw_dv * dv_dalpha;
   // Second derivative w.r.t. α for L6:
   // d²/dα² ln(w) = [d²w/dα² * w - (dw/dα)²] / w²
   // Computing d²w/dα²:
   //   dw/dα = -β * v^(β-1)*dv_dalpha,
   //   d²w/dα² = -β * [(β-1)*v^(β-2)*(dv_dalpha)² + v^(β-1)*d²v_dalpha²]
   double d2w_dalpha2 = -beta * ((beta - 1.0) * std::pow(v, beta-2.0) * (dv_dalpha * dv_dalpha)
                                   + std::pow(v, beta-1.0) * d2v_dalpha2);
   double d2L6_dalpha2 = (lambda - 1.0) * ((d2w_dalpha2 * w - (dw_dalpha * dw_dalpha)) / (w*w));
   // Derivative w.r.t. β: d/dβ ln(w). Note: d/dβ(v^β) = v^β ln(v) => d/dβ w = -v^β ln(v)
   double dw_dbeta = -v_beta * ln_v;
   // Second derivative w.r.t. β for L6:
   // d²/dβ² ln(w) = [d²w/dβ² * w - (dw/dβ)²]/w², where d²w/dβ² = -v^β (ln(v))²
   double d2w_dbeta2 = -v_beta * (ln_v * ln_v);
   double d2L6_dbeta2 = (lambda - 1.0) * ((d2w_dbeta2 * w - (dw_dbeta * dw_dbeta))/(w*w));
   // Mixed derivative L6 (α,β): d²/(dα dβ) ln(w) =
   //   = d/dβ[(dw_dalpha)/w] = (d/dβ dw_dalpha)/w - (dw_dalpha*dw_dbeta)/(w*w)
   // Approximate d/dβ dw_dalpha:
   double d_dw_dalpha_dbeta = -std::pow(v, beta-1.0) * (1.0 + beta * ln_v) * dv_dalpha;
   double d2L6_dalpha_dbeta = (lambda - 1.0) * ((d_dw_dalpha_dbeta / w) - (dw_dalpha * dw_dbeta)/(w*w));

   // Mixed derivatives with λ
   // (α,λ): d²/(dα dλ) [λ ln(w)] = d/dλ[(λ-1)(dw_dalpha/w)] = dw_dalpha/w
   double d2L6_dalpha_dlambda = dw_dalpha / w;

   // (β,λ): d²/(dβ dλ) [λ ln(w)] = d/dλ[(λ-1)(dw_dbeta/w)] = dw_dbeta/w
   double d2L6_dbeta_dlambda = dw_dbeta / w;

   // --- ACCUMULATING CONTRIBUTIONS TO THE HESSIAN MATRIX ---
   // Index: 0 = α, 1 = β, 2 = λ

   // H(α,α): sum of L2, L5, and L6 (constants already added)
   H(0,0) += d2L5_dalpha2 + d2L6_dalpha2;

   // H(α,β): mixed from L5 and L6
   H(0,1) += d2L5_dalpha_dbeta + d2L6_dalpha_dbeta;
   H(1,0) = H(0,1);

   // H(β,β): contributions from L3 and L6
   H(1,1) += d2L6_dbeta2;

   // H(α,λ): mixed derivative from L6
   H(0,2) += d2L6_dalpha_dlambda;
   H(2,0) = H(0,2);

   // H(β,λ): mixed derivative from L6
   H(1,2) += d2L6_dbeta_dlambda;
   H(2,1) = H(1,2);

 } // end of loop

 // Returns the analytic Hessian matrix of the negative log-likelihood
 return Rcpp::wrap(-H);
}



// //' @title Analytic Hessian Matrix for Exponentiated Kumaraswamy Distribution
// //'
// //' @description
// //' Computes the analytic Hessian matrix of the log-likelihood function for
// //' the Exponentiated Kumaraswamy (EKw) distribution. This function provides
// //' exact second derivatives needed for optimization and inference.
// //'
// //' @param par Numeric vector of length 3 containing the parameters
// //'        (α, β, λ) in that order. All parameters must be positive.
// //' @param data Numeric vector of observations, where all values must be
// //'        in the open interval (0,1).
// //'
// //' @return A 3×3 numeric matrix representing the Hessian of the negative
// //'         log-likelihood function. If parameters or data are invalid
// //'         (parameters ≤ 0 or data outside (0,1)), returns a matrix of
// //'         NaN values.
// //'
// //' @details
// //' The log-likelihood for the Exponentiated Kumaraswamy distribution is:
// //'
// //' \deqn{
// //' \ell(\theta) = n \ln(\lambda) + n \ln(\alpha) + n \ln(\beta)
// //' + (\alpha-1) \sum \ln(x_i)
// //' + (\beta-1) \sum \ln(1 - x_i^\alpha)
// //' + (\lambda-1) \sum \ln\{1 - (1 - x_i^\alpha)^\beta\}
// //' }
// //'
// //' The implementation computes all second derivatives analytically for each term.
// //' For computational efficiency, the following transformations are used:
// //' \itemize{
// //'   \item \deqn{A = x^α} and derivatives
// //'   \item \deqn{v = 1 - A}
// //'   \item \deqn{w = 1 - v^β}
// //' }
// //'
// //' The returned Hessian matrix has the following structure:
// //' \itemize{
// //'   \item Rows/columns 1-3 correspond to α, β, λ respectively
// //'   \item The matrix is symmetric (as expected for a Hessian)
// //'   \item The matrix represents second derivatives of the negative log-likelihood
// //' }
// //'
// //' This function is implemented in C++ for computational efficiency.
// //'
// //' @examples
// //' \dontrun{
// //' # Generate sample data from an EKw distribution
// //' set.seed(123)
// //' x <- rekw(100, 2, 3, 0.5)
// //' hist(x, breaks = 20, main = "EKw(2, 3, 0.5) Sample")
// //'
// //' # Use in optimization with Hessian-based methods
// //' result <- optim(c(0.5, 0.5, 0.5), llekw, method = "BFGS",
// //'                 hessian = TRUE, data = x)
// //'
// //' # Compare numerical and analytical derivatives
// //' num_grad <- numDeriv::grad(llekw, x = result$par, data = x)
// //' num_hess <- numDeriv::hessian(llekw, x = result$par, data = x)
// //'
// //' ana_grad <- grekw(result$par, data = x)
// //' ana_hess <- hsekw(result$par, data = x)
// //'
// //' # Check differences (should be very small)
// //' round(num_grad - ana_grad, 4)
// //' round(num_hess - ana_hess, 4)
// //'
// //' }
// //'
// //'
// //' @references
// //' Kumaraswamy, P. (1980). A generalized probability density function for double-bounded random processes.
// //' Journal of Hydrology, 46(1-2), 79-88.
// //'
// //' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized distributions.
// //' Journal of Statistical Computation and Simulation, 81(7), 883-898.
// //'
// //' @export
// // [[Rcpp::export]]
// Rcpp::NumericMatrix hsekw(const Rcpp::NumericVector& par, const Rcpp::NumericVector& data) {
//  // Parameter extraction - EKw has only 3 parameters
//  double alpha  = par[0];   // θ[0] = α
//  double beta   = par[1];   // θ[1] = β
//  double lambda = par[2];   // θ[2] = λ
//
//  // Simple parameter validation (all > 0)
//  if(alpha <= 0 || beta <= 0 || lambda <= 0) {
//    Rcpp::NumericMatrix nanH(3,3);
//    nanH.fill(R_NaN);
//    return nanH;
//  }
//
//  // Data conversion and basic validation
//  arma::vec x = Rcpp::as<arma::vec>(data);
//  if(arma::any(x <= 0) || arma::any(x >= 1)) {
//    Rcpp::NumericMatrix nanH(3,3);
//    nanH.fill(R_NaN);
//    return nanH;
//  }
//
//  int n = x.n_elem;  // sample size
//
//  // Initialize Hessian matrix H (of ℓ(θ)) as 3x3
//  arma::mat H(3,3, arma::fill::zeros);
//
//  // --- CONSTANT TERMS (do not depend on x) ---
//  // L1: n ln(λ)  => d²/dλ² = -n/λ²
//  H(2,2) += - n/(lambda*lambda);
//  // L2: n ln(α)  => d²/dα² = -n/α²
//  H(0,0) += - n/(alpha*alpha);
//  // L3: n ln(β)  => d²/dβ² = -n/β²
//  H(1,1) += - n/(beta*beta);
//
//  // --- TERMS THAT INVOLVE THE OBSERVATIONS ---
//  // Loop over each observation to accumulate contributions from:
//  // L4: (α-1) Σ ln(x_i) --> contributes only to first derivatives
//  // L5: (β-1) Σ ln(v), where v = 1 - x^α
//  // L6: (λ-1) Σ ln(w), where w = 1 - v^β
//  for (int i = 0; i < n; i++) {
//    double xi    = x(i);
//    double ln_xi = std::log(xi);
//
//    // -- Compute A = x^α and its derivatives --
//    double A = std::pow(xi, alpha);                  // A = x^α
//    double dA_dalpha = A * ln_xi;                    // dA/dα = x^α ln(x)
//    double d2A_dalpha2 = A * ln_xi * ln_xi;          // d²A/dα² = x^α (ln(x))²
//
//    // -- v = 1 - A and its derivatives --
//    double v = 1.0 - A;                              // v = 1 - x^α
//    double ln_v = std::log(v);                       // ln(v)
//    double dv_dalpha = -dA_dalpha;                   // dv/dα = -dA/dα = -x^α ln(x)
//    double d2v_dalpha2 = -d2A_dalpha2;               // d²v/dα² = -d²A/dα² = -x^α (ln(x))²
//
//    // --- L5: (β-1) ln(v) ---
//    // First derivative w.r.t. α: (β-1) * (1/v)*dv_dalpha
//    double dL5_dalpha = (beta - 1.0) * (dv_dalpha / v);
//    // Second derivative w.r.t. α: (β-1)*[(d²v/dα²*v - (dv/dα)²)/v²]
//    double d2L5_dalpha2 = (beta - 1.0) * ((d2v_dalpha2 * v - dv_dalpha * dv_dalpha) / (v*v));
//    // Mixed derivative: d²L5/(dα dβ) = d/dβ[(β-1)*(dv_dalpha/v)] = (dv_dalpha/v)
//    double d2L5_dalpha_dbeta = dv_dalpha / v;
//
//    // --- L6: (λ - 1) ln(w), where w = 1 - v^β ---
//    double v_beta = std::pow(v, beta);               // v^β
//    double w = 1.0 - v_beta;                         // w = 1 - v^β
//    double ln_w = std::log(w);                       // ln(w)
//    // Derivative of w w.r.t. v: dw/dv = -β * v^(β-1)
//    double dw_dv = -beta * std::pow(v, beta - 1.0);
//    // Chain rule: dw/dα = dw/dv * dv/dα
//    double dw_dalpha = dw_dv * dv_dalpha;
//    // First derivative w.r.t. α: d/dα ln(w) = (1/w)*dw_dalpha
//    double dL6_dalpha = (lambda - 1.0) * (dw_dalpha / w);
//    // Second derivative w.r.t. α for L6:
//    // d²/dα² ln(w) = [d²w/dα² * w - (dw/dα)²] / w²
//    // Computing d²w/dα²:
//    //   dw/dα = -β * v^(β-1)*dv_dalpha,
//    //   d²w/dα² = -β * [(β-1)*v^(β-2)*(dv_dalpha)² + v^(β-1)*d²v_dalpha²]
//    double d2w_dalpha2 = -beta * ((beta - 1.0) * std::pow(v, beta-2.0) * (dv_dalpha * dv_dalpha)
//                                    + std::pow(v, beta-1.0) * d2v_dalpha2);
//    double d2L6_dalpha2 = (lambda - 1.0) * ((d2w_dalpha2 * w - (dw_dalpha * dw_dalpha)) / (w*w));
//    // Derivative w.r.t. β: d/dβ ln(w). Note: d/dβ(v^β) = v^β ln(v) => d/dβ w = -v^β ln(v)
//    double dw_dbeta = -v_beta * ln_v;
//    double dL6_dbeta = (lambda - 1.0) * (dw_dbeta / w);
//    // Second derivative w.r.t. β for L6:
//    // d²/dβ² ln(w) = [d²w/dβ² * w - (dw/dβ)²]/w², where d²w/dβ² = -v^β (ln(v))²
//    double d2w_dbeta2 = -v_beta * (ln_v * ln_v);
//    double d2L6_dbeta2 = (lambda - 1.0) * ((d2w_dbeta2 * w - (dw_dbeta * dw_dbeta))/(w*w));
//    // Mixed derivative L6 (α,β): d²/(dα dβ) ln(w) =
//    //   = d/dβ[(dw_dalpha)/w] = (d/dβ dw_dalpha)/w - (dw_dalpha*dw_dbeta)/(w*w)
//    // Approximate d/dβ dw_dalpha:
//    double d_dw_dalpha_dbeta = -std::pow(v, beta-1.0) * (1.0 + beta * ln_v) * dv_dalpha;
//    double d2L6_dalpha_dbeta = (lambda - 1.0) * ((d_dw_dalpha_dbeta / w) - (dw_dalpha * dw_dbeta)/(w*w));
//
//    // Mixed derivatives with λ
//    // (α,λ): d²/(dα dλ) [λ ln(w)] = d/dλ[(λ-1)(dw_dalpha/w)] = dw_dalpha/w
//    double d2L6_dalpha_dlambda = dw_dalpha / w;
//
//    // (β,λ): d²/(dβ dλ) [λ ln(w)] = d/dλ[(λ-1)(dw_dbeta/w)] = dw_dbeta/w
//    double d2L6_dbeta_dlambda = dw_dbeta / w;
//
//    // --- ACCUMULATING CONTRIBUTIONS TO THE HESSIAN MATRIX ---
//    // Index: 0 = α, 1 = β, 2 = λ
//
//    // H(α,α): sum of L2, L5, and L6 (constants already added)
//    H(0,0) += d2L5_dalpha2 + d2L6_dalpha2;
//
//    // H(α,β): mixed from L5 and L6
//    H(0,1) += d2L5_dalpha_dbeta + d2L6_dalpha_dbeta;
//    H(1,0) = H(0,1);
//
//    // H(β,β): contributions from L3 and L6
//    H(1,1) += d2L6_dbeta2;
//
//    // H(α,λ): mixed derivative from L6
//    H(0,2) += d2L6_dalpha_dlambda;
//    H(2,0) = H(0,2);
//
//    // H(β,λ): mixed derivative from L6
//    H(1,2) += d2L6_dbeta_dlambda;
//    H(2,1) = H(1,2);
//
//  } // end of loop
//
//  // Returns the analytic Hessian matrix of the negative log-likelihood
//  return Rcpp::wrap(-H);
// }





/*
----------------------------------------------------------------------------
IMPORTANT NOTE:
We assume the same numeric stability functions from gkwdist.cpp or a "core" file:
- log1mexp(double)
- safe_log(double)
- safe_exp(double)
- safe_pow(double,double)
- etc.
They are NOT redefined here to avoid duplication.
*/

/*
----------------------------------------------------------------------------
BETA POWER (BP) DISTRIBUTION: BP(γ, δ, λ)
----------------------------------------------------------------------------

This arises from GKw with α=1 and β=1, leaving three parameters: (γ>0, δ≥0, λ>0).

* PDF:
f(x; γ, δ, λ) = [ λ / B(γ, δ+1) ] * x^(γλ - 1) * (1 - x^λ)^δ,   0<x<1.

* CDF:
F(x; γ, δ, λ) = I_{x^λ}(γ, δ+1) = pbeta(x^λ, γ, δ+1).

* QUANTILE:
Q(p; γ, δ, λ) = [ qbeta(p, γ, δ+1) ]^(1/λ).

* RNG:
If U ~ Beta(γ, δ+1), then X = U^(1/λ).

* NEGATIVE LOG-LIKELIHOOD:
sum( -log f(x_i) )
where
log f(x) = log(λ) - log B(γ, δ+1)
+ (γ λ -1)* log(x)
+ δ * log(1 - x^λ).

We'll define five functions:
- dmc() : PDF
- pmc() : CDF
- qmc() : quantile
- rmc() : random generator
- llmc(): negative log-likelihood

We'll also define a param-checker for (γ, δ, λ).
*/

// -----------------------------------------------------------------------------
// Parameter checker for Beta Power distribution
// BP(γ>0, δ≥0, λ>0)
// -----------------------------------------------------------------------------
inline bool check_bp_pars(double gamma, double delta, double lambda, bool strict = false) {
if (gamma <= 0.0 || delta < 0.0 || lambda <= 0.0) {
  return false;
}
if (strict) {
  const double MINP=1e-8;
  const double MAXP=1e8;
  if (gamma<MINP || lambda<MINP) return false;
  if (gamma>MAXP || delta>MAXP || lambda>MAXP) return false;
}
return true;
}

// -----------------------------------------------------------------------------
// 1) dmc: PDF of Beta Power McDonald
// -----------------------------------------------------------------------------

//' @title Density of the Beta Power Distribution
//'
//' @description
//' Computes the PDF of the Beta Power (BP or MC) distribution with parameters \eqn{γ, δ, λ}.
//'
//' @param x Vector of quantiles in (0,1).
//' @param gamma Shape parameter \eqn{γ > 0}.
//' @param delta Shape parameter \eqn{δ \ge 0}.
//' @param lambda Shape parameter \eqn{λ > 0}.
//' @param log_prob Logical; if TRUE, returns log-density.
//'
//' @details
//' The PDF is
//' \deqn{
//'   f(x) = \frac{\lambda}{B(\gamma,\delta+1)} \; x^{\gamma \lambda - 1} \; (1 - x^\lambda)^\delta, \quad 0<x<1.
//' }
//' with \eqn{B(\gamma,\delta+1)} the Beta function. This sub-family arises
//' from the Generalized Kumaraswamy (GKw) by setting \eqn{\alpha=1,\beta=1}.
//'
//' @return A vector of densities or log-densities, the same length as the broadcast of inputs.
//'
//' @export
// [[Rcpp::export]]
arma::vec dmc(
   const arma::vec& x,
   const Rcpp::NumericVector& gamma,
   const Rcpp::NumericVector& delta,
   const Rcpp::NumericVector& lambda,
   bool log_prob = false
) {
 arma::vec g_vec(gamma.begin(), gamma.size());
 arma::vec d_vec(delta.begin(), delta.size());
 arma::vec l_vec(lambda.begin(), lambda.size());

 size_t N= std::max({ x.n_elem, g_vec.n_elem, d_vec.n_elem, l_vec.n_elem });
 arma::vec out(N);

 // Pre-fill
 out.fill(log_prob ? R_NegInf : 0.0);

 for (size_t i=0; i<N; i++){
   double gg= g_vec[i % g_vec.n_elem];
   double dd= d_vec[i % d_vec.n_elem];
   double ll= l_vec[i % l_vec.n_elem];
   double xx= x[i % x.n_elem];

   if (!check_bp_pars(gg,dd,ll)) {
     // invalid => pdf=0 or logpdf=-Inf
     continue;
   }
   // domain
   if (xx<=0.0 || xx>=1.0 || !R_finite(xx)) {
     continue;
   }

   // log f(x)= log(λ) - log( B(γ, δ+1) )
   //           + (γλ -1)* log(x)
   //           + δ * log(1 - x^λ)
   double logB = R::lbeta(gg, dd+1.0);
   double logCst= std::log(ll) - logB;

   // (γ λ -1)* log(x)
   double exponent= gg*ll - 1.0;
   double lx= std::log(xx);
   double term1= exponent* lx;

   // δ * log(1 - x^λ)
   double x_pow_l= std::pow(xx, ll);
   if (x_pow_l>=1.0) {
     // => pdf=0
     continue;
   }
   double log_1_minus_xpow= std::log(1.0 - x_pow_l);
   double term2= dd * log_1_minus_xpow;

   double log_pdf= logCst + term1 + term2;
   if (log_prob) {
     out(i)= log_pdf;
   } else {
     out(i)= std::exp(log_pdf);
   }
 }

 return out;
}


// -----------------------------------------------------------------------------
// 2) pmc: CDF of Beta Power
// -----------------------------------------------------------------------------

//' @title CDF of the Beta Power Distribution
//'
//' @description
//' Computes \eqn{F(q)= P(X ≤ q) for BP(γ, δ, λ)}.
//'
//' @param q Vector of quantiles in (0,1).
//' @param gamma Shape parameter \eqn{γ > 0}.
//' @param delta Shape parameter \eqn{δ \ge 0}.
//' @param lambda Shape parameter \eqn{λ > 0}.
//' @param lower_tail Logical; if TRUE, returns F(q), else 1-F(q).
//' @param log_p Logical; if TRUE, returns log-probabilities.
//'
//' @details
//' \deqn{
//'   F(x)= I_{x^\lambda}(\gamma, \delta+1) \;=\; \mathrm{pbeta}(x^\lambda,\;\gamma,\;\delta+1).
//' }
//' The incomplete Beta approach is used for stable evaluation.
//'
//' @export
// [[Rcpp::export]]
arma::vec pmc(
   const arma::vec& q,
   const Rcpp::NumericVector& gamma,
   const Rcpp::NumericVector& delta,
   const Rcpp::NumericVector& lambda,
   bool lower_tail = true,
   bool log_p = false
) {
 arma::vec g_vec(gamma.begin(), gamma.size());
 arma::vec d_vec(delta.begin(), delta.size());
 arma::vec l_vec(lambda.begin(), lambda.size());

 size_t N= std::max({ q.n_elem, g_vec.n_elem, d_vec.n_elem, l_vec.n_elem });
 arma::vec out(N);

 for (size_t i=0; i<N; i++){
   double gg= g_vec[i % g_vec.n_elem];
   double dd= d_vec[i % d_vec.n_elem];
   double ll= l_vec[i % l_vec.n_elem];
   double xx= q[i % q.n_elem];

   if (!check_bp_pars(gg,dd,ll)) {
     out(i)= NA_REAL;
     continue;
   }

   // boundaries
   if (!R_finite(xx) || xx<=0.0) {
     double val0= (lower_tail ? 0.0 : 1.0);
     out(i)= log_p ? std::log(val0) : val0;
     continue;
   }
   if (xx>=1.0) {
     double val1= (lower_tail ? 1.0 : 0.0);
     out(i)= log_p ? std::log(val1) : val1;
     continue;
   }

   double xpow= std::pow(xx, ll);
   // pbeta(xpow, gg, dd+1, TRUE, FALSE)
   double val= R::pbeta( xpow, gg, dd+1.0, true, false );
   if (!lower_tail) {
     val= 1.0 - val;
   }
   if (log_p) {
     val= std::log(val);
   }
   out(i)= val;
 }

 return out;
}


// -----------------------------------------------------------------------------
// 3) qmc: Quantile of Beta Power
// -----------------------------------------------------------------------------

//' @title Quantile Function of the Beta Power Distribution
//'
//' @description
//' Computes \eqn{x = Q(p)}, the quantile function for \eqn{BP(γ, δ, λ)}.
//'
//' @param p Vector of probabilities in (0,1) (or in log scale if log_p=TRUE).
//' @param gamma Shape parameter \eqn{γ > 0}.
//' @param delta Shape parameter \eqn{δ \ge 0}.
//' @param lambda Shape parameter \eqn{λ > 0}.
//' @param lower_tail Logical; if TRUE, p=F(x). If FALSE, p=1-F(x).
//' @param log_p Logical; if TRUE, p is log(p).
//'
//' @details
//' The quantile function is
//' \deqn{
//'   Q(p) = \Bigl[\mathrm{qbeta}(p,\;\gamma,\;\delta+1)\Bigr]^{1/\lambda}.
//' }
//' if lower_tail=TRUE and log_p=FALSE. Adjust for tails or logs accordingly.
//'
//' @export
// [[Rcpp::export]]
arma::vec qmc(
   const arma::vec& p,
   const Rcpp::NumericVector& gamma,
   const Rcpp::NumericVector& delta,
   const Rcpp::NumericVector& lambda,
   bool lower_tail=true,
   bool log_p=false
) {
 arma::vec g_vec(gamma.begin(), gamma.size());
 arma::vec d_vec(delta.begin(), delta.size());
 arma::vec l_vec(lambda.begin(), lambda.size());

 size_t N= std::max({ p.n_elem, g_vec.n_elem, d_vec.n_elem, l_vec.n_elem });
 arma::vec out(N);

 for (size_t i=0; i<N; i++){
   double gg= g_vec[i % g_vec.n_elem];
   double dd= d_vec[i % d_vec.n_elem];
   double ll= l_vec[i % l_vec.n_elem];
   double pp= p[i % p.n_elem];

   if (!check_bp_pars(gg,dd,ll)) {
     out(i)= NA_REAL;
     continue;
   }

   // handle log_p
   if (log_p) {
     if (pp>0.0) {
       // log(p)>0 => p>1 => invalid
       out(i)= NA_REAL;
       continue;
     }
     pp= std::exp(pp);
   }
   // handle tail
   if (!lower_tail) {
     pp= 1.0 - pp;
   }

   // boundary
   if (pp<=0.0) {
     out(i)= 0.0;
     continue;
   }
   if (pp>=1.0) {
     out(i)= 1.0;
     continue;
   }

   // step1= R::qbeta(pp, gg, dd+1)
   double y= R::qbeta(pp, gg, dd+1.0, true, false);
   // step2= y^(1/λ)
   double xval;
   if (ll==1.0) {
     xval= y;
   } else {
     xval= std::pow(y, 1.0/ll);
   }
   if (!R_finite(xval) || xval<0.0) xval=0.0;
   if (xval>1.0) xval=1.0;
   out(i)= xval;
 }

 return out;
}


// -----------------------------------------------------------------------------
// 4) rmc: RNG for Beta Power
// -----------------------------------------------------------------------------

//' @title Random Generation from Beta Power Distribution
//'
//' @description
//' Generates n samples from BP(γ, δ, λ) by drawing U ~ Beta(γ, δ+1) then X= U^(1/λ).
//'
//' @param n Integer number of observations.
//' @param gamma Shape parameter \eqn{γ > 0}.
//' @param delta Shape parameter \eqn{δ \ge 0}.
//' @param lambda Shape parameter \eqn{λ > 0}.
//'
//' @return A numeric vector of length n with random draws in (0,1).
//'
//' @details
//' The distribution arises from GKw with \eqn{\alpha=1,\beta=1}. This is sometimes
//' also called the Beta-Power distribution in the literature.
//'
//' @export
// [[Rcpp::export]]
arma::vec rmc(
   int n,
   const Rcpp::NumericVector& gamma,
   const Rcpp::NumericVector& delta,
   const Rcpp::NumericVector& lambda
) {
 if (n<=0) {
   Rcpp::stop("rmc: n must be positive");
 }

 arma::vec g_vec(gamma.begin(), gamma.size());
 arma::vec d_vec(delta.begin(), delta.size());
 arma::vec l_vec(lambda.begin(), lambda.size());

 size_t k= std::max({ g_vec.n_elem, d_vec.n_elem, l_vec.n_elem });
 arma::vec out(n);

 for(int i=0; i<n; i++){
   size_t idx= i%k;
   double gg= g_vec[idx % g_vec.n_elem];
   double dd= d_vec[idx % d_vec.n_elem];
   double ll= l_vec[idx % l_vec.n_elem];

   if(!check_bp_pars(gg,dd,ll)) {
     out(i)= NA_REAL;
     Rcpp::warning("rmc: invalid parameters at index %d", i+1);
     continue;
   }

   double U= R::rbeta(gg, dd+1.0);
   double xval;
   if (ll==1.0) {
     xval= U;
   } else {
     xval= std::pow(U, 1.0/ll);
   }
   if (!R_finite(xval) || xval<0.0) xval=0.0;
   if (xval>1.0) xval=1.0;
   out(i)= xval;
 }

 return out;
}


// -----------------------------------------------------------------------------
// 5) llmc: Negative Log-Likelihood for Beta Power
// -----------------------------------------------------------------------------

//' @title Negative Log-Likelihood for Beta Power Distribution
//'
//' @description
//' Computes -log-likelihood for BP(γ, δ, λ) given data in (0,1).
//'
//' @param par NumericVector of length 3: (γ, δ, λ).
//' @param data NumericVector of observations in (0,1).
//'
//' @details
//' The PDF is
//' \deqn{
//'   f(x)= \frac{\lambda}{B(\gamma,\delta+1)} x^{\gamma\lambda -1} (1 - x^\lambda)^\delta.
//' }
//' The log-likelihood is the sum of log(f(x_i)). We return the negative of that sum.
//'
//' @return A single numeric value. \code{+Inf} if invalid parameters or data out of (0,1).
//'
//' @examples
//' \dontrun{
//' # Generate sample data from a BP distribution
//' set.seed(123)
//' x <- rmc(100, 2, 3, 0.5)
//' hist(x, breaks = 20, main = "BP(2, 3, 0.5) Sample")
//'
//' # Use in optimization with Hessian-based methods
//' result <- optim(c(0.5, 0.5, 0.5), llmc, method = "BFGS",
//'                 hessian = TRUE, data = x)
//'
//' # Compare numerical and analytical derivatives
//' num_grad <- numDeriv::grad(llmc, x = result$par, data = x)
//' num_hess <- numDeriv::hessian(llmc, x = result$par, data = x)
//'
//' ana_grad <- grmc(result$par, data = x)
//' ana_hess <- hsmc(result$par, data = x)
//'
//' # Check differences (should be very small)
//' round(num_grad - ana_grad, 4)
//' round(num_hess - ana_hess, 4)
//' }
//'
//' @export
// [[Rcpp::export]]
double llmc(const Rcpp::NumericVector& par,
           const Rcpp::NumericVector& data) {
 if (par.size()<3) {
   return R_PosInf;
 }
 double gg= par[0];
 double dd= par[1];
 double ll= par[2];

 if(!check_bp_pars(gg, dd, ll)) {
   return R_PosInf;
 }

 arma::vec x = Rcpp::as<arma::vec>(data);
 if (x.n_elem<1) {
   return R_PosInf;
 }
 if (arma::any(x<=0.0) || arma::any(x>=1.0)) {
   return R_PosInf;
 }

 int n= x.n_elem;

 // constant: n*( log(λ) - log B(γ, δ+1 ) )
 double logB= R::lbeta(gg, dd+1.0);
 double cst= n*( std::log(ll) - logB );

 // sum( (γλ -1)* log(x_i ) )
 // sum( δ * log(1 - x_i^λ ) )
 double gl= gg*ll -1.0;

 arma::vec lx= arma::log(x);
 double sum1= gl * arma::sum(lx);

 // sum( δ * log(1 - x^λ) )
 arma::vec x_pow_l = arma::pow(x, ll);

 arma::vec log1mxpw= arma::log(1.0 - x_pow_l);
 double sum2= dd * arma::sum(log1mxpw);

 double loglike= cst + sum1 + sum2;
 return -loglike;
}

//' @title Gradient Function for Beta Power Log-Likelihood
//'
//' @description
//' Calculates the gradient vector (partial derivatives) of the negative log-likelihood
//' function for the Beta Power (BP) distribution, also known as McDonald distribution.
//' This function provides the exact gradient needed for efficient optimization in
//' maximum likelihood estimation. The BP is a submodel of GKw with α = 1 and β = 1 fixed.
//'
//' @param par NumericVector of length 3 containing parameters (γ, δ, λ) in that order.
//'        All parameters must be positive.
//' @param data NumericVector of observations, where all values must be in the open interval (0,1).
//'
//' @return NumericVector of length 3 containing the gradient components (partial derivatives)
//'         of the negative log-likelihood with respect to each parameter (γ, δ, λ).
//'         Returns a vector of NaN values if any parameters or data values are invalid.
//'
//' @details
//' The gradient vector contains the following partial derivatives of the negative log-likelihood:
//'
//' \deqn{
//' \frac{\partial \ell}{\partial \gamma} = -n[\psi(\gamma+\delta+1) - \psi(\gamma)] -
//' \lambda\sum_{i=1}^{n}\log(x_i)
//' }
//'
//' \deqn{
//' \frac{\partial \ell}{\partial \delta} = -n[\psi(\gamma+\delta+1) - \psi(\delta+1)] -
//' \sum_{i=1}^{n}\log(1-x_i^{\lambda})
//' }
//'
//' \deqn{
//' \frac{\partial \ell}{\partial \lambda} = -\frac{n}{\lambda} - \gamma\sum_{i=1}^{n}\log(x_i) +
//' \delta\sum_{i=1}^{n}\frac{x_i^{\lambda}\log(x_i)}{1-x_i^{\lambda}}
//' }
//'
//' where:
//' \itemize{
//'   \item \deqn{\psi} is the digamma function (derivative of the log-gamma function)
//' }
//'
//' The implementation includes several numerical safeguards:
//' \itemize{
//'   \item Parameter and data validation with appropriate error handling
//'   \item Clamping of intermediate values to avoid numerical underflow/overflow
//'   \item Efficient vector operations using Armadillo C++ library
//' }
//'
//' The returned gradient is negated to align with minimization of negative log-likelihood
//' in optimization routines.
//'
//' @examples
//' \dontrun{
//' # Generate sample data from a BP distribution
//' set.seed(123)
//' x <- rmc(100, 2, 3, 0.5)
//' hist(x, breaks = 20, main = "BP(2, 3, 0.5) Sample")
//'
//' # Use in optimization with Hessian-based methods
//' result <- optim(c(0.5, 0.5, 0.5), llmc, method = "BFGS",
//'                 hessian = TRUE, data = x)
//'
//' # Compare numerical and analytical derivatives
//' num_grad <- numDeriv::grad(llmc, x = result$par, data = x)
//' num_hess <- numDeriv::hessian(llmc, x = result$par, data = x)
//'
//' ana_grad <- grmc(result$par, data = x)
//' ana_hess <- hsmc(result$par, data = x)
//'
//' # Check differences (should be very small)
//' round(num_grad - ana_grad, 4)
//' round(num_hess - ana_hess, 4)
//' }
//'
//' @seealso
//' \code{\link[gkwreg]{llmc}} for the negative log-likelihood function,
//' \code{\link[gkwreg]{hsmc}} for the Hessian matrix of the BP log-likelihood,
//' \code{\link[gkwreg]{dmc}} for the BP density function,
//'
//' @references
//' McDonald, J. B. (1984). Some generalized functions for the size distribution of income.
//' Econometrica, 52(3), 647-665.
//'
//' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized distributions.
//' Journal of Statistical Computation and Simulation, 81(7), 883-898.
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector grmc(const Rcpp::NumericVector& par, const Rcpp::NumericVector& data) {
 // Parameter extraction
 double gamma = par[0];   // Shape parameter γ > 0
 double delta = par[1];   // Shape parameter δ > 0
 double lambda = par[2];  // Shape parameter λ > 0

 // Parameter validation
 if (gamma <= 0 || delta <= 0 || lambda <= 0) {
   Rcpp::NumericVector grad(3, R_NaN);
   return grad;
 }

 // Data conversion and validation
 arma::vec x = Rcpp::as<arma::vec>(data);

 if (arma::any(x <= 0) || arma::any(x >= 1)) {
   Rcpp::NumericVector grad(3, R_NaN);
   return grad;
 }

 int n = x.n_elem;  // Sample size

 // Initialize gradient vector
 Rcpp::NumericVector grad(3, 0.0);

 // Small constant to avoid numerical issues
 double eps = std::numeric_limits<double>::epsilon() * 100;

 // Compute transformations and intermediate values
 arma::vec log_x = arma::log(x);                   // log(x_i)
 arma::vec x_lambda = arma::pow(x, lambda);        // x_i^λ
 arma::vec x_lambda_log_x = x_lambda % log_x;      // x_i^λ * log(x_i)

 // v_i = 1 - x_i^λ
 arma::vec v = 1.0 - x_lambda;
 v = arma::clamp(v, eps, 1.0 - eps);               // Prevent numerical issues

 arma::vec log_v = arma::log(v);                   // log(1 - x_i^λ)

 // Calculate partial derivatives for each parameter (for log-likelihood)

 // ∂ℓ/∂γ = -n[ψ(γ+δ+1) - ψ(γ)] - λΣᵢlog(xᵢ)
 double d_gamma = -n * (R::digamma(gamma + delta + 1) - R::digamma(gamma)) -
   lambda * arma::sum(log_x);

 // ∂ℓ/∂δ = -n[ψ(γ+δ+1) - ψ(δ+1)] - Σᵢlog(1-xᵢ^λ)
 double d_delta = -n * (R::digamma(gamma + delta + 1) - R::digamma(delta + 1)) -
   arma::sum(log_v);

 // ∂ℓ/∂λ = -n/λ - γΣᵢlog(xᵢ) + δΣᵢ[(xᵢ^λ*log(xᵢ))/(1-xᵢ^λ)]
 double d_lambda = -n / lambda - gamma * arma::sum(log_x) +
   delta * arma::sum(x_lambda_log_x / v);

 // Since we're optimizing negative log-likelihood, negate all derivatives
 grad[0] = -d_gamma;
 grad[1] = -d_delta;
 grad[2] = -d_lambda;

 return Rcpp::wrap(grad);
}




//' @title Hessian Matrix Function for Beta Power (McDonald) Log-Likelihood
//'
//' @description
//' Calculates the Hessian matrix (second-order partial derivatives) of the negative log-likelihood
//' function for the Beta Power (BP) distribution, also known as McDonald distribution.
//' This function provides the exact Hessian needed for efficient optimization in maximum likelihood
//' estimation and for asymptotic inference. The BP is a submodel of GKw with α = 1 and β = 1 fixed.
//'
//' @param par NumericVector of length 3 containing parameters (γ, δ, λ) in that order.
//'        All parameters must be positive.
//' @param data NumericVector of observations, where all values must be in the open interval (0,1).
//'
//' @return NumericMatrix of dimensions 3x3 containing the Hessian components (second derivatives)
//'         of the negative log-likelihood with respect to each parameter pair (γ, δ, λ).
//'         Returns a matrix of NaN values if any parameters or data values are invalid.
//'
//' @details
//' The Hessian matrix contains the following second derivatives of the negative log-likelihood:
//'
//' \deqn{
//' \frac{\partial^2 \ell}{\partial \gamma^2} = -n[\psi'(\gamma+\delta+1) - \psi'(\gamma)]
//' }
//'
//' \deqn{
//' \frac{\partial^2 \ell}{\partial \gamma \partial \delta} = -n\psi'(\gamma+\delta+1)
//' }
//'
//' \deqn{
//' \frac{\partial^2 \ell}{\partial \gamma \partial \lambda} = -\sum_{i=1}^{n}\log(x_i)
//' }
//'
//' \deqn{
//' \frac{\partial^2 \ell}{\partial \delta^2} = -n[\psi'(\gamma+\delta+1) - \psi'(\delta+1)]
//' }
//'
//' \deqn{
//' \frac{\partial^2 \ell}{\partial \delta \partial \lambda} = \sum_{i=1}^{n}\frac{x_i^{\lambda}\log(x_i)}{1-x_i^{\lambda}}
//' }
//'
//' \deqn{
//' \frac{\partial^2 \ell}{\partial \lambda^2} = \frac{n}{\lambda^2} + \gamma^2\sum_{i=1}^{n}[\log(x_i)]^2 +
//' \delta\sum_{i=1}^{n}\frac{x_i^{\lambda}[\log(x_i)]^2}{1-x_i^{\lambda}}\left(1 + \frac{x_i^{\lambda}}{1-x_i^{\lambda}}\right)
//' }
//'
//' where:
//' \itemize{
//'   \item \deqn{\psi'} is the trigamma function (second derivative of the log-gamma function)
//' }
//'
//' The implementation includes several numerical safeguards:
//' \itemize{
//'   \item Parameter and data validation with appropriate error handling
//'   \item Clamping of intermediate values to avoid numerical underflow/overflow
//'   \item Symmetry of the Hessian matrix is guaranteed by construction
//'   \item Efficient vector operations using Armadillo C++ library
//' }
//'
//' The returned Hessian is negated to align with minimization of negative log-likelihood
//' in optimization routines.
//'
//' @examples
//' \dontrun{
//' # Generate sample data from a BP distribution
//' set.seed(123)
//' x <- rmc(100, 2, 3, 0.5)
//' hist(x, breaks = 20, main = "BP(2, 3, 0.5) Sample")
//'
//' # Use in optimization with Hessian-based methods
//' result <- optim(c(0.5, 0.5, 0.5), llmc, method = "BFGS",
//'                 hessian = TRUE, data = x)
//'
//' # Compare numerical and analytical derivatives
//' num_grad <- numDeriv::grad(llmc, x = result$par, data = x)
//' num_hess <- numDeriv::hessian(llmc, x = result$par, data = x)
//'
//' ana_grad <- grmc(result$par, data = x)
//' ana_hess <- hsmc(result$par, data = x)
//'
//' # Check differences (should be very small)
//' round(num_grad - ana_grad, 4)
//' round(num_hess - ana_hess, 4)
//' }
//'
//' @seealso
//' \code{\link[gkwreg]{llmc}} for the negative log-likelihood function,
//' \code{\link[gkwreg]{grmc}} for the gradient of the BP log-likelihood,
//' \code{\link[gkwreg]{dmc}} for the BP density function,
//'
//' @references
//' McDonald, J. B. (1984). Some generalized functions for the size distribution of income.
//' Econometrica, 52(3), 647-665.
//'
//' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized distributions.
//' Journal of Statistical Computation and Simulation, 81(7), 883-898.
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix hsmc(const Rcpp::NumericVector& par, const Rcpp::NumericVector& data) {
 // Parameter extraction
 double gamma = par[0];   // Shape parameter γ > 0
 double delta = par[1];   // Shape parameter δ > 0
 double lambda = par[2];  // Shape parameter λ > 0

 // Initialize Hessian matrix
 Rcpp::NumericMatrix hess(3, 3);

 // Parameter validation
 if (gamma <= 0 || delta <= 0 || lambda <= 0) {
   hess.fill(R_NaN);
   return hess;
 }

 // Data conversion and validation
 arma::vec x = Rcpp::as<arma::vec>(data);

 if (arma::any(x <= 0) || arma::any(x >= 1)) {
   hess.fill(R_NaN);
   return hess;
 }

 int n = x.n_elem;  // Sample size

 // Small constant to avoid numerical issues
 double eps = std::numeric_limits<double>::epsilon() * 100;

 // Compute transformations and intermediate values
 arma::vec log_x = arma::log(x);                   // log(x_i)
 arma::vec log_x_squared = arma::square(log_x);    // [log(x_i)]²
 arma::vec x_lambda = arma::pow(x, lambda);        // x_i^λ
 arma::vec x_lambda_log_x = x_lambda % log_x;      // x_i^λ * log(x_i)

 // v_i = 1 - x_i^λ
 arma::vec v = 1.0 - x_lambda;
 v = arma::clamp(v, eps, 1.0 - eps);               // Prevent numerical issues

 // Additional terms for second derivatives
 arma::vec term_ratio = x_lambda / v;              // x_i^λ / (1-x_i^λ)
 arma::vec term_combined = 1.0 + term_ratio;       // 1 + x_i^λ/(1-x_i^λ)

 // Compute trigamma values
 // trigamma(x) = d²(log(Gamma(x)))/dx² = d(digamma(x))/dx
 double trigamma_gamma = R::trigamma(gamma);
 double trigamma_delta_plus_1 = R::trigamma(delta + 1.0);
 double trigamma_gamma_plus_delta_plus_1 = R::trigamma(gamma + delta + 1.0);

 // Calculate the Hessian components for negative log-likelihood

 // H[0,0] = ∂²ℓ/∂γ² = -n[ψ'(γ+δ+1) - ψ'(γ)]
 double h_gamma_gamma = -n * (trigamma_gamma_plus_delta_plus_1 - trigamma_gamma);

 // H[0,1] = H[1,0] = ∂²ℓ/∂γ∂δ = -n*ψ'(γ+δ+1)
 double h_gamma_delta = -n * trigamma_gamma_plus_delta_plus_1;

 // H[0,2] = H[2,0] = ∂²ℓ/∂γ∂λ = -Σᵢlog(xᵢ)
 double h_gamma_lambda = -arma::sum(log_x);

 // H[1,1] = ∂²ℓ/∂δ² = -n[ψ'(γ+δ+1) - ψ'(δ+1)]
 double h_delta_delta = -n * (trigamma_gamma_plus_delta_plus_1 - trigamma_delta_plus_1);

 // H[1,2] = H[2,1] = ∂²ℓ/∂δ∂λ = Σᵢ[x_i^λ*log(x_i)/(1-x_i^λ)]
 double h_delta_lambda = arma::sum(x_lambda_log_x / v);

 // H[2,2] = ∂²ℓ/∂λ² = n/λ² + δ*Σᵢ[x_i^λ*[log(x_i)]²/(1-x_i^λ)*(1 + x_i^λ/(1-x_i^λ))]
 double h_lambda_lambda = n / (lambda * lambda);

 arma::vec lambda_term = delta * x_lambda % log_x_squared % term_combined / v;
 h_lambda_lambda += arma::sum(lambda_term);

 // Fill the Hessian matrix (symmetric)
 hess(0, 0) = -h_gamma_gamma;
 hess(0, 1) = hess(1, 0) = -h_gamma_delta;
 hess(0, 2) = hess(2, 0) = -h_gamma_lambda;
 hess(1, 1) = -h_delta_delta;
 hess(1, 2) = hess(2, 1) = -h_delta_lambda;
 hess(2, 2) = -h_lambda_lambda;

 return Rcpp::wrap(hess);
}







/*
----------------------------------------------------------------------------
ASSUMPTIONS:
We rely on the same numeric stability functions from a "core" or gkwdist.cpp file:
- log1mexp(double)
- safe_log(double)
- safe_exp(double)
- safe_pow(double,double)
etc.

----------------------------------------------------------------------------
KUMARASWAMY (Kw) DISTRIBUTION
----------------------------------------------------------------------------

Parameters: alpha>0, beta>0.

* PDF:
f(x) = alpha * beta * x^(alpha -1) * (1 - x^alpha)^(beta -1),  for 0<x<1.

* CDF:
F(x)= 1 - (1 - x^alpha)^beta.

* QUANTILE:
Q(p)= {1 - [1 - p]^(1/beta)}^(1/alpha).

* RANDOM GENERATION:
If V ~ Uniform(0,1), then X= {1 - [1 - V]^(1/beta)}^(1/alpha).

* NEGATIVE LOG-LIKELIHOOD:
sum over i of -log( f(x_i) ).
log f(x_i)= log(alpha) + log(beta) + (alpha-1)*log(x_i) + (beta-1)*log(1 - x_i^alpha).

*/

// -----------------------------------------------------------------------------
// Parameter checker for Kumaraswamy
// alpha>0, beta>0
// -----------------------------------------------------------------------------
inline bool check_kw_pars(double alpha, double beta, bool strict=false) {
if (alpha <=0.0 || beta <=0.0) {
  return false;
}
if (strict) {
  const double MINP=1e-6, MAXP=1e6;
  if (alpha<MINP || beta<MINP) return false;
  if (alpha>MAXP || beta>MAXP) return false;
}
return true;
}

// -----------------------------------------------------------------------------
// 1) dkw: PDF of Kumaraswamy
// -----------------------------------------------------------------------------

//' @title PDF of the Kumaraswamy Distribution
//'
//' @description
//' Computes the PDF of Kw(\eqn{\alpha,\beta}) at values x in (0,1).
//'
//' @param x Vector of quantiles in (0,1).
//' @param alpha Shape parameter \eqn{\alpha>0}.
//' @param beta Shape parameter \eqn{\beta>0}.
//' @param log_prob If TRUE, returns log of the PDF.
//'
//' @return A vector of (log-)densities, same length as broadcast of x and parameters.
//'
//' @details
//' PDF: \eqn{f(x)= \alpha\,\beta \; x^{\alpha-1}\,(1 - x^\alpha)^{\beta-1}}, for 0<x<1.
//'
//' @export
// [[Rcpp::export]]
arma::vec dkw(
   const arma::vec& x,
   const Rcpp::NumericVector& alpha,
   const Rcpp::NumericVector& beta,
   bool log_prob=false
) {
 arma::vec a_vec(alpha.begin(), alpha.size());
 arma::vec b_vec(beta.begin(), beta.size());

 size_t N= std::max({ x.n_elem, a_vec.n_elem, b_vec.n_elem });
 arma::vec out(N);

 out.fill(log_prob ? R_NegInf : 0.0);

 for (size_t i=0; i<N; i++){
   double a= a_vec[i % a_vec.n_elem];
   double b= b_vec[i % b_vec.n_elem];
   double xx= x[i % x.n_elem];

   if (!check_kw_pars(a,b)) {
     // invalid => pdf=0 or logpdf=-Inf
     continue;
   }
   if (xx<=0.0 || xx>=1.0 || !R_finite(xx)) {
     // outside domain => 0 or -Inf
     continue;
   }

   // log f(x)= log(a)+ log(b) + (a-1)* log(x) + (b-1)* log(1- x^a)
   double la= std::log(a);
   double lb= std::log(b);

   double lx= std::log(xx);
   double xalpha= a* lx; // log(x^a)
   // log(1- x^a)= log1mexp(xalpha)
   double log_1_xalpha= log1mexp(xalpha);
   if (!R_finite(log_1_xalpha)) {
     continue;
   }

   double log_pdf= la + lb + (a-1.0)* lx + (b-1.0)* log_1_xalpha;

   if (log_prob) {
     out(i)= log_pdf;
   } else {
     out(i)= std::exp(log_pdf);
   }
 }

 return out;
}


// -----------------------------------------------------------------------------
// 2) pkw: CDF of Kumaraswamy
// -----------------------------------------------------------------------------

//' @title CDF of the Kumaraswamy Distribution
//'
//' @description
//' Computes \eqn{F(q)= P(X ≤ q)} for Kw(\eqn{\alpha,\beta}).
//'
//' @param q Vector of quantiles in (0,1).
//' @param alpha Shape parameter \eqn{\alpha>0}.
//' @param beta Shape parameter \eqn{\beta>0}.
//' @param lower_tail If TRUE (default), returns \eqn{F(q)=P(X≤q)}. If FALSE, 1-F(q).
//' @param log_p If TRUE, returns log probabilities.
//'
//' @return Vector of probabilities or log-probabilities, matching broadcast of q and parameters.
//'
//' @details
//' \eqn{F(x)= 1 - (1 - x^\alpha)^\beta}, for 0<x<1.
//'
//' @export
// [[Rcpp::export]]
arma::vec pkw(
   const arma::vec& q,
   const Rcpp::NumericVector& alpha,
   const Rcpp::NumericVector& beta,
   bool lower_tail=true,
   bool log_p=false
) {
 arma::vec a_vec(alpha.begin(), alpha.size());
 arma::vec b_vec(beta.begin(), beta.size());

 size_t N= std::max({ q.n_elem, a_vec.n_elem, b_vec.n_elem });
 arma::vec out(N);

 for (size_t i=0; i<N; i++){
   double a= a_vec[i % a_vec.n_elem];
   double b= b_vec[i % b_vec.n_elem];
   double xx= q[i % q.n_elem];

   if (!check_kw_pars(a,b)) {
     out(i)= NA_REAL;
     continue;
   }

   // boundary
   if (!R_finite(xx) || xx<=0.0) {
     double val0= (lower_tail ? 0.0 : 1.0);
     out(i)= log_p ? std::log(val0) : val0;
     continue;
   }
   if (xx>=1.0) {
     double val1= (lower_tail ? 1.0 : 0.0);
     out(i)= log_p ? std::log(val1) : val1;
     continue;
   }

   double xalpha= std::pow(xx, a);
   double tmp= 1.0 - std::pow( (1.0 - xalpha), b );
   if (tmp<=0.0) {
     double val0= (lower_tail ? 0.0 : 1.0);
     out(i)= log_p ? std::log(val0) : val0;
     continue;
   }
   if (tmp>=1.0) {
     double val1= (lower_tail ? 1.0 : 0.0);
     out(i)= log_p ? std::log(val1) : val1;
     continue;
   }

   double val= tmp;
   if (!lower_tail) {
     val= 1.0- val;
   }
   if (log_p) {
     val= std::log(val);
   }
   out(i)= val;
 }

 return out;
}

// -----------------------------------------------------------------------------
// 3) qkw: Quantile of Kumaraswamy
// -----------------------------------------------------------------------------

//' @title Quantile Function of the Kumaraswamy Distribution
//'
//' @description
//' For p in (0,1), returns x=Q(p) such that F(x)=p, where F is Kw(\eqn{\alpha,\beta}).
//'
//' @param p Vector of probabilities in (0,1) (or log scale if log_p=TRUE).
//' @param alpha Shape parameter \eqn{\alpha>0}.
//' @param beta Shape parameter \eqn{\beta>0}.
//' @param lower_tail If TRUE, p=F(x); if FALSE, p=1-F(x).
//' @param log_p If TRUE, p is log(p).
//'
//' @details
//' \eqn{Q(p)= {1 - [1 - p]^(1/\beta)}^(1/\alpha)}
//'
//' @export
// [[Rcpp::export]]
arma::vec qkw(
   const arma::vec& p,
   const Rcpp::NumericVector& alpha,
   const Rcpp::NumericVector& beta,
   bool lower_tail=true,
   bool log_p=false
) {
 arma::vec a_vec(alpha.begin(), alpha.size());
 arma::vec b_vec(beta.begin(), beta.size());

 size_t N= std::max({ p.n_elem, a_vec.n_elem, b_vec.n_elem });
 arma::vec out(N);

 for (size_t i=0; i<N; i++){
   double a= a_vec[i % a_vec.n_elem];
   double b= b_vec[i % b_vec.n_elem];
   double pp= p[i % p.n_elem];

   if (!check_kw_pars(a,b)) {
     out(i)= NA_REAL;
     continue;
   }

   // convert if log
   if (log_p) {
     if (pp>0.0) {
       // invalid => p>1
       out(i)= NA_REAL;
       continue;
     }
     pp= std::exp(pp);
   }
   if (!lower_tail) {
     pp= 1.0 - pp;
   }

   // boundary
   if (pp<=0.0) {
     out(i)= 0.0;
     continue;
   }
   if (pp>=1.0) {
     out(i)= 1.0;
     continue;
   }

   // x= {1 - [1 - p]^(1/beta)}^(1/alpha)
   double step1= 1.0 - pp;
   if (step1<0.0) step1=0.0;
   double step2= std::pow(step1, 1.0/b);
   double step3= 1.0 - step2;
   if (step3<0.0) step3=0.0;

   double xval;
   if (a==1.0) {
     xval= step3;
   } else {
     xval= std::pow(step3, 1.0/a);
     if (!R_finite(xval)|| xval<0.0) xval=0.0;
     if (xval>1.0) xval=1.0;
   }
   out(i)= xval;
 }

 return out;
}


// -----------------------------------------------------------------------------
// 4) rkw: Random Generation from Kumaraswamy
// -----------------------------------------------------------------------------

//' @title Random Generation from the Kumaraswamy Distribution
//'
//' @description
//' Generates n samples from Kw(\eqn{\alpha,\beta}), using the transformation from Uniform(0,1).
//'
//' @param n Integer number of observations.
//' @param alpha Shape parameter \eqn{\alpha>0}.
//' @param beta Shape parameter \eqn{\beta>0}.
//'
//' @return A vector of length n with samples in (0,1).
//'
//' @details
//' If V ~ Uniform(0,1), then
//' \eqn{X= \{1 - [1 - V]^{1/\beta}\}^{1/\alpha}}.
//'
//' @export
// [[Rcpp::export]]
arma::vec rkw(
   int n,
   const Rcpp::NumericVector& alpha,
   const Rcpp::NumericVector& beta
) {
 if (n<=0) {
   Rcpp::stop("rkw: n must be positive");
 }

 arma::vec a_vec(alpha.begin(), alpha.size());
 arma::vec b_vec(beta.begin(), beta.size());

 size_t k= std::max({ a_vec.n_elem, b_vec.n_elem });
 arma::vec out(n);

 for (int i=0; i<n; i++){
   size_t idx= i % k;
   double a= a_vec[idx % a_vec.n_elem];
   double b= b_vec[idx % b_vec.n_elem];

   if (!check_kw_pars(a,b)) {
     out(i)= NA_REAL;
     Rcpp::warning("rkw: invalid parameters at index %d", i+1);
     continue;
   }

   double U= R::runif(0.0,1.0);
   // X= {1 - [1 - U]^(1/beta)}^(1/alpha)
   double step1= 1.0 - U;
   if (step1<0.0) step1=0.0;
   double step2= std::pow(step1, 1.0/b);
   double step3= 1.0 - step2;
   if (step3<0.0) step3=0.0;

   double x;
   if (a==1.0) {
     x= step3;
   } else {
     x= std::pow(step3, 1.0/a);
     if (!R_finite(x)|| x<0.0) x=0.0;
     if (x>1.0) x=1.0;
   }
   out(i)= x;
 }

 return out;
}

// -----------------------------------------------------------------------------
// 5) llkw: Negative Log-Likelihood for Kumaraswamy
// -----------------------------------------------------------------------------

//' @title Negative Log-Likelihood of the Kumaraswamy Distribution
//'
//' @description
//' Computes the negative log-likelihood for Kw(\eqn{\alpha,\beta}) given data in (0,1).
//'
//' @param par NumericVector of length 2: (alpha, beta).
//' @param data NumericVector of observations in (0,1).
//'
//' @return A single numeric value. \code{+Inf} if invalid parameters or data.
//'
//' @details
//' The PDF is
//' \eqn{ f(x)= \alpha\,\beta\,x^{\alpha-1}\,(1 - x^\alpha)^{\beta-1} }.
//' Summation of log gives the log-likelihood; we return the negative of that sum.
//'
//' @examples
//' \dontrun{
//' # Generate sample data from a Kw distribution
//' set.seed(123)
//' x <- rkw(100, 2, 3)
//' hist(x, breaks = 20, main = "Kw(2, 3) Sample")
//'
//' # Use in optimization with Hessian-based methods
//' result <- optim(c(0.5, 0.5), llkw, method = "BFGS",
//'                 hessian = TRUE, data = x)
//'
//' # Compare numerical and analytical derivatives
//' num_grad <- numDeriv::grad(llkw, x = result$par, data = x)
//' num_hess <- numDeriv::hessian(llkw, x = result$par, data = x)
//'
//' ana_grad <- grkw(result$par, data = x)
//' ana_hess <- hskw(result$par, data = x)
//'
//' # Check differences (should be very small)
//' round(num_grad - ana_grad, 4)
//' round(num_hess - ana_hess, 4)
//' }
//'
//' @export
// [[Rcpp::export]]
double llkw(const Rcpp::NumericVector& par,
           const Rcpp::NumericVector& data) {
 if (par.size()<2) {
   return R_PosInf;
 }
 double a= par[0];
 double b= par[1];

 if (!check_kw_pars(a,b)) {
   return R_PosInf;
 }

 arma::vec x= Rcpp::as<arma::vec>(data);
 if (x.n_elem<1) {
   return R_PosInf;
 }
 if (arma::any(x<=0.0) || arma::any(x>=1.0)) {
   return R_PosInf;
 }

 int n= x.n_elem;
 // constant: n*( log(a)+ log(b) )
 double cst= n*( std::log(a) + std::log(b) );

 // sum( (a-1)* log(x_i ) )
 arma::vec lx= arma::log(x);
 double sum1= (a-1.0)* arma::sum(lx);

 // sum( (b-1)* log(1- x^a) )
 arma::vec xalpha= arma::pow(x,a);
 arma::vec log_1_xalpha= arma::log(1.0 - xalpha);
 double sum2= (b-1.0)* arma::sum(log_1_xalpha);

 double loglike= cst + sum1 + sum2;
 // negative
 return -loglike;
}


//' @title Gradient Function for Kumaraswamy Log-Likelihood
//'
//' @description
//' Calculates the gradient vector (partial derivatives) of the negative log-likelihood
//' function for the Kumaraswamy (Kw) distribution. This function provides the exact
//' gradient needed for efficient optimization in maximum likelihood estimation.
//' The Kw is a submodel of GKw with γ = 1, δ = 0, and λ = 1 fixed.
//'
//' @param par NumericVector of length 2 containing parameters (α, β) in that order.
//'        All parameters must be positive.
//' @param data NumericVector of observations, where all values must be in the open interval (0,1).
//'
//' @return NumericVector of length 2 containing the gradient components (partial derivatives)
//'         of the negative log-likelihood with respect to each parameter (α, β).
//'         Returns a vector of NaN values if any parameters or data values are invalid.
//'
//' @details
//' The gradient vector contains the following partial derivatives of the negative log-likelihood:
//'
//' \deqn{
//' \frac{\partial \ell}{\partial \alpha} = \frac{n}{\alpha} + \sum_{i=1}^{n}\log(x_i) -
//' \sum_{i=1}^{n}\left[\frac{(\beta-1)x_i^{\alpha}\log(x_i)}{1-x_i^{\alpha}}\right]
//' }
//'
//' \deqn{
//' \frac{\partial \ell}{\partial \beta} = \frac{n}{\beta} + \sum_{i=1}^{n}\log(1-x_i^{\alpha})
//' }
//'
//' The implementation includes several numerical safeguards:
//' \itemize{
//'   \item Parameter and data validation with appropriate error handling
//'   \item Clamping of intermediate values to avoid numerical underflow/overflow
//'   \item Efficient vector operations using Armadillo C++ library
//' }
//'
//' The returned gradient is negated to align with minimization of negative log-likelihood
//' in optimization routines.
//'
//' @examples
//' \dontrun{
//' # Generate sample data from a Kw distribution
//' set.seed(123)
//' x <- rkw(100, 2, 3)
//' hist(x, breaks = 20, main = "Kw(2, 3) Sample")
//'
//' # Use in optimization with Hessian-based methods
//' result <- optim(c(0.5, 0.5), llkw, method = "BFGS",
//'                 hessian = TRUE, data = x)
//'
//' # Compare numerical and analytical derivatives
//' num_grad <- numDeriv::grad(llkw, x = result$par, data = x)
//' num_hess <- numDeriv::hessian(llkw, x = result$par, data = x)
//'
//' ana_grad <- grkw(result$par, data = x)
//' ana_hess <- hskw(result$par, data = x)
//'
//' # Check differences (should be very small)
//' round(num_grad - ana_grad, 4)
//' round(num_hess - ana_hess, 4)
//' }
//'
//' @seealso
//' \code{\link[gkwreg]{llkw}} for the negative log-likelihood function,
//' \code{\link[gkwreg]{hskw}} for the Hessian matrix of the Kw log-likelihood,
//' \code{\link[gkwreg]{dkw}} for the Kw density function,
//'
//' @references
//' Kumaraswamy, P. (1980). A generalized probability density function for double-bounded random processes.
//' Journal of Hydrology, 46(1-2), 79-88.
//'
//' Jones, M. C. (2009). Kumaraswamy's distribution: A beta-type distribution with some tractability advantages.
//' Statistical Methodology, 6(1), 70-81.
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector grkw(const Rcpp::NumericVector& par, const Rcpp::NumericVector& data) {
 // Parameter extraction
 double alpha = par[0];   // Shape parameter α > 0
 double beta = par[1];    // Shape parameter β > 0

 // Parameter validation
 if (alpha <= 0 || beta <= 0) {
   Rcpp::NumericVector grad(2, R_NaN);
   return grad;
 }

 // Data conversion and validation
 arma::vec x = Rcpp::as<arma::vec>(data);

 if (arma::any(x <= 0) || arma::any(x >= 1)) {
   Rcpp::NumericVector grad(2, R_NaN);
   return grad;
 }

 int n = x.n_elem;  // Sample size

 // Initialize gradient vector
 Rcpp::NumericVector grad(2, 0.0);

 // Small constant to avoid numerical issues
 double eps = std::numeric_limits<double>::epsilon() * 100;

 // Compute transformations and intermediate values
 arma::vec log_x = arma::log(x);                // log(x_i)
 arma::vec x_alpha = arma::pow(x, alpha);       // x_i^α
 arma::vec x_alpha_log_x = x_alpha % log_x;     // x_i^α * log(x_i)

 // v_i = 1 - x_i^α
 arma::vec v = 1.0 - x_alpha;
 v = arma::clamp(v, eps, 1.0 - eps);            // Prevent numerical issues

 arma::vec log_v = arma::log(v);                // log(1-x_i^α)

 // Calculate partial derivatives for each parameter (for log-likelihood)

 // ∂ℓ/∂α = n/α + Σᵢlog(xᵢ) - Σᵢ[(β-1)xᵢ^α*log(xᵢ)/(1-xᵢ^α)]
 double d_alpha = n / alpha + arma::sum(log_x);

 // Calculate the term for α gradient
 arma::vec alpha_term = (beta - 1.0) * x_alpha_log_x / v;

 d_alpha -= arma::sum(alpha_term);

 // ∂ℓ/∂β = n/β + Σᵢlog(1-xᵢ^α)
 double d_beta = n / beta + arma::sum(log_v);

 // Since we're optimizing negative log-likelihood, negate all derivatives
 grad[0] = -d_alpha;
 grad[1] = -d_beta;

 return grad;
}



//' @title Hessian Matrix Function for Kumaraswamy Log-Likelihood
//'
//' @description
//' Calculates the Hessian matrix (second-order partial derivatives) of the negative log-likelihood
//' function for the Kumaraswamy (Kw) distribution. This function provides the exact Hessian needed
//' for efficient optimization in maximum likelihood estimation and for asymptotic inference.
//' The Kw is a submodel of GKw with γ = 1, δ = 0, and λ = 1 fixed.
//'
//' @param par NumericVector of length 2 containing parameters (α, β) in that order.
//'        All parameters must be positive.
//' @param data NumericVector of observations, where all values must be in the open interval (0,1).
//'
//' @return NumericMatrix of dimensions 2x2 containing the Hessian components (second derivatives)
//'         of the negative log-likelihood with respect to each parameter pair (α, β).
//'         Returns a matrix of NaN values if any parameters or data values are invalid.
//'
//' @details
//' The Hessian matrix contains the following second derivatives of the negative log-likelihood:
//'
//' \deqn{
//' \frac{\partial^2 \ell}{\partial \alpha^2} = -\frac{n}{\alpha^2} -
//' \sum_{i=1}^{n}\left[\frac{(\beta-1)x_i^{\alpha}(\log(x_i))^2}{1-x_i^{\alpha}}\left(1 + \frac{x_i^{\alpha}}{1-x_i^{\alpha}}\right)\right]
//' }
//'
//' \deqn{
//' \frac{\partial^2 \ell}{\partial \alpha \partial \beta} = -
//' \sum_{i=1}^{n}\left[\frac{x_i^{\alpha}\log(x_i)}{1-x_i^{\alpha}}\right]
//' }
//'
//' \deqn{
//' \frac{\partial^2 \ell}{\partial \beta^2} = -\frac{n}{\beta^2}
//' }
//'
//' The implementation includes several numerical safeguards:
//' \itemize{
//'   \item Parameter and data validation with appropriate error handling
//'   \item Clamping of intermediate values to avoid numerical underflow/overflow
//'   \item Symmetry of the Hessian matrix is guaranteed by construction
//'   \item Efficient vector operations using Armadillo C++ library
//' }
//'
//' The returned Hessian is negated to align with minimization of negative log-likelihood
//' in optimization routines.
//'
//' @examples
//' \dontrun{
//' # Generate sample data from a Kw distribution
//' set.seed(123)
//' x <- rkw(100, 2, 3)
//' hist(x, breaks = 20, main = "Kw(2, 3) Sample")
//'
//' # Use in optimization with Hessian-based methods
//' result <- optim(c(0.5, 0.5), llkw, method = "BFGS",
//'                 hessian = TRUE, data = x)
//'
//' # Compare numerical and analytical derivatives
//' num_grad <- numDeriv::grad(llkw, x = result$par, data = x)
//' num_hess <- numDeriv::hessian(llkw, x = result$par, data = x)
//'
//' ana_grad <- grkw(result$par, data = x)
//' ana_hess <- hskw(result$par, data = x)
//'
//' # Check differences (should be very small)
//' round(num_grad - ana_grad, 4)
//' round(num_hess - ana_hess, 4)
//' }
//'
//' @seealso
//' \code{\link[gkwreg]{llkw}} for the negative log-likelihood function,
//' \code{\link[gkwreg]{grkw}} for the gradient of the Kw log-likelihood,
//' \code{\link[gkwreg]{dkw}} for the Kw density function,
//'
//' @references
//' Kumaraswamy, P. (1980). A generalized probability density function for double-bounded random processes.
//' Journal of Hydrology, 46(1-2), 79-88.
//'
//' Jones, M. C. (2009). Kumaraswamy's distribution: A beta-type distribution with some tractability advantages.
//' Statistical Methodology, 6(1), 70-81.
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix hskw(const Rcpp::NumericVector& par, const Rcpp::NumericVector& data) {
 // Parameter extraction
 double alpha = par[0];   // Shape parameter α > 0
 double beta = par[1];    // Shape parameter β > 0

 // Initialize Hessian matrix
 Rcpp::NumericMatrix hess(2, 2);

 // Parameter validation
 if (alpha <= 0 || beta <= 0) {
   hess.fill(R_NaN);
   return hess;
 }

 // Data conversion and validation
 arma::vec x = Rcpp::as<arma::vec>(data);

 if (arma::any(x <= 0) || arma::any(x >= 1)) {
   hess.fill(R_NaN);
   return hess;
 }

 int n = x.n_elem;  // Sample size

 // Small constant to avoid numerical issues
 double eps = std::numeric_limits<double>::epsilon() * 100;

 // Compute transformations and intermediate values
 arma::vec log_x = arma::log(x);                  // log(x_i)
 arma::vec log_x_squared = arma::square(log_x);   // (log(x_i))²
 arma::vec x_alpha = arma::pow(x, alpha);         // x_i^α
 arma::vec x_alpha_log_x = x_alpha % log_x;       // x_i^α * log(x_i)

 // v_i = 1 - x_i^α
 arma::vec v = 1.0 - x_alpha;
 v = arma::clamp(v, eps, 1.0 - eps);              // Prevent numerical issues

 // Additional terms for second derivatives
 arma::vec term_ratio = x_alpha / v;              // x_i^α / (1-x_i^α)
 arma::vec term_combined = 1.0 + term_ratio;      // 1 + x_i^α/(1-x_i^α)

 // Calculate the Hessian components for negative log-likelihood

 // H[0,0] = ∂²ℓ/∂α² = -n/α² - Σᵢ[(β-1)x_i^α(log(x_i))²/(1-x_i^α)(1 + x_i^α/(1-x_i^α))]
 double h_alpha_alpha = -n / (alpha * alpha);
 arma::vec d2a_terms = (beta - 1.0) * x_alpha % log_x_squared % term_combined / v;
 h_alpha_alpha -= arma::sum(d2a_terms);

 // H[0,1] = H[1,0] = ∂²ℓ/∂α∂β = -Σᵢ[x_i^α*log(x_i)/(1-x_i^α)]
 double h_alpha_beta = -arma::sum(x_alpha_log_x / v);

 // H[1,1] = ∂²ℓ/∂β² = -n/β²
 double h_beta_beta = -n / (beta * beta);

 // Fill the Hessian matrix (symmetric)
 hess(0, 0) = -h_alpha_alpha;
 hess(0, 1) = hess(1, 0) = -h_alpha_beta;
 hess(1, 1) = -h_beta_beta;

 return hess;
}


/*
 ----------------------------------------------------------------------------
 NOTE:
 We assume the same numeric stability functions from a core or gkwdist.cpp:
 - log1mexp(double)
 - safe_log(double)
 - safe_exp(double)
 - safe_pow(double,double)
 etc.
 They are not redefined here to avoid duplication.
 */

/*
 ----------------------------------------------------------------------------
 BETA DISTRIBUTION: Beta(γ, δ)
 ----------------------------------------------------------------------------

 We use parameters gamma (γ) and delta (δ), both > 0, consistent with GKw family.
 Domain: x in (0,1).

 * PDF:
 f(x;γ,δ) = x^(γ-1) * (1-x)^δ / B(γ,δ+1),  for 0<x<1.

 * CDF:
 F(x;γ,δ) = pbeta(x, γ, δ+1).

 * QUANTILE:
 Q(p;γ,δ) = qbeta(p, γ, δ+1).

 * RNG:
 X = rbeta(γ, δ+1).

 * NEGATIVE LOG-LIKELIHOOD:
 For data x_i in (0,1),
 log f(x_i) = (γ-1)*log(x_i) + δ*log(1-x_i) - ln B(γ,δ+1).
 Summation => negative => used in MLE.

 */


// -----------------------------------------------------------------------------
// Parameter checker for Beta(gamma>0, delta>0)
// -----------------------------------------------------------------------------
inline bool check_beta_pars(double gamma, double delta, bool strict=false) {
  if (gamma <= 0.0 || delta <= 0.0) {
    return false;
  }
  if (strict) {
    const double MINP = 1e-7, MAXP = 1e7;
    if (gamma < MINP || delta < MINP) return false;
    if (gamma > MAXP || delta > MAXP) return false;
  }
  return true;
}

// -----------------------------------------------------------------------------
// 1) dbeta_: PDF for Beta distribution
// -----------------------------------------------------------------------------

//' @title Density of the Beta Distribution
//'
//' @description
//' Computes the PDF of Beta(\code{gamma, delta}) for x in (0,1).
//'
//' @param x Vector of quantiles in (0,1).
//' @param gamma First shape parameter > 0.
//' @param delta Second shape parameter > 0.
//' @param log_prob If TRUE, returns log-density.
//'
//' @details
//' The PDF is
//' \deqn{
//'   f(x) = \frac{x^{\gamma-1}\,(1-x)^{\delta}}
//'              {B(\gamma, \delta+1)}, \quad 0<x<1.
//' }
//'
//' @return A vector of densities (or log-densities).
//'
//' @export
// [[Rcpp::export]]
arma::vec dbeta_(
   const arma::vec& x,
   const Rcpp::NumericVector& gamma,
   const Rcpp::NumericVector& delta,
   bool log_prob = false
) {
 arma::vec g_vec(gamma.begin(), gamma.size());
 arma::vec d_vec(delta.begin(), delta.size());

 size_t N = std::max({x.n_elem, g_vec.n_elem, d_vec.n_elem});
 arma::vec out(N);
 out.fill(log_prob ? R_NegInf : 0.0);

 for (size_t i = 0; i < N; i++) {
   double g = g_vec[i % g_vec.n_elem];
   double d = d_vec[i % d_vec.n_elem];
   double xx = x[i % x.n_elem];

   if (!check_beta_pars(g, d)) {
     continue; // => 0 or -Inf
   }
   if (xx <= 0.0 || xx >= 1.0 || !R_finite(xx)) {
     continue;
   }

   // logBeta = R::lbeta(g, d+1)
   double lB = R::lbeta(g, d + 1.0);
   // log pdf = (g-1)*log(x) + d*log(1-x) - lB
   double lx = std::log(xx);
   double one_minus_x = 1.0 - xx;
   if (one_minus_x <= 0.0) {
     // => out of domain, effectively => 0
     continue;
   }
   double log_1_minus_x = std::log(one_minus_x);

   double log_pdf = (g - 1.0) * lx + d * log_1_minus_x - lB;

   if (log_prob) {
     out(i) = log_pdf;
   } else {
     out(i) = std::exp(log_pdf);
   }
 }

 return out;
}


// -----------------------------------------------------------------------------
// 2) pbeta_: CDF for Beta
// -----------------------------------------------------------------------------

//' @title CDF of the Beta Distribution
//'
//' @description
//' Computes F(q) = pbeta(q, gamma, delta+1) in (0,1).
//'
//' @param q Vector of quantiles in (0,1).
//' @param gamma First shape parameter > 0.
//' @param delta Second shape parameter > 0.
//' @param lower_tail If TRUE, returns F(q)=P(X≤q). If FALSE, 1-F(q).
//' @param log_p If TRUE, returns log of the probability.
//'
//' @return A vector of probabilities or log-probabilities, matching broadcast.
//'
//' @details
//' This is a wrapper calling R's \code{pbeta} internally with appropriate parameter
//' adjustment for the GKw family's Beta parameterization. We replicate the vectorized
//' approach for consistency with the GKw family style.
//'
//' @export
// [[Rcpp::export]]
arma::vec pbeta_(
   const arma::vec& q,
   const Rcpp::NumericVector& gamma,
   const Rcpp::NumericVector& delta,
   bool lower_tail = true,
   bool log_p = false
) {
 arma::vec g_vec(gamma.begin(), gamma.size());
 arma::vec d_vec(delta.begin(), delta.size());

 size_t N = std::max({q.n_elem, g_vec.n_elem, d_vec.n_elem});
 arma::vec out(N);

 for (size_t i = 0; i < N; i++) {
   double g = g_vec[i % g_vec.n_elem];
   double d = d_vec[i % d_vec.n_elem];
   double qq = q[i % q.n_elem];

   if (!check_beta_pars(g, d)) {
     out(i) = NA_REAL;
     continue;
   }

   // boundary
   if (!R_finite(qq) || qq <= 0.0) {
     double v0 = lower_tail ? 0.0 : 1.0;
     out(i) = (log_p ? std::log(v0) : v0);
     continue;
   }
   if (qq >= 1.0) {
     double v1 = lower_tail ? 1.0 : 0.0;
     out(i) = (log_p ? std::log(v1) : v1);
     continue;
   }

   // call R's pbeta with adjusted parameters for GKw-style Beta
   double val = R::pbeta(qq, g, d + 1.0, lower_tail, false); // not log
   if (log_p) {
     val = std::log(val);
   }
   out(i) = val;
 }

 return out;
}


// -----------------------------------------------------------------------------
// 3) qbeta_: Quantile function for Beta
// -----------------------------------------------------------------------------

//' @title Quantile Function of the Beta Distribution
//'
//' @description
//' For p in (0,1), returns x = Q(p) = qbeta(p, gamma, delta+1).
//'
//' @param p Vector of probabilities in (0,1) or log scale if log_p=TRUE.
//' @param gamma First shape parameter > 0.
//' @param delta Second shape parameter > 0.
//' @param lower_tail If TRUE, p=F(x). If FALSE, p=1-F(x).
//' @param log_p If TRUE, p is log(p).
//'
//' @return A vector of quantiles, same length as broadcast.
//'
//' @details
//' Wrapper around R's \code{qbeta} with parameter adjustment for GKw-style Beta.
//' Uses vectorized approach for consistency with GKw family style.
//' We handle boundaries (p=0 => Q=0, p=1 => Q=1) manually.
//'
//' @export
// [[Rcpp::export]]
arma::vec qbeta_(
   const arma::vec& p,
   const Rcpp::NumericVector& gamma,
   const Rcpp::NumericVector& delta,
   bool lower_tail = true,
   bool log_p = false
) {
 arma::vec g_vec(gamma.begin(), gamma.size());
 arma::vec d_vec(delta.begin(), delta.size());

 size_t N = std::max({p.n_elem, g_vec.n_elem, d_vec.n_elem});
 arma::vec out(N);

 for (size_t i = 0; i < N; i++) {
   double g = g_vec[i % g_vec.n_elem];
   double d = d_vec[i % d_vec.n_elem];
   double pp = p[i % p.n_elem];

   if (!check_beta_pars(g, d)) {
     out(i) = NA_REAL;
     continue;
   }

   // handle log_p
   if (log_p) {
     if (pp > 0.0) {
       // => p>1
       out(i) = NA_REAL;
       continue;
     }
     pp = std::exp(pp);
   }
   // handle lower_tail
   if (!lower_tail) {
     pp = 1.0 - pp;
   }

   // boundaries
   if (pp <= 0.0) {
     out(i) = 0.0;
     continue;
   }
   if (pp >= 1.0) {
     out(i) = 1.0;
     continue;
   }

   // Use adjusted parameters for GKw-style Beta
   double val = R::qbeta(pp, g, d + 1.0, true, false); // returns not log
   out(i) = val;
 }

 return out;
}


// -----------------------------------------------------------------------------
// 4) rbeta_: RNG for Beta distribution
// -----------------------------------------------------------------------------

//' @title Random Generation for the Beta Distribution
//'
//' @description
//' Generates n samples from Beta(\code{gamma, delta}), wrapping R's \code{rbeta} with
//' vectorized style to match GKw family approach.
//'
//' @param n Number of samples.
//' @param gamma First shape parameter > 0.
//' @param delta Second shape parameter > 0.
//'
//' @return A numeric vector of length n, each in (0,1).
//'
//' @details
//' This function generates samples from the GKw family parameterization of the Beta
//' distribution. The parameters are adjusted to call R's rbeta function correctly.
//'
//' @export
// [[Rcpp::export]]
arma::vec rbeta_(
   int n,
   const Rcpp::NumericVector& gamma,
   const Rcpp::NumericVector& delta
) {
 if (n <= 0) {
   Rcpp::stop("rbeta_: n must be positive");
 }

 arma::vec g_vec(gamma.begin(), gamma.size());
 arma::vec d_vec(delta.begin(), delta.size());

 size_t k = std::max({g_vec.n_elem, d_vec.n_elem});
 arma::vec out(n);

 for (int i = 0; i < n; i++) {
   size_t idx = i % k;
   double g = g_vec[idx % g_vec.n_elem];
   double d = d_vec[idx % d_vec.n_elem];

   if (!check_beta_pars(g, d)) {
     out(i) = NA_REAL;
     Rcpp::warning("rbeta_: invalid parameters at index %d", i+1);
     continue;
   }

   // Use adjusted parameters for GKw-style Beta
   double val = R::rbeta(g, d + 1.0);
   out(i) = val;
 }
 return out;
}


// -----------------------------------------------------------------------------
// 5) llbeta: Negative Log-Likelihood for Beta
// -----------------------------------------------------------------------------

//' @title Negative Log-Likelihood for the Beta Distribution
//'
//' @description
//' Computes \eqn{-\ell(\theta)} for Beta(\code{gamma, delta}), given data x in (0,1).
//'
//' @param par NumericVector of length 2: (gamma, delta).
//' @param data NumericVector of observations in (0,1).
//'
//' @return A single numeric value. +Inf if invalid parameters or data out of domain.
//'
//' @details
//' If \eqn{f(x)= x^{\gamma-1}(1-x)^{\delta}/B(\gamma,\delta+1)}, the log-likelihood is
//' \eqn{\sum_i [(\gamma-1)\ln x_i + \delta\ln(1-x_i)] - n \ln B(\gamma,\delta+1)}.
//' We return its negative for optimization usage.
//'
//' @examples
//' \dontrun{
//' # Generate sample data from a Beta distribution
//' set.seed(123)
//' x <- rbeta_(100, 2, 3)
//' hist(x, breaks = 20, main = "Beta(2, 3) Sample")
//'
//' # Use in optimization with Hessian-based methods
//' result <- optim(c(0.5, 0.5), llbeta, method = "BFGS",
//'                 hessian = TRUE, data = x)
//'
//' # Compare numerical and analytical derivatives
//' num_grad <- numDeriv::grad(llbeta, x = result$par, data = x)
//' num_hess <- numDeriv::hessian(llbeta, x = result$par, data = x)
//'
//' ana_grad <- grbeta(result$par, data = x)
//' ana_hess <- hsbeta(result$par, data = x)
//'
//' # Check differences (should be very small)
//' round(num_grad - ana_grad, 4)
//' round(num_hess - ana_hess, 4)
//' }
//'
//' @export
// [[Rcpp::export]]
double llbeta(const Rcpp::NumericVector& par,
             const Rcpp::NumericVector& data) {
 if (par.size() < 2) {
   return R_PosInf;
 }
 double gamma = par[0]; // gamma > 0
 double delta = par[1]; // delta > 0

 if (!check_beta_pars(gamma, delta)) {
   return R_PosInf;
 }

 arma::vec x = Rcpp::as<arma::vec>(data);
 if (x.n_elem < 1) {
   return R_PosInf;
 }
 // domain check
 if (arma::any(x <= 0.0) || arma::any(x >= 1.0)) {
   return R_PosInf;
 }

 int n = x.n_elem;
 // Use correct parameterization for GKw-style Beta
 double logB = R::lbeta(gamma, delta + 1.0);
 // constant => -n*logB
 double cst = -double(n) * logB;

 // sum((gamma-1)*log(x_i) + delta*log(1-x_i)), i=1..n
 arma::vec lx = arma::log(x);
 arma::vec l1mx = arma::log(1.0 - x);

 double sum1 = (gamma - 1.0) * arma::sum(lx);
 double sum2 = delta * arma::sum(l1mx);  // Corrected: no subtraction of 1.0

 double loglike = cst + sum1 + sum2; // that's the log-likelihood

 // We must return negative
 return -loglike;
}


//' @title Gradient Function for Beta Log-Likelihood
//'
//' @description
//' Calculates the gradient vector (partial derivatives) of the negative log-likelihood
//' function for the Beta distribution. This function provides the exact gradient needed
//' for efficient optimization in maximum likelihood estimation.
//' The Beta is a submodel of GKw with α = 1, β = 1, and λ = 1 fixed.
//'
//' @param par NumericVector of length 2 containing parameters (γ, δ) in that order.
//'        All parameters must be positive.
//' @param data NumericVector of observations, where all values must be in the open interval (0,1).
//'
//' @return NumericVector of length 2 containing the gradient components (partial derivatives)
//'         of the negative log-likelihood with respect to each parameter (γ, δ).
//'         Returns a vector of NaN values if any parameters or data values are invalid.
//'
//' @details
//' The gradient vector contains the following partial derivatives of the negative log-likelihood:
//'
//' \deqn{
//' \frac{\partial \ell}{\partial \gamma} = n[\psi(\gamma+\delta+1) - \psi(\gamma)] -
//' \sum_{i=1}^{n}\log(x_i)
//' }
//'
//' \deqn{
//' \frac{\partial \ell}{\partial \delta} = n[\psi(\gamma+\delta+1) - \psi(\delta+1)] -
//' \sum_{i=1}^{n}\log(1-x_i)
//' }
//'
//' where:
//' \itemize{
//'   \item \deqn{\psi} is the digamma function (derivative of the log-gamma function)
//' }
//'
//' Note that this implementation works with the GKw family parameterization of the Beta
//' distribution where shape1=γ and shape2=δ+1 in standard statistical notation.
//'
//' The implementation includes several numerical safeguards:
//' \itemize{
//'   \item Parameter and data validation with appropriate error handling
//'   \item Efficient vector operations using Armadillo C++ library
//' }
//'
//' The returned gradient is for the negative log-likelihood to align with minimization
//' in optimization routines.
//'
//' @examples
//' \dontrun{
//' # Generate sample data from a Beta distribution
//' set.seed(123)
//' x <- rbeta_(100, 2, 3)
//' hist(x, breaks = 20, main = "Beta(2, 3) Sample")
//'
//' # Use in optimization with Hessian-based methods
//' result <- optim(c(0.5, 0.5), llbeta, method = "BFGS",
//'                 hessian = TRUE, data = x)
//'
//' # Compare numerical and analytical derivatives
//' num_grad <- numDeriv::grad(llbeta, x = result$par, data = x)
//' ana_grad <- grbeta(result$par, data = x)
//' round(num_grad - ana_grad, 4)
//' }
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector grbeta(const Rcpp::NumericVector& par,
                          const Rcpp::NumericVector& data) {
 Rcpp::NumericVector grad(2, R_PosInf); // initialize with Inf

 if (par.size() < 2) {
   return grad;
 }

 double gamma = par[0];
 double delta = par[1];

 if (!check_beta_pars(gamma, delta)) {
   return grad;
 }

 arma::vec x = Rcpp::as<arma::vec>(data);
 if (x.n_elem < 1) {
   return grad;
 }

 // domain check
 if (arma::any(x <= 0.0) || arma::any(x >= 1.0)) {
   return grad;
 }

 int n = x.n_elem;

 // Calculate digamma terms with correct parameterization for GKw-style Beta
 double dig_g = R::digamma(gamma);
 double dig_d = R::digamma(delta + 1.0);  // Corrected: digamma(δ+1)
 double dig_gd = R::digamma(gamma + delta + 1.0);  // Corrected: digamma(γ+δ+1)

 // Sum of log terms
 arma::vec lx = arma::log(x);
 arma::vec l1mx = arma::log(1.0 - x);
 double sum_lx = arma::sum(lx);
 double sum_l1mx = arma::sum(l1mx);

 // Partial derivatives for negative log-likelihood
 grad[0] = n * (dig_g - dig_gd) - sum_lx; // wrt gamma
 grad[1] = n * (dig_d - dig_gd) - sum_l1mx; // wrt delta

 return grad; // Already negated for negative log-likelihood
}


//' @title Hessian Matrix Function for Beta Log-Likelihood
//'
//' @description
//' Calculates the Hessian matrix (second-order partial derivatives) of the negative log-likelihood
//' function for the Beta distribution. This function provides the exact Hessian needed for
//' efficient optimization in maximum likelihood estimation and for asymptotic inference.
//' The Beta is a submodel of GKw with α = 1, β = 1, and λ = 1 fixed.
//'
//' @param par NumericVector of length 2 containing parameters (γ, δ) in that order.
//'        All parameters must be positive.
//' @param data NumericVector of observations, where all values must be in the open interval (0,1).
//'
//' @return NumericMatrix of dimensions 2x2 containing the Hessian components (second derivatives)
//'         of the negative log-likelihood with respect to each parameter pair (γ, δ).
//'         Returns a matrix of NaN values if any parameters or data values are invalid.
//'
//' @details
//' The Hessian matrix contains the following second derivatives of the negative log-likelihood:
//'
//' \deqn{
//' \frac{\partial^2 \ell}{\partial \gamma^2} = n[\psi'(\gamma) - \psi'(\gamma+\delta+1)]
//' }
//'
//' \deqn{
//' \frac{\partial^2 \ell}{\partial \gamma \partial \delta} = -n\psi'(\gamma+\delta+1)
//' }
//'
//' \deqn{
//' \frac{\partial^2 \ell}{\partial \delta^2} = n[\psi'(\delta+1) - \psi'(\gamma+\delta+1)]
//' }
//'
//' where:
//' \itemize{
//'   \item \deqn{\psi'} is the trigamma function (second derivative of the log-gamma function)
//' }
//'
//' Note that this implementation works with the GKw family parameterization of the Beta
//' distribution where shape1=γ and shape2=δ+1 in standard statistical notation.
//'
//' The implementation includes several numerical safeguards:
//' \itemize{
//'   \item Parameter and data validation with appropriate error handling
//'   \item Symmetry of the Hessian matrix is guaranteed by construction
//' }
//'
//' The returned Hessian is for the negative log-likelihood to align with minimization
//' in optimization routines.
//'
//' @examples
//' \dontrun{
//' # Generate sample data from a Beta distribution
//' set.seed(123)
//' x <- rbeta_(100, 2, 3)
//' hist(x, breaks = 20, main = "Beta(2, 3) Sample")
//'
//' # Use in optimization with Hessian-based methods
//' result <- optim(c(0.5, 0.5), llbeta, method = "BFGS",
//'                 hessian = TRUE, data = x)
//'
//' # Compare numerical and analytical derivatives
//' num_hess <- numDeriv::hessian(llbeta, x = result$par, data = x)
//' ana_hess <- hsbeta(result$par, data = x)
//' round(num_hess - ana_hess, 4)
//' }
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix hsbeta(const Rcpp::NumericVector& par,
                          const Rcpp::NumericVector& data) {
 Rcpp::NumericMatrix hess(2, 2);
 // Initialize with Inf
 for (int i = 0; i < 2; i++) {
   for (int j = 0; j < 2; j++) {
     hess(i, j) = R_PosInf;
   }
 }

 if (par.size() < 2) {
   return hess;
 }

 double gamma = par[0];
 double delta = par[1];

 if (!check_beta_pars(gamma, delta)) {
   return hess;
 }

 arma::vec x = Rcpp::as<arma::vec>(data);
 if (x.n_elem < 1) {
   return hess;
 }

 // domain check
 if (arma::any(x <= 0.0) || arma::any(x >= 1.0)) {
   return hess;
 }

 int n = x.n_elem;

 // Calculate trigamma terms with correct parameterization for GKw-style Beta
 double trig_g = R::trigamma(gamma);
 double trig_d = R::trigamma(delta + 1.0);  // Corrected: trigamma(δ+1)
 double trig_gd = R::trigamma(gamma + delta + 1.0);  // Corrected: trigamma(γ+δ+1)

 // Second derivatives for negative log-likelihood
 hess(0, 0) = n * (trig_g - trig_gd);  // d²/dγ²
 hess(1, 1) = n * (trig_d - trig_gd);  // d²/dδ²
 hess(0, 1) = hess(1, 0) = -n * trig_gd;  // d²/dγdδ

 return hess; // Already for negative log-likelihood
}


// Função para calcular a Hessiana numericamente usando Armadillo
arma::mat numHess(
  NumericVector params,
  NumericVector data,
  double eps = 1e-6
) {
int n = params.size();
arma::mat hessian(n, n, arma::fill::zeros);

// Valor da função no ponto atual
double f0 = llgkw(params, data);

// Para cada par de variáveis
for (int i = 0; i < n; i++) {
  for (int j = 0; j <= i; j++) {
    // Cálculo da derivada segunda usando diferenças finitas
    if (i == j) {
      // Derivada segunda em relação à mesma variável
      NumericVector params_p = clone(params);
      NumericVector params_m = clone(params);

      params_p[i] += eps;
      params_m[i] -= eps;

      double f_p = llgkw(params_p, data);
      double f_m = llgkw(params_m, data);

      hessian(i, j) = (f_p - 2*f0 + f_m) / (eps * eps);
    } else {
      // Derivada cruzada
      NumericVector params_pp = clone(params);
      NumericVector params_pm = clone(params);
      NumericVector params_mp = clone(params);
      NumericVector params_mm = clone(params);

      params_pp[i] += eps;
      params_pp[j] += eps;

      params_pm[i] += eps;
      params_pm[j] -= eps;

      params_mp[i] -= eps;
      params_mp[j] += eps;

      params_mm[i] -= eps;
      params_mm[j] -= eps;

      double f_pp = llgkw(params_pp, data);
      double f_pm = llgkw(params_pm, data);
      double f_mp = llgkw(params_mp, data);
      double f_mm = llgkw(params_mm, data);

      hessian(i, j) = hessian(j, i) = (f_pp - f_pm - f_mp + f_mm) / (4 * eps * eps);
    }
  }
}

return hessian;
}



//' @title Newton-Raphson Optimization for Kumaraswamy Family Distributions
//'
//' @description
//' Performs maximum likelihood estimation for the parameters of any distribution in the
//' Generalized Kumaraswamy (GKw) family using a robust implementation of the Newton-Raphson
//' algorithm. This function supports all 7 submodels: GKw, BKw, KKw, EKw, Mc (McDonald),
//' Kw, and Beta.
//'
//' @details
//' The Generalized Kumaraswamy family includes the following distributions, all defined on (0,1):
//' \itemize{
//'   \item \strong{GKw} (Generalized Kumaraswamy): 5 parameters (α, β, γ, δ, λ)
//'   \item \strong{BKw} (Beta-Kumaraswamy): 4 parameters (α, β, γ, δ), with λ = 1 fixed
//'   \item \strong{KKw} (Kumaraswamy-Kumaraswamy): 4 parameters (α, β, δ, λ), with γ = 1 fixed
//'   \item \strong{EKw} (Exponentiated Kumaraswamy): 3 parameters (α, β, λ), with γ = 1, δ = 0 fixed
//'   \item \strong{Mc} (McDonald/Beta Power): 3 parameters (γ, δ, λ), with α = 1, β = 1 fixed
//'   \item \strong{Kw} (Kumaraswamy): 2 parameters (α, β), with γ = 1, δ = 0, λ = 1 fixed
//'   \item \strong{Beta}: 2 parameters (γ, δ), with α = 1, β = 1, λ = 1 fixed
//' }
//'
//' This function implements a sophisticated optimization procedure to find the maximum likelihood
//' estimates with multiple fallback strategies to handle numerical challenges:
//' 1. Cholesky decomposition (fastest, requires positive-definite Hessian)
//' 2. Standard matrix solver
//' 3. Regularized Hessian with incremental adjustment
//' 4. Pseudo-inverse for highly ill-conditioned matrices
//' 5. Gradient descent as ultimate fallback
//'
//' The function also implements backtracking line search to ensure monotonic improvement in the
//' log-likelihood, with random parameter perturbation as a recovery strategy when backtracking fails.
//'
//' @param start A numeric vector containing initial values for parameters, with length
//'        corresponding to the selected family (see Details)
//' @param data A numeric vector containing observed data. All values must be in the interval (0,1)
//' @param family A character string specifying the distribution family. One of "gkw", "bkw", "kkw",
//'        "ekw", "mc", "kw", or "beta". Default: "gkw"
//' @param tol Convergence tolerance. The algorithm stops when the gradient norm or parameter/likelihood
//'        changes are below this value. Default: 1e-6
//' @param max_iter Maximum number of iterations. Default: 100
//' @param verbose Logical flag to print detailed progress information. Default: FALSE
//' @param use_hessian Logical flag to use Hessian information for parameter updates. If FALSE,
//'        the algorithm uses gradient descent instead. Default: TRUE
//' @param step_size Initial step size for parameter updates. Default: 1.0
//' @param enforce_bounds Logical flag to enforce parameter constraints. Default: TRUE
//' @param min_param_val Minimum allowed value for parameters (except delta). Default: 1e-5
//' @param max_param_val Maximum allowed value for parameters. Default: 1e5
//' @param get_num_hess Logical flag to calculate numerical Hessian in addition to analytical Hessian.
//'        Default: FALSE
//'
//' @return A list containing the following components:
//' \describe{
//'   \item{parameters}{A numeric vector with the estimated parameters}
//'   \item{loglik}{The maximized log-likelihood value}
//'   \item{iterations}{Number of iterations performed}
//'   \item{converged}{Logical flag indicating whether the algorithm converged}
//'   \item{param_history}{Matrix of parameter values at each iteration}
//'   \item{loglik_history}{Vector of log-likelihood values at each iteration}
//'   \item{gradient}{The gradient vector at the final parameter estimates}
//'   \item{hessian}{The Hessian matrix at the final parameter estimates}
//'   \item{std_errors}{Standard errors for the estimated parameters}
//'   \item{aic}{Akaike Information Criterion: AIC = 2k - 2ln(L)}
//'   \item{bic}{Bayesian Information Criterion: BIC = k ln(n) - 2ln(L)}
//'   \item{n}{Sample size}
//'   \item{status}{Character string indicating the termination status}
//'   \item{z_values}{Z-statistics for parameter significance tests}
//'   \item{p_values}{P-values for parameter significance tests}
//'   \item{param_names}{Character vector of parameter names}
//'   \item{family}{The distribution family used in the estimation}
//'   \item{numeric_hessian}{Numerical approximation of the Hessian (only if get_num_hess=TRUE)}
//' }
//'
//' @section Warning:
//' Convergence is not guaranteed for all datasets and initial values. It's recommended to:
//' \itemize{
//'   \item Try different initial values if convergence fails
//'   \item Check the gradient norm at the final solution to verify optimality
//'   \item Examine parameter history to identify potential issues
//'   \item Use the verbose option to get detailed progress information for troubleshooting
//' }
//'
//' @examples
//' \dontrun{
//' # Generate sample data from Beta(2,5) distribution for testing
//' set.seed(123)
//' sample_data <- rbeta(200, 2, 5)
//'
//' # Fit with full GKw model
//' gkw_result <- mle_fit(c(1.5, 4.5, 1.0, 0.0, 1.0), sample_data, family = "gkw")
//' gkw_result$parameters
//'
//' # Fit with simpler Kumaraswamy model
//' kw_result <- mle_fit(c(1.5, 4.5), sample_data, family = "kw")
//' kw_result$parameters
//'
//' # Fit with Beta model
//' beta_result <- mle_fit(c(2.0, 5.0), sample_data, family = "beta")
//' beta_result$parameters
//'
//' # Compare AIC/BIC values to select the best model
//' data.frame(
//'   family = c("gkw", "kw", "beta"),
//'   AIC = c(gkw_result$aic, kw_result$aic, beta_result$aic),
//'   BIC = c(gkw_result$bic, kw_result$bic, beta_result$bic)
//' )
//' }
//'
//' @references
//' Kumaraswamy, P. (1980). A generalized probability density function for double-bounded
//' random processes. Journal of Hydrology, 46(1-2), 79-88.
//'
//' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized distributions.
//' Journal of Statistical Computation and Simulation, 81(7), 883-898.
//'
//' Fletcher, R. (1987). Practical Methods of Optimization. John Wiley & Sons.
//'
//' @author Lopes, J. E.
//'
//' @export
// [[Rcpp::export]]
List nrgkw(
   NumericVector start,
   NumericVector data,
   std::string family = "gkw",
   double tol = 1e-6,
   int max_iter = 100,
   bool verbose = false,
   bool use_hessian = true,
   double step_size = 1.0,
   bool enforce_bounds = true,
   double min_param_val = 1e-5,
   double max_param_val = 1e5,
   bool get_num_hess = false
) {
 // Final result will be a list with different components
 List result;

 // Convert family to lowercase for case-insensitive comparison
 std::string family_lower = family;
 std::transform(family_lower.begin(), family_lower.end(), family_lower.begin(), ::tolower);

 // Determine number of parameters based on family
 int n_params = 0;
 CharacterVector param_names;

 if (family_lower == "gkw") {
   n_params = 5;
   param_names = CharacterVector::create("alpha", "beta", "gamma", "delta", "lambda");
 } else if (family_lower == "bkw") {
   n_params = 4;
   param_names = CharacterVector::create("alpha", "beta", "gamma", "delta");
 } else if (family_lower == "kkw") {
   n_params = 4;
   param_names = CharacterVector::create("alpha", "beta", "delta", "lambda");
 } else if (family_lower == "ekw") {
   n_params = 3;
   param_names = CharacterVector::create("alpha", "beta", "lambda");
 } else if (family_lower == "mc" || family_lower == "mcdonald" || family_lower == "bp") {
   n_params = 3;
   param_names = CharacterVector::create("gamma", "delta", "lambda");
 } else if (family_lower == "kw") {
   n_params = 2;
   param_names = CharacterVector::create("alpha", "beta");
 } else if (family_lower == "beta") {
   n_params = 2;
   param_names = CharacterVector::create("gamma", "delta");
 } else {
   stop("Unknown family: '" + family + "'. Available options are 'gkw', 'bkw', 'kkw', 'ekw', 'mc', 'kw', 'beta'.");
 }

 // Validate initial parameters size
 if (start.size() != n_params) {
   stop("Invalid number of parameters for '" + family + "'. Expected " +
     std::to_string(n_params) + ", got " + std::to_string(start.size()));
 }

 // Check for valid data
 int n_data = data.size();
 if (n_data < n_params) {
   stop("At least " + std::to_string(n_params) + " data points are needed to estimate " +
     std::to_string(n_params) + " parameters");
 }

 // Check if all data are in the interval (0,1)
 for (int i = 0; i < n_data; i++) {
   if (data[i] <= 0.0 || data[i] >= 1.0 || !R_finite(data[i])) {
     stop("All data must be in the interval (0,1)");
   }
 }

 // Copy initial parameters and convert to standard GKw parameters where needed
 NumericVector params(5); // Always use 5 parameters internally (GKw format)

 // Set default values based on fixed parameters in specific families
 params[0] = 1.0; // α = 1 default
 params[1] = 1.0; // β = 1 default
 params[2] = 1.0; // γ = 1 default
 params[3] = 0.0; // δ = 0 default
 params[4] = 1.0; // λ = 1 default

 // Fill with provided parameters based on family
 if (family_lower == "gkw") {
   for (int j = 0; j < 5; j++) {
     params[j] = start[j];
   }
 } else if (family_lower == "bkw") {
   // α, β, γ, δ with λ = 1
   for (int j = 0; j < 4; j++) {
     params[j] = start[j];
   }
   params[4] = 1.0; // λ fixed at 1
 } else if (family_lower == "kkw") {
   // α, β, δ, λ with γ = 1
   params[0] = start[0]; // α
   params[1] = start[1]; // β
   params[2] = 1.0;             // γ fixed at 1
   params[3] = start[2]; // δ
   params[4] = start[3]; // λ
 } else if (family_lower == "ekw") {
   // α, β, λ with γ = 1, δ = 0
   params[0] = start[0]; // α
   params[1] = start[1]; // β
   params[2] = 1.0;             // γ fixed at 1
   params[3] = 0.0;             // δ fixed at 0
   params[4] = start[2]; // λ
 } else if (family_lower == "mc" || family_lower == "mcdonald" || family_lower == "bp") {
   // γ, δ, λ with α = 1, β = 1
   params[0] = 1.0;             // α fixed at 1
   params[1] = 1.0;             // β fixed at 1
   params[2] = start[0]; // γ
   params[3] = start[1]; // δ
   params[4] = start[2]; // λ
 } else if (family_lower == "kw") {
   // α, β with γ = 1, δ = 0, λ = 1
   params[0] = start[0]; // α
   params[1] = start[1]; // β
   params[2] = 1.0;             // γ fixed at 1
   params[3] = 0.0;             // δ fixed at 0
   params[4] = 1.0;             // λ fixed at 1
 } else if (family_lower == "beta") {
   // γ, δ with α = 1, β = 1, λ = 1
   params[0] = 1.0;             // α fixed at 1
   params[1] = 1.0;             // β fixed at 1
   params[2] = start[0]; // γ
   params[3] = start[1]; // δ
   params[4] = 1.0;             // λ fixed at 1
 }

 // Apply constraints to initial parameters if needed
 if (enforce_bounds) {
   for (int j = 0; j < 5; j++) {
     if (j == 3) { // delta
       params[j] = std::max(0.0, params[j]);
     } else { // other parameters must be > 0
       params[j] = std::max(min_param_val, params[j]);
     }
     params[j] = std::min(max_param_val, params[j]);
   }
 }

 // Define function pointers based on family
 std::function<double(NumericVector, NumericVector)> ll_func;
 std::function<NumericVector(NumericVector, NumericVector)> gr_func;
 std::function<NumericMatrix(NumericVector, NumericVector)> hs_func;

 // Assign appropriate functions based on family
 if (family_lower == "gkw") {
   ll_func = llgkw;
   gr_func = grgkw;
   hs_func = hsgkw;
 } else if (family_lower == "bkw") {
   ll_func = llbkw;
   gr_func = grbkw;
   hs_func = hsbkw;
 } else if (family_lower == "kkw") {
   ll_func = llkkw;
   gr_func = grkkw;
   hs_func = hskkw;
 } else if (family_lower == "ekw") {
   ll_func = llekw;
   gr_func = grekw;
   hs_func = hsekw;
 } else if (family_lower == "mc" || family_lower == "mcdonald" || family_lower == "bp") {
   ll_func = llmc;
   gr_func = grmc;
   hs_func = hsmc;
 } else if (family_lower == "kw") {
   ll_func = llkw;
   gr_func = grkw;
   hs_func = hskw;
 } else if (family_lower == "beta") {
   ll_func = llbeta;
   gr_func = grbeta;
   hs_func = hsbeta;
 }

 // Function to extract relevant parameters for specific family
 auto extract_params = [&](const NumericVector& full_params) -> NumericVector {
   NumericVector result;

   if (family_lower == "gkw") {
     result = NumericVector(5);
     for (int j = 0; j < 5; j++) result[j] = full_params[j];
   } else if (family_lower == "bkw") {
     result = NumericVector(4);
     for (int j = 0; j < 4; j++) result[j] = full_params[j];
   } else if (family_lower == "kkw") {
     result = NumericVector(4);
     result[0] = full_params[0]; // α
     result[1] = full_params[1]; // β
     result[2] = full_params[3]; // δ
     result[3] = full_params[4]; // λ
   } else if (family_lower == "ekw") {
     result = NumericVector(3);
     result[0] = full_params[0]; // α
     result[1] = full_params[1]; // β
     result[2] = full_params[4]; // λ
   } else if (family_lower == "mc" || family_lower == "mcdonald" || family_lower == "bp") {
     result = NumericVector(3);
     result[0] = full_params[2]; // γ
     result[1] = full_params[3]; // δ
     result[2] = full_params[4]; // λ
   } else if (family_lower == "kw") {
     result = NumericVector(2);
     result[0] = full_params[0]; // α
     result[1] = full_params[1]; // β
   } else if (family_lower == "beta") {
     result = NumericVector(2);
     result[0] = full_params[2]; // γ
     result[1] = full_params[3]; // δ
   }

   return result;
 };

 // Create a custom numHess function to handle specific families
 auto numHess_family = [&](NumericVector params_family, NumericVector data_family, double eps = 1e-6) {
   int n_params_family = params_family.size();
   arma::mat hessian(n_params_family, n_params_family, arma::fill::zeros);

   // Value of the function at the current point
   double f0 = ll_func(params_family, data_family);

   // For each pair of variables
   for (int i = 0; i < n_params_family; i++) {
     for (int j = 0; j <= i; j++) {
       // Calculate second derivative using finite differences
       if (i == j) {
         // Second derivative with respect to the same variable
         NumericVector params_p = clone(params_family);
         NumericVector params_m = clone(params_family);

         params_p[i] += eps;
         params_m[i] -= eps;

         double f_p = ll_func(params_p, data_family);
         double f_m = ll_func(params_m, data_family);

         hessian(i, j) = (f_p - 2*f0 + f_m) / (eps * eps);
       } else {
         // Mixed derivative
         NumericVector params_pp = clone(params_family);
         NumericVector params_pm = clone(params_family);
         NumericVector params_mp = clone(params_family);
         NumericVector params_mm = clone(params_family);

         params_pp[i] += eps;
         params_pp[j] += eps;

         params_pm[i] += eps;
         params_pm[j] -= eps;

         params_mp[i] -= eps;
         params_mp[j] += eps;

         params_mm[i] -= eps;
         params_mm[j] -= eps;

         double f_pp = ll_func(params_pp, data_family);
         double f_pm = ll_func(params_pm, data_family);
         double f_mp = ll_func(params_mp, data_family);
         double f_mm = ll_func(params_mm, data_family);

         hessian(i, j) = hessian(j, i) = (f_pp - f_pm - f_mp + f_mm) / (4 * eps * eps);
       }
     }
   }

   return hessian;
 };

 // Get family-specific parameters
 NumericVector family_params = extract_params(params);

 // Calculate initial log-likelihood
 double initial_loglik = ll_func(family_params, data);
 if (!R_finite(initial_loglik) || initial_loglik == R_PosInf) {
   stop("Initial log-likelihood is infinite or NaN. Check the initial parameters.");
 }

 // Parameter and log-likelihood history for diagnostics
 NumericMatrix param_history(max_iter + 1, n_params);
 NumericVector loglik_history(max_iter + 1);

 // Initialize history with initial values
 for (int j = 0; j < n_params; j++) {
   param_history(0, j) = family_params[j];
 }
 loglik_history[0] = initial_loglik;

 // Variables for convergence control
 bool converged = false;
 int iter = 0;
 double prev_loglik = initial_loglik;

 // Prepare to store the best result obtained
 double best_loglik = initial_loglik;
 NumericVector best_params = clone(family_params);

 // Main Newton-Raphson loop
 while (!converged && iter < max_iter) {
   iter++;

   // Calculate log-likelihood, gradient and hessian
   double current_loglik = ll_func(family_params, data);
   NumericVector gradient = gr_func(family_params, data);

   // Check if gradient has valid values
   bool valid_gradient = true;
   for (int j = 0; j < n_params; j++) {
     if (!R_finite(gradient[j])) {
       valid_gradient = false;
       break;
     }
   }

   if (!valid_gradient) {
     if (verbose) {
       Rcout << "Warning: Invalid gradient in iteration " << iter << std::endl;
     }
     result["converged"] = false;
     result["status"] = "gradient_failure";
     // Use the best parameters found so far
     family_params = best_params;
     break;
   }

   // Calculate gradient norm for stopping criterion
   double grad_norm = 0.0;
   for (int j = 0; j < n_params; j++) {
     grad_norm += gradient[j] * gradient[j];
   }
   grad_norm = std::sqrt(grad_norm);

   if (grad_norm < tol) {
     converged = true;
     if (verbose) {
       Rcout << "Convergence detected: gradient norm (" << grad_norm
             << ") < tolerance (" << tol << ")" << std::endl;
     }
     break;
   }

   // Update direction
   NumericVector update(n_params);

   if (use_hessian) {
     // Calculate the Hessian
     NumericMatrix rcpp_hessian = hs_func(family_params, data);

     // Check if Hessian has valid values
     bool valid_hessian = true;
     for (int i = 0; i < n_params; i++) {
       for (int j = 0; j < n_params; j++) {
         if (!R_finite(rcpp_hessian(i, j))) {
           valid_hessian = false;
           break;
         }
       }
       if (!valid_hessian) break;
     }

     if (!valid_hessian) {
       if (verbose) {
         Rcout << "Warning: Invalid Hessian in iteration " << iter
               << ", using only the gradient." << std::endl;
       }
       // Fallback to steepest descent if Hessian is invalid
       for (int j = 0; j < n_params; j++) {
         // Normalize gradient for step control
         update[j] = -step_size * gradient[j] / std::max(1.0, std::abs(gradient[j]));
       }
     } else {
       // Convert to arma::mat for more robust matrix operations
       arma::mat hessian = as<arma::mat>(rcpp_hessian);
       arma::vec grad_vec = as<arma::vec>(gradient);
       arma::vec neg_grad = -grad_vec;
       arma::vec update_vec;

       bool solve_success = false;

       // Try 1: Cholesky for symmetric positive definite matrices (fastest)
       try {
         update_vec = arma::solve(hessian, neg_grad, arma::solve_opts::likely_sympd);
         solve_success = true;

         if (verbose) {
           Rcout << "Cholesky decomposition successful for parameter update." << std::endl;
         }
       } catch (...) {
         if (verbose) {
           Rcout << "Warning: Cholesky decomposition failed, trying standard solver..." << std::endl;
         }

         // Try 2: Standard Armadillo solver
         try {
           update_vec = arma::solve(hessian, neg_grad);
           solve_success = true;

           if (verbose) {
             Rcout << "Standard solver successful for parameter update." << std::endl;
           }
         } catch (...) {
           if (verbose) {
             Rcout << "Warning: Standard solver failed, trying with regularization..." << std::endl;
           }

           // Try 3: Regularize the Hessian matrix
           arma::mat reg_hessian = hessian;
           double reg_factor = 1e-6;

           // Find reasonable magnitude for regularization
           double diag_max = arma::max(arma::abs(reg_hessian.diag()));
           reg_factor = std::max(reg_factor, 1e-6 * diag_max);

           // Add small value to diagonal
           reg_hessian.diag() += reg_factor;

           try {
             update_vec = arma::solve(reg_hessian, neg_grad);
             solve_success = true;

             if (verbose) {
               Rcout << "Regularized solver successful with factor: " << reg_factor << std::endl;
             }
           } catch (...) {
             // Try 4: Stronger regularization
             reg_hessian = hessian;
             reg_factor = 1e-4 * (1.0 + diag_max);
             reg_hessian.diag() += reg_factor;

             try {
               update_vec = arma::solve(reg_hessian, neg_grad);
               solve_success = true;

               if (verbose) {
                 Rcout << "Stronger regularization successful with factor: " << reg_factor << std::endl;
               }
             } catch (...) {
               // Try 5: Pseudo-inverse (very robust method)
               try {
                 arma::mat hess_pinv = arma::pinv(hessian);
                 update_vec = hess_pinv * neg_grad;
                 solve_success = true;

                 if (verbose) {
                   Rcout << "Pseudo-inverse solution successful for parameter update." << std::endl;
                 }
               } catch (...) {
                 if (verbose) {
                   Rcout << "Warning: All matrix inversion methods failed in iteration " << iter
                         << ", using only the gradient." << std::endl;
                 }

                 // If all attempts fail, use gradient descent
                 for (int j = 0; j < n_params; j++) {
                   update[j] = -step_size * gradient[j] / std::max(1.0, std::abs(gradient[j]));
                 }

                 solve_success = false;
               }
             }
           }
         }
       }

       if (solve_success) {
         // Convert solution from arma::vec to NumericVector
         update = wrap(update_vec);

         // Limit step size to avoid too large steps
         double max_update = 0.0;
         for (int j = 0; j < n_params; j++) {
           max_update = std::max(max_update, std::abs(update[j]));
         }

         // If step is too large, reduce proportionally
         const double max_step = 2.0;
         if (max_update > max_step) {
           double scale_factor = max_step / max_update;
           for (int j = 0; j < n_params; j++) {
             update[j] *= scale_factor;
           }
         }

         // Apply step_size
         for (int j = 0; j < n_params; j++) {
           update[j] *= step_size;
         }
       }
     }
   } else {
     // Use only gradient (gradient descent method)
     for (int j = 0; j < n_params; j++) {
       update[j] = -step_size * gradient[j] / std::max(1.0, std::abs(gradient[j]));
     }
   }

   // Update parameters: theta_new = theta_old + update
   NumericVector new_params(n_params);
   for (int j = 0; j < n_params; j++) {
     new_params[j] = family_params[j] + update[j];
   }

   // Enforce bounds if requested
   if (enforce_bounds) {
     for (int j = 0; j < n_params; j++) {
       bool is_delta = (family_lower == "gkw" && j == 3) ||
         (family_lower == "bkw" && j == 3) ||
         (family_lower == "kkw" && j == 2) ||
         (family_lower == "mc" && j == 1) ||
         (family_lower == "beta" && j == 1);

       // Note: for delta, we allow values down to 0
       if (is_delta) {
         new_params[j] = std::max(0.0, new_params[j]);
       } else {
         new_params[j] = std::max(min_param_val, new_params[j]);
       }
       new_params[j] = std::min(max_param_val, new_params[j]);
     }
   }

   // Calculate new objective function value
   double new_loglik = ll_func(new_params, data);

   // Line search / Backtracking if new value is not better
   bool backtracking_success = true;
   double bt_step = 1.0;
   const double bt_factor = 0.5; // reduce step by half each backtracking
   const int max_bt = 10;        // maximum backtracking iterations

   if (!R_finite(new_loglik) || new_loglik >= current_loglik) {
     backtracking_success = false;

     if (verbose) {
       Rcout << "Starting backtracking at iteration " << iter
             << ", current value: " << current_loglik
             << ", new value: " << new_loglik << std::endl;
     }

     for (int bt = 0; bt < max_bt; bt++) {
       bt_step *= bt_factor;

       // Recalculate new parameters with reduced step
       for (int j = 0; j < n_params; j++) {
         new_params[j] = family_params[j] + bt_step * update[j];
       }

       // Enforce bounds again
       if (enforce_bounds) {
         for (int j = 0; j < n_params; j++) {
           bool is_delta = (family_lower == "gkw" && j == 3) ||
             (family_lower == "bkw" && j == 3) ||
             (family_lower == "kkw" && j == 2) ||
             (family_lower == "mc" && j == 1) ||
             (family_lower == "beta" && j == 1);

           if (is_delta) {
             new_params[j] = std::max(0.0, new_params[j]);
           } else {
             new_params[j] = std::max(min_param_val, new_params[j]);
           }
           new_params[j] = std::min(max_param_val, new_params[j]);
         }
       }

       // Test new value
       new_loglik = ll_func(new_params, data);

       if (R_finite(new_loglik) && new_loglik < current_loglik) {
         backtracking_success = true;
         if (verbose) {
           Rcout << "Backtracking successful after " << (bt + 1)
                 << " attempts, new value: " << new_loglik << std::endl;
         }
         break;
       }
     }
   }

   // If we still cannot improve, evaluate the situation
   if (!backtracking_success) {
     if (verbose) {
       Rcout << "Warning: Backtracking failed at iteration " << iter << std::endl;
     }

     // If gradient is small enough, consider converged
     if (grad_norm < tol * 10) {  // Relaxed tolerance for this case
       converged = true;
       if (verbose) {
         Rcout << "Convergence detected with small gradient after backtracking failure." << std::endl;
       }
     } else {
       // If backtracking fails and we're close to max iterations,
       // check if we're in a reasonable region
       if (iter > max_iter * 0.8 && current_loglik < best_loglik * 1.1) {
         converged = true;
         if (verbose) {
           Rcout << "Forced convergence after backtracking failure near maximum iterations." << std::endl;
         }
       } else {
         // Try a small random perturbation
         NumericVector perturb(n_params);
         for (int j = 0; j < n_params; j++) {
           // Perturbation of up to 5% of current value
           double range = 0.05 * std::abs(family_params[j]);
           perturb[j] = R::runif(-range, range);
           new_params[j] = family_params[j] + perturb[j];
         }

         // Apply constraints
         if (enforce_bounds) {
           for (int j = 0; j < n_params; j++) {
             bool is_delta = (family_lower == "gkw" && j == 3) ||
               (family_lower == "bkw" && j == 3) ||
               (family_lower == "kkw" && j == 2) ||
               (family_lower == "mc" && j == 1) ||
               (family_lower == "beta" && j == 1);

             if (is_delta) {
               new_params[j] = std::max(0.0, new_params[j]);
             } else {
               new_params[j] = std::max(min_param_val, new_params[j]);
             }
             new_params[j] = std::min(max_param_val, new_params[j]);
           }
         }

         new_loglik = ll_func(new_params, data);

         if (R_finite(new_loglik) && new_loglik < current_loglik) {
           backtracking_success = true;
           if (verbose) {
             Rcout << "Recovery by random perturbation, new value: " << new_loglik << std::endl;
           }
         } else {
           // If even perturbation doesn't work, use the best result so far
           new_params = best_params;
           new_loglik = best_loglik;
           if (verbose) {
             Rcout << "Returning to the previous best result: " << -best_loglik << std::endl;
           }
         }
       }
     }
   }

   // Update parameters and history
   for (int j = 0; j < n_params; j++) {
     family_params[j] = new_params[j];
     param_history(iter, j) = family_params[j];
   }
   loglik_history[iter] = new_loglik;

   // Update the best result if this is better
   if (new_loglik < best_loglik) {
     best_loglik = new_loglik;
     for (int j = 0; j < n_params; j++) {
       best_params[j] = family_params[j];
     }
   }

   // Check convergence by parameter change
   double param_change = 0.0;
   double param_rel_change = 0.0;
   for (int j = 0; j < n_params; j++) {
     param_change += std::pow(update[j], 2);
     if (std::abs(family_params[j]) > 1e-10) {
       param_rel_change += std::pow(update[j] / family_params[j], 2);
     } else {
       param_rel_change += std::pow(update[j], 2);
     }
   }
   param_change = std::sqrt(param_change);
   param_rel_change = std::sqrt(param_rel_change / n_params);

   // Check convergence by log-likelihood change
   double loglik_change = std::abs(prev_loglik - new_loglik);
   double loglik_rel_change = loglik_change / (std::abs(prev_loglik) + 1e-10);
   prev_loglik = new_loglik;

   if (verbose) {
     Rcout << "Iteration " << iter
           << ", Log-likelihood: " << -new_loglik
           << ", Change: " << loglik_change
           << ", Rel. Change: " << loglik_rel_change
           << ", Gradient Norm: " << grad_norm
           << std::endl;

     Rcout << "Parameters:";
     for (int j = 0; j < n_params; j++) {
       Rcout << " " << family_params[j];
     }
     Rcout << std::endl;
   }

   // Convergence criteria
   if (param_change < tol || param_rel_change < tol ||
       loglik_change < tol || loglik_rel_change < tol) {
     converged = true;
     if (verbose) {
       Rcout << "Convergence detected:" << std::endl;
       if (param_change < tol) Rcout << "- Absolute parameter change < tolerance" << std::endl;
       if (param_rel_change < tol) Rcout << "- Relative parameter change < tolerance" << std::endl;
       if (loglik_change < tol) Rcout << "- Absolute log-likelihood change < tolerance" << std::endl;
       if (loglik_rel_change < tol) Rcout << "- Relative log-likelihood change < tolerance" << std::endl;
     }
   }
 }

 // If not converged, use the best parameters found
 if (!converged) {
   family_params = best_params;
   if (verbose) {
     Rcout << "Did not fully converge, using the best parameters found." << std::endl;
   }
 }

 // Prepare final result
 NumericMatrix final_param_history(iter + 1, n_params);
 NumericVector final_loglik_history(iter + 1);

 for (int i = 0; i <= iter; i++) {
   for (int j = 0; j < n_params; j++) {
     final_param_history(i, j) = param_history(i, j);
   }
   final_loglik_history[i] = loglik_history[i];
 }

 // Calculate final gradient and hessian
 NumericVector final_gradient = gr_func(family_params, data);
 NumericMatrix rcpp_hessian = hs_func(family_params, data);

 // Calculate numerical Hessian if requested
 NumericMatrix rcpp_numeric_hessian;
 if (get_num_hess) {
   arma::mat arma_numeric_hessian = numHess_family(family_params, data);
   rcpp_numeric_hessian = wrap(arma_numeric_hessian);
   result["numeric_hessian"] = rcpp_numeric_hessian;
 }

 // Calculate standard errors using Armadillo for robust matrix inversion
 NumericVector std_errors(n_params, NA_REAL);
 bool valid_se = true;

 // Convert Rcpp Hessian to Armadillo matrix
 arma::mat hessian = as<arma::mat>(rcpp_hessian);
 arma::mat cov_matrix;

 // Layered approach to calculate covariance matrix (inverse of Hessian)
 try {
   // Step 1: Try to use Cholesky decomposition (fastest, requires positive definite)
   try {
     cov_matrix = arma::inv_sympd(hessian);

     if (verbose) {
       Rcout << "Standard error calculation: Cholesky decomposition successful." << std::endl;
     }
   } catch (...) {
     // Step 2: Try standard inversion
     try {
       cov_matrix = arma::inv(hessian);

       if (verbose) {
         Rcout << "Standard error calculation: Standard inverse successful." << std::endl;
       }
     } catch (...) {
       // Step 3: Apply regularization
       arma::mat reg_hessian = hessian;
       double reg_factor = 1e-6;

       // Find reasonable magnitude for regularization
       double diag_max = arma::max(arma::abs(reg_hessian.diag()));
       reg_factor = std::max(reg_factor, 1e-6 * diag_max);

       // Add small value to diagonal
       reg_hessian.diag() += reg_factor;

       try {
         cov_matrix = arma::inv(reg_hessian);

         if (verbose) {
           Rcout << "Standard error calculation: Regularized inverse successful." << std::endl;
         }
       } catch (...) {
         // Step 4: Stronger regularization
         reg_hessian = hessian;
         reg_factor = 1e-4 * (1.0 + diag_max);
         reg_hessian.diag() += reg_factor;

         try {
           cov_matrix = arma::inv(reg_hessian);

           if (verbose) {
             Rcout << "Standard error calculation: Stronger regularized inverse successful." << std::endl;
           }
         } catch (...) {
           // Step 5: Use pseudo-inverse (more robust)
           try {
             cov_matrix = arma::pinv(hessian);

             if (verbose) {
               Rcout << "Standard error calculation: Pseudo-inverse successful." << std::endl;
             }
           } catch (...) {
             // Step 6: Try numerical Hessian if available
             if (get_num_hess) {
               arma::mat num_hessian = as<arma::mat>(rcpp_numeric_hessian);

               try {
                 cov_matrix = arma::pinv(num_hessian);

                 if (verbose) {
                   Rcout << "Standard error calculation: Numerical Hessian pseudo-inverse successful." << std::endl;
                 }
               } catch (...) {
                 valid_se = false;
               }
             } else {
               valid_se = false;
             }
           }
         }
       }
     }
   }

   // Calculate standard errors if covariance matrix is available
   if (valid_se) {
     // Extract diagonal elements and calculate square root
     arma::vec diag_cov = cov_matrix.diag();

     for (int j = 0; j < n_params; j++) {
       if (diag_cov(j) > 0) {
         std_errors[j] = std::sqrt(diag_cov(j));
       } else {
         if (verbose) {
           Rcout << "Warning: Non-positive variance detected for parameter " << j
                 << ". Standard error set to NA." << std::endl;
         }
         std_errors[j] = NA_REAL;
       }
     }
   }
 } catch (...) {
   valid_se = false;
 }

 if (!valid_se && verbose) {
   Rcout << "Warning: Could not calculate standard errors. The Hessian matrix may not be positive definite." << std::endl;
 }

 // Calculate AIC: AIC = 2k - 2ln(L) = 2k + 2*(-ln(L))
 double final_loglik = ll_func(family_params, data);
 double aic = 2 * n_params + 2 * final_loglik;

 // Calculate BIC: BIC = k ln(n) - 2ln(L) = k ln(n) + 2*(-ln(L))
 double bic = n_params * std::log(n_data) + 2 * final_loglik;

 // Fill the result
 result["parameters"] = family_params;
 result["loglik"] = -final_loglik;  // Negative because ll functions return -logL
 result["iterations"] = iter;
 result["converged"] = converged;
 result["param_history"] = final_param_history;
 result["loglik_history"] = -final_loglik_history;  // Negative for consistency
 result["gradient"] = final_gradient;
 result["hessian"] = rcpp_hessian;
 result["std_errors"] = std_errors;
 result["aic"] = aic;
 result["bic"] = bic;
 result["n"] = n_data;
 result["family"] = family;

 if (!converged && !result.containsElementNamed("status")) {
   result["status"] = "max_iterations_reached";
 } else if (converged) {
   result["status"] = "success";
 }

 // Calculate statistical significance (p-values) using normal approximation
 NumericVector z_values(n_params, NA_REAL);
 NumericVector p_values(n_params, NA_REAL);

 if (valid_se) {
   for (int j = 0; j < n_params; j++) {
     if (std_errors[j] != NA_REAL && std_errors[j] > 0) {
       z_values[j] = family_params[j] / std_errors[j];
       // Two-tailed approximation
       p_values[j] = 2.0 * R::pnorm(-std::abs(z_values[j]), 0.0, 1.0, 1, 0);
     }
   }
   result["z_values"] = z_values;
   result["p_values"] = p_values;
 }

 // Set parameter names for easier interpretation
 colnames(final_param_history) = param_names;
 result["param_names"] = param_names;

 return result;
}




// //' @title Newton-Raphson Optimization for Kumaraswamy Family Distributions
// //'
// //' @description
// //' Performs maximum likelihood estimation for the parameters of any distribution in the
// //' Generalized Kumaraswamy (GKw) family using a robust implementation of the Newton-Raphson
// //' algorithm. This function supports all 7 submodels: GKw, BKw, KKw, EKw, Mc (McDonald),
// //' Kw, and Beta.
// //'
// //' @details
// //' The Generalized Kumaraswamy family includes the following distributions, all defined on (0,1):
// //' \itemize{
// //'   \item \strong{GKw} (Generalized Kumaraswamy): 5 parameters (α, β, γ, δ, λ)
// //'   \item \strong{BKw} (Beta-Kumaraswamy): 4 parameters (α, β, γ, δ), with λ = 1 fixed
// //'   \item \strong{KKw} (Kumaraswamy-Kumaraswamy): 4 parameters (α, β, δ, λ), with γ = 1 fixed
// //'   \item \strong{EKw} (Exponentiated Kumaraswamy): 3 parameters (α, β, λ), with γ = 1, δ = 0 fixed
// //'   \item \strong{Mc} (McDonald/Beta Power): 3 parameters (γ, δ, λ), with α = 1, β = 1 fixed
// //'   \item \strong{Kw} (Kumaraswamy): 2 parameters (α, β), with γ = 1, δ = 0, λ = 1 fixed
// //'   \item \strong{Beta}: 2 parameters (γ, δ), with α = 1, β = 1, λ = 1 fixed
// //' }
// //'
// //' This function implements a sophisticated optimization procedure to find the maximum likelihood
// //' estimates with multiple fallback strategies to handle numerical challenges:
// //' 1. Cholesky decomposition (fastest, requires positive-definite Hessian)
// //' 2. Standard matrix solver
// //' 3. Regularized Hessian with incremental adjustment
// //' 4. Pseudo-inverse for highly ill-conditioned matrices
// //' 5. Gradient descent as ultimate fallback
// //'
// //' The function also implements backtracking line search to ensure monotonic improvement in the
// //' log-likelihood, with random parameter perturbation as a recovery strategy when backtracking fails.
// //'
// //' @param start A numeric vector containing initial values for parameters, with length
// //'        corresponding to the selected family (see Details)
// //' @param data A numeric vector containing observed data. All values must be in the interval (0,1)
// //' @param family A character string specifying the distribution family. One of "gkw", "bkw", "kkw",
// //'        "ekw", "mc", "kw", or "beta". Default: "gkw"
// //' @param tol Convergence tolerance. The algorithm stops when the gradient norm or parameter/likelihood
// //'        changes are below this value. Default: 1e-6
// //' @param max_iter Maximum number of iterations. Default: 100
// //' @param verbose Logical flag to print detailed progress information. Default: FALSE
// //' @param use_hessian Logical flag to use Hessian information for parameter updates. If FALSE,
// //'        the algorithm uses gradient descent instead. Default: TRUE
// //' @param step_size Initial step size for parameter updates. Default: 1.0
// //' @param enforce_bounds Logical flag to enforce parameter constraints. Default: TRUE
// //' @param min_param_val Minimum allowed value for parameters (except delta). Default: 1e-5
// //' @param max_param_val Maximum allowed value for parameters. Default: 1e5
// //' @param get_num_hess Logical flag to calculate numerical Hessian in addition to analytical Hessian.
// //'        Default: FALSE
// //'
// //' @return A list containing the following components:
// //' \describe{
// //'   \item{parameters}{A numeric vector with the estimated parameters}
// //'   \item{loglik}{The maximized log-likelihood value}
// //'   \item{iterations}{Number of iterations performed}
// //'   \item{converged}{Logical flag indicating whether the algorithm converged}
// //'   \item{param_history}{Matrix of parameter values at each iteration}
// //'   \item{loglik_history}{Vector of log-likelihood values at each iteration}
// //'   \item{gradient}{The gradient vector at the final parameter estimates}
// //'   \item{hessian}{The Hessian matrix at the final parameter estimates}
// //'   \item{std_errors}{Standard errors for the estimated parameters}
// //'   \item{aic}{Akaike Information Criterion: AIC = 2k - 2ln(L)}
// //'   \item{bic}{Bayesian Information Criterion: BIC = k ln(n) - 2ln(L)}
// //'   \item{n}{Sample size}
// //'   \item{status}{Character string indicating the termination status}
// //'   \item{z_values}{Z-statistics for parameter significance tests}
// //'   \item{p_values}{P-values for parameter significance tests}
// //'   \item{param_names}{Character vector of parameter names}
// //'   \item{family}{The distribution family used in the estimation}
// //'   \item{numeric_hessian}{Numerical approximation of the Hessian (only if get_num_hess=TRUE)}
// //' }
// //'
// //' @section Warning:
// //' Convergence is not guaranteed for all datasets and initial values. It's recommended to:
// //' \itemize{
// //'   \item Try different initial values if convergence fails
// //'   \item Check the gradient norm at the final solution to verify optimality
// //'   \item Examine parameter history to identify potential issues
// //'   \item Use the verbose option to get detailed progress information for troubleshooting
// //' }
// //'
// //' @examples
// //' \dontrun{
// //' # Generate sample data from Beta(2,5) distribution for testing
// //' set.seed(123)
// //' sample_data <- rbeta(200, 2, 5)
// //'
// //' # Fit with full GKw model
// //' gkw_result <- mle_fit(c(1.5, 4.5, 1.0, 0.0, 1.0), sample_data, family = "gkw")
// //' gkw_result$parameters
// //'
// //' # Fit with simpler Kumaraswamy model
// //' kw_result <- mle_fit(c(1.5, 4.5), sample_data, family = "kw")
// //' kw_result$parameters
// //'
// //' # Fit with Beta model
// //' beta_result <- mle_fit(c(2.0, 5.0), sample_data, family = "beta")
// //' beta_result$parameters
// //'
// //' # Compare AIC/BIC values to select the best model
// //' data.frame(
// //'   family = c("gkw", "kw", "beta"),
// //'   AIC = c(gkw_result$aic, kw_result$aic, beta_result$aic),
// //'   BIC = c(gkw_result$bic, kw_result$bic, beta_result$bic)
// //' )
// //' }
// //'
// //' @references
// //' Kumaraswamy, P. (1980). A generalized probability density function for double-bounded
// //' random processes. Journal of Hydrology, 46(1-2), 79-88.
// //'
// //' Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized distributions.
// //' Journal of Statistical Computation and Simulation, 81(7), 883-898.
// //'
// //' Fletcher, R. (1987). Practical Methods of Optimization. John Wiley & Sons.
// //'
// //' @author Lopes, J. E.
// //'
// //' @export
// // [[Rcpp::export]]
// List nrgkw(
//    NumericVector start,
//    NumericVector data,
//    std::string family = "gkw",
//    double tol = 1e-6,
//    int max_iter = 100,
//    bool verbose = false,
//    bool use_hessian = true,
//    double step_size = 1.0,
//    bool enforce_bounds = true,
//    double min_param_val = 1e-5,
//    double max_param_val = 1e5,
//    bool get_num_hess = false
// ) {
//  // Final result will be a list with different components
//  List result;
//
//  // Convert family to lowercase for case-insensitive comparison
//  std::string family_lower = family;
//  std::transform(family_lower.begin(), family_lower.end(), family_lower.begin(), ::tolower);
//
//  // Determine number of parameters based on family
//  int n_params = 0;
//  CharacterVector param_names;
//
//  if (family_lower == "gkw") {
//    n_params = 5;
//    param_names = CharacterVector::create("alpha", "beta", "gamma", "delta", "lambda");
//  } else if (family_lower == "bkw") {
//    n_params = 4;
//    param_names = CharacterVector::create("alpha", "beta", "gamma", "delta");
//  } else if (family_lower == "kkw") {
//    n_params = 4;
//    param_names = CharacterVector::create("alpha", "beta", "delta", "lambda");
//  } else if (family_lower == "ekw") {
//    n_params = 3;
//    param_names = CharacterVector::create("alpha", "beta", "lambda");
//  } else if (family_lower == "mc" || family_lower == "mcdonald" || family_lower == "bp") {
//    n_params = 3;
//    param_names = CharacterVector::create("gamma", "delta", "lambda");
//  } else if (family_lower == "kw") {
//    n_params = 2;
//    param_names = CharacterVector::create("alpha", "beta");
//  } else if (family_lower == "beta") {
//    n_params = 2;
//    param_names = CharacterVector::create("gamma", "delta");
//  } else {
//    stop("Unknown family: '" + family + "'. Available options are 'gkw', 'bkw', 'kkw', 'ekw', 'mc', 'kw', 'beta'.");
//  }
//
//  // Validate initial parameters size
//  if (start.size() != n_params) {
//    stop("Invalid number of parameters for '" + family + "'. Expected " +
//      std::to_string(n_params) + ", got " + std::to_string(start.size()));
//  }
//
//  // Check for valid data
//  int n_data = data.size();
//  if (n_data < n_params) {
//    stop("At least " + std::to_string(n_params) + " data points are needed to estimate " +
//      std::to_string(n_params) + " parameters");
//  }
//
//  // Check if all data are in the interval (0,1)
//  for (int i = 0; i < n_data; i++) {
//    if (data[i] <= 0.0 || data[i] >= 1.0 || !R_finite(data[i])) {
//      stop("All data must be in the interval (0,1)");
//    }
//  }
//
//  // Copy initial parameters and convert to standard GKw parameters where needed
//  NumericVector params(5); // Always use 5 parameters internally (GKw format)
//
//  // Set default values based on fixed parameters in specific families
//  params[0] = 1.0; // α = 1 default
//  params[1] = 1.0; // β = 1 default
//  params[2] = 1.0; // γ = 1 default
//  params[3] = 0.0; // δ = 0 default
//  params[4] = 1.0; // λ = 1 default
//
//  // Fill with provided parameters based on family
//  if (family_lower == "gkw") {
//    for (int j = 0; j < 5; j++) {
//      params[j] = start[j];
//    }
//  } else if (family_lower == "bkw") {
//    // α, β, γ, δ with λ = 1
//    for (int j = 0; j < 4; j++) {
//      params[j] = start[j];
//    }
//    params[4] = 1.0; // λ fixed at 1
//  } else if (family_lower == "kkw") {
//    // α, β, δ, λ with γ = 1
//    params[0] = start[0]; // α
//    params[1] = start[1]; // β
//    params[2] = 1.0;             // γ fixed at 1
//    params[3] = start[2]; // δ
//    params[4] = start[3]; // λ
//  } else if (family_lower == "ekw") {
//    // α, β, λ with γ = 1, δ = 0
//    params[0] = start[0]; // α
//    params[1] = start[1]; // β
//    params[2] = 1.0;             // γ fixed at 1
//    params[3] = 0.0;             // δ fixed at 0
//    params[4] = start[2]; // λ
//  } else if (family_lower == "mc" || family_lower == "mcdonald" || family_lower == "bp") {
//    // γ, δ, λ with α = 1, β = 1
//    params[0] = 1.0;             // α fixed at 1
//    params[1] = 1.0;             // β fixed at 1
//    params[2] = start[0]; // γ
//    params[3] = start[1]; // δ
//    params[4] = start[2]; // λ
//  } else if (family_lower == "kw") {
//    // α, β with γ = 1, δ = 0, λ = 1
//    params[0] = start[0]; // α
//    params[1] = start[1]; // β
//    params[2] = 1.0;             // γ fixed at 1
//    params[3] = 0.0;             // δ fixed at 0
//    params[4] = 1.0;             // λ fixed at 1
//  } else if (family_lower == "beta") {
//    // γ, δ with α = 1, β = 1, λ = 1
//    params[0] = 1.0;             // α fixed at 1
//    params[1] = 1.0;             // β fixed at 1
//    params[2] = start[0]; // γ
//    params[3] = start[1]; // δ
//    params[4] = 1.0;             // λ fixed at 1
//  }
//
//  // Apply constraints to initial parameters if needed
//  if (enforce_bounds) {
//    for (int j = 0; j < 5; j++) {
//      if (j == 3) { // delta
//        params[j] = std::max(0.0, params[j]);
//      } else { // other parameters must be > 0
//        params[j] = std::max(min_param_val, params[j]);
//      }
//      params[j] = std::min(max_param_val, params[j]);
//    }
//  }
//
//  // Define function pointers based on family
//  std::function<double(NumericVector, NumericVector)> ll_func;
//  std::function<NumericVector(NumericVector, NumericVector)> gr_func;
//  std::function<NumericMatrix(NumericVector, NumericVector)> hs_func;
//
//  // Assign appropriate functions based on family
//  if (family_lower == "gkw") {
//    ll_func = llgkw;
//    gr_func = grgkw;
//    hs_func = hsgkw;
//  } else if (family_lower == "bkw") {
//    ll_func = llbkw;
//    gr_func = grbkw;
//    hs_func = hsbkw;
//  } else if (family_lower == "kkw") {
//    ll_func = llkkw;
//    gr_func = grkkw;
//    hs_func = hskkw;
//  } else if (family_lower == "ekw") {
//    ll_func = llekw;
//    gr_func = grekw;
//    hs_func = hsekw;
//  } else if (family_lower == "mc" || family_lower == "mcdonald" || family_lower == "bp") {
//    ll_func = llmc;
//    gr_func = grmc;
//    hs_func = hsmc;
//  } else if (family_lower == "kw") {
//    ll_func = llkw;
//    gr_func = grkw;
//    hs_func = hskw;
//  } else if (family_lower == "beta") {
//    ll_func = llbeta;
//    gr_func = grbeta;
//    hs_func = hsbeta;
//  }
//
//  // Function to extract relevant parameters for specific family
//  auto extract_params = [&](const NumericVector& full_params) -> NumericVector {
//    NumericVector result;
//
//    if (family_lower == "gkw") {
//      result = NumericVector(5);
//      for (int j = 0; j < 5; j++) result[j] = full_params[j];
//    } else if (family_lower == "bkw") {
//      result = NumericVector(4);
//      for (int j = 0; j < 4; j++) result[j] = full_params[j];
//    } else if (family_lower == "kkw") {
//      result = NumericVector(4);
//      result[0] = full_params[0]; // α
//      result[1] = full_params[1]; // β
//      result[2] = full_params[3]; // δ
//      result[3] = full_params[4]; // λ
//    } else if (family_lower == "ekw") {
//      result = NumericVector(3);
//      result[0] = full_params[0]; // α
//      result[1] = full_params[1]; // β
//      result[2] = full_params[4]; // λ
//    } else if (family_lower == "mc" || family_lower == "mcdonald" || family_lower == "bp") {
//      result = NumericVector(3);
//      result[0] = full_params[2]; // γ
//      result[1] = full_params[3]; // δ
//      result[2] = full_params[4]; // λ
//    } else if (family_lower == "kw") {
//      result = NumericVector(2);
//      result[0] = full_params[0]; // α
//      result[1] = full_params[1]; // β
//    } else if (family_lower == "beta") {
//      result = NumericVector(2);
//      result[0] = full_params[2]; // γ
//      result[1] = full_params[3]; // δ
//    }
//
//    return result;
//  };
//
//  // Function to update full params given updated params for specific family
//  auto update_full_params = [&](NumericVector& full_params, const NumericVector& updated_params) {
//    if (family_lower == "gkw") {
//      for (int j = 0; j < 5; j++) full_params[j] = updated_params[j];
//    } else if (family_lower == "bkw") {
//      for (int j = 0; j < 4; j++) full_params[j] = updated_params[j];
//    } else if (family_lower == "kkw") {
//      full_params[0] = updated_params[0]; // α
//      full_params[1] = updated_params[1]; // β
//      full_params[3] = updated_params[2]; // δ
//      full_params[4] = updated_params[3]; // λ
//    } else if (family_lower == "ekw") {
//      full_params[0] = updated_params[0]; // α
//      full_params[1] = updated_params[1]; // β
//      full_params[4] = updated_params[2]; // λ
//    } else if (family_lower == "mc" || family_lower == "mcdonald" || family_lower == "bp") {
//      full_params[2] = updated_params[0]; // γ
//      full_params[3] = updated_params[1]; // δ
//      full_params[4] = updated_params[2]; // λ
//    } else if (family_lower == "kw") {
//      full_params[0] = updated_params[0]; // α
//      full_params[1] = updated_params[1]; // β
//    } else if (family_lower == "beta") {
//      full_params[2] = updated_params[0]; // γ
//      full_params[3] = updated_params[1]; // δ
//    }
//  };
//
//  // Create a custom numHess function to handle specific families
//  auto numHess_family = [&](NumericVector params_family, NumericVector data_family, double eps = 1e-6) {
//    int n_params_family = params_family.size();
//    arma::mat hessian(n_params_family, n_params_family, arma::fill::zeros);
//
//    // Value of the function at the current point
//    double f0 = ll_func(params_family, data_family);
//
//    // For each pair of variables
//    for (int i = 0; i < n_params_family; i++) {
//      for (int j = 0; j <= i; j++) {
//        // Calculate second derivative using finite differences
//        if (i == j) {
//          // Second derivative with respect to the same variable
//          NumericVector params_p = clone(params_family);
//          NumericVector params_m = clone(params_family);
//
//          params_p[i] += eps;
//          params_m[i] -= eps;
//
//          double f_p = ll_func(params_p, data_family);
//          double f_m = ll_func(params_m, data_family);
//
//          hessian(i, j) = (f_p - 2*f0 + f_m) / (eps * eps);
//        } else {
//          // Mixed derivative
//          NumericVector params_pp = clone(params_family);
//          NumericVector params_pm = clone(params_family);
//          NumericVector params_mp = clone(params_family);
//          NumericVector params_mm = clone(params_family);
//
//          params_pp[i] += eps;
//          params_pp[j] += eps;
//
//          params_pm[i] += eps;
//          params_pm[j] -= eps;
//
//          params_mp[i] -= eps;
//          params_mp[j] += eps;
//
//          params_mm[i] -= eps;
//          params_mm[j] -= eps;
//
//          double f_pp = ll_func(params_pp, data_family);
//          double f_pm = ll_func(params_pm, data_family);
//          double f_mp = ll_func(params_mp, data_family);
//          double f_mm = ll_func(params_mm, data_family);
//
//          hessian(i, j) = hessian(j, i) = (f_pp - f_pm - f_mp + f_mm) / (4 * eps * eps);
//        }
//      }
//    }
//
//    return hessian;
//  };
//
//  // Get family-specific parameters
//  NumericVector family_params = extract_params(params);
//
//  // Calculate initial log-likelihood
//  double initial_loglik = ll_func(family_params, data);
//  if (!R_finite(initial_loglik) || initial_loglik == R_PosInf) {
//    stop("Initial log-likelihood is infinite or NaN. Check the initial parameters.");
//  }
//
//  // Parameter and log-likelihood history for diagnostics
//  NumericMatrix param_history(max_iter + 1, n_params);
//  NumericVector loglik_history(max_iter + 1);
//
//  // Initialize history with initial values
//  for (int j = 0; j < n_params; j++) {
//    param_history(0, j) = family_params[j];
//  }
//  loglik_history[0] = initial_loglik;
//
//  // Variables for convergence control
//  bool converged = false;
//  int iter = 0;
//  double prev_loglik = initial_loglik;
//
//  // Prepare to store the best result obtained
//  double best_loglik = initial_loglik;
//  NumericVector best_params = clone(family_params);
//
//  // Main Newton-Raphson loop
//  while (!converged && iter < max_iter) {
//    iter++;
//
//    // Calculate log-likelihood, gradient and hessian
//    double current_loglik = ll_func(family_params, data);
//    NumericVector gradient = gr_func(family_params, data);
//
//    // Check if gradient has valid values
//    bool valid_gradient = true;
//    for (int j = 0; j < n_params; j++) {
//      if (!R_finite(gradient[j])) {
//        valid_gradient = false;
//        break;
//      }
//    }
//
//    if (!valid_gradient) {
//      if (verbose) {
//        Rcout << "Warning: Invalid gradient in iteration " << iter << std::endl;
//      }
//      result["converged"] = false;
//      result["status"] = "gradient_failure";
//      // Use the best parameters found so far
//      family_params = best_params;
//      break;
//    }
//
//    // Calculate gradient norm for stopping criterion
//    double grad_norm = 0.0;
//    for (int j = 0; j < n_params; j++) {
//      grad_norm += gradient[j] * gradient[j];
//    }
//    grad_norm = std::sqrt(grad_norm);
//
//    if (grad_norm < tol) {
//      converged = true;
//      if (verbose) {
//        Rcout << "Convergence detected: gradient norm (" << grad_norm
//              << ") < tolerance (" << tol << ")" << std::endl;
//      }
//      break;
//    }
//
//    // Update direction
//    NumericVector update(n_params);
//
//    if (use_hessian) {
//      // Calculate the Hessian
//      NumericMatrix rcpp_hessian = hs_func(family_params, data);
//
//      // Check if Hessian has valid values
//      bool valid_hessian = true;
//      for (int i = 0; i < n_params; i++) {
//        for (int j = 0; j < n_params; j++) {
//          if (!R_finite(rcpp_hessian(i, j))) {
//            valid_hessian = false;
//            break;
//          }
//        }
//        if (!valid_hessian) break;
//      }
//
//      if (!valid_hessian) {
//        if (verbose) {
//          Rcout << "Warning: Invalid Hessian in iteration " << iter
//                << ", using only the gradient." << std::endl;
//        }
//        // Fallback to steepest descent if Hessian is invalid
//        for (int j = 0; j < n_params; j++) {
//          // Normalize gradient for step control
//          update[j] = -step_size * gradient[j] / std::max(1.0, std::abs(gradient[j]));
//        }
//      } else {
//        // Convert to arma::mat for more robust matrix operations
//        arma::mat hessian = as<arma::mat>(rcpp_hessian);
//        arma::vec grad_vec = as<arma::vec>(gradient);
//        arma::vec neg_grad = -grad_vec;
//        arma::vec update_vec;
//
//        bool solve_success = false;
//
//        // Try 1: Cholesky for symmetric positive definite matrices (fastest)
//        try {
//          update_vec = arma::solve(hessian, neg_grad, arma::solve_opts::likely_sympd);
//          solve_success = true;
//
//          if (verbose) {
//            Rcout << "Cholesky decomposition successful for parameter update." << std::endl;
//          }
//        } catch (...) {
//          if (verbose) {
//            Rcout << "Warning: Cholesky decomposition failed, trying standard solver..." << std::endl;
//          }
//
//          // Try 2: Standard Armadillo solver
//          try {
//            update_vec = arma::solve(hessian, neg_grad);
//            solve_success = true;
//
//            if (verbose) {
//              Rcout << "Standard solver successful for parameter update." << std::endl;
//            }
//          } catch (...) {
//            if (verbose) {
//              Rcout << "Warning: Standard solver failed, trying with regularization..." << std::endl;
//            }
//
//            // Try 3: Regularize the Hessian matrix
//            arma::mat reg_hessian = hessian;
//            double reg_factor = 1e-6;
//
//            // Find reasonable magnitude for regularization
//            double diag_max = arma::max(arma::abs(reg_hessian.diag()));
//            reg_factor = std::max(reg_factor, 1e-6 * diag_max);
//
//            // Add small value to diagonal
//            reg_hessian.diag() += reg_factor;
//
//            try {
//              update_vec = arma::solve(reg_hessian, neg_grad);
//              solve_success = true;
//
//              if (verbose) {
//                Rcout << "Regularized solver successful with factor: " << reg_factor << std::endl;
//              }
//            } catch (...) {
//              // Try 4: Stronger regularization
//              reg_hessian = hessian;
//              reg_factor = 1e-4 * (1.0 + diag_max);
//              reg_hessian.diag() += reg_factor;
//
//              try {
//                update_vec = arma::solve(reg_hessian, neg_grad);
//                solve_success = true;
//
//                if (verbose) {
//                  Rcout << "Stronger regularization successful with factor: " << reg_factor << std::endl;
//                }
//              } catch (...) {
//                // Try 5: Pseudo-inverse (very robust method)
//                try {
//                  arma::mat hess_pinv = arma::pinv(hessian);
//                  update_vec = hess_pinv * neg_grad;
//                  solve_success = true;
//
//                  if (verbose) {
//                    Rcout << "Pseudo-inverse solution successful for parameter update." << std::endl;
//                  }
//                } catch (...) {
//                  if (verbose) {
//                    Rcout << "Warning: All matrix inversion methods failed in iteration " << iter
//                          << ", using only the gradient." << std::endl;
//                  }
//
//                  // If all attempts fail, use gradient descent
//                  for (int j = 0; j < n_params; j++) {
//                    update[j] = -step_size * gradient[j] / std::max(1.0, std::abs(gradient[j]));
//                  }
//
//                  solve_success = false;
//                }
//              }
//            }
//          }
//        }
//
//        if (solve_success) {
//          // Convert solution from arma::vec to NumericVector
//          update = wrap(update_vec);
//
//          // Limit step size to avoid too large steps
//          double max_update = 0.0;
//          for (int j = 0; j < n_params; j++) {
//            max_update = std::max(max_update, std::abs(update[j]));
//          }
//
//          // If step is too large, reduce proportionally
//          const double max_step = 2.0;
//          if (max_update > max_step) {
//            double scale_factor = max_step / max_update;
//            for (int j = 0; j < n_params; j++) {
//              update[j] *= scale_factor;
//            }
//          }
//
//          // Apply step_size
//          for (int j = 0; j < n_params; j++) {
//            update[j] *= step_size;
//          }
//        }
//      }
//    } else {
//      // Use only gradient (gradient descent method)
//      for (int j = 0; j < n_params; j++) {
//        update[j] = -step_size * gradient[j] / std::max(1.0, std::abs(gradient[j]));
//      }
//    }
//
//    // Update parameters: theta_new = theta_old + update
//    NumericVector new_params(n_params);
//    for (int j = 0; j < n_params; j++) {
//      new_params[j] = family_params[j] + update[j];
//    }
//
//    // Enforce bounds if requested
//    if (enforce_bounds) {
//      for (int j = 0; j < n_params; j++) {
//        bool is_delta = (family_lower == "gkw" && j == 3) ||
//          (family_lower == "bkw" && j == 3) ||
//          (family_lower == "kkw" && j == 2) ||
//          (family_lower == "mc" && j == 1) ||
//          (family_lower == "beta" && j == 1);
//
//        // Note: for delta, we allow values down to 0
//        if (is_delta) {
//          new_params[j] = std::max(0.0, new_params[j]);
//        } else {
//          new_params[j] = std::max(min_param_val, new_params[j]);
//        }
//        new_params[j] = std::min(max_param_val, new_params[j]);
//      }
//    }
//
//    // Calculate new objective function value
//    double new_loglik = ll_func(new_params, data);
//
//    // Line search / Backtracking if new value is not better
//    bool backtracking_success = true;
//    double bt_step = 1.0;
//    const double bt_factor = 0.5; // reduce step by half each backtracking
//    const int max_bt = 10;        // maximum backtracking iterations
//
//    if (!R_finite(new_loglik) || new_loglik >= current_loglik) {
//      backtracking_success = false;
//
//      if (verbose) {
//        Rcout << "Starting backtracking at iteration " << iter
//              << ", current value: " << current_loglik
//              << ", new value: " << new_loglik << std::endl;
//      }
//
//      for (int bt = 0; bt < max_bt; bt++) {
//        bt_step *= bt_factor;
//
//        // Recalculate new parameters with reduced step
//        for (int j = 0; j < n_params; j++) {
//          new_params[j] = family_params[j] + bt_step * update[j];
//        }
//
//        // Enforce bounds again
//        if (enforce_bounds) {
//          for (int j = 0; j < n_params; j++) {
//            bool is_delta = (family_lower == "gkw" && j == 3) ||
//              (family_lower == "bkw" && j == 3) ||
//              (family_lower == "kkw" && j == 2) ||
//              (family_lower == "mc" && j == 1) ||
//              (family_lower == "beta" && j == 1);
//
//            if (is_delta) {
//              new_params[j] = std::max(0.0, new_params[j]);
//            } else {
//              new_params[j] = std::max(min_param_val, new_params[j]);
//            }
//            new_params[j] = std::min(max_param_val, new_params[j]);
//          }
//        }
//
//        // Test new value
//        new_loglik = ll_func(new_params, data);
//
//        if (R_finite(new_loglik) && new_loglik < current_loglik) {
//          backtracking_success = true;
//          if (verbose) {
//            Rcout << "Backtracking successful after " << (bt + 1)
//                  << " attempts, new value: " << new_loglik << std::endl;
//          }
//          break;
//        }
//      }
//    }
//
//    // If we still cannot improve, evaluate the situation
//    if (!backtracking_success) {
//      if (verbose) {
//        Rcout << "Warning: Backtracking failed at iteration " << iter << std::endl;
//      }
//
//      // If gradient is small enough, consider converged
//      if (grad_norm < tol * 10) {  // Relaxed tolerance for this case
//        converged = true;
//        if (verbose) {
//          Rcout << "Convergence detected with small gradient after backtracking failure." << std::endl;
//        }
//      } else {
//        // If backtracking fails and we're close to max iterations,
//        // check if we're in a reasonable region
//        if (iter > max_iter * 0.8 && current_loglik < best_loglik * 1.1) {
//          converged = true;
//          if (verbose) {
//            Rcout << "Forced convergence after backtracking failure near maximum iterations." << std::endl;
//          }
//        } else {
//          // Try a small random perturbation
//          NumericVector perturb(n_params);
//          for (int j = 0; j < n_params; j++) {
//            // Perturbation of up to 5% of current value
//            double range = 0.05 * std::abs(family_params[j]);
//            perturb[j] = R::runif(-range, range);
//            new_params[j] = family_params[j] + perturb[j];
//          }
//
//          // Apply constraints
//          if (enforce_bounds) {
//            for (int j = 0; j < n_params; j++) {
//              bool is_delta = (family_lower == "gkw" && j == 3) ||
//                (family_lower == "bkw" && j == 3) ||
//                (family_lower == "kkw" && j == 2) ||
//                (family_lower == "mc" && j == 1) ||
//                (family_lower == "beta" && j == 1);
//
//              if (is_delta) {
//                new_params[j] = std::max(0.0, new_params[j]);
//              } else {
//                new_params[j] = std::max(min_param_val, new_params[j]);
//              }
//              new_params[j] = std::min(max_param_val, new_params[j]);
//            }
//          }
//
//          new_loglik = ll_func(new_params, data);
//
//          if (R_finite(new_loglik) && new_loglik < current_loglik) {
//            backtracking_success = true;
//            if (verbose) {
//              Rcout << "Recovery by random perturbation, new value: " << new_loglik << std::endl;
//            }
//          } else {
//            // If even perturbation doesn't work, use the best result so far
//            new_params = best_params;
//            new_loglik = best_loglik;
//            if (verbose) {
//              Rcout << "Returning to the previous best result: " << -best_loglik << std::endl;
//            }
//          }
//        }
//      }
//    }
//
//    // Update parameters and history
//    for (int j = 0; j < n_params; j++) {
//      family_params[j] = new_params[j];
//      param_history(iter, j) = family_params[j];
//    }
//    loglik_history[iter] = new_loglik;
//
//    // Update the best result if this is better
//    if (new_loglik < best_loglik) {
//      best_loglik = new_loglik;
//      for (int j = 0; j < n_params; j++) {
//        best_params[j] = family_params[j];
//      }
//    }
//
//    // Check convergence by parameter change
//    double param_change = 0.0;
//    double param_rel_change = 0.0;
//    for (int j = 0; j < n_params; j++) {
//      param_change += std::pow(update[j], 2);
//      if (std::abs(family_params[j]) > 1e-10) {
//        param_rel_change += std::pow(update[j] / family_params[j], 2);
//      } else {
//        param_rel_change += std::pow(update[j], 2);
//      }
//    }
//    param_change = std::sqrt(param_change);
//    param_rel_change = std::sqrt(param_rel_change / n_params);
//
//    // Check convergence by log-likelihood change
//    double loglik_change = std::abs(prev_loglik - new_loglik);
//    double loglik_rel_change = loglik_change / (std::abs(prev_loglik) + 1e-10);
//    prev_loglik = new_loglik;
//
//    if (verbose) {
//      Rcout << "Iteration " << iter
//            << ", Log-likelihood: " << -new_loglik
//            << ", Change: " << loglik_change
//            << ", Rel. Change: " << loglik_rel_change
//            << ", Gradient Norm: " << grad_norm
//            << std::endl;
//
//      Rcout << "Parameters:";
//      for (int j = 0; j < n_params; j++) {
//        Rcout << " " << family_params[j];
//      }
//      Rcout << std::endl;
//    }
//
//    // Convergence criteria
//    if (param_change < tol || param_rel_change < tol ||
//        loglik_change < tol || loglik_rel_change < tol) {
//      converged = true;
//      if (verbose) {
//        Rcout << "Convergence detected:" << std::endl;
//        if (param_change < tol) Rcout << "- Absolute parameter change < tolerance" << std::endl;
//        if (param_rel_change < tol) Rcout << "- Relative parameter change < tolerance" << std::endl;
//        if (loglik_change < tol) Rcout << "- Absolute log-likelihood change < tolerance" << std::endl;
//        if (loglik_rel_change < tol) Rcout << "- Relative log-likelihood change < tolerance" << std::endl;
//      }
//    }
//  }
//
//  // If not converged, use the best parameters found
//  if (!converged) {
//    family_params = best_params;
//    if (verbose) {
//      Rcout << "Did not fully converge, using the best parameters found." << std::endl;
//    }
//  }
//
//  // Prepare final result
//  NumericMatrix final_param_history(iter + 1, n_params);
//  NumericVector final_loglik_history(iter + 1);
//
//  for (int i = 0; i <= iter; i++) {
//    for (int j = 0; j < n_params; j++) {
//      final_param_history(i, j) = param_history(i, j);
//    }
//    final_loglik_history[i] = loglik_history[i];
//  }
//
//  // Calculate final gradient and hessian
//  NumericVector final_gradient = gr_func(family_params, data);
//  NumericMatrix rcpp_hessian = hs_func(family_params, data);
//
//  // Calculate numerical Hessian if requested
//  NumericMatrix rcpp_numeric_hessian;
//  if (get_num_hess) {
//    arma::mat arma_numeric_hessian = numHess_family(family_params, data);
//    rcpp_numeric_hessian = wrap(arma_numeric_hessian);
//    result["numeric_hessian"] = rcpp_numeric_hessian;
//  }
//
//  // Calculate standard errors using Armadillo for robust matrix inversion
//  NumericVector std_errors(n_params, NA_REAL);
//  bool valid_se = true;
//
//  // Convert Rcpp Hessian to Armadillo matrix
//  arma::mat hessian = as<arma::mat>(rcpp_hessian);
//  arma::mat cov_matrix;
//
//  // Layered approach to calculate covariance matrix (inverse of Hessian)
//  try {
//    // Step 1: Try to use Cholesky decomposition (fastest, requires positive definite)
//    try {
//      cov_matrix = arma::inv_sympd(hessian);
//
//      if (verbose) {
//        Rcout << "Standard error calculation: Cholesky decomposition successful." << std::endl;
//      }
//    } catch (...) {
//      // Step 2: Try standard inversion
//      try {
//        cov_matrix = arma::inv(hessian);
//
//        if (verbose) {
//          Rcout << "Standard error calculation: Standard inverse successful." << std::endl;
//        }
//      } catch (...) {
//        // Step 3: Apply regularization
//        arma::mat reg_hessian = hessian;
//        double reg_factor = 1e-6;
//
//        // Find reasonable magnitude for regularization
//        double diag_max = arma::max(arma::abs(reg_hessian.diag()));
//        reg_factor = std::max(reg_factor, 1e-6 * diag_max);
//
//        // Add small value to diagonal
//        reg_hessian.diag() += reg_factor;
//
//        try {
//          cov_matrix = arma::inv(reg_hessian);
//
//          if (verbose) {
//            Rcout << "Standard error calculation: Regularized inverse successful." << std::endl;
//          }
//        } catch (...) {
//          // Step 4: Stronger regularization
//          reg_hessian = hessian;
//          reg_factor = 1e-4 * (1.0 + diag_max);
//          reg_hessian.diag() += reg_factor;
//
//          try {
//            cov_matrix = arma::inv(reg_hessian);
//
//            if (verbose) {
//              Rcout << "Standard error calculation: Stronger regularized inverse successful." << std::endl;
//            }
//          } catch (...) {
//            // Step 5: Use pseudo-inverse (more robust)
//            try {
//              cov_matrix = arma::pinv(hessian);
//
//              if (verbose) {
//                Rcout << "Standard error calculation: Pseudo-inverse successful." << std::endl;
//              }
//            } catch (...) {
//              // Step 6: Try numerical Hessian if available
//              if (get_num_hess) {
//                arma::mat num_hessian = as<arma::mat>(rcpp_numeric_hessian);
//
//                try {
//                  cov_matrix = arma::pinv(num_hessian);
//
//                  if (verbose) {
//                    Rcout << "Standard error calculation: Numerical Hessian pseudo-inverse successful." << std::endl;
//                  }
//                } catch (...) {
//                  valid_se = false;
//                }
//              } else {
//                valid_se = false;
//              }
//            }
//          }
//        }
//      }
//    }
//
//    // Calculate standard errors if covariance matrix is available
//    if (valid_se) {
//      // Extract diagonal elements and calculate square root
//      arma::vec diag_cov = cov_matrix.diag();
//
//      for (int j = 0; j < n_params; j++) {
//        if (diag_cov(j) > 0) {
//          std_errors[j] = std::sqrt(diag_cov(j));
//        } else {
//          if (verbose) {
//            Rcout << "Warning: Non-positive variance detected for parameter " << j
//                  << ". Standard error set to NA." << std::endl;
//          }
//          std_errors[j] = NA_REAL;
//        }
//      }
//    }
//  } catch (...) {
//    valid_se = false;
//  }
//
//  if (!valid_se && verbose) {
//    Rcout << "Warning: Could not calculate standard errors. The Hessian matrix may not be positive definite." << std::endl;
//  }
//
//  // Calculate AIC: AIC = 2k - 2ln(L) = 2k + 2*(-ln(L))
//  double final_loglik = ll_func(family_params, data);
//  double aic = 2 * n_params + 2 * final_loglik;
//
//  // Calculate BIC: BIC = k ln(n) - 2ln(L) = k ln(n) + 2*(-ln(L))
//  double bic = n_params * std::log(n_data) + 2 * final_loglik;
//
//  // Fill the result
//  result["parameters"] = family_params;
//  result["loglik"] = -final_loglik;  // Negative because ll functions return -logL
//  result["iterations"] = iter;
//  result["converged"] = converged;
//  result["param_history"] = final_param_history;
//  result["loglik_history"] = -final_loglik_history;  // Negative for consistency
//  result["gradient"] = final_gradient;
//  result["hessian"] = rcpp_hessian;
//  result["std_errors"] = std_errors;
//  result["aic"] = aic;
//  result["bic"] = bic;
//  result["n"] = n_data;
//  result["family"] = family;
//
//  if (!converged && !result.containsElementNamed("status")) {
//    result["status"] = "max_iterations_reached";
//  } else if (converged) {
//    result["status"] = "success";
//  }
//
//  // Calculate statistical significance (p-values) using normal approximation
//  NumericVector z_values(n_params, NA_REAL);
//  NumericVector p_values(n_params, NA_REAL);
//
//  if (valid_se) {
//    for (int j = 0; j < n_params; j++) {
//      if (std_errors[j] != NA_REAL && std_errors[j] > 0) {
//        z_values[j] = family_params[j] / std_errors[j];
//        // Two-tailed approximation
//        p_values[j] = 2.0 * R::pnorm(-std::abs(z_values[j]), 0.0, 1.0, 1, 0);
//      }
//    }
//    result["z_values"] = z_values;
//    result["p_values"] = p_values;
//  }
//
//  // Set parameter names for easier interpretation
//  colnames(final_param_history) = param_names;
//  result["param_names"] = param_names;
//
//  return result;
// }
//
