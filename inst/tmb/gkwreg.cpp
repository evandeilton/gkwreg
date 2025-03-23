// File: gkwreg_v5.cpp
// ---------------------------------------------------------------
//  Regression Model for the Generalized Kumaraswamy (GKw) Distribution
//  using TMB, focusing on scalability with big data and many covariates.
//
//  Changes from v4:
//   - Adds dynamic scheduling in OpenMP
//   - Adds options to disable caching and/or fitted means
//   - Uses unordered_map for caching
// ---------------------------------------------------------------

#include <TMB.hpp>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <cassert>

#ifdef _OPENMP
#include <omp.h>
#else
inline int omp_get_max_threads() { return 1; }
inline int omp_get_thread_num() { return 0; }
#define omp_parallel
#define omp_critical
#endif

// ========================
// Numeric constants
// ========================
template <class Type>
struct Constants {
  static const Type eps_log;
  static const Type eps_pos;
  static const Type eps_prob;
  static const Type inf_repl;
  static const Type max_exp;
};

template <class Type>
const Type Constants<Type>::eps_log  = Type(1e-15);
template <class Type>
const Type Constants<Type>::eps_pos  = Type(1e-10);
template <class Type>
const Type Constants<Type>::eps_prob = Type(1e-12);
template <class Type>
const Type Constants<Type>::inf_repl = Type(1e10);
template <class Type>
const Type Constants<Type>::max_exp  = Type(30);

// ========================
// Safe numerical operations
// ========================
template <class Type>
Type safeLog(const Type &x) {
  // Avoid log(<=0)
  if (x <= Type(0.0)) return -Constants<Type>::inf_repl;
  return log(x + Constants<Type>::eps_log);
}

template <class Type>
Type safeExp(const Type &x) {
  // Clamp exponent to avoid overflow
  if (x >  Constants<Type>::max_exp)  return exp(Constants<Type>::max_exp);
  if (x < -Constants<Type>::max_exp)  return exp(-Constants<Type>::max_exp);
  return exp(x);
}

template <class Type>
Type safePow(const Type &base, const Type &exponent) {
  // Protect base^exponent in extreme cases
  if (base <= Constants<Type>::eps_pos) {
    return (exponent > Type(0.0)) ? Type(0.0) : Constants<Type>::inf_repl;
  }
  if (fabs(exponent) > Type(1.0)) {
    return safeExp(exponent * safeLog(base));
  }
  return pow(base, exponent);
}

// ========================
// Enforce domain
// ========================
template <class Type>
Type enforceMin(const Type &x, const Type &min_val) {
  // max(x, min_val)
  return CppAD::CondExpLt(x, min_val, min_val, x);
}

template <class Type>
Type enforceProbability(const Type &x) {
  // Clamps x into (eps_prob, 1 - eps_prob)
  Type eps   = Constants<Type>::eps_prob;
  Type upper = Type(1.0) - eps;
  Type clippedUp = CppAD::CondExpGt(x, upper, upper, x);
  return CppAD::CondExpLt(clippedUp, eps, eps, clippedUp);
}

// ========================
// Link functions
// (1=log, 2=logit, 3=probit, 4=cauchy, 5=cloglog,
//  6=identity, 7=sqrt, 8=inverse, 9=inverse-square)
// ========================
template <class Type>
Type fast_log_link(const Type &eta) {
  return safeExp(eta);
}

template <class Type>
Type fast_logit_link(const Type &eta) {
  if (eta >  Constants<Type>::max_exp) return Type(1.0) - Constants<Type>::eps_prob;
  if (eta < -Constants<Type>::max_exp) return Constants<Type>::eps_prob;
  return Type(1.0) / (Type(1.0) + safeExp(-eta));
}

template <class Type>
Type fast_probit_link(const Type &eta) {
  if (eta >  Type(5.0))  return Type(1.0) - Constants<Type>::eps_prob;
  if (eta < -Type(5.0))  return Constants<Type>::eps_prob;
  return enforceProbability(pnorm(eta));
}

template <class Type>
Type fast_cloglog_link(const Type &eta) {
  if (eta >  Constants<Type>::max_exp) return Type(1.0) - Constants<Type>::eps_prob;
  if (eta < -Constants<Type>::max_exp) return Constants<Type>::eps_prob;
  return enforceProbability(Type(1.0) - safeExp(-safeExp(eta)));
}

template <class Type>
Type apply_positive_link(const Type &eta, const int &link_type, const Type &scale_factor) {
  Type min_val = Constants<Type>::eps_pos;
  switch (link_type) {
  case 1: // log
    return enforceMin(fast_log_link(eta), min_val);
  case 2: // logit
    return enforceMin(scale_factor * fast_logit_link(eta), min_val);
  case 3: // probit
    return enforceMin(scale_factor * fast_probit_link(eta), min_val);
  case 4: // cauchy
    // cdf cauchy = 0.5 + atan(eta)/pi
    return enforceMin(scale_factor * (Type(0.5) + (atan(eta)/Type(M_PI))), min_val);
  case 5: // cloglog
    return enforceMin(scale_factor * fast_cloglog_link(eta), min_val);
  case 6: // identity (clamped)
    return enforceMin(eta, min_val);
  case 7: // sqrt
    return enforceMin(CppAD::CondExpGt(eta, Type(0.0), eta*eta, Type(0.0)), min_val);
  case 8: // inverse
    return enforceMin(Type(1.0)/enforceMin(eta, Type(1e-6)), min_val);
  case 9: // inverse-square
    return enforceMin(Type(1.0)/sqrt(enforceMin(eta, Type(1e-6))), min_val);
  default:
    // fallback to log
    return enforceMin(fast_log_link(eta), min_val);
  }
}

// ========================
// Log Beta
// ========================
template <class Type>
Type logBeta(const Type &a, const Type &b) {
  Type aa = enforceMin(a, Constants<Type>::eps_pos);
  Type bb = enforceMin(b, Constants<Type>::eps_pos);
  // For large a,b, approximate
  if (aa > Type(100.0) && bb > Type(100.0)) {
    return Type(0.5)*( log(Type(2.0)*M_PI) - log(aa+bb) )
    + (aa - Type(0.5))*log(aa)
    + (bb - Type(0.5))*log(bb)
    - (aa + bb - Type(1.0))*log(aa+bb);
  }
  // Standard
  return lgamma(aa) + lgamma(bb) - lgamma(aa + bb);
}

// ========================
// GKw log-PDF
// ========================
template <class Type>
Type log_pdf_gkw(const Type &y,
                 const Type &alpha_,
                 const Type &beta_,
                 const Type &gamma_,
                 const Type &delta_,
                 const Type &lambda_) {
  if ( y <= Constants<Type>::eps_prob || y >= Type(1.0)-Constants<Type>::eps_prob ) {
    return -Constants<Type>::inf_repl;
  }
  if ( alpha_ <= Constants<Type>::eps_pos ||
       beta_   <= Constants<Type>::eps_pos ||
       gamma_  <= Constants<Type>::eps_pos ||
       delta_  <= Constants<Type>::eps_pos ||
       lambda_ <= Constants<Type>::eps_pos ) {
    return -Constants<Type>::inf_repl;
  }

  Type log_lambda = safeLog(lambda_);
  Type log_alpha  = safeLog(alpha_);
  Type log_beta   = safeLog(beta_);
  Type logB_val   = logBeta(gamma_, delta_ + Type(1.0));

  Type ya = safePow(y, alpha_);
  Type one_minus_ya = enforceProbability(Type(1.0) - ya);
  Type log_one_minus_ya = safeLog(one_minus_ya);

  Type one_minus_ya_b = safePow(one_minus_ya, beta_);
  Type v = enforceProbability(Type(1.0) - one_minus_ya_b);
  Type log_v = safeLog(v);

  Type v_lambda = safePow(v, lambda_);
  Type u = enforceProbability(Type(1.0) - v_lambda);
  Type log_u = safeLog(u);

  Type term1 = log_lambda + log_alpha + log_beta - logB_val;
  Type term2 = (alpha_ - Type(1.0)) * safeLog(y);
  Type term3 = (beta_  - Type(1.0)) * log_one_minus_ya;
  Type term4 = (gamma_*lambda_ - Type(1.0)) * log_v;
  Type term5 = delta_ * log_u;

  Type logf = term1 + term2 + term3 + term4 + term5;
  if (!std::isfinite(asDouble(logf))) {
    return -Constants<Type>::inf_repl;
  }
  return logf;
}

// ========================
// Custom hash for vector<int>
// ========================
struct VectorIntHash {
  std::size_t operator()(const std::vector<int> &key) const {
    // Simple polynomial rolling hash
    std::size_t h = 0;
    for (auto &x : key) {
      h = 31*h + std::hash<int>()(x);
    }
    return h;
  }
};

struct VectorIntEq {
  bool operator()(const std::vector<int> &a, const std::vector<int> &b) const {
    if (a.size() != b.size()) return false;
    for (size_t i=0; i<a.size(); i++) {
      if (a[i] != b[i]) return false;
    }
    return true;
  }
};

// ========================
// make_cache_key
// ========================
template <class Type>
std::vector<int> make_cache_key(const Type &alpha,
                                const Type &beta,
                                const Type &gamma,
                                const Type &delta,
                                const Type &lambda) {
  std::vector<int> key(5);
  key[0] = static_cast<int>(asDouble(alpha)  * 100);
  key[1] = static_cast<int>(asDouble(beta)   * 100);
  key[2] = static_cast<int>(asDouble(gamma)  * 100);
  key[3] = static_cast<int>(asDouble(delta)  * 100);
  key[4] = static_cast<int>(asDouble(lambda) * 100);
  return key;
}

// ========================
// Numeric approximation of mean (quadrature 30 pts)
// ========================
template <class Type>
Type calc_mean_gkw(const Type &alpha,
                   const Type &beta,
                   const Type &gamma,
                   const Type &delta,
                   const Type &lambda) {
  const int n_points = 30;

  static const Type points[30] = {
    Type(0.0052995325041789), Type(0.0277124884633837), Type(0.0671843988060841),
    Type(0.1222977958224985), Type(0.1910618777986781), Type(0.2709916111713514),
    Type(0.3591982246103705), Type(0.4524937450811611), Type(0.5475062549188389),
    Type(0.6408017753896295), Type(0.7290083888286486), Type(0.8089381222013219),
    Type(0.8777022041775015), Type(0.9328156011939159), Type(0.9722875115366163),
    Type(0.9947004674958211), Type(0.0016634282895682), Type(0.0088218260005356),
    Type(0.0216951734782546), Type(0.0401505773180499), Type(0.0640198684962854),
    Type(0.0929719876996177), Type(0.1266873881927669), Type(0.1648092571058728),
    Type(0.2069463985939003), Type(0.2526493772311021), Type(0.3014937918994291),
    Type(0.3529709288365058), Type(0.4065775351876358), Type(0.4618179845446256)
  };

  static const Type weights[30] = {
    Type(0.0135762297058770), Type(0.0311267619693239), Type(0.0475792558412463),
    Type(0.0623144856277781), Type(0.0747979944082893), Type(0.0845782596975012),
    Type(0.0913017075224617), Type(0.0947253052275342), Type(0.0947253052275342),
    Type(0.0913017075224617), Type(0.0845782596975012), Type(0.0747979944082893),
    Type(0.0623144856277781), Type(0.0475792558412463), Type(0.0311267619693239),
    Type(0.0135762297058770), Type(0.0042582355019693), Type(0.0098975679009239),
    Type(0.0153793884993804), Type(0.0207860520784162), Type(0.0260583032078977),
    Type(0.0311490754242281), Type(0.0360154830389962), Type(0.0406283004740704),
    Type(0.0449535797324026), Type(0.0489611395857007), Type(0.0526254269148138),
    Type(0.0559249517732422), Type(0.0588415244791467), Type(0.0613600687415760)
  };

  Type mean(0.0), total_weight(0.0);

  for (int i = 0; i < n_points; i++) {
    Type y_i = points[i];
    if (y_i <= Constants<Type>::eps_prob || y_i >= Type(1.0)-Constants<Type>::eps_prob) {
      continue;
    }
    Type logpdf_val = log_pdf_gkw(y_i, alpha, beta, gamma, delta, lambda);
    Type pdf_val    = (logpdf_val > -Constants<Type>::max_exp) ? safeExp(logpdf_val) : Type(0.0);
    mean         += y_i * weights[i] * pdf_val;
    total_weight += weights[i]       * pdf_val;
  }

  if (total_weight > Constants<Type>::eps_pos) {
    mean /= total_weight;
  }
  return enforceProbability(mean);
}

// ========================
// objective_function
// ========================
template<class Type>
Type objective_function<Type>::operator() () {

  // == 1) DATA ==
  DATA_VECTOR(y);
  DATA_MATRIX(X1);
  DATA_MATRIX(X2);
  DATA_MATRIX(X3);
  DATA_MATRIX(X4);
  DATA_MATRIX(X5);

  DATA_INTEGER(link_type1);
  DATA_INTEGER(link_type2);
  DATA_INTEGER(link_type3);
  DATA_INTEGER(link_type4);
  DATA_INTEGER(link_type5);

  DATA_SCALAR(scale1);
  DATA_SCALAR(scale2);
  DATA_SCALAR(scale3);
  DATA_SCALAR(scale4);
  DATA_SCALAR(scale5);

  // Parâmetros opcionais de controle de performance
  DATA_INTEGER(useMeanCache);   // 0 ou 1
  DATA_INTEGER(calcFitted);     // 0 ou 1
  DATA_INTEGER(userChunkSize);  // Tamanho do escalonamento dinâmico

  // == 2) PARAMETERS ==
  PARAMETER_VECTOR(beta1);
  PARAMETER_VECTOR(beta2);
  PARAMETER_VECTOR(beta3);
  PARAMETER_VECTOR(beta4);
  PARAMETER_VECTOR(beta5);

  // == 3) VALIDATION ==
  int n = y.size();
  if (n <= 0) return Constants<Type>::inf_repl;
  // Dimensões
  if ( (X1.rows()!=n) || (X2.rows()!=n) || (X3.rows()!=n) ||
       (X4.rows()!=n) || (X5.rows()!=n) ) {
    return Constants<Type>::inf_repl;
  }
  if ( (X1.cols()!=beta1.size()) || (X2.cols()!=beta2.size()) ||
       (X3.cols()!=beta3.size()) || (X4.cols()!=beta4.size()) ||
       (X5.cols()!=beta5.size()) ) {
    return Constants<Type>::inf_repl;
  }

  // == 4) Preditores lineares e parâmetros transformados ==
  vector<Type> eta1 = X1 * beta1;
  vector<Type> eta2 = X2 * beta2;
  vector<Type> eta3 = X3 * beta3;
  vector<Type> eta4 = X4 * beta4;
  vector<Type> eta5 = X5 * beta5;

  // Alocação para parâmetros
  vector<Type> alphaVec(n), betaVec(n), gammaVec(n), deltaVec(n), lambdaVec(n);

  for (int i=0; i<n; i++) {
    alphaVec(i)  = apply_positive_link(eta1(i), link_type1, scale1);
    betaVec(i)   = apply_positive_link(eta2(i), link_type2, scale2);
    gammaVec(i)  = apply_positive_link(eta3(i), link_type3, scale3);
    deltaVec(i)  = apply_positive_link(eta4(i), link_type4, scale4);
    lambdaVec(i) = apply_positive_link(eta5(i), link_type5, scale5);
  }

  // == Vetor fitted? ==
  vector<Type> fitted;
  if (calcFitted == 1) {
    fitted.resize(n);
  }

  // == Caching das médias usando unordered_map ==
  typedef std::unordered_map< std::vector<int>, double, VectorIntHash, VectorIntEq > MeanMap;
  MeanMap meanCache;
  meanCache.reserve(std::min(n, 10000)); // heurística

  // == 5) NLL via paralelismo com schedule(dynamic) ==
  Type nll(0.0);

  // Se o tipo é double, podemos paralelizar
  // (caso contrário, TMB está no modo autodiff e não é seguro)
  if (isDouble<Type>::value) {
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, userChunkSize) reduction(+:nll)
#endif
    for (int i=0; i<n; i++) {
      Type alpha_i  = alphaVec(i);
      Type beta_i   = betaVec(i);
      Type gamma_i  = gammaVec(i);
      Type delta_i  = deltaVec(i);
      Type lambda_i = lambdaVec(i);

      // penalidade se parâmetros não são positivos
      if ( alpha_i<Constants<Type>::eps_pos || beta_i<Constants<Type>::eps_pos ||
           gamma_i<Constants<Type>::eps_pos || delta_i<Constants<Type>::eps_pos ||
           lambda_i<Constants<Type>::eps_pos ) {
        nll += Constants<Type>::inf_repl;
        continue;
      }

      Type logf = log_pdf_gkw(y(i), alpha_i, beta_i, gamma_i, delta_i, lambda_i);
      if (!std::isfinite(asDouble(logf))) {
        nll += Constants<Type>::inf_repl;
      } else {
        nll -= logf;
      }

      // Se calcFitted==1, precisamos do valor médio
      if (calcFitted == 1) {
        // Checa se vamos cachear (useMeanCache==1)
        if (useMeanCache == 1) {
          std::vector<int> key = make_cache_key(alpha_i, beta_i, gamma_i, delta_i, lambda_i);
          // Sessão crítica p/ acesso ao map
#ifdef _OPENMP
#pragma omp critical
#endif
{
  auto it = meanCache.find(key);
  if (it != meanCache.end()) {
    fitted(i) = Type(it->second);
  } else {
    Type mval = calc_mean_gkw(alpha_i, beta_i, gamma_i, delta_i, lambda_i);
    meanCache[key] = asDouble(mval);
    fitted(i) = mval;
  }
}
        } else {
          // sem cache, calcula diretamente
          fitted(i) = calc_mean_gkw(alpha_i, beta_i, gamma_i, delta_i, lambda_i);
        }
      }
    } // end for
  } else {
    // Modo autodiff: loop sequencial
    for (int i=0; i<n; i++) {
      Type alpha_i  = alphaVec(i);
      Type beta_i   = betaVec(i);
      Type gamma_i  = gammaVec(i);
      Type delta_i  = deltaVec(i);
      Type lambda_i = lambdaVec(i);

      if ( alpha_i<Constants<Type>::eps_pos || beta_i<Constants<Type>::eps_pos ||
           gamma_i<Constants<Type>::eps_pos || delta_i<Constants<Type>::eps_pos ||
           lambda_i<Constants<Type>::eps_pos ) {
        nll += Constants<Type>::inf_repl;
        continue;
      }

      Type logf = log_pdf_gkw(y(i), alpha_i, beta_i, gamma_i, delta_i, lambda_i);
      if (!std::isfinite(asDouble(logf))) {
        nll += Constants<Type>::inf_repl;
      } else {
        nll -= logf;
      }

      // Se calcFitted==1, computa fitted sem cache
      if (calcFitted == 1) {
        fitted(i) = calc_mean_gkw(alpha_i, beta_i, gamma_i, delta_i, lambda_i);
      }
    }
  }

  // == 6) Métricas ==
  int k = beta1.size() + beta2.size() + beta3.size() + beta4.size() + beta5.size();
  Type deviance = Type(2.0)*nll;
  Type aic      = deviance + Type(2.0)*Type(k);
  Type bic      = deviance + Type(k)*log(Type(n));

  // == 7) REPORT ==
  ADREPORT(beta1);
  ADREPORT(beta2);
  ADREPORT(beta3);
  ADREPORT(beta4);
  ADREPORT(beta5);

  Type alpha_mean  = alphaVec.sum()/Type(n);
  Type beta_mean   = betaVec.sum()/Type(n);
  Type gamma_mean  = gammaVec.sum()/Type(n);
  Type delta_mean  = deltaVec.sum()/Type(n);
  Type lambda_mean = lambdaVec.sum()/Type(n);

  REPORT(alpha_mean);
  REPORT(beta_mean);
  REPORT(gamma_mean);
  REPORT(delta_mean);
  REPORT(lambda_mean);

  REPORT(nll);
  REPORT(deviance);
  REPORT(aic);
  REPORT(bic);

  REPORT(alphaVec);
  REPORT(betaVec);
  REPORT(gammaVec);
  REPORT(deltaVec);
  REPORT(lambdaVec);

  if (calcFitted == 1) {
    REPORT(fitted);
  }

  // Tamanho do cache
  if (isDouble<Type>::value && (useMeanCache == 1)) {
    int cache_size = meanCache.size();
    REPORT(cache_size);
  }

  return nll;
}















// // File: gkwreg_v4_optimized.cpp
// // ---------------------------------------------------------------
// //  High Performance Regression Model for Generalized Kumaraswamy (GKw)
// //  Distribution using TMB (Template Model Builder).
// //
// //  Optimized version with residual calculation removed for better performance.
// //  Residuals are calculated separately using RcppArmadillo.
// // ---------------------------------------------------------------
//
// #include <TMB.hpp>
// #include <vector>
//
// #ifdef _OPENMP
// #include <omp.h>
// #else
// inline int omp_get_max_threads() { return 1; }
// inline int omp_get_thread_num() { return 0; }
// #define omp_parallel
// #define omp_critical
// #endif
//
// // ========================
// // Numeric constants
// // ========================
// template <class Type>
// struct Constants {
//   static const Type eps_log;
//   static const Type eps_pos;
//   static const Type eps_prob;
//   static const Type inf_repl;
//   static const Type max_exp;
//   static const Type discretize_step;
// };
//
// template <class Type>
// const Type Constants<Type>::eps_log = Type(1e-15);
// template <class Type>
// const Type Constants<Type>::eps_pos = Type(1e-10);
// template <class Type>
// const Type Constants<Type>::eps_prob = Type(1e-12);
// template <class Type>
// const Type Constants<Type>::inf_repl = Type(1e10);
// template <class Type>
// const Type Constants<Type>::max_exp = Type(30);
// template <class Type>
// const Type Constants<Type>::discretize_step = Type(0.01);
//
// // ========================
// // Safe numerical operations
// // ========================
// template <class Type>
// Type safeLog(const Type &x) {
//   if (x <= Type(0.0)) return -Constants<Type>::inf_repl;
//   return log(x + Constants<Type>::eps_log);
// }
//
// template <class Type>
// Type safeExp(const Type &x) {
//   if (x > Constants<Type>::max_exp) return exp(Constants<Type>::max_exp);
//   if (x < -Constants<Type>::max_exp) return exp(-Constants<Type>::max_exp);
//   return exp(x);
// }
//
// template <class Type>
// Type safePow(const Type &base, const Type &exponent) {
//   if (base <= Constants<Type>::eps_pos) {
//     return (exponent > Type(0.0)) ? Type(0.0) : Constants<Type>::inf_repl;
//   }
//   if (fabs(exponent) > Type(1.0)) {
//     return safeExp(exponent * safeLog(base));
//   }
//   return pow(base, exponent);
// }
//
// template <class Type>
// Type logSumExp(const Type &a, const Type &b) {
//   if (!std::isfinite(asDouble(a))) return b;
//   if (!std::isfinite(asDouble(b))) return a;
//   Type m = std::max(a, b);
//   return m + log(exp(a - m) + exp(b - m));
// }
//
// // ========================
// // Enforcement for domain
// // ========================
// template <class Type>
// Type enforceMin(const Type &x, const Type &min_val) {
//   return CppAD::CondExpLt(x, min_val, min_val, x);
// }
//
// template <class Type>
// Type enforceProbability(const Type &x) {
//   Type eps = Constants<Type>::eps_prob;
//   Type upper = Type(1.0) - eps;
//   Type clippedUp = CppAD::CondExpGt(x, upper, upper, x);
//   return CppAD::CondExpLt(clippedUp, eps, eps, clippedUp);
// }
//
// // ========================
// // Discretize function
// // ========================
// template <class Type>
// Type discretize(const Type &val) {
//   return floor(val / Constants<Type>::discretize_step) * Constants<Type>::discretize_step;
// }
//
// // ========================
// // Fast optimized link functions
// // ========================
// template <class Type>
// Type fast_log_link(const Type &eta) {
//   if (eta < -Constants<Type>::max_exp) return Constants<Type>::eps_pos;
//   if (eta > Constants<Type>::max_exp) return exp(Constants<Type>::max_exp);
//   return exp(eta);
// }
//
// template <class Type>
// Type fast_logit_link(const Type &eta) {
//   if (eta > Constants<Type>::max_exp) return Type(1.0) - Constants<Type>::eps_prob;
//   if (eta < -Constants<Type>::max_exp) return Constants<Type>::eps_prob;
//   return Type(1.0) / (Type(1.0) + exp(-eta));
// }
//
// template <class Type>
// Type fast_probit_link(const Type &eta) {
//   if (eta > Type(5.0)) return Type(1.0) - Constants<Type>::eps_prob;
//   if (eta < Type(-5.0)) return Constants<Type>::eps_prob;
//   return enforceProbability(pnorm(eta));
// }
//
// template <class Type>
// Type fast_cloglog_link(const Type &eta) {
//   if (eta > Constants<Type>::max_exp) return Type(1.0) - Constants<Type>::eps_prob;
//   if (eta < -Constants<Type>::max_exp) return Constants<Type>::eps_prob;
//   return enforceProbability(Type(1.0) - safeExp(-safeExp(eta)));
// }
//
// // ========================
// // Optimized link application for R->(0,∞)
// // ========================
// template <class Type>
// Type apply_positive_link(const Type &eta, const int &link_type, const Type &scale_factor) {
//   Type min_val = Constants<Type>::eps_pos;
//   switch(link_type) {
//   case 1:
//     return enforceMin(fast_log_link(eta), min_val);
//   case 2:
//     return enforceMin(scale_factor * fast_logit_link(eta), min_val);
//   case 3:
//     return enforceMin(scale_factor * fast_probit_link(eta), min_val);
//   case 4:
//     return enforceMin(scale_factor * (Type(0.5) + (atan(eta) / Type(M_PI))), min_val);
//   case 5:
//     return enforceMin(scale_factor * fast_cloglog_link(eta), min_val);
//   case 6:
//     return enforceMin(eta, min_val);
//   case 7:
//     return enforceMin(CppAD::CondExpGt(eta, Type(0.0), eta * eta, Type(0.0)), min_val);
//   case 8:
//     return enforceMin(Type(1.0) / enforceMin(eta, Type(1e-6)), min_val);
//   case 9:
//     return enforceMin(Type(1.0) / sqrt(enforceMin(eta, Type(1e-6))), min_val);
//   default:
//     return enforceMin(fast_log_link(eta), min_val);
//   }
// }
//
// // ========================
// // Enhanced Beta function utilities
// // ========================
// template <class Type>
// Type logBeta(const Type &a, const Type &b) {
//   Type aa = enforceMin(a, Constants<Type>::eps_pos);
//   Type bb = enforceMin(b, Constants<Type>::eps_pos);
//   if (aa > Type(100.0) && bb > Type(100.0)) {
//     return Type(0.5) * (log(Type(2.0) * M_PI) - log(aa + bb)) +
//       (aa - Type(0.5)) * log(aa) +
//       (bb - Type(0.5)) * log(bb) -
//       (aa + bb - Type(1.0)) * log(aa + bb);
//   }
//   return lgamma(aa) + lgamma(bb) - lgamma(aa + bb);
// }
//
// // ========================
// // Optimized PDF calculation (log scale)
// // ========================
// template <class Type>
// Type log_pdf_gkw(const Type &y,
//                  const Type &alpha_,
//                  const Type &beta_,
//                  const Type &gamma_,
//                  const Type &delta_,
//                  const Type &lambda_) {
//   if (y <= Constants<Type>::eps_prob || y >= Type(1.0) - Constants<Type>::eps_prob) {
//     return -Constants<Type>::inf_repl;
//   }
//   if (alpha_ <= Constants<Type>::eps_pos ||
//       beta_ <= Constants<Type>::eps_pos ||
//       gamma_ <= Constants<Type>::eps_pos ||
//       delta_ <= Constants<Type>::eps_pos ||
//       lambda_ <= Constants<Type>::eps_pos) {
//     return -Constants<Type>::inf_repl;
//   }
//
//   Type log_y = safeLog(y);
//   Type log_lambda = safeLog(lambda_);
//   Type log_alpha = safeLog(alpha_);
//   Type log_beta = safeLog(beta_);
//   Type logB = logBeta(gamma_, delta_ + Type(1.0));
//   Type ya = safePow(y, alpha_);
//   Type one_minus_ya = enforceProbability(Type(1.0) - ya);
//   Type log_one_minus_ya = safeLog(one_minus_ya);
//   Type one_minus_ya_b = safePow(one_minus_ya, beta_);
//   Type v = enforceProbability(Type(1.0) - one_minus_ya_b);
//   Type log_v = safeLog(v);
//   Type v_lambda = safePow(v, lambda_);
//   Type u = enforceProbability(Type(1.0) - v_lambda);
//   Type log_u = safeLog(u);
//
//   Type term1 = log_lambda + log_alpha + log_beta - logB;
//   Type term2 = (alpha_ - Type(1.0)) * log_y;
//   Type term3 = (beta_ - Type(1.0)) * log_one_minus_ya;
//   Type term4 = (gamma_ * lambda_ - Type(1.0)) * log_v;
//   Type term5 = delta_ * log_u;
//
//   Type logf = term1 + term2 + term3 + term4 + term5;
//   if (!std::isfinite(asDouble(logf))) {
//     return -Constants<Type>::inf_repl;
//   }
//
//   return logf;
// }
//
// // ========================
// // Parameter caching helpers
// // ========================
// template <class Type>
// std::vector<int> make_cache_key(const Type &alpha, const Type &beta,
//                                 const Type &gamma, const Type &delta,
//                                 const Type &lambda) {
//   std::vector<int> key(5);
//   key[0] = static_cast<int>(asDouble(alpha) * 100);
//   key[1] = static_cast<int>(asDouble(beta) * 100);
//   key[2] = static_cast<int>(asDouble(gamma) * 100);
//   key[3] = static_cast<int>(asDouble(delta) * 100);
//   key[4] = static_cast<int>(asDouble(lambda) * 100);
//   return key;
// }
//
// // ========================
// // Optimized mean calculation without cache interno
// // ========================
// template <class Type>
// Type calc_mean_gkw(const Type &alpha,
//                    const Type &beta,
//                    const Type &gamma,
//                    const Type &delta,
//                    const Type &lambda) {
//   const int n_points = 30;
//
//   static const Type points[30] = {
//     Type(0.0052995325041789), Type(0.0277124884633837), Type(0.0671843988060841),
//     Type(0.1222977958224985), Type(0.1910618777986781), Type(0.2709916111713514),
//     Type(0.3591982246103705), Type(0.4524937450811611), Type(0.5475062549188389),
//     Type(0.6408017753896295), Type(0.7290083888286486), Type(0.8089381222013219),
//     Type(0.8777022041775015), Type(0.9328156011939159), Type(0.9722875115366163),
//     Type(0.9947004674958211), Type(0.0016634282895682), Type(0.0088218260005356),
//     Type(0.0216951734782546), Type(0.0401505773180499), Type(0.0640198684962854),
//     Type(0.0929719876996177), Type(0.1266873881927669), Type(0.1648092571058728),
//     Type(0.2069463985939003), Type(0.2526493772311021), Type(0.3014937918994291),
//     Type(0.3529709288365058), Type(0.4065775351876358), Type(0.4618179845446256)
//   };
//
//   static const Type weights[30] = {
//     Type(0.0135762297058770), Type(0.0311267619693239), Type(0.0475792558412463),
//     Type(0.0623144856277781), Type(0.0747979944082893), Type(0.0845782596975012),
//     Type(0.0913017075224617), Type(0.0947253052275342), Type(0.0947253052275342),
//     Type(0.0913017075224617), Type(0.0845782596975012), Type(0.0747979944082893),
//     Type(0.0623144856277781), Type(0.0475792558412463), Type(0.0311267619693239),
//     Type(0.0135762297058770), Type(0.0042582355019693), Type(0.0098975679009239),
//     Type(0.0153793884993804), Type(0.0207860520784162), Type(0.0260583032078977),
//     Type(0.0311490754242281), Type(0.0360154830389962), Type(0.0406283004740704),
//     Type(0.0449535797324026), Type(0.0489611395857007), Type(0.0526254269148138),
//     Type(0.0559249517732422), Type(0.0588415244791467), Type(0.0613600687415760)
//   };
//
//   Type mean = Type(0.0);
//   Type total_weight = Type(0.0);
//
//   for (int i = 0; i < n_points; i++) {
//     Type y_i = points[i];
//     if (y_i < Constants<Type>::eps_prob || y_i > Type(1.0) - Constants<Type>::eps_prob)
//       continue;
//
//     Type logpdf_val = log_pdf_gkw(y_i, alpha, beta, gamma, delta, lambda);
//     Type pdf_val = (logpdf_val > -Constants<Type>::max_exp) ? safeExp(logpdf_val) : Type(0.0);
//     mean += y_i * weights[i] * pdf_val;
//     total_weight += weights[i] * pdf_val;
//   }
//
//   if (total_weight > Constants<Type>::eps_pos) {
//     mean /= total_weight;
//   }
//
//   return enforceProbability(mean);
// }
//
// // ========================
// // Chunking for parallel processing
// // ========================
// template <class Type>
// struct ChunkInfo {
//   int start;
//   int end;
//   Type partial_nll;
// };
//
// // ========================
// // TMB objective
// // ========================
// template<class Type>
// Type objective_function<Type>::operator() () {
//
//   // 1) DATA
//   DATA_VECTOR(y);
//   DATA_MATRIX(X1);
//   DATA_MATRIX(X2);
//   DATA_MATRIX(X3);
//   DATA_MATRIX(X4);
//   DATA_MATRIX(X5);
//   DATA_INTEGER(link_type1);
//   DATA_INTEGER(link_type2);
//   DATA_INTEGER(link_type3);
//   DATA_INTEGER(link_type4);
//   DATA_INTEGER(link_type5);
//   DATA_SCALAR(scale1);
//   DATA_SCALAR(scale2);
//   DATA_SCALAR(scale3);
//   DATA_SCALAR(scale4);
//   DATA_SCALAR(scale5);
//
//   // 2) PARAMETERS
//   PARAMETER_VECTOR(beta1);
//   PARAMETER_VECTOR(beta2);
//   PARAMETER_VECTOR(beta3);
//   PARAMETER_VECTOR(beta4);
//   PARAMETER_VECTOR(beta5);
//
//   // 3) VALIDATION
//   int n = y.size();
//   if(n <= 0) return Constants<Type>::inf_repl;
//   if((X1.rows() != n) || (X2.rows() != n) || (X3.rows() != n) ||
//      (X4.rows() != n) || (X5.rows() != n))
//     return Constants<Type>::inf_repl;
//   if((X1.cols() != beta1.size()) || (X2.cols() != beta2.size()) ||
//      (X3.cols() != beta3.size()) || (X4.cols() != beta4.size()) ||
//      (X5.cols() != beta5.size()))
//     return Constants<Type>::inf_repl;
//
//   // 4) LINEAR PREDICTORS & TRANSFORMAÇÃO DOS PARÂMETROS
//   vector<Type> eta1 = X1 * beta1;
//   vector<Type> eta2 = X2 * beta2;
//   vector<Type> eta3 = X3 * beta3;
//   vector<Type> eta4 = X4 * beta4;
//   vector<Type> eta5 = X5 * beta5;
//
//   vector<Type> alphaVec(n), betaVec(n), gammaVec(n), deltaVec(n), lambdaVec(n);
//   for (int i = 0; i < n; i++) {
//     alphaVec(i)  = apply_positive_link(eta1(i), link_type1, scale1);
//     betaVec(i)   = apply_positive_link(eta2(i), link_type2, scale2);
//     gammaVec(i)  = apply_positive_link(eta3(i), link_type3, scale3);
//     deltaVec(i)  = apply_positive_link(eta4(i), link_type4, scale4);
//     lambdaVec(i) = apply_positive_link(eta5(i), link_type5, scale5);
//   }
//
//   vector<Type> fitted(n);
//
//   // Cache para o cálculo da média
//   std::map<std::vector<int>, double> meanCache;
//
//   // 5) PARALLEL CHUNKING PARA NEGATIVE LOG-LIKELIHOOD
//   int n_threads = omp_get_max_threads();
//   int chunk_size = std::max(10, n / (n_threads * 2));
//   int n_chunks = (n + chunk_size - 1) / chunk_size;
//   std::vector<ChunkInfo<Type>> chunks(n_chunks);
//   for (int c = 0; c < n_chunks; c++) {
//     chunks[c].start = c * chunk_size;
//     chunks[c].end = std::min((c+1) * chunk_size, n);
//     chunks[c].partial_nll = Type(0.0);
//   }
//
// #ifdef _OPENMP
// #pragma omp parallel for
// #endif
//   for (int c = 0; c < n_chunks; c++) {
//     int start = chunks[c].start;
//     int end = chunks[c].end;
//     Type partial_nll = Type(0.0);
//
//     for (int i = start; i < end; i++) {
//       Type alpha_i  = alphaVec(i);
//       Type beta_i   = betaVec(i);
//       Type gamma_i  = gammaVec(i);
//       Type delta_i  = deltaVec(i);
//       Type lambda_i = lambdaVec(i);
//
//       if (alpha_i < Constants<Type>::eps_pos || beta_i < Constants<Type>::eps_pos ||
//           gamma_i < Constants<Type>::eps_pos || delta_i < Constants<Type>::eps_pos ||
//           lambda_i < Constants<Type>::eps_pos) {
//         partial_nll += Constants<Type>::inf_repl;
//         continue;
//       }
//
//       Type logf = log_pdf_gkw(y(i), alpha_i, beta_i, gamma_i, delta_i, lambda_i);
//       if (!std::isfinite(asDouble(logf))) {
//         partial_nll += Constants<Type>::inf_repl;
//       } else {
//         partial_nll -= logf;
//       }
//
//       // Geração da chave para cache
//       std::vector<int> key = make_cache_key(alpha_i, beta_i, gamma_i, delta_i, lambda_i);
//
// #ifdef _OPENMP
// #pragma omp critical
// #endif
// {
//   if (isDouble<Type>::value) {
//     auto it = meanCache.find(key);
//     if (it != meanCache.end()) {
//       fitted(i) = it->second;
//     } else {
//       Type mean_val = calc_mean_gkw(alpha_i, beta_i, gamma_i, delta_i, lambda_i);
//       meanCache[key] = asDouble(mean_val);
//       fitted(i) = mean_val;
//     }
//   } else {
//     fitted(i) = calc_mean_gkw(alpha_i, beta_i, gamma_i, delta_i, lambda_i);
//   }
// }
//     }
//     chunks[c].partial_nll = partial_nll;
//   }
//
//   Type nll = Type(0.0);
//   for (int c = 0; c < n_chunks; c++) {
//     nll += chunks[c].partial_nll;
//   }
//
//   // 6) CALCULA MÉTRICAS
//   int k = beta1.size() + beta2.size() + beta3.size() + beta4.size() + beta5.size();
//   Type deviance = Type(2.0) * nll;
//   Type aic = deviance + Type(2.0) * Type(k);
//   Type bic = deviance + Type(k) * log(Type(n));
//
//   // 7) REPORT - MINIMIZED FOR PERFORMANCE
//   ADREPORT(beta1);
//   ADREPORT(beta2);
//   ADREPORT(beta3);
//   ADREPORT(beta4);
//   ADREPORT(beta5);
//
//   Type alpha_mean = alphaVec.sum() / Type(n);
//   Type beta_mean = betaVec.sum() / Type(n);
//   Type gamma_mean = gammaVec.sum() / Type(n);
//   Type delta_mean = deltaVec.sum() / Type(n);
//   Type lambda_mean = lambdaVec.sum() / Type(n);
//
//   REPORT(alpha_mean);
//   REPORT(beta_mean);
//   REPORT(gamma_mean);
//   REPORT(delta_mean);
//   REPORT(lambda_mean);
//   REPORT(nll);
//   REPORT(deviance);
//   REPORT(aic);
//   REPORT(bic);
//
//   // Export parameter vectors for residual calculation in R
//   REPORT(alphaVec);
//   REPORT(betaVec);
//   REPORT(gammaVec);
//   REPORT(deltaVec);
//   REPORT(lambdaVec);
//   REPORT(fitted);
//
//   // Export cache size for performance monitoring
//   if (isDouble<Type>::value) {
//     int cache_size = meanCache.size();
//     REPORT(cache_size);
//   }
//
//   return nll;
// }






// // File: gkwreg_v3.cpp
// // ---------------------------------------------------------------
// //  High Performance Regression Model for Generalized Kumaraswamy (GKw)
// //  Distribution using TMB (Template Model Builder).
// //
// //  [Comentários e informações do cabeçalho permanecem inalterados]
// // ---------------------------------------------------------------
//
// #include <TMB.hpp>
// #include <vector>
//
// #ifdef _OPENMP
// #include <omp.h>
// #else
// inline int omp_get_max_threads() { return 1; }
// inline int omp_get_thread_num() { return 0; }
// #define omp_parallel
// #define omp_critical
// #endif
//
// // ========================
// // Numeric constants
// // ========================
// template <class Type>
// struct Constants {
//   static const Type eps_log;
//   static const Type eps_pos;
//   static const Type eps_prob;
//   static const Type inf_repl;
//   static const Type max_exp;
//   static const Type discretize_step;
// };
//
// template <class Type>
// const Type Constants<Type>::eps_log = Type(1e-15);
// template <class Type>
// const Type Constants<Type>::eps_pos = Type(1e-10);
// template <class Type>
// const Type Constants<Type>::eps_prob = Type(1e-12);
// template <class Type>
// const Type Constants<Type>::inf_repl = Type(1e10);
// template <class Type>
// const Type Constants<Type>::max_exp = Type(30);
// template <class Type>
// const Type Constants<Type>::discretize_step = Type(0.01);
//
// // ========================
// // Safe numerical operations
// // ========================
// template <class Type>
// Type safeLog(const Type &x) {
//   if (x <= Type(0.0)) return -Constants<Type>::inf_repl;
//   return log(x + Constants<Type>::eps_log);
// }
//
// template <class Type>
// Type safeExp(const Type &x) {
//   if (x > Constants<Type>::max_exp) return exp(Constants<Type>::max_exp);
//   if (x < -Constants<Type>::max_exp) return exp(-Constants<Type>::max_exp);
//   return exp(x);
// }
//
// template <class Type>
// Type safePow(const Type &base, const Type &exponent) {
//   if (base <= Constants<Type>::eps_pos) {
//     return (exponent > Type(0.0)) ? Type(0.0) : Constants<Type>::inf_repl;
//   }
//   if (fabs(exponent) > Type(1.0)) {
//     return safeExp(exponent * safeLog(base));
//   }
//   return pow(base, exponent);
// }
//
// template <class Type>
// Type logSumExp(const Type &a, const Type &b) {
//   if (!std::isfinite(asDouble(a))) return b;
//   if (!std::isfinite(asDouble(b))) return a;
//   Type m = std::max(a, b);
//   return m + log(exp(a - m) + exp(b - m));
// }
//
// // ========================
// // Enforcement for domain
// // ========================
// template <class Type>
// Type enforceMin(const Type &x, const Type &min_val) {
//   return CppAD::CondExpLt(x, min_val, min_val, x);
// }
//
// template <class Type>
// Type enforceProbability(const Type &x) {
//   Type eps = Constants<Type>::eps_prob;
//   Type upper = Type(1.0) - eps;
//   Type clippedUp = CppAD::CondExpGt(x, upper, upper, x);
//   return CppAD::CondExpLt(clippedUp, eps, eps, clippedUp);
// }
//
// // ========================
// // Discretize function (única definição)
// // ========================
// template <class Type>
// Type discretize(const Type &val) {
//   return floor(val / Constants<Type>::discretize_step) * Constants<Type>::discretize_step;
// }
//
// // ========================
// // Fast optimized link functions
// // ========================
// template <class Type>
// Type fast_log_link(const Type &eta) {
//   if (eta < -Constants<Type>::max_exp) return Constants<Type>::eps_pos;
//   if (eta > Constants<Type>::max_exp) return exp(Constants<Type>::max_exp);
//   return exp(eta);
// }
//
// template <class Type>
// Type fast_logit_link(const Type &eta) {
//   if (eta > Constants<Type>::max_exp) return Type(1.0) - Constants<Type>::eps_prob;
//   if (eta < -Constants<Type>::max_exp) return Constants<Type>::eps_prob;
//   return Type(1.0) / (Type(1.0) + exp(-eta));
// }
//
// template <class Type>
// Type fast_probit_link(const Type &eta) {
//   if (eta > Type(5.0)) return Type(1.0) - Constants<Type>::eps_prob;
//   if (eta < Type(-5.0)) return Constants<Type>::eps_prob;
//   return enforceProbability(pnorm(eta));
// }
//
// template <class Type>
// Type fast_cloglog_link(const Type &eta) {
//   if (eta > Constants<Type>::max_exp) return Type(1.0) - Constants<Type>::eps_prob;
//   if (eta < -Constants<Type>::max_exp) return Constants<Type>::eps_prob;
//   return enforceProbability(Type(1.0) - safeExp(-safeExp(eta)));
// }
//
// // ========================
// // Optimized link application for R->(0,∞)
// // ========================
// template <class Type>
// Type apply_positive_link(const Type &eta, const int &link_type, const Type &scale_factor) {
//   Type min_val = Constants<Type>::eps_pos;
//   switch(link_type) {
//   case 1:
//     return enforceMin(fast_log_link(eta), min_val);
//   case 2:
//     return enforceMin(scale_factor * fast_logit_link(eta), min_val);
//   case 3:
//     return enforceMin(scale_factor * fast_probit_link(eta), min_val);
//   case 4:
//     return enforceMin(scale_factor * (Type(0.5) + (atan(eta) / Type(M_PI))), min_val);
//   case 5:
//     return enforceMin(scale_factor * fast_cloglog_link(eta), min_val);
//   case 6:
//     return enforceMin(eta, min_val);
//   case 7:
//     return enforceMin(CppAD::CondExpGt(eta, Type(0.0), eta * eta, Type(0.0)), min_val);
//   case 8:
//     return enforceMin(Type(1.0) / enforceMin(eta, Type(1e-6)), min_val);
//   case 9:
//     return enforceMin(Type(1.0) / sqrt(enforceMin(eta, Type(1e-6))), min_val);
//   default:
//     return enforceMin(fast_log_link(eta), min_val);
//   }
// }
//
// // ========================
// // Enhanced Beta function utilities
// // ========================
// template <class Type>
// Type logBeta(const Type &a, const Type &b) {
//   Type aa = enforceMin(a, Constants<Type>::eps_pos);
//   Type bb = enforceMin(b, Constants<Type>::eps_pos);
//   if (aa > Type(100.0) && bb > Type(100.0)) {
//     return Type(0.5) * (log(Type(2.0) * M_PI) - log(aa + bb)) +
//       (aa - Type(0.5)) * log(aa) +
//       (bb - Type(0.5)) * log(bb) -
//       (aa + bb - Type(1.0)) * log(aa + bb);
//   }
//   return lgamma(aa) + lgamma(bb) - lgamma(aa + bb);
// }
//
// // ========================
// // Optimized PDF calculation (log scale)
// // ========================
// template <class Type>
// Type log_pdf_gkw(const Type &y,
//                  const Type &alpha_,
//                  const Type &beta_,
//                  const Type &gamma_,
//                  const Type &delta_,
//                  const Type &lambda_) {
//   if (y <= Constants<Type>::eps_prob || y >= Type(1.0) - Constants<Type>::eps_prob) {
//     return -Constants<Type>::inf_repl;
//   }
//   if (alpha_ <= Constants<Type>::eps_pos ||
//       beta_ <= Constants<Type>::eps_pos ||
//       gamma_ <= Constants<Type>::eps_pos ||
//       delta_ <= Constants<Type>::eps_pos ||
//       lambda_ <= Constants<Type>::eps_pos) {
//     return -Constants<Type>::inf_repl;
//   }
//
//   Type log_y = safeLog(y);
//   Type log_lambda = safeLog(lambda_);
//   Type log_alpha = safeLog(alpha_);
//   Type log_beta = safeLog(beta_);
//   Type logB = logBeta(gamma_, delta_ + Type(1.0));
//   Type ya = safePow(y, alpha_);
//   Type one_minus_ya = enforceProbability(Type(1.0) - ya);
//   Type log_one_minus_ya = safeLog(one_minus_ya);
//   Type one_minus_ya_b = safePow(one_minus_ya, beta_);
//   Type v = enforceProbability(Type(1.0) - one_minus_ya_b);
//   Type log_v = safeLog(v);
//   Type v_lambda = safePow(v, lambda_);
//   Type u = enforceProbability(Type(1.0) - v_lambda);
//   Type log_u = safeLog(u);
//
//   Type term1 = log_lambda + log_alpha + log_beta - logB;
//   Type term2 = (alpha_ - Type(1.0)) * log_y;
//   Type term3 = (beta_ - Type(1.0)) * log_one_minus_ya;
//   Type term4 = (gamma_ * lambda_ - Type(1.0)) * log_v;
//   Type term5 = delta_ * log_u;
//
//   Type logf = term1 + term2 + term3 + term4 + term5;
//   if (!std::isfinite(asDouble(logf))) {
//     return -Constants<Type>::inf_repl;
//   }
//
//   return logf;
// }
//
// // ========================
// // Parameter caching helpers
// // ========================
// template <class Type>
// std::vector<int> make_cache_key(const Type &alpha, const Type &beta,
//                                 const Type &gamma, const Type &delta,
//                                 const Type &lambda) {
//   std::vector<int> key(5);
//   key[0] = static_cast<int>(asDouble(alpha) * 100);
//   key[1] = static_cast<int>(asDouble(beta) * 100);
//   key[2] = static_cast<int>(asDouble(gamma) * 100);
//   key[3] = static_cast<int>(asDouble(delta) * 100);
//   key[4] = static_cast<int>(asDouble(lambda) * 100);
//   return key;
// }
//
// // ========================
// // Optimized mean calculation without cache interno
// // ========================
// template <class Type>
// Type calc_mean_gkw(const Type &alpha,
//                    const Type &beta,
//                    const Type &gamma,
//                    const Type &delta,
//                    const Type &lambda) {
//   const int n_points = 30;
//
//   static const Type points[30] = {
//     Type(0.0052995325041789), Type(0.0277124884633837), Type(0.0671843988060841),
//     Type(0.1222977958224985), Type(0.1910618777986781), Type(0.2709916111713514),
//     Type(0.3591982246103705), Type(0.4524937450811611), Type(0.5475062549188389),
//     Type(0.6408017753896295), Type(0.7290083888286486), Type(0.8089381222013219),
//     Type(0.8777022041775015), Type(0.9328156011939159), Type(0.9722875115366163),
//     Type(0.9947004674958211), Type(0.0016634282895682), Type(0.0088218260005356),
//     Type(0.0216951734782546), Type(0.0401505773180499), Type(0.0640198684962854),
//     Type(0.0929719876996177), Type(0.1266873881927669), Type(0.1648092571058728),
//     Type(0.2069463985939003), Type(0.2526493772311021), Type(0.3014937918994291),
//     Type(0.3529709288365058), Type(0.4065775351876358), Type(0.4618179845446256)
//   };
//
//   static const Type weights[30] = {
//     Type(0.0135762297058770), Type(0.0311267619693239), Type(0.0475792558412463),
//     Type(0.0623144856277781), Type(0.0747979944082893), Type(0.0845782596975012),
//     Type(0.0913017075224617), Type(0.0947253052275342), Type(0.0947253052275342),
//     Type(0.0913017075224617), Type(0.0845782596975012), Type(0.0747979944082893),
//     Type(0.0623144856277781), Type(0.0475792558412463), Type(0.0311267619693239),
//     Type(0.0135762297058770), Type(0.0042582355019693), Type(0.0098975679009239),
//     Type(0.0153793884993804), Type(0.0207860520784162), Type(0.0260583032078977),
//     Type(0.0311490754242281), Type(0.0360154830389962), Type(0.0406283004740704),
//     Type(0.0449535797324026), Type(0.0489611395857007), Type(0.0526254269148138),
//     Type(0.0559249517732422), Type(0.0588415244791467), Type(0.0613600687415760)
//   };
//
//   Type mean = Type(0.0);
//   Type total_weight = Type(0.0);
//
//   for (int i = 0; i < n_points; i++) {
//     Type y_i = points[i];
//     if (y_i < Constants<Type>::eps_prob || y_i > Type(1.0) - Constants<Type>::eps_prob)
//       continue;
//
//     Type logpdf_val = log_pdf_gkw(y_i, alpha, beta, gamma, delta, lambda);
//     Type pdf_val = (logpdf_val > -Constants<Type>::max_exp) ? safeExp(logpdf_val) : Type(0.0);
//     mean += y_i * weights[i] * pdf_val;
//     total_weight += weights[i] * pdf_val;
//   }
//
//   if (total_weight > Constants<Type>::eps_pos) {
//     mean /= total_weight;
//   }
//
//   return enforceProbability(mean);
// }
//
// // ========================
// // Chunking for parallel processing
// // ========================
// template <class Type>
// struct ChunkInfo {
//   int start;
//   int end;
//   Type partial_nll;
// };
//
// // ========================
// // TMB objective
// // ========================
// template<class Type>
// Type objective_function<Type>::operator() () {
//
//   // 1) DATA
//   DATA_VECTOR(y);
//   DATA_MATRIX(X1);
//   DATA_MATRIX(X2);
//   DATA_MATRIX(X3);
//   DATA_MATRIX(X4);
//   DATA_MATRIX(X5);
//   DATA_INTEGER(link_type1);
//   DATA_INTEGER(link_type2);
//   DATA_INTEGER(link_type3);
//   DATA_INTEGER(link_type4);
//   DATA_INTEGER(link_type5);
//   DATA_SCALAR(scale1);
//   DATA_SCALAR(scale2);
//   DATA_SCALAR(scale3);
//   DATA_SCALAR(scale4);
//   DATA_SCALAR(scale5);
//
//   // 2) PARAMETERS
//   PARAMETER_VECTOR(beta1);
//   PARAMETER_VECTOR(beta2);
//   PARAMETER_VECTOR(beta3);
//   PARAMETER_VECTOR(beta4);
//   PARAMETER_VECTOR(beta5);
//
//   // 3) VALIDATION
//   int n = y.size();
//   if(n <= 0) return Constants<Type>::inf_repl;
//   if((X1.rows() != n) || (X2.rows() != n) || (X3.rows() != n) ||
//      (X4.rows() != n) || (X5.rows() != n))
//     return Constants<Type>::inf_repl;
//   if((X1.cols() != beta1.size()) || (X2.cols() != beta2.size()) ||
//      (X3.cols() != beta3.size()) || (X4.cols() != beta4.size()) ||
//      (X5.cols() != beta5.size()))
//     return Constants<Type>::inf_repl;
//
//   // 4) LINEAR PREDICTORS & TRANSFORMAÇÃO DOS PARÂMETROS
//   vector<Type> eta1 = X1 * beta1;
//   vector<Type> eta2 = X2 * beta2;
//   vector<Type> eta3 = X3 * beta3;
//   vector<Type> eta4 = X4 * beta4;
//   vector<Type> eta5 = X5 * beta5;
//
//   vector<Type> alphaVec(n), betaVec(n), gammaVec(n), deltaVec(n), lambdaVec(n);
//   for (int i = 0; i < n; i++) {
//     alphaVec(i)  = apply_positive_link(eta1(i), link_type1, scale1);
//     betaVec(i)   = apply_positive_link(eta2(i), link_type2, scale2);
//     gammaVec(i)  = apply_positive_link(eta3(i), link_type3, scale3);
//     deltaVec(i)  = apply_positive_link(eta4(i), link_type4, scale4);
//     lambdaVec(i) = apply_positive_link(eta5(i), link_type5, scale5);
//   }
//
//   vector<Type> fitted(n);
//   vector<Type> residuals(n);
//
//   // Cache para o cálculo da média
//   std::map<std::vector<int>, double> meanCache;
//
//   // 5) PARALLEL CHUNKING PARA NEGATIVE LOG-LIKELIHOOD
//   int n_threads = omp_get_max_threads();
//   int chunk_size = std::max(10, n / (n_threads * 2));
//   int n_chunks = (n + chunk_size - 1) / chunk_size;
//   std::vector<ChunkInfo<Type>> chunks(n_chunks);
//   for (int c = 0; c < n_chunks; c++) {
//     chunks[c].start = c * chunk_size;
//     chunks[c].end = std::min((c+1) * chunk_size, n);
//     chunks[c].partial_nll = Type(0.0);
//   }
//
// #ifdef _OPENMP
// #pragma omp parallel for
// #endif
//   for (int c = 0; c < n_chunks; c++) {
//     int start = chunks[c].start;
//     int end = chunks[c].end;
//     Type partial_nll = Type(0.0);
//
//     for (int i = start; i < end; i++) {
//       Type alpha_i  = alphaVec(i);
//       Type beta_i   = betaVec(i);
//       Type gamma_i  = gammaVec(i);
//       Type delta_i  = deltaVec(i);
//       Type lambda_i = lambdaVec(i);
//
//       if (alpha_i < Constants<Type>::eps_pos || beta_i < Constants<Type>::eps_pos ||
//           gamma_i < Constants<Type>::eps_pos || delta_i < Constants<Type>::eps_pos ||
//           lambda_i < Constants<Type>::eps_pos) {
//         partial_nll += Constants<Type>::inf_repl;
//         continue;
//       }
//
//       Type logf = log_pdf_gkw(y(i), alpha_i, beta_i, gamma_i, delta_i, lambda_i);
//       if (!std::isfinite(asDouble(logf))) {
//         partial_nll += Constants<Type>::inf_repl;
//       } else {
//         partial_nll -= logf;
//       }
//
//       // Geração da chave para cache (sem cast desnecessário)
//       std::vector<int> key = make_cache_key(alpha_i, beta_i, gamma_i, delta_i, lambda_i);
//
// #ifdef _OPENMP
// #pragma omp critical
// #endif
// {
//   if (isDouble<Type>::value) {
//     auto it = meanCache.find(key);
//     if (it != meanCache.end()) {
//       fitted(i) = it->second;
//     } else {
//       Type mean_val = calc_mean_gkw(alpha_i, beta_i, gamma_i, delta_i, lambda_i);
//       meanCache[key] = asDouble(mean_val);
//       fitted(i) = mean_val;
//     }
//   } else {
//     fitted(i) = calc_mean_gkw(alpha_i, beta_i, gamma_i, delta_i, lambda_i);
//   }
// }
//
// residuals(i) = y(i) - fitted(i);
//     }
//     chunks[c].partial_nll = partial_nll;
//   }
//
//   Type nll = Type(0.0);
//   for (int c = 0; c < n_chunks; c++) {
//     nll += chunks[c].partial_nll;
//   }
//
//   // 6) CALCULA MÉTRICAS
//   int k = beta1.size() + beta2.size() + beta3.size() + beta4.size() + beta5.size();
//   Type deviance = Type(2.0) * nll;
//   Type aic = deviance + Type(2.0) * Type(k);
//   Type bic = deviance + Type(k) * log(Type(n));
//   Type sum_sq_residuals = (residuals * residuals).sum();
//   Type rmse = sqrt(sum_sq_residuals / Type(n));
//   Type y_mean = y.sum() / Type(n);
//   Type sst = Type(0.0);
//   for (int i = 0; i < n; i++) {
//     Type diff = y(i) - y_mean;
//     sst += diff * diff;
//   }
//   Type efron_r2 = Type(1.0) - (sum_sq_residuals / sst);
//   Type mean_absolute_error = residuals.abs().sum() / Type(n);
//
//   // 7) REPORT
//   ADREPORT(beta1);
//   ADREPORT(beta2);
//   ADREPORT(beta3);
//   ADREPORT(beta4);
//   ADREPORT(beta5);
//   ADREPORT(rmse);
//   ADREPORT(efron_r2);
//   ADREPORT(mean_absolute_error);
//
//   Type alpha_mean = alphaVec.sum() / Type(n);
//   Type beta_mean = betaVec.sum() / Type(n);
//   Type gamma_mean = gammaVec.sum() / Type(n);
//   Type delta_mean = deltaVec.sum() / Type(n);
//   Type lambda_mean = lambdaVec.sum() / Type(n);
//
//   REPORT(alpha_mean);
//   REPORT(beta_mean);
//   REPORT(gamma_mean);
//   REPORT(delta_mean);
//   REPORT(lambda_mean);
//   REPORT(nll);
//   REPORT(deviance);
//   REPORT(aic);
//   REPORT(bic);
//   REPORT(rmse);
//   REPORT(efron_r2);
//
//   ADREPORT(fitted);
//   ADREPORT(residuals);
//
//   if (isDouble<Type>::value) {
//     int cache_size = meanCache.size();
//     REPORT(cache_size);
//   }
//
//   return nll;
// }







// // File: gkwreg_v2.cpp
// // ---------------------------------------------------------------
// //  Optimized Regression Model for Generalized Kumaraswamy (GKw)
// //  Distribution using TMB (Template Model Builder).
// //
// //  Distribution PDF (for y in (0,1)):
// //    f(y; alpha, beta, gamma, delta, lambda) =
// //      (lambda * alpha * beta / B(gamma, delta+1)) * y^(alpha - 1) * (1 - y^alpha)^(beta - 1) *
// //      [1 - (1 - y^alpha)^beta]^(gamma * lambda - 1) * {1 - [1 - (1 - y^alpha)^beta]^lambda}^delta
// //
// //  Five parameters (alpha, beta, gamma, delta, lambda) are each modeled
// //  via linear predictors (eta) + link functions for positivity.
// //
// //  In this version (v2):
// //   - The expensive integral for the mean is removed from the main loop
// //     to avoid high computational cost in large datasets.
// //   - Additional domain checks and reports have been streamlined.
// //   - Comments and code are all in English, with inline formula references.
// //
// // ---------------------------------------------------------------
//
// #include <TMB.hpp>
//
// // ========================
// // Numeric constants
// // ========================
// template <class Type>
// struct Constants {
//   static const Type eps_log;    // small constant for log
//   static const Type eps_pos;    // small constant for positivity constraints
//   static const Type eps_prob;   // small constant for probability boundaries
//   static const Type inf_repl;   // large penalty for invalid values
//   static const Type max_exp;    // max exponent to avoid overflow in exp
// };
//
// template <class Type>
// const Type Constants<Type>::eps_log = Type(1e-15);
//
// template <class Type>
// const Type Constants<Type>::eps_pos = Type(1e-10);
//
// template <class Type>
// const Type Constants<Type>::eps_prob = Type(1e-12);
//
// template <class Type>
// const Type Constants<Type>::inf_repl = Type(1e10);
//
// template <class Type>
// const Type Constants<Type>::max_exp = Type(30);  // Increased from 20 to allow bigger exponent range
//
// // ========================
// // Safe log and safe exp
// // ========================
// template <class Type>
// Type safeLog(const Type &x) {
//   // log(x + eps_log) helps avoid log(0)
//   return log(x + Constants<Type>::eps_log);
// }
//
// template <class Type>
// Type safeExp(const Type &x) {
//   // clips x to [-max_exp, max_exp]
//   Type val = CppAD::CondExpGt(x,  Constants<Type>::max_exp,  Constants<Type>::max_exp,
//                               CppAD::CondExpLt(x, -Constants<Type>::max_exp, -Constants<Type>::max_exp, x));
//   return exp(val);
// }
//
// // ========================
// // Enforcement for domain
// // ========================
// template <class Type>
// Type enforceMin(const Type &x, const Type &min_val) {
//   // ensures x >= min_val
//   return CppAD::CondExpLt(x, min_val, min_val, x);
// }
//
// template <class Type>
// Type enforceProbability(const Type &x) {
//   // ensures 0 < x < 1 by clipping
//   Type eps = Constants<Type>::eps_prob;
//   Type upper = Type(1.0) - eps;
//   Type clippedUp = CppAD::CondExpGt(x, upper, upper, x);
//   return CppAD::CondExpLt(clippedUp, eps, eps, clippedUp);
// }
//
// // ========================
// // Link functions: R -> (0,1)
// // ========================
// template <class Type>
// Type logit_link(const Type &x) {
//   // 1 / (1 + exp(-x))
//   return Type(1.0) / (Type(1.0) + safeExp(-x));
// }
//
// template <class Type>
// Type probit_link(const Type &x) {
//   // pnorm(x) clipped to (0,1)
//   return enforceProbability(pnorm(x));
// }
//
// template <class Type>
// Type cauchy_link(const Type &x) {
//   // 0.5 + (1/pi)*atan(x)
//   return enforceProbability(Type(0.5) + (atan(x) / Type(M_PI)));
// }
//
// template <class Type>
// Type cloglog_link(const Type &x) {
//   // 1 - exp(-exp(x))
//   return enforceProbability(Type(1.0) - safeExp(-safeExp(x)));
// }
//
// // ========================
// // Link application for R->(0,∞)
// // including scale factors for links that yield (0,1)
// // ========================
// template <class Type>
// Type apply_positive_link(const Type &eta, const int &link_type, const Type &scale_factor) {
//   // minimal positive value
//   Type min_val = Constants<Type>::eps_pos;
//
//   switch(link_type) {
//   case 1: // log link: maps R -> (0,∞)
//     return safeExp(eta);
//
//   case 2: // logit link -> (0,1), then scaled to (0,∞)
//     return enforceMin(scale_factor * logit_link(eta), min_val);
//
//   case 3: // probit link -> (0,1), then scaled to (0,∞)
//     return enforceMin(scale_factor * probit_link(eta), min_val);
//
//   case 4: // cauchy link -> (0,1), then scaled
//     return enforceMin(scale_factor * cauchy_link(eta), min_val);
//
//   case 5: // cloglog link -> (0,1), then scaled
//     return enforceMin(scale_factor * cloglog_link(eta), min_val);
//
//   case 6: // identity link (ensure positivity)
//     return enforceMin(eta, min_val);
//
//   case 7: // sqrt link: g^(-1)(eta) = eta^2, ensure non-negative
//     {
//       // if eta <= 0, we clamp to 0, then clamp to min_val
//       Type sq = CppAD::CondExpGt(eta, Type(0.0), eta * eta, Type(0.0));
//       return enforceMin(sq, min_val);
//     }
//
//   case 8: // inverse link: g^(-1)(eta) = 1/eta
//     {
//       // ensure denominator is positive
//       Type safe_den = enforceMin(eta, Type(1e-6));
//       return enforceMin(Type(1.0) / safe_den, min_val);
//     }
//
//   case 9: // inverse-square link: g^(-1)(eta) = 1 / sqrt(eta)
//     {
//       // ensure denominator is positive
//       Type safe_den = enforceMin(eta, Type(1e-6));
//       return enforceMin(Type(1.0) / sqrt(safe_den), min_val);
//     }
//
//   default:
//     // fallback: log link
//     return safeExp(eta);
//   }
// }
//
// // ========================
// // Beta function utilities
// // ========================
// template <class Type>
// Type logBeta(const Type &a, const Type &b) {
//   // B(a,b) = Gamma(a)*Gamma(b)/Gamma(a+b)
//   // we use safe min to avoid log(0)
//   Type aa = enforceMin(a, Constants<Type>::eps_pos);
//   Type bb = enforceMin(b, Constants<Type>::eps_pos);
//   return lgamma(aa) + lgamma(bb) - lgamma(aa + bb);
// }
//
// // ========================
// // PDF of GKw (log scale)
// // log f(y; alpha, beta, gamma, delta, lambda)
// // ========================
// template <class Type>
// Type log_pdf_gkw(const Type &y,
//                  const Type &alpha_,
//                  const Type &beta_,
//                  const Type &gamma_,
//                  const Type &delta_,
//                  const Type &lambda_) {
//   /*
//    * PDF formula in log-scale:
//    *   log f(y; α, β, γ, δ, λ) =
//    *     log(λ) + log(α) + log(β) - log[B(γ, δ+1)]
//    *     + (α - 1)*log(y) + (β - 1)*log(1 - y^α)
//    *     + (γλ - 1)*log[1 - (1 - y^α)^β]
//    *     + δ*log{1 - [1 - (1 - y^α)^β]^λ}
//    */
//
//   // check y in (0,1)
//   if( (y <= Type(0.0)) || (y >= Type(1.0)) ) {
//     // invalid domain -> logf is -∞ in practice
//     return -Constants<Type>::inf_repl;
//   }
//
//   // safely compute powers and logs
//   Type ya = pow(y, alpha_);                 // y^α
//   Type one_minus_ya = enforceProbability(Type(1.0) - ya);  // 1 - y^α
//   Type one_minus_ya_b = pow(one_minus_ya, beta_);           // (1 - y^α)^β
//   one_minus_ya_b = enforceProbability(one_minus_ya_b);
//   Type v = enforceProbability(Type(1.0) - one_minus_ya_b);  // 1 - (1 - y^α)^β
//   Type v_lambda = pow(v, lambda_);                          // v^λ
//   v_lambda = enforceProbability(v_lambda);
//   Type u = enforceProbability(Type(1.0) - v_lambda);         // 1 - v^λ
//
//   // log-beta(gamma, delta+1)
//   Type logB = logBeta(gamma_, delta_ + Type(1.0));
//
//   // assemble logf
//   Type logf = safeLog(lambda_) + safeLog(alpha_) + safeLog(beta_) - logB
//   + (alpha_ - Type(1.0))*safeLog(y)
//     + (beta_  - Type(1.0))*safeLog(one_minus_ya)
//     + (gamma_*lambda_ - Type(1.0))*safeLog(v)
//     + (delta_)*safeLog(u);
//
//     return logf;
// }
//
// // ========================
// // Calculate mean of GKw distribution using numerical integration
// // ========================
// template <class Type>
// Type calc_mean_gkw(const Type &alpha,
//                    const Type &beta,
//                    const Type &gamma,
//                    const Type &delta,
//                    const Type &lambda) {
//   // Use Gauss-Legendre quadrature for numerical integration with 30 points
//   // This provides higher accuracy than the previous 15-point quadrature
//
//   // Number of quadrature points (30 points for higher accuracy)
//   const int n_points = 30;
//
//   // Gauss-Legendre quadrature points for the interval [0,1]
//   Type points[30] = {
//     Type(0.0052995325041789), Type(0.0277124884633837), Type(0.0671843988060841),
//     Type(0.1222977958224985), Type(0.1910618777986781), Type(0.2709916111713514),
//     Type(0.3591982246103705), Type(0.4524937450811611), Type(0.5475062549188389),
//     Type(0.6408017753896295), Type(0.7290083888286486), Type(0.8089381222013219),
//     Type(0.8777022041775015), Type(0.9328156011939159), Type(0.9722875115366163),
//     Type(0.9947004674958211), Type(0.0016634282895682), Type(0.0088218260005356),
//     Type(0.0216951734782546), Type(0.0401505773180499), Type(0.0640198684962854),
//     Type(0.0929719876996177), Type(0.1266873881927669), Type(0.1648092571058728),
//     Type(0.2069463985939003), Type(0.2526493772311021), Type(0.3014937918994291),
//     Type(0.3529709288365058), Type(0.4065775351876358), Type(0.4618179845446256)
//   };
//
//   // Gauss-Legendre quadrature weights for the interval [0,1]
//   Type weights[30] = {
//     Type(0.0135762297058770), Type(0.0311267619693239), Type(0.0475792558412463),
//     Type(0.0623144856277781), Type(0.0747979944082893), Type(0.0845782596975012),
//     Type(0.0913017075224617), Type(0.0947253052275342), Type(0.0947253052275342),
//     Type(0.0913017075224617), Type(0.0845782596975012), Type(0.0747979944082893),
//     Type(0.0623144856277781), Type(0.0475792558412463), Type(0.0311267619693239),
//     Type(0.0135762297058770), Type(0.0042582355019693), Type(0.0098975679009239),
//     Type(0.0153793884993804), Type(0.0207860520784162), Type(0.0260583032078977),
//     Type(0.0311490754242281), Type(0.0360154830389962), Type(0.0406283004740704),
//     Type(0.0449535797324026), Type(0.0489611395857007), Type(0.0526254269148138),
//     Type(0.0559249517732422), Type(0.0588415244791467), Type(0.0613600687415760)
//   };
//
//   Type mean = Type(0.0);
//
//   // For each quadrature point, calculate density and accumulate weighted contribution
//   for (int i = 0; i < n_points; i++) {
//     Type y_i = points[i];
//
//     // Skip calculation if y_i is too close to 0 or 1 to avoid numerical issues
//     if (y_i < Constants<Type>::eps_prob || y_i > Type(1.0) - Constants<Type>::eps_prob) {
//       continue;
//     }
//
//     // Calculate PDF value at y_i (normalized so we get proper weights)
//     Type pdf_val = safeExp(log_pdf_gkw(y_i, alpha, beta, gamma, delta, lambda));
//
//     // Add weighted contribution to mean: y_i * weight_i * pdf_val
//     mean += y_i * weights[i] * pdf_val;
//   }
//
//   // We need to multiply by the normalizing constant
//   // In theory, PDF should already be normalized, but for numerical stability
//   // we normalize again by the total integrated value
//   Type total_weight = Type(0.0);
//   for (int i = 0; i < n_points; i++) {
//     Type y_i = points[i];
//
//     if (y_i < Constants<Type>::eps_prob || y_i > Type(1.0) - Constants<Type>::eps_prob) {
//       continue;
//     }
//
//     Type pdf_val = safeExp(log_pdf_gkw(y_i, alpha, beta, gamma, delta, lambda));
//     total_weight += weights[i] * pdf_val;
//   }
//
//   // Normalize by total weight
//   if (total_weight > Constants<Type>::eps_pos) {
//     mean /= total_weight;
//   }
//
//   // Ensure the result is in (0,1)
//   return enforceProbability(mean);
// }
//
// // ========================
// // TMB objective
// // ========================
// template<class Type>
// Type objective_function<Type>::operator() () {
//
//   // ------------------------------------------------------
//   // 1) DATA
//   // ------------------------------------------------------
//   DATA_VECTOR(y); // response, must be in (0,1)
//
//   // design matrices for each parameter
//   DATA_MATRIX(X1); // alpha
//   DATA_MATRIX(X2); // beta
//   DATA_MATRIX(X3); // gamma
//   DATA_MATRIX(X4); // delta
//   DATA_MATRIX(X5); // lambda
//
//   // link function types for each parameter (1=log,2=logit,3=probit,4=cauchy,5=cloglog,6=identity,7=sqrt,8=inverse,9=inverse-square)
//   DATA_INTEGER(link_type1); // alpha
//   DATA_INTEGER(link_type2); // beta
//   DATA_INTEGER(link_type3); // gamma
//   DATA_INTEGER(link_type4); // delta
//   DATA_INTEGER(link_type5); // lambda
//
//   // scale factors for links that map (R->(0,1)->(0,∞))
//   DATA_SCALAR(scale1);
//   DATA_SCALAR(scale2);
//   DATA_SCALAR(scale3);
//   DATA_SCALAR(scale4);
//   DATA_SCALAR(scale5);
//
//   // ------------------------------------------------------
//   // 2) PARAMETERS
//   // ------------------------------------------------------
//   PARAMETER_VECTOR(beta1); // coefficients for alpha
//   PARAMETER_VECTOR(beta2); // coefficients for beta
//   PARAMETER_VECTOR(beta3); // coefficients for gamma
//   PARAMETER_VECTOR(beta4); // coefficients for delta
//   PARAMETER_VECTOR(beta5); // coefficients for lambda
//
//   // ------------------------------------------------------
//   // 3) VALIDATION
//   // ------------------------------------------------------
//   int n = y.size();
//   if (n <= 0) {
//     // no data -> infinite penalty
//     return Constants<Type>::inf_repl;
//   }
//
//   // check design matrices match n
//   if (   (X1.rows() != n) || (X2.rows() != n)
//            || (X3.rows() != n) || (X4.rows() != n)
//            || (X5.rows() != n) ) {
//            return Constants<Type>::inf_repl;
//   }
//
//   // check param vector sizes match # of columns
//   if ( (X1.cols() != beta1.size()) || (X2.cols() != beta2.size()) ||
//        (X3.cols() != beta3.size()) || (X4.cols() != beta4.size()) ||
//        (X5.cols() != beta5.size()) ) {
//     return Constants<Type>::inf_repl;
//   }
//
//   // ------------------------------------------------------
//   // 4) LINEAR PREDICTORS & PARAM TRANSFORM
//   // ------------------------------------------------------
//   // linear predictors
//   vector<Type> eta1 = X1 * beta1;  // alpha
//   vector<Type> eta2 = X2 * beta2;  // beta
//   vector<Type> eta3 = X3 * beta3;  // gamma
//   vector<Type> eta4 = X4 * beta4;  // delta
//   vector<Type> eta5 = X5 * beta5;  // lambda
//
//   // create vectors of distribution parameters
//   vector<Type> alphaVec(n), betaVec(n), gammaVec(n), deltaVec(n), lambdaVec(n);
//
//   // for fitted values and residuals
//   vector<Type> fitted(n);
//   vector<Type> residuals(n);
//
//   // Vector to cache mean calculations
//   vector<Type> meanCache(n);
//   bool meansCached = false;
//
//   // Apply link transformations first
//   for(int i=0; i<n; i++){
//     alphaVec(i)  = apply_positive_link(eta1(i), link_type1, scale1);
//     betaVec(i)   = apply_positive_link(eta2(i), link_type2, scale2);
//     gammaVec(i)  = apply_positive_link(eta3(i), link_type3, scale3);
//     deltaVec(i)  = apply_positive_link(eta4(i), link_type4, scale4);
//     lambdaVec(i) = apply_positive_link(eta5(i), link_type5, scale5);
//   }
//
//   // Precompute lookup table for fitted values - doing this approach to improve performance
//   // Only for unique parameter combinations to avoid redundant calculations
//   if (!meansCached) {
//     parallel_for(0, n, [&](int i) {
//       // Calculate fitted value (mean of GKw distribution with these parameters)
//       fitted(i) = calc_mean_gkw(alphaVec(i), betaVec(i), gammaVec(i), deltaVec(i), lambdaVec(i));
//
//       // Calculate residual (observed - fitted)
//       residuals(i) = y(i) - fitted(i);
//     });
//     meansCached = true;
//   }
//
//   // ------------------------------------------------------
//   // 5) NEGATIVE LOG-LIKELIHOOD
//   // ------------------------------------------------------
//   Type nll = 0.0;
//
//   // Loop over observations
//   for(int i=0; i<n; i++){
//     // extract parameters for observation i
//     Type alpha_i  = alphaVec(i);
//     Type beta_i   = betaVec(i);
//     Type gamma_i  = gammaVec(i);
//     Type delta_i  = deltaVec(i);
//     Type lambda_i = lambdaVec(i);
//
//     // if any param is below eps_pos, penalize
//     if (   alpha_i  < Constants<Type>::eps_pos
//              || beta_i   < Constants<Type>::eps_pos
//              || gamma_i  < Constants<Type>::eps_pos
//              || delta_i  < Constants<Type>::eps_pos
//              || lambda_i < Constants<Type>::eps_pos ) {
//              nll += Constants<Type>::inf_repl;
//       continue;
//     }
//
//     // compute log-pdf
//     Type logf = log_pdf_gkw(y(i), alpha_i, beta_i, gamma_i, delta_i, lambda_i);
//
//     // if logf is -inf_repl => invalid
//     if(!R_FINITE(asDouble(logf))){
//       nll += Constants<Type>::inf_repl;
//     } else {
//       nll -= logf; // accumulate negative log-likelihood
//     }
//   }
//
//   // ------------------------------------------------------
//   // 6) REPORT & ADREPORT
//   // ------------------------------------------------------
//   // We removed the ADREPORT statements as requested
//
//   // Calculate metrics for model fit evaluation
//   Type deviance = Type(2.0) * nll;
//   // number of regression coefficients
//   int k = beta1.size() + beta2.size() + beta3.size() + beta4.size() + beta5.size();
//   Type aic = deviance + Type(2.0) * Type(k);
//
//   // Calculate BIC: BIC = -2*log-likelihood + k*log(n)
//   Type bic = deviance + Type(k) * log(Type(n));
//
//   // Calculate additional goodness-of-fit metrics
//
//   // Root Mean Square Error (RMSE)
//   Type sum_sq_residuals = Type(0.0);
//   for(int i=0; i<n; i++) {
//     sum_sq_residuals += residuals(i) * residuals(i);
//   }
//   Type rmse = sqrt(sum_sq_residuals / Type(n));
//
//   // Calculate mean of observed values for R² calculation
//   Type y_mean = Type(0.0);
//   for(int i=0; i<n; i++) {
//     y_mean += y(i);
//   }
//   y_mean /= Type(n);
//
//   // Calculate total sum of squares (SST) for R² calculation
//   Type sst = Type(0.0);
//   for(int i=0; i<n; i++) {
//     Type diff = y(i) - y_mean;
//     sst += diff * diff;
//   }
//
//   // Efron's Pseudo R-Squared: 1 - [Σ(yi-πˆi)²]/[Σ(yi-ȳ)²]
//   Type efron_r2 = Type(1.0) - (sum_sq_residuals / sst);
//
//   // McFadden's Pseudo R-Squared: 1 - [ln LL(Mˆfull)]/[ln LL(Mˆintercept)]
//   // We need to calculate the log-likelihood of an intercept-only model
//
//   // First, calculate the mean of y as estimate for intercept-only model
//   // The intercept-only model would predict the same value for all observations
//   Type intercept_ll = Type(0.0);
//
//   // For GKw distribution, we'd need to find the parameters that maximize
//   // likelihood when all observations get same parameters (intercept-only)
//   // This is an approximation using the log-likelihood of predicting the mean
//   for(int i=0; i<n; i++) {
//     intercept_ll -= log_pdf_gkw(y(i), Type(1.0), Type(1.0), Type(1.0), Type(1.0), Type(1.0));
//   }
//
//   // McFadden's Pseudo R-Squared
//   Type mcfadden_r2 = Type(1.0) - (nll / intercept_ll);
//
//   // McFadden's Adjusted Pseudo R-Squared: 1 - [ln LL(Mˆfull)-K]/[ln LL(Mˆintercept)]
//   Type mcfadden_r2_adj = Type(1.0) - ((nll - Type(k)) / intercept_ll);
//
//   // Report metrics
//   REPORT(nll);
//   REPORT(deviance);
//   REPORT(aic);
//   REPORT(bic);
//   REPORT(fitted);  // Renamed from fitted_values to fitted
//   REPORT(residuals);
//   REPORT(rmse);
//   REPORT(efron_r2);
//   REPORT(mcfadden_r2);
//   REPORT(mcfadden_r2_adj);
//
//   // Return negative log-likelihood
//   return nll;
// }


// // File: gkwreg.cpp
// // Generalized Kumaraswamy (GKw) Distribution Regression Model with TMB
// //
// // This implements a regression model for the five-parameter GKw distribution:
// // f(y; α, β, γ, δ, λ) = (λαβ/B(γ,δ+1)) * y^(α-1) * (1-y^α)^(β-1) *
// //                      [1-(1-y^α)^β]^(γλ-1) * {1-[1-(1-y^α)^β]^λ}^δ
// //
// // All five parameters (α, β, γ, δ, λ) can be modeled with covariates
// // using flexible link functions.
// //
// // Available link functions:
// // 1. log (default): R -> (0,∞), g^(-1)(η) = exp(η)
// // 2. logit: R -> (0,1) -> scale to (0,∞)
// // 3. probit: R -> (0,1) -> scale to (0,∞)
// // 4. cauchy: R -> (0,1) -> scale to (0,∞)
// // 5. cloglog: R -> (0,1) -> scale to (0,∞)
// // 6. identity: R -> R, constrained to (0,∞)
// // 7. sqrt: R -> (0,∞), g^(-1)(η) = η²
// // 8. inverse: R -> (0,∞), g^(-1)(η) = 1/η
// // 9. inverse-square: R -> (0,∞), g^(-1)(η) = 1/√η
//
// #include <TMB.hpp>
//
// // ========================
// // Auxiliary functions
// // ========================
//
// // Constants for numerical stability
// template <class Type>
// struct Constants {
//   static const Type eps_log;   // Small constant for log
//   static const Type eps_pos;   // Small constant for positivity constraint
//   static const Type eps_prob;  // Small constant for probability boundary
//   static const Type inf_repl;  // Replacement for infinity in penalties
//   static const Type max_exp;   // Maximum value for exp to avoid overflow
//   static const Type max_iter;  // Maximum iterations for numerical calculation
//   static const Type tol;       // Tolerance for numerical calculation
// };
//
// template <class Type>
// const Type Constants<Type>::eps_log = Type(1.0) / Type(1000000000000000.0); // 1e-15
//
// template <class Type>
// const Type Constants<Type>::eps_pos = Type(1.0) / Type(10000000000.0); // 1e-10
//
// template <class Type>
// const Type Constants<Type>::eps_prob = Type(1.0) / Type(1000000000000.0); // 1e-12
//
// template <class Type>
// const Type Constants<Type>::inf_repl = Type(10000000000.0); // 1e10
//
// template <class Type>
// const Type Constants<Type>::max_exp = Type(20.0);
//
// template <class Type>
// const Type Constants<Type>::max_iter = Type(1000.0);
//
// template <class Type>
// const Type Constants<Type>::tol = Type(1.0) / Type(100000000.0); // 1e-8
//
// // Function to avoid numerical issues with log(x)
// template <class Type>
// Type safeLog(const Type &x) {
//   return log(x + Constants<Type>::eps_log);
// }
//
// // Function to safely compute exp(x) without overflow
// template <class Type>
// Type safeExp(const Type &x) {
//   // Clip x to avoid exp overflow
//   Type val = x;
//   if (x > Constants<Type>::max_exp) val = Constants<Type>::max_exp;
//   if (x < -Constants<Type>::max_exp) val = -Constants<Type>::max_exp;
//   return exp(val);
// }
//
// // Function to safely ensure x > min_val
// template <class Type>
// Type enforceMin(const Type &x, const Type &min_val) {
//   return CppAD::CondExpGt(x, min_val, x, min_val);
// }
//
// // Function to safely ensure x is within probability range (0,1)
// template <class Type>
// Type enforceProbability(const Type &x) {
//   Type eps = Constants<Type>::eps_prob;
//   Type val = CppAD::CondExpLt(x, Type(1.0) - eps, x, Type(1.0) - eps);
//   return CppAD::CondExpGt(val, eps, val, eps);
// }
//
// // ========================
// // Link functions
// // ========================
//
// // LOGIT link function: R -> (0,1)
// template <class Type>
// Type logit_link(const Type &x) {
//   // Using safe exp to avoid overflow
//   return Type(1.0) / (Type(1.0) + safeExp(-x));
// }
//
// // PROBIT link function: R -> (0,1)
// template <class Type>
// Type probit_link(const Type &x) {
//   return enforceProbability(pnorm(x));
// }
//
// // CAUCHY link function: R -> (0,1)
// template <class Type>
// Type cauchy_link(const Type &x) {
//   return enforceProbability(Type(0.5) + atan(x) / Type(M_PI));
// }
//
// // COMPLEMENTARY LOG-LOG link function: R -> (0,1)
// template <class Type>
// Type cloglog_link(const Type &x) {
//   return enforceProbability(Type(1.0) - safeExp(-safeExp(x)));
// }
//
// // Apply link function for positive parameters
// template <class Type>
// Type apply_positive_link(const Type &x, const int &link_type, const Type &scale_factor) {
//   Type min_val = Constants<Type>::eps_pos;
//
//   switch(link_type) {
//   case 1: // LOG link (default): maps R -> (0,∞)
//     return safeExp(x);
//   case 2: // LOGIT: maps R -> (0,1) -> scale to (0,∞)
//     return enforceMin(scale_factor * logit_link(x), min_val);
//   case 3: // PROBIT: maps R -> (0,1) -> scale to (0,∞)
//     return enforceMin(scale_factor * probit_link(x), min_val);
//   case 4: // CAUCHY: maps R -> (0,1) -> scale to (0,∞)
//     return enforceMin(scale_factor * cauchy_link(x), min_val);
//   case 5: // CLOGLOG: maps R -> (0,1) -> scale to (0,∞)
//     return enforceMin(scale_factor * cloglog_link(x), min_val);
//   case 6: // IDENTITY: ensure positivity
//     return enforceMin(x, min_val);
//   case 7: // SQRT: maps R -> (0,∞), g^(-1)(η) = η²
//     // Ensure x is non-negative before squaring
//     return enforceMin(CppAD::CondExpGt(x, Type(0.0), x * x, Type(0.0)), min_val);
//   case 8: // INVERSE: maps R -> (0,∞), g^(-1)(η) = 1/η
//     // Ensure denominator is positive and not too close to zero
//     return enforceMin(Type(1.0) / enforceMin(x, Type(1.0) / Type(1000000.0)), min_val); // 1e-6
//   case 9: // INVERSE-SQUARE: maps R -> (0,∞), g^(-1)(η) = 1/√η
//     // Ensure denominator is positive before taking square root
//     return enforceMin(Type(1.0) / sqrt(enforceMin(x, Type(1.0) / Type(1000000.0))), min_val); // 1e-6
//   default: // Default to LOG
//     return safeExp(x);
//   }
// }
//
// // ========================
// // Beta function utilities
// // ========================
//
// // Compute log of the Beta function: log[B(a,b)]
// template <class Type>
// Type logBeta(const Type &a, const Type &b) {
//   // Ensure parameters are positive
//   Type safe_a = enforceMin(a, Constants<Type>::eps_pos);
//   Type safe_b = enforceMin(b, Constants<Type>::eps_pos);
//
//   // Use lgamma for numerical stability
//   return lgamma(safe_a) + lgamma(safe_b) - lgamma(safe_a + safe_b);
// }
//
// // ========================
// // GKw Mean calculation
// // ========================
//
// // Calculate the mean of GKw distribution
// template <class Type>
// Type gkw_mean(const Type &alpha, const Type &beta, const Type &gamma,
//               const Type &delta, const Type &lambda) {
//   // Ensure all parameters are positive
//   Type safe_alpha = enforceMin(alpha, Constants<Type>::eps_pos);
//   Type safe_beta = enforceMin(beta, Constants<Type>::eps_pos);
//   Type safe_gamma = enforceMin(gamma, Constants<Type>::eps_pos);
//   Type safe_delta = enforceMin(delta, Constants<Type>::eps_pos);
//   Type safe_lambda = enforceMin(lambda, Constants<Type>::eps_pos);
//
//   // For numerical stability, we approximate the mean using a direct method
//   // The mean is calculated using numerical integration with 100 points
//   // Equally spaced in the (0,1) interval
//   Type mean = 0.0;
//   int n_points = 100;
//
//   for (int i = 1; i <= n_points; i++) {
//     // Calculate y value for this point (avoid exact 0 or 1)
//     Type yi = Type(i) / Type(n_points + 1);
//
//     // Calculate pdf at this point
//     Type ya = pow(yi, safe_alpha);
//     Type one_minus_ya = Type(1.0) - ya;
//     Type one_minus_ya_b = pow(one_minus_ya, safe_beta);
//     Type v = Type(1.0) - one_minus_ya_b;
//     Type v_lambda = pow(v, safe_lambda);
//     Type u = Type(1.0) - v_lambda;
//
//     // Calculate log-beta function: log[B(γ, δ+1)]
//     Type log_beta_func = logBeta(safe_gamma, safe_delta + Type(1.0));
//
//     // Compute pdf (density)
//     Type const_part = safe_lambda * safe_alpha * safe_beta / exp(log_beta_func);
//     Type pdf = const_part * pow(yi, safe_alpha - Type(1.0)) *
//       pow(one_minus_ya, safe_beta - Type(1.0)) *
//       pow(v, safe_gamma * safe_lambda - Type(1.0)) *
//       pow(u, safe_delta);
//
//     // Add contribution to mean
//     mean += yi * pdf;
//   }
//
//   // Normalize by total probability
//   mean = mean / Type(n_points);
//
//   // Ensure mean is in (0,1)
//   mean = enforceProbability(mean);
//
//   return mean;
// }
//
// // ========================
// // Main TMB template
// // ========================
// template<class Type>
// Type objective_function<Type>::operator() () {
//
//   // ========================
//   // Data section
//   // ========================
//   DATA_VECTOR(y);        // response vector (0,1)
//
//   // Covariate matrices for each parameter
//   DATA_MATRIX(X1);       // design matrix for alpha
//   DATA_MATRIX(X2);       // design matrix for beta
//   DATA_MATRIX(X3);       // design matrix for gamma
//   DATA_MATRIX(X4);       // design matrix for delta
//   DATA_MATRIX(X5);       // design matrix for lambda
//
//   // Link function types for each parameter
//   // 1=log (default), 2=logit, 3=probit, 4=cauchy, 5=cloglog, 6=identity
//   // 7=sqrt, 8=inverse, 9=inverse-square
//   DATA_INTEGER(link_type1);  // link type for alpha
//   DATA_INTEGER(link_type2);  // link type for beta
//   DATA_INTEGER(link_type3);  // link type for gamma
//   DATA_INTEGER(link_type4);  // link type for delta
//   DATA_INTEGER(link_type5);  // link type for lambda
//
//   // Scale factors for links that map to (0,1)
//   DATA_SCALAR(scale1);   // default=10
//   DATA_SCALAR(scale2);   // default=10
//   DATA_SCALAR(scale3);   // default=10
//   DATA_SCALAR(scale4);   // default=10
//   DATA_SCALAR(scale5);   // default=10
//
//   // ========================
//   // Parameter section
//   // ========================
//   PARAMETER_VECTOR(beta1);  // coefficients for alpha
//   PARAMETER_VECTOR(beta2);  // coefficients for beta
//   PARAMETER_VECTOR(beta3);  // coefficients for gamma
//   PARAMETER_VECTOR(beta4);  // coefficients for delta
//   PARAMETER_VECTOR(beta5);  // coefficients for lambda
//
//   // ========================
//   // Procedure section
//   // ========================
//   int n = y.size();  // sample size
//
//   // Input validation
//   if (n == 0) {
//     REPORT(n);
//     return Type(Constants<Type>::inf_repl);  // No data provided
//   }
//
//   if (X1.rows() != n || X2.rows() != n || X3.rows() != n || X4.rows() != n || X5.rows() != n) {
//     REPORT( static_cast<int>(X1.rows()) );
//     REPORT( static_cast<int>(X2.rows()) );
//     REPORT( static_cast<int>(X3.rows()) );
//     REPORT( static_cast<int>(X4.rows()) );
//     REPORT( static_cast<int>(X5.rows()) );
//     return Type(Constants<Type>::inf_repl);  // Design matrices don't match n
//   }
//
//   // Fixed: Cast to int to avoid overloaded function call ambiguity
//   int beta1_size = int(beta1.size());
//   int beta2_size = int(beta2.size());
//   int beta3_size = int(beta3.size());
//   int beta4_size = int(beta4.size());
//   int beta5_size = int(beta5.size());
//
//   int X1_cols = int(X1.cols());
//   int X2_cols = int(X2.cols());
//   int X3_cols = int(X3.cols());
//   int X4_cols = int(X4.cols());
//   int X5_cols = int(X5.cols());
//
//   if (X1_cols != beta1_size || X2_cols != beta2_size ||
//       X3_cols != beta3_size || X4_cols != beta4_size ||
//       X5_cols != beta5_size) {
//     REPORT(beta1_size);
//     REPORT(beta2_size);
//     REPORT(beta3_size);
//     REPORT(beta4_size);
//     REPORT(beta5_size);
//     return Type(Constants<Type>::inf_repl);  // Parameter vectors don't match design matrices
//   }
//
//   // Linear predictors
//   vector<Type> eta1 = X1 * beta1;  // for alpha
//   vector<Type> eta2 = X2 * beta2;  // for beta
//   vector<Type> eta3 = X3 * beta3;  // for gamma
//   vector<Type> eta4 = X4 * beta4;  // for delta
//   vector<Type> eta5 = X5 * beta5;  // for lambda
//
//   // Transformed parameters (all must be positive)
//   vector<Type> alpha(n), beta(n), gamma(n), delta(n), lambda(n);
//
//   // Apply appropriate link functions with scaling
//   for(int i=0; i<n; i++) {
//     alpha(i) = apply_positive_link(eta1(i), link_type1, scale1);
//     beta(i) = apply_positive_link(eta2(i), link_type2, scale2);
//     gamma(i) = apply_positive_link(eta3(i), link_type3, scale3);
//     delta(i) = apply_positive_link(eta4(i), link_type4, scale4);
//     lambda(i) = apply_positive_link(eta5(i), link_type5, scale5);
//   }
//
//   // Initialize negative log-likelihood
//   Type nll = 0.0;
//
//   // Flag indicating if any invalid values were encountered
//   int has_invalid_values = 0;
//
//   // Vector for fitted values
//   vector<Type> mu(n);
//   vector<Type> pearson_residuals(n);
//
//   // Evaluate negative log-likelihood for each observation
//   for(int i=0; i<n; i++) {
//     Type yi = y(i);
//
//     // Check if response is in valid range (0,1)
//     if((yi <= Type(0.0)) || (yi >= Type(1.0))) {
//       has_invalid_values = 1;
//       nll += Constants<Type>::inf_repl;  // Add penalty for out-of-range values
//       mu(i) = Type(0.5);  // Set a default value for mu if response is invalid
//       pearson_residuals(i) = Type(0.0);  // Set a default value for residuals
//       continue;
//     }
//
//     // Extract parameters for observation i
//     Type alpha_i = alpha(i);
//     Type beta_i = beta(i);
//     Type gamma_i = gamma(i);
//     Type delta_i = delta(i);
//     Type lambda_i = lambda(i);
//
//     // Additional check for extreme parameter values that could cause numerical problems
//     if (alpha_i < Constants<Type>::eps_pos || beta_i < Constants<Type>::eps_pos ||
//         gamma_i < Constants<Type>::eps_pos || delta_i < Constants<Type>::eps_pos ||
//         lambda_i < Constants<Type>::eps_pos) {
//       has_invalid_values = 1;
//       nll += Constants<Type>::inf_repl;  // Add penalty for invalid parameter values
//       mu(i) = Type(0.5);  // Set a default value for mu if parameters are invalid
//       pearson_residuals(i) = Type(0.0);  // Set a default value for residuals
//       continue;
//     }
//
//     /*
//      * Calculate log-pdf for the GKw distribution:
//      * log f(y; α, β, γ, δ, λ) =
//      *    log(λ) + log(α) + log(β) - log[B(γ, δ+1)]
//      *    + (α-1)log(y) + (β-1)log(1-y^α)
//      *    + (γλ-1)log[1-(1-y^α)^β]
//      *    + δlog{1-[1-(1-y^α)^β]^λ}
//      */
//
//     // Calculate auxiliary terms with numerical safeguards
//     Type ya = pow(yi, alpha_i);              // y^α
//     Type one_minus_ya = Type(1.0) - ya;        // 1-y^α
//
//     // Ensure one_minus_ya is within (0,1) to avoid domain errors in log and pow
//     one_minus_ya = enforceProbability(one_minus_ya);
//
//     Type one_minus_ya_b = pow(one_minus_ya, beta_i);  // (1-y^α)^β
//
//     // Ensure one_minus_ya_b is within (0,1)
//     one_minus_ya_b = enforceProbability(one_minus_ya_b);
//
//     Type v = Type(1.0) - one_minus_ya_b;     // 1-(1-y^α)^β
//
//     // Ensure v is within (0,1) to prevent domain errors
//     v = enforceProbability(v);
//
//     Type v_lambda = pow(v, lambda_i);        // [1-(1-y^α)^β]^λ
//
//     // Ensure v_lambda is within (0,1)
//     v_lambda = enforceProbability(v_lambda);
//
//     Type u = Type(1.0) - v_lambda;           // 1-[1-(1-y^α)^β]^λ
//
//     // Ensure u is within (0,1)
//     u = enforceProbability(u);
//
//     // Calculate log-beta function: log[B(γ, δ+1)]
//     Type log_beta_func = logBeta(gamma_i, delta_i + Type(1.0));
//
//     // Compute log-pdf with safe log operations
//     Type logf = safeLog(lambda_i) + safeLog(alpha_i) + safeLog(beta_i) - log_beta_func
//     + (alpha_i - Type(1.0)) * safeLog(yi)
//       + (beta_i - Type(1.0)) * safeLog(one_minus_ya)
//       + (gamma_i * lambda_i - Type(1.0)) * safeLog(v)
//       + delta_i * safeLog(u);
//
//       // Add to negative log-likelihood, checking for NaN or Inf
//       if (!R_FINITE(asDouble(logf))) {
//         has_invalid_values = 1;
//         nll += Constants<Type>::inf_repl;
//         mu(i) = Type(0.5);  // Set a default value for mu if likelihood is invalid
//         pearson_residuals(i) = Type(0.0);  // Set a default value for residuals
//       } else {
//         nll -= logf;
//
//         // Calculate fitted value (mean) for this observation
//         mu(i) = gkw_mean(alpha_i, beta_i, gamma_i, delta_i, lambda_i);
//
//         // Calculate Pearson residual: (y - μ) / √(Var(Y))
//         // For GKw, we'll use a simple approximation for the variance
//         // as calculating the exact variance is computationally intensive
//         Type approx_var = mu(i) * (Type(1.0) - mu(i)) / Type(5.0);
//         pearson_residuals(i) = (yi - mu(i)) / sqrt(approx_var);
//       }
//   }
//
//   // Report if any invalid values were encountered
//   REPORT(has_invalid_values);
//
//   // ========================
//   // Report section
//   // ========================
//
//   // Report transformed parameters
//   ADREPORT(alpha);
//   ADREPORT(beta);
//   ADREPORT(gamma);
//   ADREPORT(delta);
//   ADREPORT(lambda);
//
//   // Report fitted values and residuals
//   ADREPORT(mu);
//   ADREPORT(pearson_residuals);
//
//   // Calculate and report parameter means (useful for diagnostics)
//   Type mean_alpha = alpha.sum() / Type(n);
//   Type mean_beta = beta.sum() / Type(n);
//   Type mean_gamma = gamma.sum() / Type(n);
//   Type mean_delta = delta.sum() / Type(n);
//   Type mean_lambda = lambda.sum() / Type(n);
//
//   REPORT(mean_alpha);
//   REPORT(mean_beta);
//   REPORT(mean_gamma);
//   REPORT(mean_delta);
//   REPORT(mean_lambda);
//
//   // Calculate parameter standard deviations
//   Type sd_alpha = 0.0;
//   Type sd_beta = 0.0;
//   Type sd_gamma = 0.0;
//   Type sd_delta = 0.0;
//   Type sd_lambda = 0.0;
//
//   for(int i=0; i<n; i++) {
//     sd_alpha += pow(alpha(i) - mean_alpha, 2);
//     sd_beta += pow(beta(i) - mean_beta, 2);
//     sd_gamma += pow(gamma(i) - mean_gamma, 2);
//     sd_delta += pow(delta(i) - mean_delta, 2);
//     sd_lambda += pow(lambda(i) - mean_lambda, 2);
//   }
//
//   sd_alpha = sqrt(sd_alpha / Type(n));
//   sd_beta = sqrt(sd_beta / Type(n));
//   sd_gamma = sqrt(sd_gamma / Type(n));
//   sd_delta = sqrt(sd_delta / Type(n));
//   sd_lambda = sqrt(sd_lambda / Type(n));
//
//   REPORT(sd_alpha);
//   REPORT(sd_beta);
//   REPORT(sd_gamma);
//   REPORT(sd_delta);
//   REPORT(sd_lambda);
//
//   // Calculate and report model fit statistics
//   Type deviance = Type(2.0) * nll;
//   Type aic = deviance + Type(2.0) * Type(beta1_size + beta2_size + beta3_size + beta4_size + beta5_size);
//
//   REPORT(deviance);
//   REPORT(aic);
//
//   // Report goodness-of-fit statistics
//   Type mean_pearson_resid = pearson_residuals.sum() / Type(n);
//   Type sum_sq_pearson_resid = 0.0;
//
//   for(int i=0; i<n; i++) {
//     sum_sq_pearson_resid += pow(pearson_residuals(i), 2);
//   }
//
//   Type rmse = sqrt(sum_sq_pearson_resid / Type(n));
//
//   REPORT(mean_pearson_resid);
//   REPORT(rmse);
//
//   // Report parameter ranges for diagnostics
//   Type min_alpha = alpha(0);
//   Type max_alpha = alpha(0);
//   Type min_beta = beta(0);
//   Type max_beta = beta(0);
//   Type min_gamma = gamma(0);
//   Type max_gamma = gamma(0);
//   Type min_delta = delta(0);
//   Type max_delta = delta(0);
//   Type min_lambda = lambda(0);
//   Type max_lambda = lambda(0);
//
//   for(int i=1; i<n; i++) {
//     min_alpha = CppAD::CondExpLt(alpha(i), min_alpha, alpha(i), min_alpha);
//     max_alpha = CppAD::CondExpGt(alpha(i), max_alpha, alpha(i), max_alpha);
//     min_beta = CppAD::CondExpLt(beta(i), min_beta, beta(i), min_beta);
//     max_beta = CppAD::CondExpGt(beta(i), max_beta, beta(i), max_beta);
//     min_gamma = CppAD::CondExpLt(gamma(i), min_gamma, gamma(i), min_gamma);
//     max_gamma = CppAD::CondExpGt(gamma(i), max_gamma, gamma(i), max_gamma);
//     min_delta = CppAD::CondExpLt(delta(i), min_delta, delta(i), min_delta);
//     max_delta = CppAD::CondExpGt(delta(i), max_delta, delta(i), max_delta);
//     min_lambda = CppAD::CondExpLt(lambda(i), min_lambda, lambda(i), min_lambda);
//     max_lambda = CppAD::CondExpGt(lambda(i), max_lambda, lambda(i), max_lambda);
//   }
//
//   REPORT(min_alpha);
//   REPORT(max_alpha);
//   REPORT(min_beta);
//   REPORT(max_beta);
//   REPORT(min_gamma);
//   REPORT(max_gamma);
//   REPORT(min_delta);
//   REPORT(max_delta);
//   REPORT(min_lambda);
//   REPORT(max_lambda);
//
//   // Report the negative log-likelihood
//   REPORT(nll);
//
//   // Return negative log-likelihood
//   return nll;
// }
