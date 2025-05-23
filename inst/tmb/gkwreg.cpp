// File: gkwreg_v6.cpp
// ---------------------------------------------------------------
//  Regression Model for the Generalized Kumaraswamy (GKw) Distribution
//  using TMB, focusing on scalability with big data and many covariates.
//
//  Changes from v5:
//   - Removed all OpenMP dependencies for better portability
//   - Maintains all optimizations including caching with unordered_map
//   - Preserves computation scheduling flexibility via user chunk size
// ---------------------------------------------------------------

#include <TMB.hpp>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <cassert>

// Provide stub functions to replace OpenMP functions
// inline int omp_get_max_threads() { return 1; }
// inline int omp_get_thread_num() { return 0; }

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
// Process observation block helper function
// ========================
template <class Type>
void process_observation_block(
    Type& nll,
    vector<Type>& fitted,
    std::unordered_map<std::vector<int>, double, VectorIntHash, VectorIntEq>& meanCache,
    const vector<Type>& y,
    const vector<Type>& alphaVec,
    const vector<Type>& betaVec,
    const vector<Type>& gammaVec,
    const vector<Type>& deltaVec,
    const vector<Type>& lambdaVec,
    int start_idx,
    int end_idx,
    int calcFitted,
    int useMeanCache) {

  for (int i = start_idx; i < end_idx; i++) {
    Type alpha_i  = alphaVec(i);
    Type beta_i   = betaVec(i);
    Type gamma_i  = gammaVec(i);
    Type delta_i  = deltaVec(i);
    Type lambda_i = lambdaVec(i);

    // penalidade se parâmetros não são positivos
    if ( alpha_i < Constants<Type>::eps_pos || beta_i < Constants<Type>::eps_pos ||
         gamma_i < Constants<Type>::eps_pos || delta_i < Constants<Type>::eps_pos ||
         lambda_i < Constants<Type>::eps_pos ) {
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

        auto it = meanCache.find(key);
        if (it != meanCache.end()) {
          fitted(i) = Type(it->second);
        } else {
          Type mval = calc_mean_gkw(alpha_i, beta_i, gamma_i, delta_i, lambda_i);
          meanCache[key] = asDouble(mval);
          fitted(i) = mval;
        }
      } else {
        // sem cache, calcula diretamente
        fitted(i) = calc_mean_gkw(alpha_i, beta_i, gamma_i, delta_i, lambda_i);
      }
    }
  }
}

// ========================
// Sequential chunk processing to mimic OpenMP scheduling
// ========================
template <class Type>
void process_observations_in_chunks(
    Type& nll,
    vector<Type>& fitted,
    std::unordered_map<std::vector<int>, double, VectorIntHash, VectorIntEq>& meanCache,
    const vector<Type>& y,
    const vector<Type>& alphaVec,
    const vector<Type>& betaVec,
    const vector<Type>& gammaVec,
    const vector<Type>& deltaVec,
    const vector<Type>& lambdaVec,
    int n,
    int calcFitted,
    int useMeanCache,
    int userChunkSize) {

  // Default chunk size - this emulates the dynamic scheduling pattern
  int chunkSize = (userChunkSize > 0) ? userChunkSize : std::max(1, n / 20);

  for (int start = 0; start < n; start += chunkSize) {
    int end = std::min(start + chunkSize, n);
    process_observation_block(
      nll, fitted, meanCache, y,
      alphaVec, betaVec, gammaVec, deltaVec, lambdaVec,
      start, end, calcFitted, useMeanCache
    );
  }
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
  DATA_INTEGER(userChunkSize);  // Tamanho do escalonamento - mantido para compatibilidade

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
  typedef std::unordered_map<std::vector<int>, double, VectorIntHash, VectorIntEq> MeanMap;
  MeanMap meanCache;
  meanCache.reserve(std::min(n, 10000)); // heurística

  // == 5) NLL via processamento sequencial com chunks ==
  Type nll(0.0);

  // Se o tipo é double, podemos usar otimizações específicas
  // (caso contrário, TMB está no modo autodiff)
  if (isDouble<Type>::value) {
    process_observations_in_chunks(
      nll, fitted, meanCache, y,
      alphaVec, betaVec, gammaVec, deltaVec, lambdaVec,
      n, calcFitted, useMeanCache, userChunkSize
    );
  } else {
    // Modo autodiff: loop sequencial simples
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
