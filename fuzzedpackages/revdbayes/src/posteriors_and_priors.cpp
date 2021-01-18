// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::interfaces(r, cpp)]]

// Miscellaneous functions

// [[Rcpp::export]]
bool any_nonpos(const Rcpp::NumericVector& x) {
  return Rcpp::is_true(Rcpp::any(x <= 0)) ;
}

// Generalized Pareto log-likelihood

// [[Rcpp::export]]
double cpp_gp_loglik(const Rcpp::NumericVector& x, const Rcpp::List& ss) {
  double xm = ss["xm"] ;
  if (x[0] <= 0 || x[1] <= -x[0] / xm)
    return R_NegInf ;
  double loglik ;
  Rcpp::NumericVector gpd_data = ss["data"] ;
  Rcpp::NumericVector sdat = gpd_data / x[0] ;
  Rcpp::NumericVector zz = 1 + x[1] * sdat ;
  int m = ss["m"] ;
  if (std::abs(x[1]) > 1e-6) {
    loglik = -m * log(x[0]) - (1 + 1 / x[1]) * sum(log(zz)) ;
  } else {
    double sum_gp = ss["sum_gp"] ;
    double t1, t2, sdatj ;
    double total = 0.0;
    for(int j = 0; j < m; ++j) {
      sdatj = sdat[j] ;
      for(int i = 1; i < 5; ++i) {
        t1 = pow(sdatj, i) ;
        t2 = (i * sdatj - i - 1) ;
        total += pow(-1.0, i) * t1 * t2 * pow(x[1], i) / i / (i + 1) ;
      }
    }
    loglik = -m * log(x[0]) - sum_gp / x[0] - total ;
  }
  return loglik ;
}

// Generalized Extreme Value log-likelihood

// [[Rcpp::export]]
double cpp_gev_loglik(const Rcpp::NumericVector& x, const Rcpp::List& ss) {
  if (x[1] <= 0)
    return R_NegInf ;
  Rcpp::NumericVector data = ss["data"] ;
  Rcpp::NumericVector sdat = (data - x[0]) / x[1] ;
  Rcpp::NumericVector zz = 1 + x[2] * sdat ;
  if (any_nonpos(zz))
    return R_NegInf ;
  int m = ss["m"] ;
  double val = -m * log(x[1]) ;
  if (std::abs(x[2]) > 1e-6) {
    val = val - (1 + 1 / x[2]) * sum(log(zz)) - sum(pow(zz, (-1 / x[2]))) ;
  } else {
    double sum_gev = ss["sum_gev"] ;
    double t1, t2, sdatj, temp ;
    double t0 = (sum_gev - m * x[0]) / x[1] ;
    double tot = 0.0;
    double tot2 = 0.0;
    for(int j = 0; j < m; ++j) {
      sdatj = sdat[j] ;
      temp = 0.0 ;
      for(int i = 1; i < 5; ++i) {
        t1 = pow(sdatj, i) ;
        t2 = (i * sdatj - i - 1) ;
        tot += pow(-1.0, i) * t1 * t2 * pow(x[2], i) / i / (i + 1) ;
        temp += pow(-1.0, i) * pow(sdatj, i + 1) * pow(x[2], i) / (i + 1) ;
      }
      tot2 += exp(-sdatj - temp) ;
    }
    val = val - t0 - tot - tot2 ;
  }
  return val ;
}

// Order statistics log-likelihood

// [[Rcpp::export]]
double cpp_os_loglik(const Rcpp::NumericVector& x, const Rcpp::List& ss) {
  if (x[1] <= 0)
    return R_NegInf ;
  Rcpp::NumericVector data = ss["data"] ;
  Rcpp::NumericVector sdat = (data - x[0]) / x[1] ;
  Rcpp::NumericVector zz = 1 + x[2] * sdat ;
  if (any_nonpos(zz))
    return R_NegInf ;
  Rcpp::NumericVector min_data = ss["min_data"] ;
  Rcpp::NumericVector smindat = (min_data - x[0]) / x[1] ;
  Rcpp::NumericVector zz_min = 1 + x[2] * smindat ;
  int nos = ss["nos"] ;
  double val = -nos * log(x[1]) ;
  if (std::abs(x[2]) > 1e-6) {
    val = val - (1 + 1 / x[2]) * sum(log(zz)) - sum(pow(zz_min, (-1 / x[2]))) ;
  } else {
    double sum_os = ss["sum_os"] ;
    double t1, t2, sdatj, smindatj, temp ;
    double t0 = (sum_os - nos * x[0]) / x[1] ;
    double tot = 0.0;
    double tot2 = 0.0;
    for(int j = 0; j < nos; ++j) {
      sdatj = sdat[j] ;
      for(int i = 1; i < 5; ++i) {
        t1 = pow(sdatj, i) ;
        t2 = (i * sdatj - i - 1) ;
        tot += pow(-1.0, i) * t1 * t2 * pow(x[2], i) / i / (i + 1) ;
      }
    }
    int nmax = ss["nmax"] ;
    for(int j = 0; j < nmax; ++j) {
      temp = 0.0 ;
      for(int i = 1; i < 5; ++i) {
        smindatj = smindat[j] ;
        temp += pow(-1.0, i) * pow(smindatj, i + 1) * pow(x[2], i) / (i + 1) ;
      }
      tot2 += exp(-smindatj - temp) ;
    }
      val = val - t0 - tot - tot2 ;
  }
  return val ;
}

// Point process log-likelihood

// [[Rcpp::export]]
double cpp_pp_loglik(const Rcpp::NumericVector& x, const Rcpp::List& ss) {
  if (x[1] <= 0)
    return R_NegInf ;
  double thresh = ss["thresh"] ;
  double udat = (thresh - x[0]) / x[1] ;
  double zz_u = 1 + x[2] * udat ;
  if (zz_u <= 0)
    return R_NegInf ;
  Rcpp::NumericVector data = ss["data"] ;
  Rcpp::NumericVector sdat = (data - x[0]) / x[1] ;
  Rcpp::NumericVector zz = 1 + x[2] * sdat ;
  if (any_nonpos(zz))
    return R_NegInf ;
  double n_exc = ss["n_exc"] ;
  double noy = ss["noy"] ;
  double val = -n_exc * log(x[1]) ;
  if (std::abs(x[2]) > 1e-6) {
    val = val - (1 + 1 / x[2]) * sum(log(zz)) - noy * pow(zz_u, -1 / x[2]) ;
  } else {
    double sum_pp = ss["sum_pp"] ;
    double t1, t2, sdatj ;
    double t0 = (sum_pp - n_exc * x[0]) / x[1] ;
    double tot = 0.0;
    double tot2 = 0.0;
    for(int j = 0; j < n_exc; ++j) {
      sdatj = sdat[j] ;
      for(int i = 1; i < 5; ++i) {
        t1 = pow(sdatj, i) ;
        t2 = (i * sdatj - i - 1) ;
        tot += pow(-1.0, i) * t1 * t2 * pow(x[2], i) / i / (i + 1) ;
      }
    }
    for(int i = 1; i < 5; ++i) {
      tot2 += pow(-1.0, i) * pow(udat, i + 1) * pow(x[2], i) / (i + 1) ;
    }
    val = val - t0 - tot - noy * exp(-udat - tot2) ;
  }
  return val ;
}

// Generalized Pareto log-priors

// [[Rcpp::export]]
double cpp_gp_norm(const Rcpp::NumericVector& x, const Rcpp::List& ppars) {
  double min_xi = ppars["min_xi"] ;
  double max_xi = ppars["max_xi"] ;
  if (x[0] <= 0 || x[1] < min_xi || x[1] > max_xi)
    return R_NegInf ;
  Rcpp::NumericVector mean = ppars["mean"] ;
  Rcpp::NumericVector icov = ppars["icov"] ;
  double c0 = log(x[0]) - mean[0] ;
  double c1 = x[1] - mean[1] ;
  double ld = icov[0] * pow(c0, 2.0) + 2 * icov[1] * c0 * c1 +
    icov[2] * pow(c1, 2.0) ;
  return (-ld / 2 - log(x[0])) ;
}

// [[Rcpp::export]]
double cpp_gp_mdi(const Rcpp::NumericVector& x, const Rcpp::List& ppars) {
  double min_xi = ppars["min_xi"] ;
  double max_xi = ppars["max_xi"] ;
  if (x[0] <= 0 || x[1] < min_xi || x[1] > max_xi)
    return R_NegInf ;
  double a = ppars["a"] ;
  return -log(x[0]) - a * x[1] ;
}

// [[Rcpp::export]]
double cpp_gp_flat(const Rcpp::NumericVector& x, const Rcpp::List& ppars) {
  double min_xi = ppars["min_xi"] ;
  double max_xi = ppars["max_xi"] ;
  if (x[0] <= 0 || x[1] < min_xi || x[1] > max_xi)
    return R_NegInf ;
  return -log(x[0]) ;
}

// [[Rcpp::export]]
double cpp_gp_flatflat(const Rcpp::NumericVector& x, const Rcpp::List& ppars) {
  double min_xi = ppars["min_xi"] ;
  double max_xi = ppars["max_xi"] ;
  if (x[0] <= 0 || x[1] < min_xi || x[1] > max_xi)
    return R_NegInf ;
  return 0.0 ;
}

// [[Rcpp::export]]
double cpp_gp_jeffreys(const Rcpp::NumericVector& x, const Rcpp::List& ppars) {
  double min_xi = ppars["min_xi"] ;
  double max_xi = ppars["max_xi"] ;
  if (x[0] <= 0 || x[1] < min_xi || x[1] > max_xi)
    return R_NegInf ;
  return -log(x[0]) - log(1.0 + x[1]) - log(1.0 + 2.0 * x[1]) / 2.0 ;
}

// [[Rcpp::export]]
double cpp_gp_beta(const Rcpp::NumericVector& x, const Rcpp::List& ppars) {
  double min_xi = ppars["min_xi"] ;
  double max_xi = ppars["max_xi"] ;
  if (x[0] <= 0 || x[1] < min_xi || x[1] > max_xi)
    return R_NegInf ;
  Rcpp::NumericVector pq = ppars["pq"] ;
  return -log(x[0]) + (pq[0] - 1.0) * log(x[1] - min_xi) +
         (pq[1] - 1.0) * log(max_xi - x[1]) ;
}

// Generalized Extreme Value log-priors

// [[Rcpp::export]]
double cpp_gev_norm(const Rcpp::NumericVector& x, const Rcpp::List& ppars) {
  double min_xi = ppars["min_xi"] ;
  double max_xi = ppars["max_xi"] ;
  if (x[1] <= 0 || x[2] < min_xi || x[2] > max_xi)
    return R_NegInf ;
  Rcpp::NumericVector mean = ppars["mean"] ;
  Rcpp::NumericVector icov = ppars["icov"] ;
  double c0 = x[0] - mean[0] ;
  double c1 = log(x[1]) - mean[1] ;
  double c2 = x[2] - mean[2] ;
  double ld = icov[0] * pow(c0, 2.0) + 2 * icov[1] * c0 * c1 +
    2 * icov[2] * c0 * c2 + icov[3] * pow(c1, 2.0) + 2 * icov[4] * c1 * c2 +
    icov[5] * pow(c2, 2.0) ;
  return (-ld / 2 - log(x[1])) ;
}

// [[Rcpp::export]]
double cpp_gev_loglognorm(const Rcpp::NumericVector& x,
                          const Rcpp::List& ppars) {
  double min_xi = ppars["min_xi"] ;
  double max_xi = ppars["max_xi"] ;
  if (x[0] <= 0 || x[1] <= 0 || x[2] < min_xi || x[2] > max_xi)
    return R_NegInf ;
  Rcpp::NumericVector mean = ppars["mean"] ;
  Rcpp::NumericVector icov = ppars["icov"] ;
  double c0 = log(x[0]) - mean[0] ;
  double c1 = log(x[1]) - mean[1] ;
  double c2 = x[2] - mean[2] ;
  double ld = icov[0] * pow(c0, 2.0) + 2 * icov[1] * c0 * c1 +
    2 * icov[2] * c0 * c2 + icov[3] * pow(c1, 2.0) + 2 * icov[4] * c1 * c2 +
    icov[5] * pow(c2, 2.0) ;
  return (-ld / 2 - log(x[1]) - log(x[0])) ;
}

// [[Rcpp::export]]
double cpp_gev_mdi(const Rcpp::NumericVector& x, const Rcpp::List& ppars) {
  double min_xi = ppars["min_xi"] ;
  double max_xi = ppars["max_xi"] ;
  if (x[1] <= 0 || x[2] < min_xi || x[2] > max_xi)
    return R_NegInf ;
  double a = ppars["a"] ;
  return -log(x[1]) - a * x[2] ;
}

// [[Rcpp::export]]
double cpp_gev_flat(const Rcpp::NumericVector& x, const Rcpp::List& ppars) {
  double min_xi = ppars["min_xi"] ;
  double max_xi = ppars["max_xi"] ;
  if (x[1] <= 0 || x[2] < min_xi || x[2] > max_xi)
    return R_NegInf ;
  return -log(x[1]) ;
}

// [[Rcpp::export]]
double cpp_gev_flatflat(const Rcpp::NumericVector& x,
                        const Rcpp::List& ppars) {
  double min_xi = ppars["min_xi"] ;
  double max_xi = ppars["max_xi"] ;
  if (x[1] <= 0 || x[2] < min_xi || x[2] > max_xi)
    return R_NegInf ;
  return 0.0 ;
}

// [[Rcpp::export]]
double cpp_gev_beta(const Rcpp::NumericVector& x, const Rcpp::List& ppars) {
  double min_xi = ppars["min_xi"] ;
  double max_xi = ppars["max_xi"] ;
  if (x[1] <= 0 || x[2] < min_xi || x[2] > max_xi)
    return R_NegInf ;
  Rcpp::NumericVector pq = ppars["pq"] ;
  return -log(x[1]) + (pq[0] - 1.0) * log(x[2] - min_xi) +
    (pq[1] - 1.0) * log(max_xi - x[2]) ;
}

// GEV functions

// [[Rcpp::export]]
Rcpp::NumericVector lgdgev_cpp(const Rcpp::NumericVector& x, const double& loc,
                             const double& scale, const double& shape) {
  if (scale <= 0.0) {
    stop("invalid scale: scale must be positive.") ;
  }
  Rcpp::NumericVector xs = (x - loc) / scale ;
  Rcpp::NumericVector d = 1 + shape * xs ;
  for(int i = 0; i < x.size(); ++i) {
    if (d[i] < 0.0) {
      d[i] = R_NegInf ;
    } else {
      if (std::abs(shape) > 1e-6) {
        d[i] = -(1 + 1 / shape) * log(d[i]) - pow(d[i], -1 / shape) ;
      } else {
        d[i] = -xs[i] + shape * xs[i] * (xs[i] - 2) / 2 -
          exp(-xs[i] + shape * pow(xs[i], 2.0) / 2) ;
      }
    }
  }
  return d - log(scale) ;
}

// [[Rcpp::export]]
Rcpp::NumericVector pgev_cpp(const Rcpp::NumericVector& q, const double& loc,
                             const double& scale, const double& shape) {
  if (scale <= 0.0) {
    stop("invalid scale: scale must be positive.") ;
  }
  Rcpp::NumericVector qs = (q - loc) / scale ;
  Rcpp::NumericVector p = 1 + shape * qs ;
  for(int i = 0; i < q.size(); ++i) {
    if ((std::abs(shape) > 1e-6) || (p[i] < 0.0)) {
      p[i] = exp(-pow(std::max(p[i], 0.0), -1.0 / shape)) ;
    } else {
      p[i] = exp(-exp(-qs[i] + shape * pow(qs[i], 2.0) / 2)) ;
    }
  }
  return p ;
}

// [[Rcpp::export]]
Rcpp::NumericVector qgev_cpp(const Rcpp::NumericVector& p, const double& loc,
                             const double& scale, const double& shape) {
  if (scale <= 0.0) {
    stop("invalid scale: scale must be positive.") ;
  }
  int nq = p.size() ;
  Rcpp::NumericVector q(nq) ;
  Rcpp::NumericVector xp = -log(p) ;
  for(int i = 0; i < nq; ++i) {
    if (std::abs(shape) > 1e-6) {
      q[i] = -(pow(xp[i], -shape) - 1) / shape;
    } else {
      q[i] = log(xp[i]) * (1.0 - shape / 2.0) ;
    }
  }
  return loc - scale * q  ;
}

// [[Rcpp::export]]
double cpp_gev_prob(const Rcpp::NumericVector& x, const Rcpp::List& ppars) {
  // Extract GEV parameter values.
  double mu = x[0] ;
  double sigma = x[1] ;
  double xi = x[2] ;
  // prior is zero if scale parameter is non-positive
  if (sigma <= 0.0)
    return R_NegInf ;
  Rcpp::NumericVector quant = ppars["quant"] ;
  Rcpp::NumericVector y = (quant - mu) / sigma ;
  Rcpp::NumericVector lin = 1 + xi * y ;
  // prior is zero if any component of lin is not positive.
  if (Rcpp::is_true(Rcpp::any(lin <= 0)))
    return R_NegInf ;
  // If abs(xi) < xi_tol then xi is treated as being close to zero.
  double xi_tol = 1.e-6 ;
  Rcpp::NumericVector h(3) ;
  if (std::abs(xi) > xi_tol) {
    h = y / xi - lin * log(lin) / pow(xi, 2.0) ;
  } else {
    for(int i = 0; i < 3; ++i) {
      double tot = 0.0 ;
      for(int j = 0; j < 5; ++j) {
        tot += pow(-1.0, j + 1.0) * pow(y[i], j + 2.0) * pow(xi, j) /
          ((j + 1) * (j + 2)) ;
      }
      h[i] = tot ;
    }
  }
  arma::mat mat = ones(3,3) ;
  mat.col(1) = as<arma::vec>(y) ;
  mat.col(2) = as<arma::vec>(h) ;
  double det_mat = std::abs(arma::det(mat)) ;
  double log_det_mat = log(det_mat) ;
  // log-density of GEV at q
  Rcpp::NumericVector log_g = lgdgev_cpp(quant, mu, sigma, xi) ;
  Rcpp::NumericVector pq = pgev_cpp(quant, mu, sigma, xi) ;
  // log-Jacobian
  double log_J = log(sigma) + sum(log_g) + log_det_mat ;
  // Combine the Jacobian with the Dirichlet prior
  Rcpp::NumericVector alpha = ppars["alpha"] ;
  // Calculate the log Dirichlet prior
  Rcpp::NumericVector diff_pq =
    Rcpp::NumericVector::create(pq[0], pq[1] - pq[0], pq[2] - pq[1],
                                1.0 - pq[2]) ;
  if (Rcpp::is_true(Rcpp::any(diff_pq <= 0)))
    return R_NegInf ;
  // log-prior, up to an additive constant
  double log_dir_prior = sum((alpha - 1.0) * log(diff_pq)) ;
  double val = log_J + log_dir_prior ;
  return val ;
}

// [[Rcpp::export]]
double cpp_gev_quant(const Rcpp::NumericVector& x, const Rcpp::List& ppars) {
  // Extract GEV parameter values.
  double mu = x[0] ;
  double sigma = x[1] ;
  double xi = x[2] ;
  // prior is zero if scale parameter is non-positive
  if (sigma <= 0.0)
    return R_NegInf ;
  Rcpp::NumericVector prob = ppars["prob"] ;
// Calculate quantiles. Note: prob contains exceedance probabilities
  Rcpp::NumericVector quant = qgev_cpp(1.0 - prob, mu, sigma, xi) ;
// If the combination of (mu, sigma, xi) is such that any of the quantiles
// (quant[1] is the smallest) are non-positive then return -Inf.
  if (quant[0] <= 0) {
    return R_NegInf ;
  }
  Rcpp::NumericVector y = (quant - mu) / sigma ;
  Rcpp::NumericVector lin = 1 + xi * y ;
  // prior is zero if any component of lin is not positive.
  if (Rcpp::is_true(Rcpp::any(lin <= 0)))
    return R_NegInf ;
  // If abs(xi) < xi_tol then xi is treated as being close to zero.
  double xi_tol = 1.e-6 ;
  Rcpp::NumericVector h(3) ;
  if (std::abs(xi) > xi_tol) {
    h = y / xi - lin * log(lin) / pow(xi, 2.0) ;
  } else {
    for(int i = 0; i < 3; ++i) {
      double tot = 0.0 ;
      for(int j = 0; j < 5; ++j) {
        tot += pow(-1.0, j + 1.0) * pow(y[i], j + 2.0) * pow(xi, j) /
          ((j + 1) * (j + 2)) ;
      }
      h[i] = tot ;
    }
  }
  arma::mat mat = ones(3,3) ;
  mat.col(1) = as<arma::vec>(y) ;
  mat.col(2) = as<arma::vec>(h) ;
  double det_mat = std::abs(arma::det(mat)) ;
  double log_det_mat = log(det_mat) ;
  // log-Jacobian
  double log_J = log(sigma) + log_det_mat ;
  // Combine the Jacobian with the gamma prior
  // Calculate the log gamma prior
  // differences between the quant values
    Rcpp::NumericVector diff_quant =
    Rcpp::NumericVector::create(quant[0], quant[1] - quant[0],
                                quant[2] - quant[1]) ;
  // log-prior, up to an additive constant
  Rcpp::NumericVector shape = ppars["shape"] ;
  Rcpp::NumericVector scale = ppars["scale"] ;
  double log_gamma_prior = sum((shape - 1) * log(diff_quant) -
                               diff_quant / scale) ;
  double val = log_J + log_gamma_prior ;
  return val ;
}

// Model-specific log-posteriors for a user-defined prior

// [[Rcpp::export]]
double gp_user_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  SEXP prior_ptr = pars["prior"] ;
  // Unwrap pointer to the log-prior function.
  typedef double (*priorPtr)(const Rcpp::NumericVector& x,
                  const Rcpp::List& ppars) ;
  Rcpp::XPtr<priorPtr> xpfun(prior_ptr) ;
  priorPtr priorfun = *xpfun ;
  return cpp_gp_loglik(x, pars) + priorfun(x, pars) ;
}

// [[Rcpp::export]]
double gev_user_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  SEXP prior_ptr = pars["prior"] ;
  // Unwrap pointer to the log-prior function.
  typedef double (*priorPtr)(const Rcpp::NumericVector& x,
                  const Rcpp::List& ppars) ;
  Rcpp::XPtr<priorPtr> xpfun(prior_ptr) ;
  priorPtr priorfun = *xpfun ;
  return cpp_gev_loglik(x, pars) + priorfun(x, pars) ;
}

// [[Rcpp::export]]
double os_user_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  SEXP prior_ptr = pars["prior"] ;
  // Unwrap pointer to the log-prior function.
  typedef double (*priorPtr)(const Rcpp::NumericVector& x,
                  const Rcpp::List& ppars) ;
  Rcpp::XPtr<priorPtr> xpfun(prior_ptr) ;
  priorPtr priorfun = *xpfun ;
  return cpp_os_loglik(x, pars) + priorfun(x, pars) ;
}

// [[Rcpp::export]]
double pp_user_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  SEXP prior_ptr = pars["prior"] ;
  // Unwrap pointer to the log-prior function.
  typedef double (*priorPtr)(const Rcpp::NumericVector& x,
                  const Rcpp::List& ppars) ;
  Rcpp::XPtr<priorPtr> xpfun(prior_ptr) ;
  priorPtr priorfun = *xpfun ;
  return cpp_pp_loglik(x, pars) + priorfun(x, pars) ;
}

// GP Posteriors for specific in-built priors

// [[Rcpp::export]]
double gp_mdi_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_gp_loglik(x, pars) + cpp_gp_mdi(x, pars) ;
}

// [[Rcpp::export]]
double gp_norm_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_gp_loglik(x, pars) + cpp_gp_norm(x, pars) ;
}

// [[Rcpp::export]]
double gp_flat_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_gp_loglik(x, pars) + cpp_gp_flat(x, pars) ;
}

// [[Rcpp::export]]
double gp_flatflat_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_gp_loglik(x, pars) + cpp_gp_flatflat(x, pars) ;
}

// [[Rcpp::export]]
double gp_jeffreys_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_gp_loglik(x, pars) + cpp_gp_jeffreys(x, pars) ;
}

// [[Rcpp::export]]
double gp_beta_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_gp_loglik(x, pars) + cpp_gp_beta(x, pars) ;
}

// [[Rcpp::export]]
double gev_mdi_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_gev_loglik(x, pars) + cpp_gev_mdi(x, pars) ;
}

// GEV Posteriors for specific in-built priors

// [[Rcpp::export]]
double gev_norm_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_gev_loglik(x, pars) + cpp_gev_norm(x, pars) ;
}

// [[Rcpp::export]]
double gev_loglognorm_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_gev_loglik(x, pars) + cpp_gev_loglognorm(x, pars) ;
}

// [[Rcpp::export]]
double gev_flat_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_gev_loglik(x, pars) + cpp_gev_flat(x, pars) ;
}

// [[Rcpp::export]]
double gev_flatflat_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_gev_loglik(x, pars) + cpp_gev_flatflat(x, pars) ;
}

// [[Rcpp::export]]
double gev_beta_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_gev_loglik(x, pars) + cpp_gev_beta(x, pars) ;
}

// [[Rcpp::export]]
double gev_prob_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_gev_loglik(x, pars) + cpp_gev_prob(x, pars) ;
}

// [[Rcpp::export]]
double gev_quant_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_gev_loglik(x, pars) + cpp_gev_quant(x, pars) ;
}

// PP Posteriors for specific in-built priors

// [[Rcpp::export]]
double pp_mdi_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_pp_loglik(x, pars) + cpp_gev_mdi(x, pars) ;
}

// [[Rcpp::export]]
double pp_norm_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_pp_loglik(x, pars) + cpp_gev_norm(x, pars) ;
}

// [[Rcpp::export]]
double pp_loglognorm_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_pp_loglik(x, pars) + cpp_gev_loglognorm(x, pars) ;
}

// [[Rcpp::export]]
double pp_flat_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_pp_loglik(x, pars) + cpp_gev_flat(x, pars) ;
}

// [[Rcpp::export]]
double pp_flatflat_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_pp_loglik(x, pars) + cpp_gev_flatflat(x, pars) ;
}

// [[Rcpp::export]]
double pp_beta_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_pp_loglik(x, pars) + cpp_gev_beta(x, pars) ;
}

// [[Rcpp::export]]
double pp_prob_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_pp_loglik(x, pars) + cpp_gev_prob(x, pars) ;
}

// [[Rcpp::export]]
double pp_quant_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_pp_loglik(x, pars) + cpp_gev_quant(x, pars) ;
}

// OS Posteriors for specific in-built priors

// [[Rcpp::export]]
double os_mdi_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_os_loglik(x, pars) + cpp_gev_mdi(x, pars) ;
}

// [[Rcpp::export]]
double os_norm_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_os_loglik(x, pars) + cpp_gev_norm(x, pars) ;
}

// [[Rcpp::export]]
double os_loglognorm_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_os_loglik(x, pars) + cpp_gev_loglognorm(x, pars) ;
}

// [[Rcpp::export]]
double os_flat_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_os_loglik(x, pars) + cpp_gev_flat(x, pars) ;
}

// [[Rcpp::export]]
double os_flatflat_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_os_loglik(x, pars) + cpp_gev_flatflat(x, pars) ;
}

// [[Rcpp::export]]
double os_beta_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_os_loglik(x, pars) + cpp_gev_beta(x, pars) ;
}

// [[Rcpp::export]]
double os_prob_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_os_loglik(x, pars) + cpp_gev_prob(x, pars) ;
}

// [[Rcpp::export]]
double os_quant_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  return cpp_os_loglik(x, pars) + cpp_gev_quant(x, pars) ;
}

// [[Rcpp::export]]
SEXP gp_logpost_xptr(std::string fstr) {
  typedef double (*logpostPtr)(const Rcpp::NumericVector& x,
                  const Rcpp::List& pars) ;
  if (fstr == "gp_mdi")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gp_mdi_logpost))) ;
  else if (fstr == "gp_norm")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gp_norm_logpost))) ;
  else if (fstr == "gp_flat")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gp_flat_logpost))) ;
  else if (fstr == "gp_flatflat")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gp_flatflat_logpost))) ;
  else if (fstr == "gp_jeffreys")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gp_jeffreys_logpost))) ;
  else if (fstr == "gp_beta")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gp_beta_logpost))) ;
  else if (fstr == "gp_user")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gp_user_logpost))) ;
  else
    return(Rcpp::XPtr<logpostPtr>(R_NilValue)) ;
}

// [[Rcpp::export]]
SEXP gev_logpost_xptr(std::string fstr) {
  typedef double (*logpostPtr)(const Rcpp::NumericVector& x,
                  const Rcpp::List& pars) ;
  if (fstr == "gev_mdi")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gev_mdi_logpost))) ;
  else if (fstr == "gev_norm")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gev_norm_logpost))) ;
  else if (fstr == "gev_loglognorm")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gev_loglognorm_logpost))) ;
  else if (fstr == "gev_flat")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gev_flat_logpost))) ;
  else if (fstr == "gev_flatflat")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gev_flatflat_logpost))) ;
  else if (fstr == "gev_beta")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gev_beta_logpost))) ;
  else if (fstr == "gev_prob")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gev_prob_logpost))) ;
  else if (fstr == "gev_quant")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gev_quant_logpost))) ;
  else if (fstr == "gev_user")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gev_user_logpost))) ;
  else
    return(Rcpp::XPtr<logpostPtr>(R_NilValue)) ;
}

// [[Rcpp::export]]
SEXP os_logpost_xptr(std::string fstr) {
  typedef double (*logpostPtr)(const Rcpp::NumericVector& x,
                  const Rcpp::List& pars) ;
  if (fstr == "gev_mdi")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&os_mdi_logpost))) ;
  else if (fstr == "gev_norm")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&os_norm_logpost))) ;
  else if (fstr == "gev_loglognorm")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&os_loglognorm_logpost))) ;
  else if (fstr == "gev_flat")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&os_flat_logpost))) ;
  else if (fstr == "gev_flatflat")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&os_flatflat_logpost))) ;
  else if (fstr == "gev_beta")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&os_beta_logpost))) ;
  else if (fstr == "gev_prob")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&os_prob_logpost))) ;
  else if (fstr == "gev_quant")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&os_quant_logpost))) ;
  else if (fstr == "gev_user")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&os_user_logpost))) ;
  else
    return(Rcpp::XPtr<logpostPtr>(R_NilValue)) ;
}

// [[Rcpp::export]]
SEXP pp_logpost_xptr(std::string fstr) {
  typedef double (*logpostPtr)(const Rcpp::NumericVector& x,
                  const Rcpp::List& pars) ;
  if (fstr == "gev_mdi")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&pp_mdi_logpost))) ;
  else if (fstr == "gev_norm")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&pp_norm_logpost))) ;
  else if (fstr == "gev_loglognorm")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&pp_loglognorm_logpost))) ;
  else if (fstr == "gev_flat")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&pp_flat_logpost))) ;
  else if (fstr == "gev_flatflat")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&pp_flatflat_logpost))) ;
  else if (fstr == "gev_beta")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&pp_beta_logpost))) ;
  else if (fstr == "gev_prob")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&pp_prob_logpost))) ;
  else if (fstr == "gev_quant")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&pp_quant_logpost))) ;
  else if (fstr == "gev_user")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&pp_user_logpost))) ;
  else
    return(Rcpp::XPtr<logpostPtr>(R_NilValue)) ;
}

// Generalized Pareto phi_to_theta

// [[Rcpp::export]]
Rcpp::NumericVector gp_phi_to_theta(const Rcpp::NumericVector& phi,
                                    const Rcpp::List& user_args) {
  double xm = user_args["xm"] ;
  Rcpp::NumericVector val(2);
  val[0] = phi[0] ;
  val[1] = phi[1] - phi[0] / xm ;
  return val ;
}

// Generalized Extreme Value and (Order Statistics model) phi_to_theta

// [[Rcpp::export]]
Rcpp::NumericVector gev_phi_to_theta(const Rcpp::NumericVector& phi,
                                     const Rcpp::List& user_args) {
  double x1 = user_args["x1"] ;
  double xm = user_args["xm"] ;
  double sr = sqrt(xm - x1) ;
  Rcpp::NumericVector val(3);
  val[0] = phi[0] ;
  val[2] = (phi[2] - phi[1]) / sr ;
  val[1] = ((xm - phi[0]) * phi[1] + (phi[0] - x1) * phi[2]) / sr ;
  return val ;
}

// Point process model phi_to_theta

// [[Rcpp::export]]
Rcpp::NumericVector pp_phi_to_theta(const Rcpp::NumericVector& phi,
                                    const Rcpp::List& user_args) {
  double thresh = user_args["thresh"] ;
  double xm = user_args["xm"] ;
  double sr = sqrt(xm - thresh) ;
  Rcpp::NumericVector val(3);
  val[0] = phi[0] ;
  val[2] = (phi[2] - phi[1]) / sr ;
  val[1] = ((xm - phi[0]) * phi[1] + (phi[0] - thresh) * phi[2]) / sr ;
  return val ;
}

// K-gaps model phi_to_theta

// [[Rcpp::export]]
Rcpp::NumericVector kgaps_phi_to_theta(const Rcpp::NumericVector& phi,
                                       const Rcpp::List& user_args) {
  Rcpp::NumericVector ephi = exp(phi) ;
  return ephi / (1 + ephi) ;
}

// Create external pointers to phi_to_theta functions

// [[Rcpp::export]]
SEXP phi_to_theta_xptr(std::string fstr) {
  typedef Rcpp::NumericVector (*p2tPtr)(const Rcpp::NumericVector& phi,
                               const Rcpp::List& user_args) ;
  if (fstr == "gp")
    return(Rcpp::XPtr<p2tPtr>(new p2tPtr(&gp_phi_to_theta))) ;
  else if (fstr == "gev")
    return(Rcpp::XPtr<p2tPtr>(new p2tPtr(&gev_phi_to_theta))) ;
  else if (fstr == "os")
    return(Rcpp::XPtr<p2tPtr>(new p2tPtr(&gev_phi_to_theta))) ;
  else if (fstr == "pp")
    return(Rcpp::XPtr<p2tPtr>(new p2tPtr(&pp_phi_to_theta))) ;
  else if (fstr == "kgaps")
    return(Rcpp::XPtr<p2tPtr>(new p2tPtr(&kgaps_phi_to_theta))) ;
  else
    return(Rcpp::XPtr<p2tPtr>(R_NilValue)) ;
}

// GP posteriors for in-built priors, after transformation from theta to phi.

// [[Rcpp::export]]
double gp_mdi_logpost_phi(const Rcpp::NumericVector& phi,
                          const Rcpp::List& pars){
  Rcpp::NumericVector theta = gp_phi_to_theta(phi, pars) ;
  return gp_mdi_logpost(theta, pars) ;
}

// [[Rcpp::export]]
double gp_norm_logpost_phi(const Rcpp::NumericVector& phi,
                          const Rcpp::List& pars){
  Rcpp::NumericVector theta = gp_phi_to_theta(phi, pars) ;
  return gp_norm_logpost(theta, pars) ;
}

// [[Rcpp::export]]
double gp_flat_logpost_phi(const Rcpp::NumericVector& phi,
                          const Rcpp::List& pars){
  Rcpp::NumericVector theta = gp_phi_to_theta(phi, pars) ;
  return gp_flat_logpost(theta, pars) ;
}

// [[Rcpp::export]]
double gp_flatflat_logpost_phi(const Rcpp::NumericVector& phi,
                           const Rcpp::List& pars){
  Rcpp::NumericVector theta = gp_phi_to_theta(phi, pars) ;
  return gp_flatflat_logpost(theta, pars) ;
}

// [[Rcpp::export]]
double gp_jeffreys_logpost_phi(const Rcpp::NumericVector& phi,
                           const Rcpp::List& pars){
  Rcpp::NumericVector theta = gp_phi_to_theta(phi, pars) ;
  return gp_jeffreys_logpost(theta, pars) ;
}

// [[Rcpp::export]]
double gp_beta_logpost_phi(const Rcpp::NumericVector& phi,
                           const Rcpp::List& pars){
  Rcpp::NumericVector theta = gp_phi_to_theta(phi, pars) ;
  return gp_beta_logpost(theta, pars) ;
}

// [[Rcpp::export]]
double gp_user_logpost_phi(const Rcpp::NumericVector& phi,
                           const Rcpp::List& pars){
  Rcpp::NumericVector theta = gp_phi_to_theta(phi, pars) ;
  return gp_user_logpost(theta, pars) ;
}

// GEV posteriors for in-built priors, after transformation from theta to phi.

// [[Rcpp::export]]
double gev_mdi_logpost_phi(const Rcpp::NumericVector& phi,
                           const Rcpp::List& pars){
  Rcpp::NumericVector theta = gev_phi_to_theta(phi, pars) ;
  return gev_mdi_logpost(theta, pars) ;
}

// [[Rcpp::export]]
double gev_norm_logpost_phi(const Rcpp::NumericVector& phi,
                            const Rcpp::List& pars){
  Rcpp::NumericVector theta = gev_phi_to_theta(phi, pars) ;
  return gev_norm_logpost(theta, pars) ;
}

// [[Rcpp::export]]
double gev_loglognorm_logpost_phi(const Rcpp::NumericVector& phi,
                                  const Rcpp::List& pars){
  Rcpp::NumericVector theta = gev_phi_to_theta(phi, pars) ;
  return gev_loglognorm_logpost(theta, pars) ;
}

// [[Rcpp::export]]
double gev_flat_logpost_phi(const Rcpp::NumericVector& phi,
                            const Rcpp::List& pars){
  Rcpp::NumericVector theta = gev_phi_to_theta(phi, pars) ;
  return gev_flat_logpost(theta, pars) ;
}

// [[Rcpp::export]]
double gev_flatflat_logpost_phi(const Rcpp::NumericVector& phi,
                                const Rcpp::List& pars){
  Rcpp::NumericVector theta = gev_phi_to_theta(phi, pars) ;
  return gev_flatflat_logpost(theta, pars) ;
}

// [[Rcpp::export]]
double gev_beta_logpost_phi(const Rcpp::NumericVector& phi,
                            const Rcpp::List& pars){
  Rcpp::NumericVector theta = gev_phi_to_theta(phi, pars) ;
  return gev_beta_logpost(theta, pars) ;
}

// [[Rcpp::export]]
double gev_prob_logpost_phi(const Rcpp::NumericVector& phi,
                            const Rcpp::List& pars){
  Rcpp::NumericVector theta = gev_phi_to_theta(phi, pars) ;
  return gev_prob_logpost(theta, pars) ;
}

// [[Rcpp::export]]
double gev_quant_logpost_phi(const Rcpp::NumericVector& phi,
                             const Rcpp::List& pars){
  Rcpp::NumericVector theta = gev_phi_to_theta(phi, pars) ;
  return gev_quant_logpost(theta, pars) ;
}

// [[Rcpp::export]]
double gev_user_logpost_phi(const Rcpp::NumericVector& phi,
                            const Rcpp::List& pars){
  Rcpp::NumericVector theta = gev_phi_to_theta(phi, pars) ;
  return gev_user_logpost(theta, pars) ;
}

// OS posteriors for in-built priors, after transformation from theta to phi.

// [[Rcpp::export]]
double os_mdi_logpost_phi(const Rcpp::NumericVector& phi,
                          const Rcpp::List& pars){
  Rcpp::NumericVector theta = gev_phi_to_theta(phi, pars) ;
  return os_mdi_logpost(theta, pars) ;
}

// [[Rcpp::export]]
double os_norm_logpost_phi(const Rcpp::NumericVector& phi,
                           const Rcpp::List& pars){
  Rcpp::NumericVector theta = gev_phi_to_theta(phi, pars) ;
  return os_norm_logpost(theta, pars) ;
}

// [[Rcpp::export]]
double os_loglognorm_logpost_phi(const Rcpp::NumericVector& phi,
                                 const Rcpp::List& pars){
  Rcpp::NumericVector theta = gev_phi_to_theta(phi, pars) ;
  return os_loglognorm_logpost(theta, pars) ;
}

// [[Rcpp::export]]
double os_flat_logpost_phi(const Rcpp::NumericVector& phi,
                           const Rcpp::List& pars){
  Rcpp::NumericVector theta = gev_phi_to_theta(phi, pars) ;
  return os_flat_logpost(theta, pars) ;
}

// [[Rcpp::export]]
double os_flatflat_logpost_phi(const Rcpp::NumericVector& phi,
                               const Rcpp::List& pars){
  Rcpp::NumericVector theta = gev_phi_to_theta(phi, pars) ;
  return os_flatflat_logpost(theta, pars) ;
}

// [[Rcpp::export]]
double os_beta_logpost_phi(const Rcpp::NumericVector& phi,
                           const Rcpp::List& pars){
  Rcpp::NumericVector theta = gev_phi_to_theta(phi, pars) ;
  return os_beta_logpost(theta, pars) ;
}

// [[Rcpp::export]]
double os_prob_logpost_phi(const Rcpp::NumericVector& phi,
                            const Rcpp::List& pars){
  Rcpp::NumericVector theta = gev_phi_to_theta(phi, pars) ;
  return os_prob_logpost(theta, pars) ;
}

// [[Rcpp::export]]
double os_quant_logpost_phi(const Rcpp::NumericVector& phi,
                             const Rcpp::List& pars){
  Rcpp::NumericVector theta = gev_phi_to_theta(phi, pars) ;
  return os_quant_logpost(theta, pars) ;
}

// [[Rcpp::export]]
double os_user_logpost_phi(const Rcpp::NumericVector& phi,
                           const Rcpp::List& pars){
  Rcpp::NumericVector theta = gev_phi_to_theta(phi, pars) ;
  return os_user_logpost(theta, pars) ;
}

// PP posteriors for in-built priors, after transformation from theta to phi.

// [[Rcpp::export]]
double pp_mdi_logpost_phi(const Rcpp::NumericVector& phi,
                          const Rcpp::List& pars){
  Rcpp::NumericVector theta = pp_phi_to_theta(phi, pars) ;
  return pp_mdi_logpost(theta, pars) ;
}

// [[Rcpp::export]]
double pp_norm_logpost_phi(const Rcpp::NumericVector& phi,
                           const Rcpp::List& pars){
  Rcpp::NumericVector theta = pp_phi_to_theta(phi, pars) ;
  return pp_norm_logpost(theta, pars) ;
}

// [[Rcpp::export]]
double pp_loglognorm_logpost_phi(const Rcpp::NumericVector& phi,
                                 const Rcpp::List& pars){
  Rcpp::NumericVector theta = pp_phi_to_theta(phi, pars) ;
  return pp_loglognorm_logpost(theta, pars) ;
}

// [[Rcpp::export]]
double pp_flat_logpost_phi(const Rcpp::NumericVector& phi,
                           const Rcpp::List& pars){
  Rcpp::NumericVector theta = pp_phi_to_theta(phi, pars) ;
  return pp_flat_logpost(theta, pars) ;
}

// [[Rcpp::export]]
double pp_flatflat_logpost_phi(const Rcpp::NumericVector& phi,
                               const Rcpp::List& pars){
  Rcpp::NumericVector theta = pp_phi_to_theta(phi, pars) ;
  return pp_flatflat_logpost(theta, pars) ;
}

// [[Rcpp::export]]
double pp_beta_logpost_phi(const Rcpp::NumericVector& phi,
                           const Rcpp::List& pars){
  Rcpp::NumericVector theta = pp_phi_to_theta(phi, pars) ;
  return pp_beta_logpost(theta, pars) ;
}

// [[Rcpp::export]]
double pp_prob_logpost_phi(const Rcpp::NumericVector& phi,
                           const Rcpp::List& pars){
  Rcpp::NumericVector theta = pp_phi_to_theta(phi, pars) ;
  return pp_prob_logpost(theta, pars) ;
}

// [[Rcpp::export]]
double pp_quant_logpost_phi(const Rcpp::NumericVector& phi,
                            const Rcpp::List& pars){
  Rcpp::NumericVector theta = pp_phi_to_theta(phi, pars) ;
  return pp_quant_logpost(theta, pars) ;
}

// [[Rcpp::export]]
double pp_user_logpost_phi(const Rcpp::NumericVector& phi,
                           const Rcpp::List& pars){
  Rcpp::NumericVector theta = pp_phi_to_theta(phi, pars) ;
  return pp_user_logpost(theta, pars) ;
}

// [[Rcpp::export]]
SEXP gp_logpost_phi_xptr(std::string fstr) {
  typedef double (*logpostPtr)(const Rcpp::NumericVector& x,
                  const Rcpp::List& pars) ;
  if (fstr == "gp_mdi")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gp_mdi_logpost_phi))) ;
  else if (fstr == "gp_norm")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gp_norm_logpost_phi))) ;
  else if (fstr == "gp_flat")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gp_flat_logpost_phi))) ;
  else if (fstr == "gp_flatflat")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gp_flatflat_logpost_phi))) ;
  else if (fstr == "gp_jeffreys")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gp_jeffreys_logpost_phi))) ;
  else if (fstr == "gp_beta")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gp_beta_logpost_phi))) ;
  else if (fstr == "user")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gp_user_logpost_phi))) ;
  else
    return(Rcpp::XPtr<logpostPtr>(R_NilValue)) ;
}

// [[Rcpp::export]]
SEXP gev_logpost_phi_xptr(std::string fstr) {
  typedef double (*logpostPtr)(const Rcpp::NumericVector& x,
                  const Rcpp::List& pars) ;
  if (fstr == "gev_mdi")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gev_mdi_logpost_phi))) ;
  else if (fstr == "gev_norm")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gev_norm_logpost_phi))) ;
  else if (fstr == "gev_loglognorm")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gev_loglognorm_logpost_phi))) ;
  else if (fstr == "gev_flat")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gev_flat_logpost_phi))) ;
  else if (fstr == "gev_flatflat")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gev_flatflat_logpost_phi))) ;
  else if (fstr == "gev_beta")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gev_beta_logpost_phi))) ;
  else if (fstr == "gev_prob")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gev_prob_logpost_phi))) ;
  else if (fstr == "gev_quant")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gev_quant_logpost_phi))) ;
  else if (fstr == "user")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&gev_user_logpost_phi))) ;
  else
    return(Rcpp::XPtr<logpostPtr>(R_NilValue)) ;
}

// [[Rcpp::export]]
SEXP pp_logpost_phi_xptr(std::string fstr) {
  typedef double (*logpostPtr)(const Rcpp::NumericVector& x,
                  const Rcpp::List& pars) ;
  if (fstr == "gev_mdi")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&pp_mdi_logpost_phi))) ;
  else if (fstr == "gev_norm")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&pp_norm_logpost_phi))) ;
  else if (fstr == "gev_loglognorm")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&pp_loglognorm_logpost_phi))) ;
  else if (fstr == "gev_flat")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&pp_flat_logpost_phi))) ;
  else if (fstr == "gev_flatflat")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&pp_flatflat_logpost_phi))) ;
  else if (fstr == "gev_beta")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&pp_beta_logpost_phi))) ;
  else if (fstr == "gev_prob")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&pp_prob_logpost_phi))) ;
  else if (fstr == "gev_quant")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&pp_quant_logpost_phi))) ;
  else if (fstr == "user")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&pp_user_logpost_phi))) ;
  else
    return(Rcpp::XPtr<logpostPtr>(R_NilValue)) ;
}

// [[Rcpp::export]]
SEXP os_logpost_phi_xptr(std::string fstr) {
  typedef double (*logpostPtr)(const Rcpp::NumericVector& x,
                  const Rcpp::List& pars) ;
  if (fstr == "gev_mdi")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&os_mdi_logpost_phi))) ;
  else if (fstr == "gev_norm")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&os_norm_logpost_phi))) ;
  else if (fstr == "gev_loglognorm")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&os_loglognorm_logpost_phi))) ;
  else if (fstr == "gev_flat")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&os_flat_logpost_phi))) ;
  else if (fstr == "gev_flatflat")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&os_flatflat_logpost_phi))) ;
  else if (fstr == "gev_beta")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&os_beta_logpost_phi))) ;
  else if (fstr == "gev_prob")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&os_prob_logpost_phi))) ;
  else if (fstr == "gev_quant")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&os_quant_logpost_phi))) ;
  else if (fstr == "user")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&os_user_logpost_phi))) ;
  else
    return(Rcpp::XPtr<logpostPtr>(R_NilValue)) ;
}

// K-gaps

// [[Rcpp::export]]
double kgaps_logpost(const Rcpp::NumericVector& x, const Rcpp::List& pars) {
  if (x[0] < 0 || x[0] > 1)
    return R_NegInf ;
  int N0 = pars["N0"] ;
  int N1 = pars["N1"] ;
  double sum_qs = pars["sum_qs"] ;
  double loglik = 0.0 ;
  if (N1 > 0)
    loglik = loglik + 2 * N1 * log(x[0]) - sum_qs * x[0] ;
  if (N0 > 0)
    loglik = loglik + N0 * log(1 - x[0]) ;
  double alpha = pars["alpha"] ;
  double beta = pars["beta"] ;
  double logprior = (alpha - 1) * log(x[0]) + (beta - 1) * log(1 - x[0]) ;
  return loglik + logprior ;
}

// [[Rcpp::export]]
SEXP kgaps_logpost_xptr(std::string fstr) {
  typedef double (*logpostPtr)(const Rcpp::NumericVector& x,
                  const Rcpp::List& pars) ;
  if (fstr == "kgaps")
    return(Rcpp::XPtr<logpostPtr>(new logpostPtr(&kgaps_logpost))) ;
  else
    return(Rcpp::XPtr<logpostPtr>(R_NilValue)) ;
}

// K-gaps model log_j

// [[Rcpp::export]]
double kgaps_log_j(const Rcpp::NumericVector& theta,
                   const Rcpp::List& user_args) {
  return -log(theta[0]) - log(1 - theta[0]) ;
}

// Create external pointers to log_j functions

// [[Rcpp::export]]
SEXP log_j_xptr(std::string fstr) {
  typedef double (*p2tPtr)(const Rcpp::NumericVector& theta,
                  const Rcpp::List& user_args) ;
  if (fstr == "kgaps")
    return(Rcpp::XPtr<p2tPtr>(new p2tPtr(&kgaps_log_j))) ;
  else
    return(Rcpp::XPtr<p2tPtr>(R_NilValue)) ;
}
