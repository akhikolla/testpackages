#include <RcppArmadillo.h>
#include <test_stats.h>
#include <string>
#include <functional>

using namespace Rcpp;

std::function<double (const NumericVector&)>
  get_ts_fun(std::string stat, double k = 0) {

  if (stat == "RW") {
    return RW_Cpp;
  }
  if (stat == "BH2") {
    return BH2_Cpp;
  }
  if (stat == "K2") {
    return K2_Cpp;
  }
  if (stat == "K2U") {
    return K2U_Cpp;
  }
  if (stat == "KS") {
    return KS_Cpp;
  }
  if (stat == "SGN") {
    return SGN_Cpp;
  }
  if (stat == "WCX") {
    return WCX_Cpp;
  }
  if (stat == "CM") {
    return CM_Cpp;
  }
  if (stat == "MI") {
    return MI_Cpp;
  }
  if (stat == "MGG") {
    return MGG_Cpp;
  }
  if (stat == "B1") {
    return B1_Cpp;
  }
  if (stat == "FM") {
    return FM_Cpp;
  }
  if (stat == "BHI") {
    return BHI_Cpp;
  }
  if (stat == "BHK") {
    return BHK_Cpp;
  }
  if (stat == "MOI") {
    return std::bind(MOI_Cpp, std::placeholders::_1, k);
  }
  if (stat == "MOK") {
    return std::bind(MOK_Cpp, std::placeholders::_1, k);
  }
  if (stat == "NAI") {
    return std::bind(NAI_Cpp, std::placeholders::_1, k);
  }
  if (stat == "NAK") {
    return std::bind(NAK_Cpp, std::placeholders::_1, k);
  }
  if (stat == "NAC1") {
    return std::bind(NAC1_Cpp, std::placeholders::_1, k);
  }
  if (stat == "NAC2") {
    return std::bind(NAC2_Cpp, std::placeholders::_1, k);
  }
  if (stat == "BHC1") {
    return std::bind(BHC1_Cpp, std::placeholders::_1, k);
  }
  if (stat == "BHC2") {
    return std::bind(BHC2_Cpp, std::placeholders::_1, k);
  }
  if (stat == "HM") {
    return std::bind(HM_Cpp, std::placeholders::_1, k);
  }
  return NULL;
}

// [[Rcpp::export]]
NumericVector randomize_sign(const NumericVector& X, double mu) {
  int n = X.size();
  NumericVector res(X - mu);
  LogicalVector negative = runif(n, -1, 1) < 0;
  for (int i = 0; i < n; i++) {
    if(negative[i]) {
      res[i] = -res[i];
    }
  }
  return res;
}

NumericVector sample_with_replacement(NumericVector x, int n) {
  return x[floor(runif(n, 0, x.size()))];
}

NumericVector reflect_sample(const NumericVector& X, double mu, int n) {
  NumericVector reflected(2*n);
  for (int i = 0; i < n; i++) {
    reflected[i] = X[i];
  }
  for (int i = 0; i < n; i++) {
    reflected[i + n] = 2*mu - X[i];
  }
  return reflected;
}

// [[Rcpp::export]]
NumericVector reflected_boot(const NumericVector& X, double mu) {
  int n = X.size();
  return sample_with_replacement(reflect_sample(X, mu, n), n);
}

std::function<NumericVector (const NumericVector&, double)>
  get_null_fun(std::string null_method) {

    if (null_method == "sign") {
      return randomize_sign;
    }
    if (null_method == "reflect") {
      return reflected_boot;
    }
    return NULL;
  }


double trimmed_mean(const NumericVector& X, double alpha = 0) {
  if (alpha == 0) {
    return mean(X);
  }

  int n = X.size();

  NumericVector Xs(X);
  std::sort(Xs.begin(), Xs.end());

  int trim, lwr, upr;

  if (alpha < 0.5) {
    trim = floor(n * alpha);
    lwr = trim;
    upr = n - trim - 1;
  } else {
    trim = floor(n * 0.5);
    lwr = n - trim - 1;
    upr = trim;
  }

  return mean(Xs[Range(lwr, upr)]);
}

// [[Rcpp::export]]
NumericVector boot_sample(const NumericVector& X, double mu_param,
                          int B, std::string null_method,
                          std::string stat, double k = 0,
                          bool known_mean = false) {
  auto ts_fun = get_ts_fun(stat, k);
  auto null_sample_fun = get_null_fun(null_method);

  double mu = known_mean ? mu_param : trimmed_mean(X, mu_param);

  double mu_sym;
  NumericVector X_sym;

  NumericVector boot_sample(B);

  for (int i = 0; i < B; i++) {
    X_sym = null_sample_fun(X, mu);
    mu_sym = known_mean ? mu_param : trimmed_mean(X_sym, mu_param);
    boot_sample[i] = ts_fun(X_sym - mu_sym);
  }

  return boot_sample;
}

// [[Rcpp::export]]
NumericVector mn_boot_sample(const NumericVector& X, double mu_param,
                             int B, std::string stat, double k = 0,
                             double q = 8.0/9, bool known_mean = false) {
  auto ts_fun = get_ts_fun(stat, k);

  IntegerVector m_pre(21);
  int n = X.size();
  for (int i = 0; i < 21; i++) {
    m_pre[i] = round(n * pow(q, i));
  }
  IntegerVector m_filter = m_pre[m_pre > 4];
  IntegerVector m = unique(m_filter).sort(true);
  int m_size = m.size();

  NumericVector best_boot_sample(B);
  NumericVector last_boot_sample(B);
  NumericVector curr_boot_sample(B);
  NumericVector temp_diff(B);
  double best_dist = -1;
  double dist;

  NumericVector X_boot;
  double mu_boot;
  for (int i = 0; i < m_size; i++) {
    std::copy(curr_boot_sample.begin(), curr_boot_sample.end(),
              last_boot_sample.begin());

    for (int b = 0; b < B; b++) {
      X_boot = sample_with_replacement(X, m[i]);
      mu_boot = known_mean ? mu_param : trimmed_mean(X_boot, mu_param);
      curr_boot_sample[b] = ts_fun(X_boot - mu_boot);
    }

    if (i != 0) {
      temp_diff = last_boot_sample - curr_boot_sample;
      dist = sum(temp_diff * temp_diff);

      if (best_dist < 0 || dist < best_dist) {
        best_dist = dist;
        std::copy(curr_boot_sample.begin(), curr_boot_sample.end(),
                  best_boot_sample.begin());
      }
    }
  }

  return best_boot_sample;
}

// Get regression residuals given model matrix and response vector
// The code is a stripped version of Rcpp gallery example:
// http://gallery.rcpp.org/articles/fast-linear-model-with-armadillo/
NumericVector lm_resid(const arma::mat& X, NumericVector& yr) {
  arma::colvec y(yr.begin(), yr.size(), false);

  arma::colvec coef = arma::solve(X, y);    // fit model y ~ X
  arma::colvec res  = y - X*coef;           // residuals

  return NumericVector(res.begin(), res.end());
}

// [[Rcpp::export]]
NumericVector boot_sample_lm(const arma::mat& model_matrix,
                             const NumericVector& fitted,
                             const NumericVector& residuals,
                             int B, std::string null_method,
                             std::string stat, double k = 0) {
  auto ts_fun = get_ts_fun(stat, k);
  auto null_sample_fun = get_null_fun(null_method);

  NumericVector boot_resid;
  NumericVector new_resid;
  NumericVector boot_y;

  NumericVector boot_sample(B);

  for (int i = 0; i < B; i++) {
    boot_resid = null_sample_fun(residuals, 0);
    boot_y = fitted + boot_resid;
    new_resid = lm_resid(model_matrix, boot_y);
    boot_sample[i] = ts_fun(new_resid);
  }

  return boot_sample;
}

// [[Rcpp::export]]
NumericVector simulate_garch(const NumericVector& resid,
                             const NumericVector& y,
                             const NumericVector& cfit,
                             double omega,
                             const NumericVector& alpha,
                             const NumericVector& beta) {
  int q = alpha.size(), p = beta.size(), m = std::max(p, q);
  int n = resid.size();

  NumericVector booty(n+m);
  NumericVector bootc(n+m);
  NumericVector yrec(q);
  NumericVector crec(p);

  for (int i = 0; i < m; i++) {
    booty[i] = y[i];
    bootc[i] = cfit[i];
  }

  for (int i = 0; i < n; i++) {
    yrec = booty[Range(i + m - q, i + m - 1)];
    crec = bootc[Range(i + m - p, i + m - 1)];
    bootc[i + m] = sqrt(omega + sum(yrec * yrec * alpha) + sum(crec * crec * beta));
    booty[i + m] = bootc[i + m] * resid[i];
  }

  return booty[Range(m, n + m - 1)];
}
