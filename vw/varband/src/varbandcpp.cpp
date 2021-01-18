#include <cmath>
// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]
#define TOL 1.0e-6
#define SMALL 1.0e-15
#define MAXIT 1.0e+3

////' Soft-Threshold operator on a r-dim vector with the first element unpenalized
////'
////' @param beta An r-dim vector.
////' @param lam tuning paramter, i.e., the amount of penalization.
////' @param r row index, i.e., the length of beta.
// [[Rcpp::export]]
void soft_threshold(const arma::vec& beta, const double lam, const int r, arma::vec& result){
  // the last element is the same as beta(unpenalized)
  result(r - 1) = beta(r - 1);
  // This performs the soft-threshold for beta
  // for the first r - 1 elements
  for (int l = 0; l < (r - 1); l++){
    // element-wise soft-threshold
    if (beta(l) > lam)
      result(l) = beta(l) - lam;
    else if (beta(l) < -lam)
      result(l) = beta(l) + lam;
    else
      result(l) = 0;
  }
}

////' Close-form update of beta in Algorithm 1
////'
////' This function solves (6) in the paper with a closed form solution
////'
////' @param S An r-by-r submatrix of sample covariance matrix.
////' @param S_inv inverse of (2S_{-r,-r} + rho I)
////' @param r row index
////' @param rho parameter rho used in ADMM
////' @param u dual variable in ADMM
////' @param gamma priaml variable in ADMM
// [[Rcpp::export]]
void close_update(const arma::mat& S, const arma::mat& S_inv, const int r,
                       const double rho, const arma::vec& u, const arma::vec& gamma, arma::vec& res){
  // This performs the closed updated for beta
  // see Section 3 in the paper
  arma::vec vec_tmp = S_inv * S.col(r-1).head(r-1);
  double A = 4 * dot(vec_tmp, S.col(r-1).head(r-1)) - 2 * S(r-1, r-1) - rho;
  double B = 2 * dot(vec_tmp, u.head(r-1) - rho * gamma.head(r-1)) - u(r-1) + rho * gamma(r-1);
  // beta_r
  res(r-1) = (-std::sqrt(B * B - 8 * A) - B) / (2 * A);
  // beta_[-r]
  res.head(r-1) = -2 * res(r-1) * vec_tmp - S_inv * (u.head(r-1) - rho * gamma.head(r-1));
}

// [[Rcpp::export]]
void inverse_update(const arma::mat& S, double rho, arma::mat& S_inv){
  S_inv = 2*S;
  S_inv.diag() += rho;
  S_inv = inv_sympd(S_inv);
}

////' Evaluate the proximal operator of the hierarchical group lasso with simple weights
////'
////' This function solves (7) in the paper for unweighted version(w = 1)
////' by solving its dual by performing Newton's method on at most
////' r-1 univariate functions.
////' See Algorithm 2 in the paper
////'
////' @param y An r-dimensional vector.
////' @param tau lambda/rho
// [[Rcpp::export]]
void elliproj_u(const arma::vec& y, const double tau, arma::vec& pp){
  // This function performs the ellipsoid projection
  // of the unweighted estimator, which is very easy
  // See algorithm 2 in the paper
  int r = y.n_elem;
  // pp is the z vector in the paper
  pp = y;
  double tmpnorm = 0;
  for(int l = 0; l < (r - 1); l++){
    tmpnorm = arma::norm(pp.head(l + 1), 2);
    if(tmpnorm <= tau)
      pp.head(l+1).zeros();
    else
      pp.head(l+1) = (1 - tau / tmpnorm) * pp.head(l+1);
  }
}

// [[Rcpp::export]]
double rootfind(const arma::vec& pp, const arma::vec& ww, double tau, int l){
  // perform rooting finding
  // using a combinationg of Newton's method and bisection
  // see Numerical Recipes (3rd Edition, 2007) pp.460-461
  double xup = std::sqrt(arma::accu(arma::square(ww % pp.head(l+1)))) / tau;
  double xlo = ((xup-ww(l) * ww(l)) > 0.0 ? (xup-ww(l) * ww(l)) : 0.0);
  double fh = 1 - tau / std::sqrt(arma::accu(arma::square(pp.head(l+1) / (ww + xup / ww))));
  double fl = 1 - tau / std::sqrt(arma::accu(arma::square(pp.head(l+1) / (ww + xlo / ww))));

  if(std::fabs(fh) <= SMALL)
    return xup;
  if(std::fabs(fl) <= SMALL)
    return xlo;
  double xh, xl;
  if (fl < 0.0) {
    xl = xlo;
    xh = xup;
  }
  else {
    xh = xlo;
    xl = xup;
  }

  double rts = 0.5 * (xlo + xup);
  double dxold = std::fabs(xup - xlo);
  double dx = dxold;
  double f = 1 - tau / std::sqrt(arma::accu(arma::square(pp.head(l+1) / (ww + rts / ww))));
  double df = -tau * arma::accu(arma::square(ww % pp.head(l+1)) / arma::pow(arma::square(ww)+rts, 3)) / std::pow(arma::accu(arma::square(pp.head(l+1) / (ww + rts / ww))), 1.5);

  for (int j = 0; j < MAXIT; j++) {
    // if Newton's out of range
    // or if not decreasing fast enough
    if ((((rts - xh) * df - f) * ((rts - xl) * df - f) > 0.0)
          || (std::fabs(2.0 * f) > std::fabs(dxold * df))) {
      // bisection
      dxold = dx;
      dx = 0.5 * (xh - xl);
      rts = xl + dx;
      if (xl == rts)
        return rts;
    }
    else {
      // Newton
      dxold = dx;
      dx = f / df;
      double temp = rts;
      rts -= dx;
      if (temp == rts)
        return rts;
    }

    // convergence criterion
    if (std::fabs(dx) < TOL)
      return rts;

    // if not converged, next iteration
    // new function and direvative function values evaluations
    f = 1 - tau / std::sqrt(arma::accu(arma::square(pp.head(l+1) / (ww + rts / ww))));
    df = -tau * arma::accu(arma::square(ww % pp.head(l+1))/pow(arma::square(ww) + rts, 3))/std::pow(arma::accu(arma::square(pp.head(l+1) / (ww + rts / ww))), 1.5);
    if (f < 0.0)
      xl = rts;
    else
      xh = rts;
  }
  Rcpp::Rcout << "root finding fails to converge" << std::endl;
  return rts;
}

////' Evaluate the proximal operator of the hierarchical group lasso with general weights
////'
////' This function solves (7) in the paper for general weight w
////' by solving its dual by performing Newton's method on at most
////' r-1 univariate functions.
////' See Algorithm 1 and Theorem 1 in the online supplemenatry.
////'
////' @param y An r-dimensional vector.
////' @param tau lambda/rho
// [[Rcpp::export]]
void elliproj_w(const arma::vec& y, const double tau, arma::vec& pp){
  // This function performs the ellipsoid projection
  // See supplementary material
  int r = y.n_elem;

  arma::vec nu(r-1);
  nu.zeros();
  // pp is the z vector in the paper
  pp = y;

  for(int l = 0; l < (r - 1); l++){
    // ww[m] = w_{lm}
    arma::vec ww = arma::linspace<arma::vec>(l+1, 1, l+1);
    ww = 1 / arma::square(ww);
    // check if it lies in the ellipsoid
    if (arma::accu(arma::square(pp.head(l + 1) / ww)) <= tau * tau){
      nu(l) = 0;
      pp.head(l + 1).zeros();
    }
    else{
      // project onto the elliposid
      nu(l) = rootfind(pp, ww, tau, l);
      pp.head(l + 1) = pp.head(l + 1) * nu(l) / ( arma::square(ww) + nu(l) );
    }
  }
}

////' Compute one row of varband estimate with hierarchical group lasso penalty for a fixed tuning parameter
////'
////' This function solve the following r-th row estimation problem \deqn{min_{beta_r>0} -2 log beta_r + 1/n ||X beta||^2 + lambda P(beta)}
////' using an ADMM algorithm with changing rho.
////'
////' See algorithm 1 in the paper.
////'
////' @param S An r-by-r submatrix of sample covariance matrix.
////' @param init_row The initial estimate of the row.
////' @param lambda Non-negative tuning parameter. Controls sparsity level.
////' @param w Logical. Should we use weighted version of the penalty or not? If \code{TRUE}, we use general weight. If \code{FALSE}, use unweighted penalty. Default is \code{FALSE}.
////' @param tol Tolerance for convergence.
////' @param itermax Maximum number of iterations of ADMM to perform.
// [[Rcpp::export]]
arma::vec rowadmm(const arma::mat& S, const arma::vec& init_row,
                  const double lambda,
                  const bool w = false, double tol = 1.0e-4,
                  const int itermax = 1e+6){
  // This function solve the following row estimation problem
  // \min_{\beta_r>0} -2 log \beta_r + 1/n ||X\beta||^2
  // + \lambda P(\beta)
  // using an ADMM algorithm with changing rho
  int r = S.n_cols;
  // could use a lower tolerance for unweighted version
  if (!w)
    tol = 1.0e-8;

  // Default parameter in ADMM
  double tolabs = tol;
  double tolrel = tol;
  // Changing rho
  double rho = 2.0;
  double mu = 10.0;
  double inc = 2.0;
  double dec = 2.0;

  double pres = 0.0;
  double dres = 0.0;
  double peps = 0.0;
  double deps = 0.0;

  // Initialize the result
  arma::vec beta(init_row);
  arma::vec gamma(init_row);
  arma::vec beta_new(r);
  arma::vec gamma_new(r);

  // dual variable
  arma::vec u(r);
  u.zeros();

  arma::mat S_inv(r-1, r-1);
  inverse_update(S.submat(0, 0, r-2, r-2), rho, S_inv);
  for(int i = 0; i < itermax; i++) {
    // Primal&Dual Updates
    close_update(S, S_inv, r, rho, u, gamma, beta_new);
    if (w)
      elliproj_w(beta_new + u/rho, lambda/rho, gamma_new);
    else
      elliproj_u(beta_new + u/rho, lambda/rho, gamma_new);

    u = u + rho*(beta_new - gamma_new);

    // check convergence See pp 22 Boyd(2011)
    // primal residual
    pres = norm(beta_new - gamma_new, 2);
    // dual residual
    dres = rho * norm(gamma_new - gamma, 2);
    // primal tolerance
    peps = tolabs * std::sqrt(r) + tolrel * std::max(arma::norm(beta_new, 2), arma::norm(gamma_new, 2));
    // dual tolerance
    deps = tolabs * std::sqrt(r) + tolrel * arma::norm(u, 2);

    if(pres <= peps && dres <= deps)
      return gamma_new;
    else{
      // if not, update estimates and rho
      beta = beta_new;
      gamma = gamma_new;

      // Update rho if needed and corresponding S_inv
      if(pres > mu * dres){
        rho *= inc;
        inverse_update(S.submat(0, 0, r-2, r-2), rho, S_inv);
      }
      else if(dres > mu * pres){
        rho /= dec;
        inverse_update(S.submat(0, 0, r-2, r-2), rho, S_inv);
      }
    }
  }
  Rcpp::Rcout << "ADMM fails to converge" << std::endl;
  return gamma_new;
}

////' Compute one row of varband estimate with l1 penalty for a fixed tuning parameter
////'
////'This function solve the following r-th row estimation problem \deqn{min_{beta_r>0} -2 log beta_r + 1/n ||X beta||^2 + lambda |beta|_1}
////' using an ADMM algorithm with changing rho.
////'
////' See algorithm 1 in the paper.
////'
////' @param S An r-by-r submatrix of sample covariance matrix.
////' @param init_row The initial estimate of the row.
////' @param lambda Non-negative tuning parameter. Controls sparsity level.
////' @param tol Tolerance for convergence.
////' @param itermax Maximum number of iterations of ADMM to perform.
// [[Rcpp::export]]
arma::vec rowadmm_lasso(const arma::mat& S, const arma::vec& init_row,
                  const double lambda,
                  double tol = 1.0e-4,
                  const int itermax = 1e+6){
  // This function solve the following row estimation problem
  // \min_{\beta_r>0} -2 log \beta_r + 1/n ||X\beta||^2
  // + \lambda |\beta|_1
  // using an ADMM algorithm with changing rho
  int r = S.n_cols;
  // could use a lower tolerance for unweighted version

  // Default parameter in ADMM
  double tolabs = tol;
  double tolrel = tol;
  // Changing rho
  double rho = 2.0;
  double mu = 10.0;
  double inc = 2.0;
  double dec = 2.0;

  double pres = 0.0;
  double dres = 0.0;
  double peps = 0.0;
  double deps = 0.0;

  // Initialize the result
  arma::vec beta(init_row);
  arma::vec gamma(init_row);
  arma::vec beta_new(r);
  arma::vec gamma_new(r);

  // dual variable
  arma::vec u(r);
  u.zeros();

  arma::mat S_inv(r-1, r-1);
  inverse_update(S.submat(0, 0, r-2, r-2), rho, S_inv);
  for(int i = 0; i < itermax; i++) {
    // Primal&Dual Updates
    close_update(S, S_inv, r, rho, u, gamma, beta_new);
    soft_threshold(beta_new + u/rho, lambda/rho, r, gamma_new);

    u = u + rho*(beta_new - gamma_new);

    // check convergence See pp 22 Boyd(2011)
    // primal residual
    pres = norm(beta_new - gamma_new, 2);
    // dual residual
    dres = rho * norm(gamma_new - gamma, 2);
    // primal tolerance
    peps = tolabs * std::sqrt(r) + tolrel * std::max(arma::norm(beta_new, 2), arma::norm(gamma_new, 2));
    // dual tolerance
    deps = tolabs * std::sqrt(r) + tolrel * arma::norm(u, 2);

    if(pres <= peps && dres <= deps)
      return gamma_new;
    else{
      // if not, update estimates and rho
      beta = beta_new;
      gamma = gamma_new;

      // Update rho if needed and corresponding S_inv
      if(pres > mu * dres){
        rho *= inc;
        inverse_update(S.submat(0, 0, r-2, r-2), rho, S_inv);
      }
      else if(dres > mu * pres){
        rho /= dec;
        inverse_update(S.submat(0, 0, r-2, r-2), rho, S_inv);
      }
    }
  }
  Rcpp::Rcout << "ADMM fails to converge" << std::endl;
  return gamma_new;
}

//' Compute the varband estimate for a fixed tuning parameter value with different penalty options.
//'
//' Solves the main optimization problem in Yu & Bien (2016):
//' \deqn{min_L -2 \sum_{r=1}^p L_{rr} + tr(SLL^T) + lam * \sum_{r=2}^p P_r(L_{r.})}{min_L -2 sum_{r=1}^p L_{rr} + tr(SLL^T) + lam * sum_{r=2}^p P_r(L_{r.})}
//' where \deqn{P_r(L_{r.}) = \sum_{\ell = 2}^{r-1} \left(\sum_{m=1}^\ell w_{\ell m}^2 L_{rm}^2\right)^{1/2}}{P_r(L_r.) = sum_{l=2}^{r-1} (sum_m=1^l w^2_lm L^2_rm)^{1/2}}
//' or \deqn{P_r(L_{r.}) = \sum_{\ell = 1}^{r-1} |L_{r\ell}|}
//'
//' The function decomposes into p independent row problems,
//' each of which is solved by an ADMM algorithm.
//' see paper for more explanation.
//' @param S The sample covariance matrix
//' @param lambda Non-negative tuning parameter. Controls sparsity level.
//' @param w Logical. Should we use weighted version of the penalty or not? If \code{TRUE}, we use general weight. If \code{FALSE}, use unweighted penalty. Default is \code{FALSE}.
//' @param lasso Logical. Should we use l1 penalty instead of hierarchical group lasso penalty? Note that by using l1 penalty, we lose the banded structure in the resulting estimate. Default is \code{FALSE}.
//' @param init Initial estimate of L. Default is a closed-form diagonal estimate of L.
//' @return Returns the variable banding estimate of L, where L^TL = Omega.
//'
//' @examples
//' set.seed(123)
//' n <- 50
//' true <- varband_gen(p = 50, block = 5)
//' x <- sample_gen(L = true, n = n)
//' S <- crossprod(scale(x, center = TRUE, scale = FALSE)) / n
//' init <- diag(1/sqrt(diag(S)))
//' # unweighted estimate
//' L_unweighted <- varband(S, lambda = 0.1, init, w = FALSE)
//' # weighted estimate
//' L_weighted <- varband(S, lambda = 0.1, init, w = TRUE)
//' # lasso estimate
//' L_lasso <- varband(S, lambda = 0.1, init, w = TRUE, lasso = TRUE)
//' @seealso \code{\link{varband_path}} \code{\link{varband_cv}}
//'
//' @export
// [[Rcpp::export]]
arma::mat varband(arma::mat S, double lambda, arma::mat init, bool w = false,
                  bool lasso = false){
  int p = S.n_rows;
  arma::mat L(p, p);
  L.zeros();
  L(0, 0) = 1/(std::sqrt(S(0, 0)));
  init = init.t();
  if (lasso){
    for (int r = 1; r < p; r++)
      L.col(r).head(r+1) = rowadmm_lasso(S.submat(0, 0, r, r), init.col(r).head(r+1), lambda);
  }
  else{
    for (int r = 1; r < p; r++)
      L.col(r).head(r+1) = rowadmm(S.submat(0, 0, r, r), init.col(r).head(r+1), lambda, w);
  }
  return L.t();
}
