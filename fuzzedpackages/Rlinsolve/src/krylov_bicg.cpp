#define ARMA_DONT_PRINT_ERRORS
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//------------------------------------------------------------------------------
// Krylov Iterative Solver 2-1. Biconjugate Gradient
//------------------------------------------------------------------------------
//' @keywords internal
//' @noRd
// [[Rcpp::export(linsolve.bicg.single)]]
Rcpp::List single_bicg(const arma::mat& A, const arma::colvec& b, arma::colvec& xinit,
                       const double reltol, const int maxiter, const arma::mat& M){
  // 1-1. parameter settings
  int n = A.n_rows;
  int iter=0;
  int flag=0; // 0 found; 1 no cvgt; -1 breakdown

  // 1-2. Preiteration
  double bnrm2 = norm(b);
  if (bnrm2==0){
    bnrm2=1.0;
  }

  arma::colvec x = xinit;
  if (norm(b-A*xinit)<reltol){
    arma::colvec xtmp(n,fill::randn);
    x = xtmp;
  }
  arma::colvec r = b-A*x;
  double error = norm(r)/bnrm2;
  arma::colvec r_tld = r;

  // 1-3. main iteration
  double rho = 0;
  double rho_1 = 0;
  double alpha = 0;
  double beta = 0;

  arma::colvec p(n,fill::zeros);
  arma::colvec p_tld(n,fill::zeros);
  arma::colvec z(n,fill::zeros);
  arma::colvec z_tld(n,fill::zeros);
  arma::colvec q(n,fill::zeros);
  arma::colvec q_tld(n,fill::zeros);

  arma::vec errors(maxiter,fill::zeros);
  errors(0) = error;
  while (iter<maxiter){
    z = solve(M,r);
    z_tld = solve(M.t(),r_tld);
    rho   = dot(z,r_tld);
    if (rho==0.0){
      break;
    }

    if (iter>0){ // direction vectors
      beta  = rho/rho_1;
      p     = z+beta*p;
      p_tld = z_tld+beta*p_tld;
    } else {
      p     = z;
      p_tld = z_tld;
    }

    q = A*p; // compute residual pair
    q_tld = A.t()*p_tld;
    alpha = rho/dot(p_tld,q);

    x = x+alpha*p;
    r = r-alpha*q;
    r_tld = r_tld - alpha*q_tld;

    error = norm(r)/bnrm2;
    if (iter<(maxiter-1)){
      errors(iter+1) = error;
    }
    if (error<reltol){
      break;
    }
    rho_1 = rho;
    iter += 1;
  }

  // 1-4. post processing
  if (error<=reltol){
    flag = 0;
  } else if (rho==0.0){
    flag = -1;
  } else {
    flag = 1;
  }

  // 1-4. return outputs
  List res;
  res["x"] = x;
  res["iter"] = iter;
  if (iter>=maxiter){
    res["errors"] = errors;
  } else {
    arma::vec newerrors = errors.subvec(0,iter);
    res["errors"] = newerrors;
  }
  res["flag"] = flag;
  return res;
}


//------------------------------------------------------------------------------
// Krylov Iterative Solver 2-2. Biconjugate Gradient : Sparse Way
//------------------------------------------------------------------------------
//' @keywords internal
//' @noRd
// [[Rcpp::export(linsolve.bicg.single.sparse)]]
Rcpp::List single_bicg_sparse(const arma::sp_mat A, const arma::sp_mat b, arma::colvec& xinit,
                       const double reltol, const int maxiter, const arma::sp_mat M){
  // 1-1. parameter settings
  int n = A.n_rows;
  int iter=0;
  int flag=0; // 0 found; 1 no cvgt; -1 breakdown

  // 1-2. Preiteration
  double bnrm2 = norm(b);
  if (bnrm2==0){
    bnrm2=1.0;
  }
  arma::colvec x = xinit;
  if (norm(b-A*xinit)<reltol){
    arma::colvec xtmp(n,fill::randn);
    x = xtmp;
  }
  arma::colvec r = b-A*x;
  double error = norm(r)/bnrm2;
  arma::colvec r_tld = r;

  // 1-3. main iteration
  double rho = 0;
  double rho_1 = 0;
  double alpha = 0;
  double beta = 0;

  arma::colvec p(n,fill::zeros);
  arma::colvec p_tld(n,fill::zeros);
  arma::colvec z(n,fill::zeros);
  arma::colvec z_tld(n,fill::zeros);
  arma::colvec q(n,fill::zeros);
  arma::colvec q_tld(n,fill::zeros);

  arma::vec errors(maxiter,fill::zeros);
  errors(0) = error;
  while (iter<maxiter){
    z     = spsolve(M,r,"lapack");
    z_tld = spsolve(M.t(),r_tld,"lapack");
    rho   = dot(z,r_tld);
    if (rho==0.0){
      break;
    }

    if (iter>0){ // direction vectors
      beta  = rho/rho_1;
      p     = z+beta*p;
      p_tld = z_tld+beta*p_tld;
    } else {
      p     = z;
      p_tld = z_tld;
    }

    q = A*p; // compute residual pair
    q_tld = A.t()*p_tld;
    alpha = rho/dot(p_tld,q);

    x = x+alpha*p;
    r = r-alpha*q;
    r_tld = r_tld - alpha*q_tld;

    error = norm(r)/bnrm2;
    if (iter<(maxiter-1)){
      errors(iter+1) = error;
    }
    if (error<reltol){
      break;
    }
    rho_1 = rho;
    iter += 1;
  }

  // 1-4. post processing
  if (error<=reltol){
    flag = 0;
  } else if (rho==0.0){
    flag = -1;
  } else {
    flag = 1;
  }

  // 1-4. return outputs
  List res;
  res["x"] = x;
  res["iter"] = iter;
  if (iter>=maxiter){
    res["errors"] = errors;
  } else {
    arma::vec newerrors = errors.subvec(0,iter);
    res["errors"] = newerrors;
  }
  res["flag"] = flag;
  return res;
}
