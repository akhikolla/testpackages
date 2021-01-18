#define ARMA_DONT_PRINT_ERRORS
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//------------------------------------------------------------------------------
// Krylov Iterative Solver 4-1. Conjugate Gradient Squared
//------------------------------------------------------------------------------
//' @keywords internal
//' @noRd
// [[Rcpp::export(linsolve.cgs.single)]]
Rcpp::List single_cgs(const arma::mat& A, const arma::colvec& b, arma::colvec& xinit,
                     const double reltol, const int maxiter, const arma::mat& M){
  // 1-1. parameter settings
  const int n = A.n_cols;
  int iter=0;
  int flag=0; // 0 found; 1 no convergence

  // 1-2. pre-iteration
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
  arma::colvec r_tld = r;

  // 1-3. main iteration
  double rho = 0;
  double rho_1 = 0;
  double alpha = 0;
  double beta = 0;
  
  arma::colvec u(n,fill::zeros);
  arma::colvec p(n,fill::zeros);
  arma::colvec q(n,fill::zeros);
  arma::colvec p_hat(n,fill::zeros);
  arma::colvec v_hat(n,fill::zeros);
  arma::colvec u_hat(n,fill::zeros);

  double error = norm(r)/bnrm2;
  arma::vec errors(maxiter,fill::zeros);
  errors(0) = error;
  while (iter<maxiter){
    rho = dot(r_tld,r);
    if (rho==0.0){
      break;
    }

    if (iter>0){ // direction vectors
      beta = rho/rho_1;
      u = r+beta*q;
      p = u+beta*(q+beta*p);
    } else {
      u = r;
      p = u;
    }

    p_hat = solve(M,p);
    v_hat = A*p_hat; // adjusting scalars
    alpha = rho/dot(r_tld,v_hat);
    q     = u-alpha*v_hat;
    u_hat = solve(M,u+q);

    x = x+alpha*u_hat;
    r = r-alpha*A*u_hat;
    error = norm(r)/bnrm2;
    if (iter<(maxiter-1)){
      errors(iter+1) = error;
    }
    if (error <= reltol){
      break;
    }
    rho_1 = rho;
    iter += 1;
  }

  if (error <= reltol){
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
// Krylov Iterative Solver 4-2. Conjugate Gradient Squared : Sparse Way
//------------------------------------------------------------------------------
//' @keywords internal
//' @noRd
// [[Rcpp::export(linsolve.cgs.single.sparse)]]
Rcpp::List single_cgs_sparse(const arma::sp_mat A, const arma::sp_mat b, arma::colvec& xinit,
                             const double reltol, const int maxiter, const arma::sp_mat M){
  // 1-1. parameter settings
  const int n = A.n_cols;
  int iter=0;
  int flag=0; // 0 found; 1 no convergence

  // 1-2. pre-iteration
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
  arma::colvec r_tld = r;


  // 1-3. main iteration
  double rho = 0;
  double rho_1 = 0;
  double alpha = 0;
  double beta = 0;
  
  arma::colvec u(n,fill::zeros);
  arma::colvec p(n,fill::zeros);
  arma::colvec q(n,fill::zeros);
  arma::colvec p_hat(n,fill::zeros);
  arma::colvec v_hat(n,fill::zeros);
  arma::colvec u_hat(n,fill::zeros);

  double error = norm(r)/bnrm2;
  arma::vec errors(maxiter,fill::zeros);
  errors(0) = error;
  while (iter<maxiter){
    rho = dot(r_tld,r);
    if (rho==0.0){
      break;
    }

    if (iter>0){ // direction vectors
      beta = rho/rho_1;
      u = r+beta*q;
      p = u+beta*(q+beta*p);
    } else {
      u = r;
      p = u;
    }

    p_hat = spsolve(M,p,"lapack");
    v_hat = A*p_hat; // adjusting scalars
    alpha = rho/dot(r_tld,v_hat);
    q     = u-alpha*v_hat;
    u_hat = spsolve(M,u+q,"lapack");

    x = x+alpha*u_hat;
    r = r-alpha*A*u_hat;
    error = norm(r)/bnrm2;
    if (iter<(maxiter-1)){
      errors(iter+1) = error;
    }
    if (error <= reltol){
      break;
    }
    rho_1 = rho;
    iter += 1;
  }

  if (error <= reltol){
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
