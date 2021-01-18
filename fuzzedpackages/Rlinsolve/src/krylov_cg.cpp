#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//------------------------------------------------------------------------------
// Krylov Iterative Solver 1-1. Conjugate Gradient
//------------------------------------------------------------------------------
//' @keywords internal
//' @noRd
// [[Rcpp::export(linsolve.cg.single)]]
Rcpp::List single_cg(const arma::mat& A, const arma::colvec& b, arma::colvec& xinit,
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
  double error = norm(r)/bnrm2;

  // 1-3. main iteration
  double rho = 0;
  double rho_1 = 0;
  double alpha = 0;
  double beta = 0;
  
  arma::colvec z(n,fill::zeros);
  arma::colvec p(n,fill::zeros);
  arma::colvec q(n,fill::zeros);

  arma::vec errors(maxiter,fill::zeros);
  errors(0) = error;
  while (iter<maxiter){
    z   = solve(M,r);
    rho = dot(r,z);

    if (iter>0){ // direction vector
      beta = rho/rho_1;
      p    = z+beta*p;
    } else {
      p    = z;
    }

    q     = A*p;
    alpha = rho/dot(p,q);
    x     = x+alpha*p; // update approximation vector

    r     = r-alpha*q; // compute residual
    error = norm(r)/bnrm2;
    if (iter<(maxiter-1)){
      errors(iter+1) = error;
    }
    if (error<=reltol){// check convergence
      break;
    }
    rho_1 = rho;
    iter += 1;
  }

  if (error>reltol){
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
// Krylov Iterative Solver 1-2. Conjugate Gradient : Sparse Way
//------------------------------------------------------------------------------
//' @keywords internal
//' @noRd
// [[Rcpp::export(linsolve.cg.single.sparse)]]
Rcpp::List single_cg_sparse(const arma::sp_mat A, const arma::sp_mat b, arma::colvec& xinit,
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
  double error = norm(r)/bnrm2;

  // 1-3. main iteration
  double rho = 0;
  double rho_1 = 0;
  double alpha = 0;
  double beta = 0;
  arma::colvec z(n,fill::zeros);
  arma::colvec p(n,fill::zeros);
  arma::colvec q(n,fill::zeros);

  arma::vec errors(maxiter,fill::zeros);
  errors(0) = error;
  while (iter<maxiter){
    z   = spsolve(M,r,"lapack");
    rho = dot(r,z);

    if (iter>0){ // direction vector
      beta = rho/rho_1;
      p    = z+beta*p;
    } else {
      p    = z;
    }

    q     = A*p;
    alpha = rho/dot(p,q);
    x     = x+alpha*p; // update approximation vector

    r     = r-alpha*q; // compute residual
    error = norm(r)/bnrm2;
    if (iter<(maxiter-1)){
      errors(iter+1) = error;
    }
    if (error<=reltol){// check convergence
      break;
    }
    rho_1 = rho;
    iter += 1;
  }

  if (error>reltol){
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

