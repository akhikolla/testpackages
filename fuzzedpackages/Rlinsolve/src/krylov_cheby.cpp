#define ARMA_DONT_PRINT_ERRORS
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//------------------------------------------------------------------------------
// Krylov Iterative Solver 6-1. Chebyshev
//------------------------------------------------------------------------------
//' @keywords internal
//' @noRd
// [[Rcpp::export(linsolve.cheby.single)]]
Rcpp::List single_cheby(const arma::mat& A, const arma::colvec& b, arma::colvec& xinit,
                        const double reltol, const int maxiter, arma::mat& M,
                        const double eigmax, const double eigmin){
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

  double c = (eigmax-eigmin)/2.0;
  double d = (eigmax+eigmin)/2.0;

  // 1-3. main iteration
  arma::vec errors(maxiter,fill::zeros);
  errors(0) = error;

  double alpha = 0.0;
  double beta  = 0.0;

  arma::colvec z(n,fill::zeros);
  arma::colvec p(n,fill::zeros);
  while (iter<maxiter){ // begin iteration
    z = solve(M,r);
    if (iter>0){
      beta = pow(c*alpha/2.0,2);
      alpha= 1.0/(d-beta);
      p    = z+beta*p;
    } else {
      p = z;
      alpha = 2.0/d;
    }

    x = x + alpha*p; // update approximation
    r = r-alpha*A*p; // update residual
    error = norm(r)/bnrm2;
    if (error<=reltol){ // check convergence
      break;
    }
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
// Krylov Iterative Solver 6-2. Chebyshev : Sparse Way
//------------------------------------------------------------------------------
//' @keywords internal
//' @noRd
// [[Rcpp::export(linsolve.cheby.single.sparse)]]
Rcpp::List single_cheby_sparse(const arma::sp_mat A, const arma::sp_mat b, arma::colvec& xinit,
                               const double reltol, const int maxiter, arma::sp_mat M,
                               const double eigmax, const double eigmin){
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

  double c = (eigmax-eigmin)/2.0;
  double d = (eigmax+eigmin)/2.0;

  // 1-3. main iteration
  arma::vec errors(maxiter,fill::zeros);
  errors(0) = error;

  double alpha = 0.0;
  double beta  = 0.0;

  arma::colvec z(n,fill::zeros);
  arma::colvec p(n,fill::zeros);
  while (iter<maxiter){ // begin iteration
    z = spsolve(M,r,"lapack");
    if (iter>0){
      beta = pow(c*alpha/2.0,2);
      alpha= 1.0/(d-beta);
      p    = z+beta*p;
    } else {
      p = z;
      alpha = 2.0/d;
    }

    x = x + alpha*p; // update approximation
    r = r-alpha*A*p; // update residual
    error = norm(r)/bnrm2;
    if (error<=reltol){ // check convergence
      break;
    }
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
