#define ARMA_DONT_PRINT_ERRORS
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//------------------------------------------------------------------------------
// Basic Iterative Solver 3-1. SOR
//------------------------------------------------------------------------------
//' @keywords internal
//' @noRd
// [[Rcpp::export(linsolve.sor.single)]]
Rcpp::List single_sor(const arma::mat& A, const arma::colvec& b, arma::colvec& xinit,
                      const double reltol, const int maxiter, const double w){
  // 1. seperate terms
  vec d = diagvec(A);
  mat D = diagmat(d);
  mat L = trimatl(A,-1); // omit diagonal for lower
  mat U = trimatu(A, 1); // omit diagonal for upper

  // 2. get parameters and set ready
  int n = A.n_rows;
  arma::colvec xold = xinit;
  if (norm(b-A*xinit)<reltol){
    arma::colvec xtmp(n,fill::randn);
    xold = xtmp;
  }
  arma::colvec xnew(n,1,fill::zeros);
  double xinc = 0.0;
  arma::vec errors(maxiter,fill::zeros);

  // 2. run iteration
  double bnrm2 = norm(b-A*xold);
  int iter=0;
  arma::mat Lstar = D+w*L;
  arma::mat Rstar = (w*U+(w-1)*D);
  arma::colvec wb = w*b;
  arma::colvec rhs(n,1,fill::zeros);
  while (iter<maxiter){
    rhs  = wb-Rstar*xold;
    xnew = solve(trimatl(Lstar),rhs);
    xinc = norm(b-A*xnew)/bnrm2;
    errors(iter) = xinc;
    xold = xnew;
    if (xinc<reltol){
      break;
    }
    iter += 1;
  }

  // 3. return outputs
  List res;
  res["x"] = xold;
  res["iter"] = iter;
  if (iter>=maxiter){
    res["errors"] = errors;
  } else {
    arma::vec newerrors = errors.subvec(0,iter);
    res["errors"] = newerrors;
  }
  return res;
}


//------------------------------------------------------------------------------
// Basic Iterative Solver 3-2. SOR : Sparse Way
//------------------------------------------------------------------------------
//' @keywords internal
//' @noRd
// [[Rcpp::export(linsolve.sor.single.sparse)]]
Rcpp::List single_sor_sparse(const arma::sp_mat A, const arma::sp_mat b, arma::colvec& xinit,
                             const double reltol, const int maxiter, const double w){
  // 1. seperate terms
  int n = A.n_rows;
  arma::sp_mat D = diagmat(A);
  arma::sp_mat L(n,n);
  for (int i=1;i<n;i++){
    for (int j=0;j<i;j++){
      L(i,j) = A(i,j);
    }
  }
  arma::sp_mat U(n,n);
  for (int i=0;i<(n-1);i++){
    for (int j=(i+1);j<n;j++){
      U(i,j) = A(i,j);
    }
  }

  // 2. get parameters and set ready
  arma::colvec xold = xinit;
  if (norm(b-A*xinit)<reltol){
    arma::colvec xtmp(n,fill::randn);
    xold = xtmp;
  }
  arma::colvec xnew(n,1,fill::zeros);
  double xinc = 0.0;
  arma::vec errors(maxiter,fill::zeros);

  // 2. run iteration
  double bnrm2 = norm(b-A*xold);
  int iter=0;
  arma::sp_mat Lstar = D+w*L;
  arma::sp_mat Rstar = (w*U+(w-1)*D);
  arma::sp_mat wb = w*b;
  while (iter<maxiter){
    xnew = spsolve(trimatl(Lstar),wb-Rstar*xold,"lapack");
    xinc = norm(b-A*xnew)/bnrm2;
    errors(iter) = xinc;
    xold = xnew;
    if (xinc<reltol){
      break;
    }
    iter += 1;
  }

  // 3. return outputs
  List res;
  res["x"] = xold;
  res["iter"] = iter;
  if (iter>=maxiter){
    res["errors"] = errors;
  } else {
    arma::vec newerrors = errors.subvec(0,iter);
    res["errors"] = newerrors;
  }
  return res;
}
