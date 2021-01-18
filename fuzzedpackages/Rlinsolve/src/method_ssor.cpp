#define ARMA_DONT_PRINT_ERRORS
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//------------------------------------------------------------------------------
// Basic Iterative Solver 4-1. SSOR
//------------------------------------------------------------------------------
//' @keywords internal
//' @noRd
// [[Rcpp::export(linsolve.ssor.single)]]
Rcpp::List single_ssor(const arma::mat& A, const arma::colvec& b, arma::colvec& xinit,
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
  arma::colvec xhalf(n,1,fill::zeros);
  arma::colvec xnew(n,1,fill::zeros);
  double xinc = 0.0;
  arma::vec errors(maxiter,fill::zeros);

  // 3. preliminary computations
  arma::mat Lstar = D+w*L;
  arma::mat Ustar = D+w*U;
  arma::mat S1    = -w*L+(1-w)*D;
  arma::mat S2    = -w*U+(1-w)*D;

  arma::colvec c1(n,1,fill::zeros);
  arma::colvec c2(n,1,fill::zeros);
  c1 = solve(trimatl(Lstar),b);
  c2 = solve(trimatu(Ustar),D*c1);
  c2 *= w*(2-w);

  // 4. run iteration
  double bnrm2 = norm(b-A*xold);
  int iter=0;
  while (iter<maxiter){
    xhalf = solve(trimatl(Lstar),S2*xold);
    xnew  = solve(trimatu(Ustar),S1*xhalf);
    xnew += c2;
    xinc = norm(b-A*xnew)/bnrm2;
    errors(iter) = xinc;
    xold = xnew;
    if (xinc<reltol){
      break;
    }
    iter += 1;
  }

  // 5. return outputs
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
// Basic Iterative Solver 4-2. SSOR : Sparse Way
//------------------------------------------------------------------------------
//' @keywords internal
//' @noRd
// [[Rcpp::export(linsolve.ssor.single.sparse)]]
Rcpp::List single_ssor_sparse(const arma::sp_mat A, const arma::sp_mat b, arma::colvec& xinit,
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
  arma::colvec xhalf(n,1,fill::zeros);
  arma::colvec xnew(n,1,fill::zeros);
  double xinc = 0.0;
  arma::vec errors(maxiter,fill::zeros);

  // 3. preliminary computations
  arma::sp_mat Lstar = D+w*L;
  arma::sp_mat Ustar = D+w*U;
  arma::sp_mat S1    = -w*L+(1-w)*D;
  arma::sp_mat S2    = -w*U+(1-w)*D;
  arma::colvec rhsb  = b-xold+xold;

  arma::colvec c1(n,1,fill::zeros);
  arma::colvec c2(n,1,fill::zeros);
  c1 = spsolve(Lstar,rhsb,"lapack");
  c2 = spsolve(trimatu(Ustar),D*c1,"lapack");
  c2 *= w*(2-w);

  // 4. run iteration
  double bnrm2 = norm(b-A*xold);
  int iter=0;
  while (iter<maxiter){
    xhalf = spsolve(trimatl(Lstar),S2*xold,"lapack");
    xnew  = spsolve(trimatu(Ustar),S1*xhalf,"lapack");
    xnew += c2;
    xinc = norm(b-A*xnew)/bnrm2;
    errors(iter) = xinc;
    xold = xnew;
    if (xinc<reltol){
      break;
    }
    iter += 1;
  }

  // 5. return outputs
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
