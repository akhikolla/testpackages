#define ARMA_DONT_PRINT_ERRORS
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//------------------------------------------------------------------------------
// Basic Iterative Solver 1-1. Jacobi Method
//------------------------------------------------------------------------------
//' @keywords internal
//' @noRd
// [[Rcpp::export(linsolve.jacobi.single)]]
Rcpp::List single_jacobi(const arma::mat& A, const arma::colvec& b, arma::colvec& xinit,
                         const double reltol, const int maxiter, const double weight){
  // 1. seperate terms
  int n = A.n_rows;
  arma::mat Dinv(n,n,fill::zeros);
  for (int it=0;it<n;it++){
    if (A(it,it)==0){
      Dinv(it,it) = 0;
    } else {
      Dinv(it,it)=1/A(it,it);
    }
  }
  arma::mat R = A-diagmat(A);

  // 2. get parameters and set ready
  arma::colvec xold = xinit;
  if (norm(b-A*xinit)<reltol){
    arma::colvec xtmp(n,fill::randn);
    xold = xtmp;
  }
  arma::colvec xnew(n,1,fill::zeros);
  double xinc = 0.0;
  arma::vec errors(maxiter,fill::zeros);

  double bnrm2 = norm(b-A*xold);
  // 2. run iteration
  int iter=0;
  while (iter<maxiter){
    xnew = weight*Dinv*(b-R*xold)+(1-weight)*xold;
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
// Basic Iterative Solver 1-2. Jacobi Method : Sparse Way
//------------------------------------------------------------------------------
//' @keywords internal
//' @noRd
// [[Rcpp::export(linsolve.jacobi.single.sparse)]]
Rcpp::List single_jacobi_sparse(const arma::sp_mat A, const arma::sp_mat b, arma::colvec& xinit,
                                const double reltol, const int maxiter, const double weight){
  // 1. seperate terms
  const int n = A.n_rows;
  arma::sp_mat Dinv(n,n);
  for (int i=0;i<n;i++){
    if (A(i,i)==0){
      Dinv(i,i) = 0;
    } else {
      Dinv(i,i) = 1/A(i,i);
    }

  }
  arma::sp_mat R = A-diagmat(A);

  // 2. get parameters and set ready
  arma::colvec xold = xinit;
  if (norm(b-A*xinit)<reltol){
    arma::colvec xtmp(n,fill::randn);
    xold = xtmp;
  }
  arma::colvec xnew(n,1,fill::zeros);
  double xinc = 0.0;
  arma::vec errors(maxiter,fill::zeros);

  double bnrm2 = norm(b-A*xold);
  // 2. run iteration
  // Dinv*(b-R*xold) = term1 - term2*xold
  // w*Dinv*(b-R*xold)+(1-w)*xold
  // = w*term1-w*term2+(1-w)*xold
  int iter=0;
  arma::sp_mat term1 = Dinv*b;
  arma::sp_mat term2 = Dinv*R;
  while (iter<maxiter){
    xnew = weight*Dinv*(b-R*xold)+(1-weight)*xold;
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



