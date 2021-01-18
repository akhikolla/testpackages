#define ARMA_DONT_PRINT_ERRORS
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//------------------------------------------------------------------------------
// Basic Iterative Solver 2-1. Gauss-Seidel Method
//------------------------------------------------------------------------------
//' @keywords internal
//' @noRd
// [[Rcpp::export(linsolve.gs.single)]]
Rcpp::List single_gs(const arma::mat& A, const arma::colvec& b, arma::colvec& xinit,
                     const double reltol, const int maxiter, const int dflagval){
  // 1. seperate terms
  vec d = diagvec(A);
  mat D = diagmat(d);
  mat E = -trimatl(A,-1); // omit diagonal for lower
  mat F = -trimatu(A, 1); // omit diagonal for upper

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

  double bnrm2 = norm(b-A*xold);
  // 2. run iteration
  int iter=0;
  if (dflagval==1){
    mat Lstar = D-E;
    while (iter<maxiter){
      xnew = solve(trimatl(Lstar),b+F*xold);
      xinc = norm(b-A*xnew)/bnrm2;
      errors(iter) = xinc;
      xold = xnew;
      if (xinc<reltol){
        break;
      }
      iter += 1;
    }
  } else if (dflagval==3){
    mat Lstar = D-E;
    mat Ustar = D-F;
    colvec xhalf(n,1,fill::zeros);
    while (iter<maxiter){
      xhalf = solve(trimatl(Lstar),b+F*xold);
      xnew  = solve(trimatu(Ustar),b+E*xhalf);
      xinc = norm(b-A*xnew)/bnrm2;
      errors(iter) = xinc;
      xold = xnew;
      if (xinc<reltol){
        break;
      }
      iter += 1;
    }
  } else {
    mat Ustar = D-F;
    while (iter<maxiter){
      xnew = solve(trimatu(Ustar),b+E*xold);
      xinc = norm(b-A*xnew)/bnrm2;
      errors(iter) = xinc;
      xold = xnew;
      if (xinc<reltol){
        break;
      }
      iter += 1;
    }
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
// Basic Iterative Solver 2-2. Gauss-Seidel Method : Sparse Way
//------------------------------------------------------------------------------
//' @keywords internal
//' @noRd
// [[Rcpp::export(linsolve.gs.single.sparse)]]
Rcpp::List single_gs_sparse(const arma::sp_mat A, const arma::sp_mat b, arma::colvec& xinit,
                            const double reltol, const int maxiter, const int dflagval){
  // 1. seperate terms
  int n = A.n_rows;
  arma::sp_mat D = diagmat(A);
  arma::sp_mat E(n,n);
  for (int i=1;i<n;i++){
    for (int j=0;j<i;j++){
      E(i,j) = -A(i,j);
    }
  }
  arma::sp_mat F(n,n);
  for (int i=0;i<(n-1);i++){
    for (int j=(i+1);j<n;j++){
      F(i,j) = -A(i,j);
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

  double bnrm2 = norm(b-A*xold);
  // 2. run iteration
  int iter=0;
  if (dflagval==1){
    arma::sp_mat Lstar = D-E;
    while (iter<maxiter){
      xnew = spsolve(trimatl(Lstar),b+F*xold,"lapack");
      xinc = norm(b-A*xnew)/bnrm2;
      errors(iter) = xinc;
      xold = xnew;
      if (xinc<reltol){
        break;
      }
      iter += 1;
    }
  } else if (dflagval==3){
    arma::sp_mat Lstar = D-E;
    arma::sp_mat Ustar = D-F;
    colvec xhalf(n,1,fill::zeros);
    while (iter<maxiter){
      xhalf = spsolve(trimatl(Lstar),b+F*xold,"lapack");
      xnew  = spsolve(trimatu(Ustar),b+E*xhalf,"lapack");
      xinc = norm(b-A*xnew)/bnrm2;
      errors(iter) = xinc;
      xold = xnew;
      if (xinc<reltol){
        break;
      }
      iter += 1;
    }
  } else {
    arma::sp_mat Ustar = D-F;
    while (iter<maxiter){
      xnew = spsolve(trimatu(Ustar),b+E*xold,"lapack");
      xinc = norm(b-A*xnew)/bnrm2;
      errors(iter) = xinc;
      xold = xnew;
      if (xinc<reltol){
        break;
      }
      iter += 1;
    }
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
