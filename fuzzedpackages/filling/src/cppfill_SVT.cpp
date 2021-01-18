#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

arma::mat cpp_nSVD_thrs(arma::vec s, const double lambda, const int n, const int p){
  arma::mat output(n,p,fill::zeros);
  int m = 0;
  if (n < p){
    m = n;
  } else {
    m = p;
  }
  for (int i=0;i<m;i++){
    if (s(i)>lambda){
      output(i,i) = s(i)-lambda;
    } else {
      output(i,i) = 0.0;
    }
  }
  return(output);
}

// [[Rcpp::export]]
arma::mat cpp_nSVD(arma::mat& X, arma::mat& idmat, arma::mat& Minit,
                   const double lambda, const double tol, const int maxiter){
  // 1. settings
  const int n = X.n_rows;
  const int p = X.n_cols;

  // 2. matrices
  arma::mat Xhat(n,p,fill::zeros);
  arma::mat Mold = Minit;
  arma::mat Mnew(n,p,fill::zeros);

  // 3. indexes
  arma::uvec idx_miss = find(idmat > 0.5);
  arma::uvec idx_obss = find(idmat <=0.5);

  // 4. main computation
  int iter = 0;
  double increment = 1000.0;

  arma::mat U;
  arma::vec s;
  arma::mat V;

  while (increment > tol){
    // 4-1. replace missing entries, oh yeah.
    Xhat(idx_obss) = X(idx_obss);
    Xhat(idx_miss) = Mold(idx_miss);

    // 4-2. soft-thresholded SVD of Xhat
    svd(U, s, V, Xhat);
    arma::mat thrS = cpp_nSVD_thrs(s, lambda, n, p);
    Mnew = U*thrS*V.t();

    // 4-3. update
    increment = arma::norm(Mnew-Mold,"fro");
    Mold = Mnew;
    iter = iter+1;
    if (iter>=maxiter){
      break;
    }
  }

  // 5. return results
  return(Mold);
}
