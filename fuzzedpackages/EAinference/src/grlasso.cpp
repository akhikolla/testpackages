#include <RcppArmadillo.h>
#include "Utdbeta.h"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec grlasso(NumericMatrix X, NumericVector Y, NumericMatrix XY,
                     NumericVector weights, arma::rowvec group,
                     double lbd,
                     NumericVector Gamma, NumericVector initBeta, double eps)
{
  //int n = Y.size();
  //int p = X.ncol();
  arma::vec Beta = initBeta;
  arma::vec oldBeta;
  int niter = 0;
  int J = max(group);
  /*
  NumericMatrix XY(n, p);
  for (int i=0; i<n; i++) {
  XY(i,_) = Y(i) * X(i,_);
  }
  */
  do{
    oldBeta = Beta;
    for (int j=0; j<J; j++) {
      arma::vec U = Utdbeta(X,Y,XY,Beta);
      arma::uvec ind = arma::find(group == (j+1)); // Find indices
      arma::vec TEMP = U.elem(ind) + Gamma[j] * Beta.elem(ind); // selected indices

      NumericVector TEMP2 = as<NumericVector>(wrap(TEMP));
      arma::vec TEMP3 = std::max(0.0, 1 - lbd * weights[j] / sqrt(sum(pow(TEMP2, 2)))) / Gamma[j] * TEMP;

      for (int i=0; i< as<int>(wrap(ind.n_elem)); i++) { Beta[ind(i)] = TEMP3[i]; }
    }
    niter += 1;
  } while (sum(pow(Beta-oldBeta, 2)) >  eps*eps && niter < 500);
  return Beta;
}
