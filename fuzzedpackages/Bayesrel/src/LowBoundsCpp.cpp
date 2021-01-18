

#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]



//[[Rcpp::export]]
double alphaArma(const arma::mat& X) {
  double k = X.n_cols;
    double out = k / (k - 1.0) * (1.0 - arma::trace(X) / arma::accu(X));
    return out;
}

//[[Rcpp::export]]
double l2Arma(const arma::mat& X) {
    double k  = X.n_cols,
           s  = 0.0,
           s2 = 0.0;
    for (unsigned int i = 0; i < X.n_cols - 1; i++)
        for (unsigned int j = i + 1; j < X.n_rows; j++)
        {
            s  += 2.0 * X(i, j);
            s2 += 2.0 * X(i, j) * X(i, j);
        }
    double s3 = s + arma::accu(X.diag());
    double out = (s + std::sqrt(k / (k - 1) * s2)) / s3;
    return out;
}

//[[Rcpp::export]]
double l6Arma(const arma::mat& X) {
  // correlation matrix from covariance matrix:
  arma::vec sds = 1/arma::sqrt(X.diag());
  arma::mat Xcor = arma::diagmat(sds) * X * arma::diagmat(sds);
  Xcor.diag().ones();
  arma::mat XCorInv = arma::inv_sympd(Xcor);
  arma::vec smc = 1 - 1 / XCorInv.diag();
  double out = 1 - arma::accu(1 - smc) / arma::accu(Xcor);
  return out;
}




