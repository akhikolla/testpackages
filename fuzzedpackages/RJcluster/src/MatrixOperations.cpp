#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat tcrossprod_c(const arma::mat& X, const arma::mat& Y)
{
  arma::mat XXt_return = X * Y.t();
  return(XXt_return);
}

// [[Rcpp::export]]
arma::mat scale_c(const arma::mat& X, bool medians = false)
{
  int p = X.n_cols;
  int n = X.n_rows;

  // get medians or mean
  arma::mat scaleValue = arma::mean(X, 0);
  if (medians)
  {
    scaleValue = arma::median(X, 0);
  }

  arma::mat colMeans(n,p);
  for (int i = 0; i < n; i++)
  {
    colMeans.row(i) = scaleValue;
  }

  // center data
  arma::mat centered = X - colMeans; // this may need to be changed

  // get norms and weight
  arma::mat weights = arma::sqrt(arma::sum(arma::pow(centered, 2) / n, 0));

  // get final Xtilde
  arma::mat Xtilde = centered * arma::diagmat(pow(weights, -1));

  return(Xtilde);
}

