#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat cderivate_mat(arma::mat X, arma::rowvec v)
{
  arma::mat dydx(arma::size(X), arma::fill::zeros);
  dydx.shed_col(0);
  
  arma::rowvec diff_v = arma::diff(v);
  
  for(arma::uword i = 0; i < X.n_rows; ++i)
  {
    dydx.row(i) = arma::diff(X.row(i)) - diff_v;
  }
  
  return dydx;
}

// [[Rcpp::export]]
arma::cube cderivate_cube(arma::cube myCube, arma::vec v)
{
  arma::cube dydx(arma::size(myCube), arma::fill::zeros);
  dydx.shed_slice(0);
  
  arma::vec diff_v = arma::diff(v);
  
  for(arma::uword i = 0; i < myCube.n_rows; ++i)
  {
    for(arma::uword j = 0; j < myCube.n_cols; ++j)
    {
      dydx.tube(i,j) = arma::diff(vectorise(myCube.tube(i,j))) - diff_v;
    }
  }
  
  return dydx;
}
