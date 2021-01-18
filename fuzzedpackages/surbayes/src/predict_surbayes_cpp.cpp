// [[Rcpp::depends(RcppArmadillo)]]
#include<RcppArmadillo.h>
using namespace Rcpp;


//' Get one sample from predictive posterior of SUR
//' 
//' C++ implementation to obtain one sample from predictive posterior
//' density
//'
//' @param mu vector of means
//' @param Sigma covariance matrix shared among all observations
//' @param n number of observations
//' @param J number of endpoints
// [[Rcpp::export]]
arma::mat predict_surbayes_helper(
    arma::vec const& mu,
    arma::mat const& Sigma,
    int const& n,
    int const& J
) {
  // Sample z ~ N(0, 1)
  arma::mat z = arma::mat( n, J, arma::fill::randn );
  
  // Take cholesky decomposition of Sigma
  arma::mat U = arma::chol(Sigma);
  
  // Reshape (n*J x 1) vector mu to be (n x J) matrix
  arma::mat M = reshape( mu, n, J );
  
  // return mu + Z * U;
  return M + z * U;
}








//' Sample from predictive posterior density C++ helper
//' 
//' C++ implementation to obtain a matrix of samples from predictive posterior density
//'
//' @param Mu matrix of means
//' @param Sigmalist list of covariance matrices
//' @param n number of observations
//' @param J number of endpoints
//' @param nsims Number of simulations (number of rows in Mu)
// [[Rcpp::export]]
arma::cube predict_surbayes_cpp(
    arma::mat const& Mu,
    List const& Sigmalist,
    int const& n,
    int const& J,
    int const& nsims
) {
  // Initialize result which is a n x J x nsims cube
  arma::cube res = arma::cube( n, J, nsims, arma::fill::zeros );
  for( int i = 0; i < nsims; i++ ) {
    arma::vec mu = Mu.col(i);
    arma::mat Sigma = Sigmalist[i];
    res.slice(i) = predict_surbayes_helper( mu, Sigma, n, J);
  }
  return res;
}














