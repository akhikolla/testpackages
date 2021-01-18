#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

using namespace Rcpp;
using namespace arma;

//' A function for sampling from conditional multivariate normal distributions with mean A^{-1}b and covariance matrix A^{-1}.
//'
//' @param a \code{a} A scalar for the Gaussian full conditional distribution precision.
//' @param b \code{b} A \eqn{d} \code{vector} for the Gaussian full conditional distribution mean.
//'
//' @examples
//' set.seed(111)
//' a <- 4
//' b <- rnorm(1)
//' sample <- rmvn_arma_scalar(a, b)
//'
//' @export
//[[Rcpp::export]]
double rmvn_arma_scalar(const double & a, const double & b) {
    double a_inv = 1.0 / a;
    double z = R::rnorm(0, 1);
    return(b * a_inv + z * sqrt(a_inv));
}
