#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
void Test() {

}


/*** R
# returns true if the cholesky or inverse is succesful
library(mvtnorm)
set.seed(54)
X <- mvtnorm::pmvnorm(upper = 0, mean = , diag(4))
Test(X)

log(det(X))

*/
