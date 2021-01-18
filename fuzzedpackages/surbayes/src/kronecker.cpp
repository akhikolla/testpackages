

#include "post.h"
// [[Rcpp::depends(RcppArmadillo)]]



//' Fast kronecker product with response vector
//' 
//' This is a c++ implementation of the fast kronecker product with response vector
//'
//' @param Sigma covariance matrix
//' @param Y matrix of response variables (Y1, ..., YJ)
//' @param n number of observations
//' @param J number of endpoints
//' @return Returns a vector with result of \code{ kron(Sigma, diag(n)) \% y }
//' @keywords internal
// [[Rcpp::export]]
arma::vec fastKronEye_Y (
    arma::mat const& Sigma,
    arma::mat const& Y,
    int const& n,
    int const& J
) {
  arma::mat const& Sigmarep = repelem(Sigma, n, 1);   // repeat each row of Sigma n times
  arma::mat const& Yrep = repmat(Y, J, 1);            // repeat each row of Y J times
  arma::vec res = sum(Sigmarep % Yrep, 1);            // rowSums(Sigmarep % Yrep) % = element-wise product
  return(res);
}



//' Fast kronecker product of crossproduct matrix
//' 
//' This is a c++ implementation of the fast kronecker product
//' t(X) %*% kron(Sigma, I) %*% X. It avoids computing the kronecker product
//'
//' @param XtX a matrix that is crossprod((X1, ..., XJ)) in R
//' @param Sigma JxJ covariance matrix
//' @param pvec J-dimensional vector giving number of observations for each endpoint
//' @param n number of observations
//' @param J number of endpoints
//' @keywords internal
//' @return \code{matrix} result of \eqn{X' (\Sigma \otimes I_n) X}
// [[Rcpp::export]]
arma::mat fastKronEye_crossprod (
    arma::mat const& XtX,
    arma::mat const& Sigma,
    arma::vec const& pvec,
    int const& n,
    int const& J
) {
  // Initialize two matrices
  //   Sigmabig is a block matrix whose ij block is Sigma[i,j] %*% ones(p_i, p_j).
  //   This is because since X is blkdiag, 
  //      \code{t(X) %*% kron(Sigma, I) %*% X = Sigmabig % [ t(X) %*% X]} where % = element-wise product
  arma::mat Sigmabig;
  arma::mat res;
  for(int j = 0; j < J; j++) {
    int pj = pvec(j);                    // number of covars in jth regression
    arma::mat tempj;                     // temporary matrix holder
    for(int k = 0; k < J; k++) {
      int pk = pvec(k);                   // number of covars in kth regression
      tempj = join_rows(tempj, Sigma(j, k) * arma::ones(pj, pk)); // join_rows = cbind; This is making Sigma[j,k] * Ones(p_j, p_k) matrix
    }
    Sigmabig = join_cols(Sigmabig, tempj);   // join_cols = rbind. We did the above block-column-wise so now we need to bind the rows together for each j
  }
  res = Sigmabig % XtX;                     // % = element-wise product 
  return(res);
}

