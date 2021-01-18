#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
//' @title Column sums for sparseMatrix
//'
//' @description
//' Form column sums for sparseMatrix.
//'
//' @param x A sparse matrix, i.e., inheriting from \code{\link[Matrix]{sparseMatrix}}.
//'
//' @details
//' This function is designed to be used for internal \code{RcppArmadillo} functions. Nevertheless it could be applied in R.
//' It loops on the non-zero entries of the \code{\link[Matrix]{sparseMatrix}}. For general uses, the function
//' \code{\link[Matrix]{colSums}} should be prefered.
//'
//' @return column sums of x.
//' 
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//' 
//' @seealso
//' \code{\link[Matrix]{colSums}}, \code{\link[Matrix]{rowSums}}.
//'
arma::vec colSumsiter(const arma::sp_mat& x) {

  int N = x.n_cols;

  // initialize memory for result
  arma::vec result(N);
  result.fill(0.0);

  arma::sp_mat::const_iterator end = x.end();
  //loop only on non-zero entry to sum up the row
  for(arma::sp_mat::const_iterator it = x.begin(); it != end; ++it){
    result[it.col()] = result[it.col()] + *it;
  }
  return result;
}

// [[Rcpp::depends(RcppArmadillo)]]
//' @title Row sums on sparse matrix.
//'
//' @description
//' Form row sums for sparseMatrix.
//'
//' @param x A sparse matrix, i.e., inheriting from \code{\link[Matrix]{sparseMatrix}}.
//'
//' @details
//' This function is designed to be used for internal \code{RcppArmadillo} functions. Nevertheless it could be applied in R.
//' It loops on the non-zero entries of the \code{\link[Matrix]{sparseMatrix}}. For general uses, the function \code{\link[Matrix]{rowSums}} should
//' be prefered.
//'
//' @return row sums of x. 
//' 
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//' 
//' @seealso
//' \code{\link[Matrix]{colSums}}, \code{\link[Matrix]{rowSums}}.
//' 
//' 
arma::vec rowSumsiter(const arma::sp_mat& x) {
  return colSumsiter(x.t());
}


/*** R
sp <- rsparsematrix(5000,5000,density = 0.4)
system.time(test1 <- Matrix::colSums(sp))
system.time(test2 <- colSumsiter(sp))

rm(list = ls())
X <- as.matrix(cbind(runif(30),runif(30)))
pik <- rep(1/5,30)
A <- wpik(X,pik,bound = 1,tore = TRUE,shift = TRUE,toreBound = 0)
# A <- A+ t(A)
image(A)
## No longer working as function not exported
# WaveSampling:::colSumsiter(A)
colSums(A)
## No longer working as function not exported
# rowSumsiter(A)
rowSums(A)

## No longer working as function not exported

# Was a tiny bit slower but we could use this in a c++ code

# sp <- rsparsematrix(5000,500,density = 0.4)
# system.time(test1 <- Matrix::colSums(sp))
# system.time(test2 <- colSumsiter(sp))
# 
# system.time(svd(sp))
# system.time(svds(sp,NSvals = 1,which = "S"))
# 
# system.time(rowSums(sp))
# system.time(rowSumsiter(sp))


*/
