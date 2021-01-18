// Matt Galloway

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;


//' @title Soft threshold (elementwise) (c++)
//' @description Soft thresholding function.
//'
//' @param s scalar
//' @param tau scalar
//' @return scalar
//' @keywords internal
//'

double softc(const double &s, const double &tau) {

  // soft-thresholding
  double d = 0;
  if (s > 0 && s > tau) return(s - tau);
  else if (s < 0 && -s > tau) return(s + tau);
  else return(d);

}



////-----------------------------------------------------



//' @title Soft threshold (matrix) (c++)
//' @description Elementwise soft thresholding function for matrices. Requires `softc`.
//'
//' @param s matrix
//' @param Tau scalar
//' @keywords internal
//'

void softmatrixc(arma::mat &S, const arma::mat &Tau) {

  // loop over all elements
  for (int i = 0; i < S.n_rows; ++i){
    for (int j = 0; j < S.n_cols; ++j){

      // soft threshold each element
      S(i, j) = softc(S(i, j), Tau(i, j));

    }
  }

}



////-----------------------------------------------------



//' @title Nonzeros (c++)
//' @description This function counts the number of nonzero elements in a matrix
//' @param X matrix
//' @keywords internal
//'

int numzeros(arma::mat &X) {
  
  // loop over all elements of X and count nonzeros
  int num = 0;
  arma::mat::iterator it = X.begin();
  arma::mat::iterator it_end = X.end();

  for (; it != it_end; ++it){
    if (*it != 0){
      num++;
    }
  }

  return(num);
}
