// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include "utils.h"
using namespace Rcpp;


// [[Rcpp::export]]
arma::mat get_Z_mat(arma::vec ZOneDim, int m, int n){

  arma::mat Zmat = arma::zeros<arma::mat>(m, n);

  for (int j = 0; j < n; j++) {
    Zmat(ZOneDim(j) - 1, j) = 1;
  }
  return(Zmat);
}

const double log2pi = std::log(2.0 * M_PI);

// [[Rcpp::export]]
arma::vec dmvnrm_arma(arma::mat x,
                      arma::rowvec mean,
                      arma::mat sigma,
                      bool logd) {
  int n = x.n_rows;
  int xdim = x.n_cols;
  arma::vec out(n);
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;

  for (int i=0; i < n; i++) {
    arma::vec z = rooti * arma::trans( x.row(i) - mean) ;
    out(i)      = constants - 0.5 * arma::sum(z%z) + rootisum;
  }

  if (logd == false) {
    out = exp(out);
  }
  return(out);
}

// [[Rcpp::export]]
double calculate_Ratio(double logDeno, arma::vec logNume){

  int n = logNume.n_elem;
  double maxNume = arma::max(logNume);
  double transDeno = logDeno - maxNume;

  arma::vec repMaxNume = rep(maxNume,n);

  arma::vec transNume = logNume - repMaxNume;
  double ratio = exp(transDeno)/sum(exp(transNume));

  return(ratio);
}





