#include <RcppArmadillo.h>

double disE(arma::rowvec x, arma::rowvec y){
  double de = std::sqrt( arma::sum( arma::pow( x - y, 2) ) );
  return de;
}

RcppExport SEXP sil_inter(SEXP x_, SEXP clusRes_, SEXP clusNo_, SEXP obsNo_, SEXP breakp_){

  arma::mat x = Rcpp::as<arma::mat>(x_);
  arma::vec clusRes = Rcpp::as<arma::vec>(clusRes_);
  int clusNo = Rcpp::as<int>(clusNo_);
  int obsNo = Rcpp::as<int>(obsNo_);
  arma::vec breakp = Rcpp::as<arma::vec>(breakp_);

  arma::vec a(obsNo);
  arma::vec b(obsNo);
  arma::mat d(obsNo, clusNo);
  arma::mat d0(obsNo, clusNo);
  arma::vec sil(obsNo);
  arma::mat disMat(obsNo, obsNo);
  arma::mat tmp;
  disMat.zeros();

  for( int i = 0; i < obsNo - 1; i++){
    for( int j = i + 1; j < obsNo; j++){
      disMat(j, i) = disMat(i, j) = disE( x.row(i), x.row(j) );
    }
  }

  for( int i = 0; i < obsNo; i++){
    for( int j = 0; j < clusNo; j++){
      if(j == clusRes[i] - 1){
		d0(i, j) = a[i] = arma::sum( disMat(i, arma::span( breakp[j], breakp[j+1] - 1) ) ) / (breakp[j+1] - breakp[j] - 1);
		d(i, j) = std::numeric_limits<double>::infinity();
      } else {
		d0(i, j) = d(i, j) = arma::mean( disMat(i, arma::span( breakp[j], breakp[j+1] - 1) ) );
      }
    }
  }

  for( int i = 0; i < obsNo; i++){
    b[i] = arma::min( d.row(i) );
    sil[i] = (b[i] - a[i]) / std::max(a[i], b[i]);
  }

   return Rcpp::List::create( Rcpp::Named("sil",sil), Rcpp::Named("d0", d0) );
}

