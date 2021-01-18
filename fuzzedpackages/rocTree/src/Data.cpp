#include "Data.h"

Data::Data(
  Rcpp::IntegerMatrix mat1Z,
  Rcpp::NumericMatrix mat1f,
  Rcpp::List mat2Z,
  Rcpp::IntegerMatrix range0)
{
  this->mat1Z = mat1Z;
  this->mat1f = mat1f;
  // List mat2(10);
  // for(int i = 0; i< 10; i++)
  // {
  //   mat2[i] = as<arma::umat>(mat2Z[i]);
  // }
  this->mat2Z = mat2Z;
  this->range0 = range0;
}

// [[Rcpp::export]]
arma::field<arma::mat> make_mat2(arma::vec& tk0,
				 arma::vec& Y0,
				 arma::uvec& id0,
				 arma::mat& z0) {
  int K = tk0.size();
  arma::field<arma::mat> out(K);
  for(int k = 0; k < K; k++) {
    arma::uvec ind = arma::find( Y0 >= tk0[k] );
    arma::uvec idAR = id0(ind);
    arma::mat subX = z0.rows(ind);
    arma::uvec ind2 = arma::find(Rcpp::as<arma::uvec>(Rcpp::duplicated(Rcpp::IntegerVector(idAR.begin(), idAR.end()) )) == 0);
    out[k] = subX.rows(ind2);
  }
  return out;
}

// [[Rcpp::export]]
arma::field<arma::mat> make_mat2_t(arma::vec& tk0,
				   arma::vec& Y0,
				   arma::uvec& id0,
				   arma::mat& z0) {
  int K = tk0.size();
  arma::field<arma::mat> out(K);
  for(int k = 0; k < K; k++) {
    arma::uvec ind = arma::find( Y0 >= tk0[k] );
    arma::uvec idAR = id0(ind);
    arma::mat subX = z0.rows(ind);
    arma::uvec ind2 = arma::find(Rcpp::as<arma::uvec>(Rcpp::duplicated(Rcpp::IntegerVector(idAR.begin(), idAR.end()) )) == 0);
    out[k] = subX.rows(ind2).t();
  }
  return out;
}


Rcpp::IntegerMatrix Data::get_mat1Z() const {
  return mat1Z;
}

Rcpp::NumericMatrix Data::get_mat1f() const {
  return mat1f;
}

Rcpp::IntegerMatrix Data::get_range0() const {
  return range0;
}

Rcpp::List Data::get_mat2ZK() const {
  return mat2Z;
}

Rcpp::IntegerMatrix Data::get_mat2Zk( int k ) const {
  return mat2Z[k];
}

arma::uword Data::get_len_mat2k( int k ) const {
  arma::umat m = mat2Z[k];
  return m.n_rows;
}
