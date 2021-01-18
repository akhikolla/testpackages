#ifndef Data_h
#define Data_h


#include <RcppArmadillo.h>
#include <armadillo>

// This class is to store the input data for building the trees.

class Data {
 public:
  Data(Rcpp::IntegerMatrix mat1Z,
       Rcpp::NumericMatrix mat1f,
       Rcpp::List mat2Z,
       Rcpp::IntegerMatrix range0);

  void sort();
  //mat get_mat1f() const;
  //need a few get_XXX functions

  SEXP make_mat2(arma::vec& tk0, arma::vec& Y0, arma::uvec& id0, arma::mat& z0);
  SEXP make_mat2_t(arma::vec& tk0, arma::vec& Y0, arma::uvec& id0, arma::mat& z0);

  arma::uword get_len_mat2k( int k ) const;
  Rcpp::IntegerMatrix get_mat1Z() const;
  Rcpp::NumericMatrix get_mat1f() const;
  Rcpp::IntegerMatrix get_mat2Zk( int k ) const;
  Rcpp::List get_mat2ZK() const;
  Rcpp::IntegerMatrix get_range0() const;

private:
  Rcpp::IntegerMatrix mat1Z;
  Rcpp::NumericMatrix mat1f;
  Rcpp::List mat2Z;
  Rcpp::IntegerMatrix range0;
  // umat mat1Z;
  // mat mat1f;
  // field<umat> mat2Z;
  // umat range0;

  //std::vector<std::vector<std::vector<double>>> mat2k;
};

#endif /* Data_h */
