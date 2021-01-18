#ifndef Data2_h
#define Data2_h


#include <RcppArmadillo.h>
#include <armadillo>

/* This class is to store the input data for the second set of data, which is
 * used for prediction (calulating the weights after a tree is grown). It may contain data
 * from the same set of subjects as the training sample. It may also contain data
 * from another independent sample (honest tree).
 */

class Data2 {
public:

  Data2(Rcpp::List& zt0, Rcpp::IntegerMatrix& zy0,
        Rcpp::NumericVector& y0, Rcpp::IntegerVector& e0) : zt(zt0), zy(zy0), Y(y0), E(e0) { };

  arma::field<arma::umat> get_zt() const;
  arma::umat get_zy() const;
  arma::uvec get_e() const;
  arma::vec get_Y() const;
private:
  // zt is a matrix of Z(t): each row is one subject, each column is one time point t
  // We only need uncensored time points ?
  // We allow ties in t
  // zy is a vector of Z(Y): each element is one subject
  Rcpp::List zt;
  Rcpp::IntegerMatrix zy;
  Rcpp::NumericVector Y;
  Rcpp::IntegerVector E;

};

#endif /* Data2_h */
