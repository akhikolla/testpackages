#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

// using namespace Rcpp;

/*

  Fast Fourier Transform for Rcpp complex vectors

  This function performs the fft by transforming it to a RcppArmadillo
  vector ann then use the function in this library.

*/

Rcpp::ComplexVector fft_rcpp(const Rcpp::NumericVector & real,
                             const Rcpp::NumericVector & imag)
{

  arma::cx_vec out(Rcpp::as<arma::vec>(real), Rcpp::as<arma::vec>(imag));

  out =  arma::fft(out);

  return Rcpp::wrap(out);
}

Rcpp::ComplexVector fft_rcpp_2(const Rcpp::ComplexVector & v)
{
  return Rcpp::wrap(arma::fft(Rcpp::as<arma::cx_vec>(v)));
}

Rcpp::ComplexVector fft_rcpp_3(const Rcpp::ComplexVector & v)
{
  arma::cx_vec out = Rcpp::as<arma::cx_vec>(v);
  return Rcpp::wrap(arma::fft(out));
}

arma::cx_vec fft_rcpp_4(const arma::cx_vec & v)
{
  return arma::fft(v);
}
