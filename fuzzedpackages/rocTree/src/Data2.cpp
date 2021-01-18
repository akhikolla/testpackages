#include "Data2.h"
#include "globals.h"


arma::field<arma::umat> Data2::get_zt() const
{
  arma::field<arma::umat> czt(zt.size());
  for(int i = 0; i < zt.size(); i++) {
    czt(i) = Rcpp::as<arma::umat>(zt[i]);
  }
  return czt;
}

arma::umat Data2::get_zy() const
{
  arma::umat mzy = Rcpp::as<arma::umat>(zy);
  return mzy;
}

arma::uvec Data2::get_e() const
{
  arma::uvec ve = Rcpp::as<arma::uvec>(E);
  return ve;
}

arma::vec Data2::get_Y() const
{
  arma::vec vy = Rcpp::as<arma::vec>(Y);
  return vy;
}
