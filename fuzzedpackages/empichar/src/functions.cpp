#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

using namespace arma;

/* -------------------------------------------------------------------

   EMPIRICAL CHARACTERISTIC FUNCTIONS

   We will compute the real and imaginary parts of the empirical
   characteristic function as well as its modulus and the empirical
   characteristic function itself. All of this for univariate and
   multivariate cases.

   The purpose of doing it separately it is just for avoiding
   unnecesary computations or use of complex numbers.

--------------------------------------------------------------------*/

// ------------------------------------------------------------------
//           Real part of ecf
// ------------------------------------------------------------------

//[[Rcpp::export]]
arma::vec ecf_re_cpp(const arma::mat & t,
		     const arma::mat & smp)
{
  //  Display error dimensions are different in sample and
  //  eval. points.
  if(t.n_cols != smp.n_cols)
  {
    Rcpp::stop("t and smp must have the same number of columns");
  }

  return mean(cos(t * trans(smp)), 1);
}

// ------------------------------------------------------------------
//           Imaginary part of ecf
// ------------------------------------------------------------------

//[[Rcpp::export]]
arma::vec ecf_im_cpp(const arma::mat & t,
		     const arma::mat & smp)
{
  //  Display error dimensions are different in sample and
  //  eval. points.
  if(t.n_cols != smp.n_cols)
  {
    Rcpp::stop("t and smp must have the same number of columns");
  }

  return mean(sin(t * trans(smp)), 1);
}

// ------------------------------------------------------------------
//           Modulus of ecf
// ------------------------------------------------------------------

//[[Rcpp::export]]
arma::vec ecf_mod_cpp(const arma::mat & t, const arma::mat & smp)
{
  arma::vec real, imag;
  arma::mat arg;

  //  Display error dimensions are different in sample and
  //  eval. points.
  if(t.n_cols != smp.n_cols)
  {
    Rcpp::stop("t and smp must have the same number of columns");
  }

  arg = t * trans(smp);
  real = mean(cos(arg), 1);
  imag = mean(sin(arg), 1);

  return sqrt(real % real + imag % imag);
}

/* --------------------------------------------------------------- */
/*                ecf                                              */
/* --------------------------------------------------------------- */

//[[Rcpp::export]]
arma::cx_vec ecf_cpp(const arma::mat & t, const arma::mat & smp)
{
  // t is mxd
  // smp is nxd
  // m is the number of points and d the dimension.
  arma::vec real, imag;
  arma::mat arg;

  //  Display error dimensions are different in sample and
  //  eval. points.
  if(t.n_cols != smp.n_cols)
  {
    Rcpp::stop("t and smp must have the same number of columns");
  }

  arg = t * trans(smp);
  real = mean(cos(arg), 1);
  imag = mean(sin(arg), 1);

  return arma::cx_vec(real, imag);
}
