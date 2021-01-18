#ifndef _RANDOM_R_HH
#define _RANDOM_R_HH

#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;


class Rnd
{
public:
  Rnd();
  ~Rnd();

  double uniform();
  double normal( void );
  double gamma( const double a, const double b );
  double invGamma( const double a, const double b );
  int poisson( const double a );
  int poissont0( const double a );
  double beta( const double a, const double b );
  void normal( vec& x );
  vec mnormal(const vec& m, const mat& V); 
  double normal_lt( const double a, const double var );
  double normal_lt( const double a, const double mean, const double var );
  double normal_rt( const double a, const double mean, const double var );  
  double t( const double mu, const double sigma2, const double nu );  
};

#endif

