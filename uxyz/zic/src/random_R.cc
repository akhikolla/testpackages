#include "random_R.h"

Rnd::Rnd()
{
  GetRNGstate();
}

Rnd::~Rnd() 
{ 
  PutRNGstate(); 
};

double Rnd::uniform() 
{ 
  return unif_rand(); 
};

double Rnd::normal() 
{ 
  return norm_rand(); 
};

double Rnd::gamma( const double a, const double b ) 
{ 
  return rgamma( 1, a, 1.0 / b )[0]; 
};

double Rnd::invGamma( const double a, const double b ) 
{
  // generates from inverse gamma distribution with pdf = b^a / Gamma(a) * x^(-a-1) * exp(-b/x)
  // E = b / (a-1) for a > 1
  // V = E^2 / (a-1) for a > 2

  return 1.0 / gamma( a, b );
};

int Rnd::poisson( const double a ) 
{ 
  return rpois( 1, a )[0]; 
};

int Rnd::poissont0( const double a )
{
  int x;
  do {
    x = poisson( a );
  } while( x == 0 );
  return x;
}

double Rnd::beta( const double a, const double b )
{
  return rbeta( 1, a, b )[0];
}

void Rnd::normal( vec& x )
{
  for(unsigned int i=0; i<x.n_rows; i++)
    x(i) = normal();
}

vec Rnd::mnormal( const vec& m, const mat& V ) 
{
  vec d( m.n_rows );
  normal( d );
  mat L = chol( V );
  return m + trans(L) * d;
};

double Rnd::normal_lt( const double a, const double var )
{
  // as in the GSL scientific library; see also Knuth, v2, 3rd ed, pp 123-128

  double s = a / sqrt( var );
  
  if( s < 1.0 )
    {
      double x;
      do
	{
	  x = normal();
	} 
      while( x < s );
      return x * sqrt( var );
    }

  else
    {
      double u, v, x;
      do
	{
	  u = uniform();
	  do
	    {
	      v = uniform();
	    }
	  while( v == 0.0 );
	  x = sqrt (s * s - 2 * log (v));
	}
      while( x * u > s );
      return x * sqrt( var );
    }
};

double Rnd::normal_lt( const double a, const double mean, const double var )
{
  return mean + normal_lt( a - mean, var );
};

double Rnd::normal_rt( const double a, const double mean, const double var )
{
  return -normal_lt( -a, -mean, var );
};

double Rnd::t( const double mu, const double sigma2, const double nu )
{
  // Draws from Student's t distribution with density 
  //    p = Gamma[(nu+1)/2] / Gamma(nu/2) / sqrt(pi*sigma2*nu) * [1+(x-mu)^2/(nu*sigma2)]^[-(nu+1)/2]
  // Parameters: location parameter mu, scale parameter sigma2 > 0, degrees of freedom nu > 0
  // Moments: E = mu, var = sigma2 * nu / (nu-2)

  return mu + sqrt(sigma2) * rt(1, nu)[0];
}
