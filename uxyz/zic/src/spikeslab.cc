#include "spikeslab.h"

SpikeSlabParameter::SpikeSlabParameter( const int Ak, const SpikeSlabPrior& Aprior )
  : k(Ak), beta(zeros<vec>(k+1)), tau(ones<vec>(k)), kappa(ones<vec>(k)), omega(0.5), pVinv(zeros<mat>(k+1,k+1)), prior(Aprior)
{
  pVinv(0,0) = 1.0 / prior.Va;
  for( int i=1; i<=k; i++ )
    pVinv(i,i) = 1.0 / prior.Vb;
};

void SpikeSlabParameter::update( const vec& y, const mat& X, const double sigma2, Rnd& rnd )
{
  if( prior.svs )
    {
      for( int i=0; i<k; i++ )
	pVinv(i+1,i+1) = 1.0 / ( tau(i) * kappa(i) );
    }
  
  const mat V = inv( X.t() * X / sigma2 + pVinv );
  const vec m = V * ( X.t() * y / sigma2 );
  beta = rnd.mnormal( m, V );  

  if( prior.svs )
    {
      kappa_update( rnd );
      tau_update( rnd );
      omega_update( rnd );
    }
}

void SpikeSlabParameter::kappa_update( Rnd& rnd )
{
  for( int i=0; i<k; i++ )
    kappa(i) = rnd.invGamma( prior.g + 0.5, prior.h + 0.5 * beta(i+1) * beta(i+1) / tau(i) );
}

void SpikeSlabParameter::tau_update( Rnd& rnd )
{
  for( int i=0; i<k; i++ )
    {
      double bsq = beta(i+1) * beta(i+1);
      double p0 = ( 1.0 - omega ) / sqrt( prior.nu ) * exp( -0.5 * bsq / ( prior.nu * kappa(i) ) );
      double p1 =         omega                      * exp( -0.5 * bsq /              kappa(i)   );

      if ( p1 / ( p0+p1 ) > rnd.uniform() )
	tau(i) = 1.0;
      else
	tau(i) = prior.nu;
    }
}

void SpikeSlabParameter::omega_update( Rnd& rnd )
{
  int n_in = accu( tau > 0.99 );
  omega = rnd.beta( n_in + prior.r, k - n_in + prior.s );
}

vec SpikeSlabParameter::operator()( void )
{
  return beta;
}

vec SpikeSlabParameter::get_tau( void )
{
  return tau;
}

double SpikeSlabParameter::get_omega( void )
{
  return omega;
}
