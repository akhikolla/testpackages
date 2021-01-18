#ifndef _SPIKESLAB_H
#define _SPIKESLAB_H

#include <RcppArmadillo.h>
#include "random_R.h"

using namespace arma;

struct SpikeSlabPrior
{
  double Va, Vb, g, h, nu, r, s, e, f;
  bool svs;
};

class SpikeSlabParameter
{
public:
  SpikeSlabParameter( const int Ak, const SpikeSlabPrior& Aprior );

  void update( const vec& y, const mat& X, const double sigma2, Rnd& rnd );  
  void kappa_update( Rnd& );
  void tau_update( Rnd& );
  void omega_update( Rnd& );

  vec operator()( void );
  vec get_tau( void );
  double get_omega( void );

  friend class ZicModel;
  
private:
  const int k;        // dimension
  vec beta;           // 1st level parameter
  vec tau, kappa;     // lower level parameter
  double omega;       // lower level parameter
  mat pVinv;          // inverse prior variance
  
  SpikeSlabPrior prior;
};

#endif
