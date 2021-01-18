#include "zicmodel.h"

ZicModel::ZicModel( const ivec& Ay, const mat& AX, 
                    const SpikeSlabPrior& Abetaprior, const SpikeSlabPrior& Adeltaprior, const double Atune )
  : y(Ay), X(AX), n(y.n_rows), n0(accu(y==0)), k(X.n_cols-1), eta(zeros<vec>(n)), dstar(zeros<vec>(n)),
    beta(X.n_cols-1,Abetaprior), delta(X.n_cols-1,Adeltaprior), sigma2(1.0), tune(Atune), betaprior(Abetaprior), deltaprior(Adeltaprior)
{
}; 
 
void ZicModel::sample( const int nburnin, const int nmcmc, const int nthin, 
                       vec& alphamem, mat& betamem, vec& gammamem, mat& deltamem, vec& sigma2mem, 
                       imat& Ibetamem, imat& Ideltamem, vec& omegabetamem, vec& omegadeltamem, vec& accmem )
{
  Rcout << "Gibbs sampler is running.\n";    
  Rcout << "Progress: [                    ]\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b";
  
  int save_idx = 0;
  for( int i=1; i<=nburnin+nmcmc; i++ )
    {
      if( i%((nburnin+nmcmc)/20) == 0 )
	Rcout << "#";

      beta.update( eta, X, sigma2, rnd );
      sigma2_update();
      delta.update( dstar, X, 1.0, rnd );
      latvar_update();     
      
      if( ( i>nburnin ) && ( i%nthin==0 ) )
	{
	  alphamem(save_idx) = beta()(0);
	  betamem.row(save_idx) = beta().subvec(1,k).t(); 
	  gammamem(save_idx) = delta()(0);
	  deltamem.row(save_idx) = delta().subvec(1,k).t(); 
	  sigma2mem(save_idx) = sigma2;
	  Ibetamem.row(save_idx) = conv_to<ivec>::from( beta.get_tau() > 0.99 ).t(); 
	  Ideltamem.row(save_idx) = conv_to<ivec>::from( delta.get_tau() > 0.99 ).t(); 
	  omegabetamem(save_idx) = beta.get_omega();
	  omegadeltamem(save_idx) = delta.get_omega();
	  accmem(save_idx) = acc;

	  save_idx++;
	}
    }
  
  Rcout << "]\nGibbs sampler is finished.\n";
}

void ZicModel::latvar_update()   // checked 
{
  const vec etaold = eta;
  const vec Xbeta = X * beta();
  const vec Xdelta = X * delta();

  for( int i=0; i<n0; i++ )
    {
      // update eta
      double p = Rcpp::stats::pnorm_0( Xdelta(i), true, false );   
      double eta_prob = eta(i) + tune * rnd.t( 0.0, sigma2, 5.0 );
      double logp_prob = log( 1.0 - p + p * exp(-exp(eta_prob)) ) - 0.5 * (eta_prob-Xbeta(i)) * (eta_prob-Xbeta(i)) / sigma2;
      double logp_old  = log( 1.0 - p + p * exp(-exp(eta(i)  )) ) - 0.5 * (eta(i)  -Xbeta(i)) * (eta(i)  -Xbeta(i)) / sigma2;
      if( exp( logp_prob - logp_old ) > rnd.uniform() )
	eta(i) = eta_prob;
    
      // update dstar
      if( rnd.uniform() > (1.0 - p) / ( 1.0 - p + p * exp(-exp(eta(i))) ) )
	dstar(i) = rnd.normal_lt( 0.0, Xdelta(i), 1.0 );
      else
	dstar(i) = rnd.normal_rt( 0.0, Xdelta(i), 1.0 );
    }

  for( int i=n0; i<n; i++ )
    {
      // update eta
      double eta_prob = eta(i) + tune * rnd.t( 0.0, sigma2, 5.0 );
      double logp_prob = -exp(eta_prob) + eta_prob * y(i) - 0.5 * (eta_prob-Xbeta(i)) * (eta_prob-Xbeta(i)) / sigma2;
      double logp_old  = -exp(eta(i)  ) + eta(i)   * y(i) - 0.5 * (eta(i)  -Xbeta(i)) * (eta(i)  -Xbeta(i)) / sigma2;
      if( exp( logp_prob - logp_old ) > rnd.uniform() )
	eta(i) = eta_prob;
 
      // update dstar
      dstar(i) = rnd.normal_lt( 0.0, Xdelta(i), 1.0 );
    }

  acc = 1.0 - mean( conv_to<vec>::from( eta == etaold ) );
}

void ZicModel::sigma2_update()   // checked
{
  const vec resid = eta - X * beta();
  const double ssr = dot( resid, resid );
  sigma2 = rnd.invGamma( betaprior.e + 0.5 * n, betaprior.f + 0.5 * ssr );
}


