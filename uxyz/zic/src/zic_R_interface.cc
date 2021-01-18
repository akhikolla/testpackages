#include <RcppArmadillo.h>
#include "zicmodel.h"

extern "C" SEXP zic_sample( SEXP Sy, SEXP SX, 
			    SEXP Sa0, SEXP Sb0, SEXP Sgbeta0,  SEXP Shbeta0, SEXP Snubeta0, 
			    SEXP Srbeta0, SEXP Ssbeta0, SEXP Se0, SEXP Sf0, 
			    SEXP Sc0, SEXP Sd0, SEXP Sgdelta0, SEXP Shdelta0, SEXP Snudelta0, 
			    SEXP Srdelta0, SEXP Ssdelta0,
			    SEXP Ssvs0, SEXP Snburnin, SEXP Snmcmc, SEXP Snthin, SEXP Stune ) 
{
  try 
    {
       // data matrices and dimensions
       Rcpp::IntegerVector yR( Sy );  
       Rcpp::NumericMatrix XR( SX );  
       const int n = XR.nrow();
       const int k = XR.ncol() - 1;
       const ivec y( yR.begin(), yR.size(), false );
       const mat X( XR.begin(), n, k+1, false );   

       // prior
       SpikeSlabPrior betaprior, deltaprior;
       
       betaprior.Va = Rcpp::as<double>( Sa0 );
       betaprior.Vb = Rcpp::as<double>( Sb0 );
       betaprior.g  = Rcpp::as<double>( Sgbeta0 );
       betaprior.h  = Rcpp::as<double>( Shbeta0 );
       betaprior.nu = Rcpp::as<double>( Snubeta0 );
       betaprior.r  = Rcpp::as<double>( Srbeta0 );
       betaprior.s  = Rcpp::as<double>( Ssbeta0 );
       betaprior.e   = Rcpp::as<double>( Se0 );
       betaprior.f   = Rcpp::as<double>( Sf0 );
       betaprior.svs = Rcpp::as<bool>( Ssvs0 );
       
       deltaprior.Va   = Rcpp::as<double>( Sc0 );
       deltaprior.Vb   = Rcpp::as<double>( Sd0 );
       deltaprior.g    = Rcpp::as<double>( Sgdelta0 );
       deltaprior.h    = Rcpp::as<double>( Shdelta0 );
       deltaprior.nu   = Rcpp::as<double>( Snudelta0 );
       deltaprior.r    = Rcpp::as<double>( Srdelta0 );
       deltaprior.s    = Rcpp::as<double>( Ssdelta0 );
       deltaprior.e    = -9.9;
       deltaprior.f    = -9.9;
       deltaprior.svs  = Rcpp::as<bool>( Ssvs0 );
       
       // mcmc iterations
       const int nburnin = Rcpp::as<int>( Snburnin );
       const int nmcmc = Rcpp::as<int>( Snmcmc );
       const int nthin = Rcpp::as<int>( Snthin );

       // tuning parameter
       const double tune = Rcpp::as<double>( Stune );
              
       // memories
       vec alphamem( nmcmc/nthin ), gammamem( nmcmc/nthin ), sigma2mem( nmcmc/nthin ), 
           omegabetamem( nmcmc/nthin ), omegadeltamem( nmcmc/nthin ), accmem( nmcmc/nthin );
       mat betamem( nmcmc/nthin, k ), deltamem( nmcmc/nthin, k ); 
       imat Ibetamem( nmcmc/nthin, k ), Ideltamem( nmcmc/nthin, k );
            
       // initialize model and sample
       ZicModel m( y, X, betaprior, deltaprior, tune );
       m.sample( nburnin, nmcmc, nthin, alphamem, betamem, gammamem, deltamem, sigma2mem, Ibetamem, Ideltamem, 
                 omegabetamem, omegadeltamem, accmem );
       
       // return memories
       Rcpp::List memlist;

       if( betaprior.svs )
	 {
	   return Rcpp::List::create( Rcpp::Named( "alpha" )       = alphamem,
	                              Rcpp::Named( "beta" )        = betamem,
	                              Rcpp::Named( "gamma" )       = gammamem,
		                      Rcpp::Named( "delta" )       = deltamem,
                                      Rcpp::Named( "sigma2" )      = sigma2mem,
                                      Rcpp::Named( "I.beta" )      = Ibetamem,
		                      Rcpp::Named( "I.delta" )     = Ideltamem,
		                      Rcpp::Named( "omega.beta" )  = omegabetamem,
		                      Rcpp::Named( "omega.delta" ) = omegadeltamem,
				      Rcpp::Named( "acc" )         = accmem         );
	 }
       else
         {
	   return Rcpp::List::create( Rcpp::Named( "alpha" )       = alphamem,
	                              Rcpp::Named( "beta" )        = betamem,
	                              Rcpp::Named( "gamma" )       = gammamem,
		                      Rcpp::Named( "delta" )       = deltamem,
                                      Rcpp::Named( "sigma2" )      = sigma2mem,
				      Rcpp::Named( "acc" )         = accmem         );
	 }
     } 
  catch( std::exception &ex ) 
    {
      forward_exception_to_r( ex );
    } 
  catch( ... ) 
    { 
      ::Rf_error( "C++ exception (unknown reason)" ); 
    }
  return R_NilValue; 
}





