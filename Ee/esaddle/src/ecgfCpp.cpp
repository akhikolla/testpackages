#include "esaddle.h"

/*
 * Fast computation of Empirical Cumulant Generationg Function ECGF
*/
  
SEXP ecgfCpp(SEXP lambda_, SEXP X_, SEXP mix_, SEXP grad_, SEXP kum1_, SEXP kum2_)
{
    using namespace arma;
    
    try{
      
      vec lambda = Rcpp::as<vec>(lambda_);
      mat X = Rcpp::as<mat>(X_);
      double mix = Rcpp::as<double>(mix_);
      int grad = Rcpp::as<int>(grad_);
      vec kum1 = Rcpp::as<vec>(kum1_);
      mat kum2 = Rcpp::as<mat>(kum2_);
      
      uint32_t n = X.n_rows;
      uint32_t d = X.n_cols;
       
      if( d != kum1.n_elem)       Rcpp::stop("X.n_cols != kum1.n_elem");
      if( d != lambda.n_elem )    Rcpp::stop("X.n_cols != lambda.n_elem");
      if( d != kum2.n_cols )      Rcpp::stop("X.n_cols != kum2.n_cols");
      if( d != kum2.n_rows )      Rcpp::stop("X.n_cols != kum2.n_rows");
      
      vec lx = X * lambda;
      
      double alpha = max( lx );       // constant for preventing overflow
      vec elx = exp( lx - alpha );    // exp(lambda'x_i - alpha) vector
       
      double sum_elx = sum( elx );
      mat elx_by_X(n, d); 
          
      for(unsigned int ii = 0; ii < d; ii++) elx_by_X.col(ii) = X.col(ii) % elx;
      
      double tmp_K = log( sum_elx / n ) + alpha;
      double K =  mix * tmp_K + (1 - mix) * ( dot(kum1, lambda) + 0.5 * dot(lambda, kum2 * lambda) );
      
      // We want only the value of K
      if(grad == 0) return( Rcpp::List::create( Rcpp::Named("K") = K,
                                                Rcpp::Named("tmp_K") = tmp_K, 
                                                Rcpp::Named("elx") = elx) );
        
      vec tmp_dK = conv_to<vec>::from( sum( elx_by_X, 0 ) / sum_elx );
      vec dK = mix * tmp_dK + ( 1-mix ) * ( kum1 + kum2 * lambda );
        
      // We want K and dK
      if(grad == 1) return( Rcpp::List::create( Rcpp::Named("K") = K, Rcpp::Named("dK") = Rcpp::wrap( dK.memptr(), dK.memptr() + dK.n_elem ),
                                                Rcpp::Named("tmp_K") = tmp_K, Rcpp::Named("tmp_dK") = tmp_dK, Rcpp::Named("elx") = elx ) );
        
      mat tmp_d2K = X.t() * elx_by_X / sum_elx - tmp_dK * tmp_dK.t();
      mat d2K =  mix * tmp_d2K  + (1-mix) * kum2;
      
      // We wanto also k, dK and d2K      
      if(grad == 2){
        
        return( Rcpp::List::create( Rcpp::Named("K") = K, Rcpp::Named("dK") = Rcpp::wrap( dK.memptr(), dK.memptr() + dK.n_elem ), Rcpp::Named("d2K") = d2K,
                                    Rcpp::Named("tmp_K") = tmp_K, Rcpp::Named("tmp_dK") = Rcpp::wrap( tmp_dK.memptr(), tmp_dK.memptr() + tmp_dK.n_elem ), 
                                    Rcpp::Named("tmp_d2K") = tmp_d2K, Rcpp::Named("elx") = elx) );
        
      } {
        Rcpp::stop("\"grad\" should be an interger equal to 0, 1 or 2");
      }
      
    } catch( std::exception& __ex__){
      forward_exception_to_r(__ex__);
    } catch(...){
      ::Rf_error( "c++ exception (unknown reason)" );
    }
    return Rcpp::wrap(NA_REAL);
}
