#include "helpers.h"
#include "hmg.h"

using namespace Rcpp;
using namespace arma;






// [[Rcpp::export]]
Rcpp::List fitGrove(  arma::mat D, 
                      arma::mat X, 
                      arma::vec p,
                      arma::vec tau_par,
                      arma::vec eta_par,
                      arma::vec gamma_par,
                      arma::vec init_state,
                      double nu,
                      double sigma0,
                      double alpha,
                      double beta, 
                      int n_samp = 500,
                      int transition_mode = 1) // 1 is markov transition, 0 is independent transition, default to 1
{
  HMG H(  D, X, p, tau_par, 
          eta_par, gamma_par,  
          init_state, nu, sigma0, alpha, beta, transition_mode);
  
  std::vector<NumericMatrix> S = H.get_post_states();
  std::vector<NumericMatrix> Sprior = H.get_prior_states();
  
  arma::vec prior_null = H.get_prior_null(init_state);
  double ML = H.get_marginal_likelihood(); 

  std::vector<cube> Samples = H.post_sample_coeff(n_samp);
  cube SamplesCube = H.get_post_sample_coeff(Samples);
  vec post_null = H.get_posterior_null();
           
  List data = Rcpp::List::create(  
      Rcpp::Named( "design_matrix" ) = X,
      Rcpp::Named( "p" ) = p
    ) ; 
    
  List samples = Rcpp::List::create(  
      Rcpp::Named( "mean" ) = SamplesCube,
      Rcpp::Named( "S" ) = S,
      Rcpp::Named( "Sprior" ) = Sprior
    ) ; 
  
  return Rcpp::List::create(  
      Rcpp::Named( "samples" ) = samples,
      Rcpp::Named( "marginal_likelihood" ) = ML,
      Rcpp::Named( "prior_null" ) = prior_null,
      Rcpp::Named( "post_null" ) = post_null,
      Rcpp::Named( "data" ) = data
    ) ;    
};






// [[Rcpp::export]]
Rcpp::List fitGroveML(arma::mat D, 
                      arma::mat X, 
                      arma::vec p,
                      arma::vec tau_par,
                      arma::vec eta_par,
                      arma::vec gamma_par,
                      arma::vec init_state,
                      double nu,
                      double sigma0,
                      double alpha,
                      double beta,
                      int transition_mode = 1)
{
  HMG H(  D, X, p, tau_par, 
          eta_par, gamma_par,  
          init_state, nu, sigma0, alpha, beta, transition_mode);
  
  std::vector<NumericMatrix> S = H.get_post_states();
  arma::vec prior_null = H.get_prior_null(init_state);
  double ML = H.get_marginal_likelihood(); 

    
  return Rcpp::List::create(  
      Rcpp::Named( "marginal_likelihood" ) = ML
    ) ;    
};







