#ifndef HMG_H
#define HMG_H

#include "RcppArmadillo.h"

using namespace Rcpp;
using namespace arma;

const double log2pi = std::log(2.0 * M_PI);



class HMG {
  
  private:
    mat X, FullLambda;
    double a_0, b_0;
    int J, n_tot, tot_states, n_factors;
    vec p, tau_par, eta_par, gamma_par;
    double alpha, beta;
    std::vector<mat> w;   
    std::vector<mat> Phi, Psi, PostStates, PriorStates;
    std::vector<cube> PostTrans;
    std::vector<mat> PriorTrans;
    // std::vector<cube> PriorTrans;
    vec initial_state; 
    int transition_mode;
    
  public:  
    int n_sample;
  
  // constructor
  HMG(  mat D, 
        mat X, 
        vec dims, 
        vec taus_par, 
        vec eta_par, 
        vec gamma_par, 
        vec init_state,
        double nu, 
        double sigma0, 
        double alpha, 
        double beta,
        int transition_mode );
  
  // organize initial state of the product space
  vec init_init_state(vec input_init_state);
  
  // organize data on a tree
  std::vector<mat> init_data(mat D);
  
  // compute marginal likelihood of every node for every state
  std::vector<mat> init_marg();
  
  // vector of active factors from integer
  uvec active_columns(int s);  
  
  // extract submatrix according to active factors
  mat DesignMatrix(int s);
  
  // compute the full precision matrix
  mat FullPrecision();
  
  // extract submatrix for computing marginal likelihood 
  mat Lambda(int j, int s);
  
  // evaluate marginal likelihood at given node and given state
  double MargLike(int j, int k, int s);  
  
  // compute prior transition probability
  double prior_trans_elem(int j,  int s, int t);  
  
  // compute marginal likelihood recursively
  std::vector<mat> update_marg(); 
  
  // compute posterior transition probabilities  
  std::vector<cube> post_trans();
  
  // compute prior transition probability matrices
  std::vector<mat> prior_trans();
  //std::vector<cube> prior_trans();
 
  
  // compute posterior transition probability for a given node
  double post_trans_elem(int j, int k, int s, int t);
  
  // compute posterior state probabilities
  std::vector<mat> post_state();
  
  // compute prior state probabilities
  std::vector<mat> prior_state();
  
  // get marginal likelihood
  double get_marginal_likelihood();
  
  // get posterior null for each factor
  vec get_posterior_null();
  
  // get posterior null for each node and each factor
  double posterior_null(int j, int k, int s, int g);
  
  // get prior null probability for each factor
  vec get_prior_null(vec init_state);
  
  // get marginal posterior probability of hidden states
  std::vector<NumericMatrix> get_post_states();
  
  // get marginal prior probability of hidden states
  std::vector<NumericMatrix> get_prior_states();
  
  std::vector<NumericMatrix> Sample_States(int n_samples);
  
  std::vector<mat> Count_Sample_States(std::vector<NumericMatrix> StatesSample);
  
  mat post_sample_single_coeff(int j, int k, int s, int n);  
  
  std::vector<cube>  post_sample_coeff(int n_samp);
  
  cube get_post_sample_coeff(std::vector<cube> input);








  
  
};



#endif
