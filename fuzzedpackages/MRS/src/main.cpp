#include "RcppArmadillo.h"
#include "helpers.h"
#include "recursion.h"

using namespace Rcpp;
using namespace arma;
using namespace std;


// [[Rcpp::export]]
Rcpp::List fitMRScpp( arma::mat X, 
                      arma::vec G, 
                      int n_groups, 
                      arma::vec init_state,
                      arma::mat Omega,  
                      int K = 5, 
                      double alpha = 0.5, 
                      double beta = 1.0, 
                      double gamma = 0.3,
                      double delta = 0.3,
                      double eta = 0.3,
                      bool return_global_null = true,
                      bool return_tree = true,
                      int min_n_node = 0
                  )
{

  /* **************************************** */
  /* center the data and transform in binary  */
  /* **************************************** */
  
  // number of dimensions
  int p = X.n_cols;
  // total number of observations
  int n_tot = X.n_rows;
  // total number of observations in each subgroup (used in dANOVA)
  Col<int> n_subgroups(n_groups); n_subgroups.fill(1);     
  // sub-gorups label for each observation (used in dANOVA)
  vec H(n_tot); H.fill(1);
  // prior on parameter nu (used in dANOVA)
  arma::vec nu_vec(1); nu_vec << 1 ;   
  vec a = 1.0 / (Omega.col(1) - Omega.col(0));
  vec b = - Omega.col(0) % a;
  // Initialize matrix of observations normalized to the p-dimensional hypercube
  Mat<unsigned int> X_binary(n_tot,p);    
  // Map each observation to the p-dimensional hypercube
  // Transform each observation in p-dimensional binary vector (up to resolution K+1)
  for(int i = 0; i < n_tot; i++ )
  {
    for(int d = 0; d < p; d++ )
      X_binary(i,d) = convert_to_inverse_base_2(a(d)*X(i,d)+b(d), K+1 );
  }
    
  /* ***************** */
  /* Compute posterior */
  /* ***************** */
  class_tree my_tree( X_binary, 
                      G, 
                      H,
                      init_state, 
                      n_groups, 
                      n_subgroups,
                      K, 
                      nu_vec,
                      alpha, 
                      beta, 
                      gamma, 
                      delta,
                      eta, 
                      return_global_null, 
                      return_tree,
                      min_n_node);
  my_tree.update();
  
  double loglike = my_tree.get_marginal_loglikelihood();

  double post_glob_null = NA_REAL;  
  double prior_glob_null = my_tree.get_prior_global_null();
  if(return_global_null == true)
    post_glob_null = my_tree.get_posterior_global_null();
    
  Rcpp::List data = Rcpp::List::create(
    Rcpp::Named( "X") = X,
    Rcpp::Named( "G") = G,
    Rcpp::Named( "Groups") = n_groups,
    Rcpp::Named( "Omega" ) = Omega
  );
    
  Rcpp::List other_info = Rcpp::List::create(
    Rcpp::Named("K") = K,
    Rcpp::Named("alpha") = alpha,
    Rcpp::Named("eta") = eta,
    Rcpp::Named("beta") = beta,
    Rcpp::Named("gamma") = gamma,
    Rcpp::Named("delta") = delta
  );
  
  if(return_tree == true)
  {
    my_tree.compute_varphi_post();  
    my_tree.representative_tree();     
    vector<unsigned short> levels = my_tree.get_level_nodes();
    vector<double> alt_probs = my_tree.get_alt_prob_nodes();
    vector<vec> effect_sizes = my_tree.get_effect_size_nodes(); 
    vector<int> directions = my_tree.get_direction_nodes();
    vector< vector<double> > sides = my_tree.get_sides_nodes(a, b);   
    vector< Col< unsigned > > data_points =  my_tree.get_data_points_nodes();
    vector<unsigned short> node_idx = my_tree.get_idx_nodes();
    
    my_tree.clear();
        
    Rcpp::List rep_tree = Rcpp::List::create( 
      Rcpp::Named( "Levels") = levels,
      Rcpp::Named( "AltProbs") = alt_probs,
      Rcpp::Named( "EffectSizes") = effect_sizes,
      Rcpp::Named( "Directions") = directions,
      Rcpp::Named( "Regions") = sides,
      Rcpp::Named( "DataPoints") = data_points,
      Rcpp::Named( "Ids") = node_idx 
    );
    
    return Rcpp::List::create(  
      Rcpp::Named( "PostGlobNull") = post_glob_null,
      Rcpp::Named( "PriorGlobNull") = prior_glob_null,
      Rcpp::Named( "LogLikelihood") = loglike,
      Rcpp::Named( "RepresentativeTree") = rep_tree,
      Rcpp::Named( "Data" ) = data,
      Rcpp::Named( "OtherInfo" ) = other_info
    ) ;   
  }
  else
  {
    my_tree.clear();
    Rcpp::List rep_tree = Rcpp::List::create();
        
    return Rcpp::List::create(  
      Rcpp::Named( "PostGlobNull") = post_glob_null,
      Rcpp::Named( "PriorGlobNull") = prior_glob_null,
      Rcpp::Named( "LogLikelihood") = loglike,
      Rcpp::Named( "RepresentativeTree") = rep_tree,
      Rcpp::Named( "Data" ) = data,
      Rcpp::Named( "OtherInfo" ) = other_info
    ) ;   

  }

    
}





// [[Rcpp::export]]
Rcpp::List fitMRSNESTEDcpp( arma::mat X, 
                            arma::vec G, 
                            arma::vec H, 
                            int n_groups, 
                            arma::Col<int> n_subgroups,
                            arma::vec init_state,
                            arma::mat Omega, 
                            arma::vec nu_vec,
                            int K = 5,                             
                            double alpha = 0.5, 
                            double beta = 1.0, 
                            double gamma = 0.07, 
                            double delta = 0.4,
                            double eta = 0,
                            bool return_global_null = true,
                            bool return_tree = true
                          )
{

  /* **************************************** */
  /* center the data and transform in binary  */
  /* **************************************** */
  
  // number of dimensions
  int p = X.n_cols;
  // total number of observations
  int n_tot = X.n_rows;
  vec a = 1.0 / (Omega.col(1) - Omega.col(0));
  vec b = - Omega.col(0) % a;
  // Initialize matrix of observations normalized to the p-dimensional hypercube
  Mat<unsigned int> X_binary(n_tot,p);    
  // Map each observation to the p-dimensional hypercube
  // Transform each observation in p-dimensional binary vector (up to resolution K+1)
  for(int i = 0; i < n_tot; i++ )
  {
    for(int d = 0; d < p; d++ )
      X_binary(i,d) = convert_to_inverse_base_2(a(d)*X(i,d)+b(d), K+1 );
  }
    
  /* ***************** */
  /* Compute posterior */
  /* ***************** */
  class_tree my_tree( X_binary, 
                      G, 
                      H,
                      init_state, 
                      n_groups, 
                      n_subgroups,
                      K, 
                      nu_vec,
                      alpha, 
                      beta, 
                      gamma, 
                      delta,
                      eta, 
                      return_global_null, 
                      return_tree);
  my_tree.update();
  
  double loglike = my_tree.get_marginal_loglikelihood();

  double post_glob_null = NA_REAL;  
  double prior_glob_null = my_tree.get_prior_global_null();
  if(return_global_null == true)
    post_glob_null = my_tree.get_posterior_global_null();
    
  Rcpp::List data = Rcpp::List::create(
    Rcpp::Named( "X") = X,
    Rcpp::Named( "G") = G,
    Rcpp::Named( "H") = H,
    Rcpp::Named( "Groups") = n_groups,
    Rcpp::Named( "Replicates") = n_subgroups,
    Rcpp::Named( "Omega" ) = Omega
  );
    
  Rcpp::List other_info = Rcpp::List::create(
    Rcpp::Named("K") = K,
    Rcpp::Named("alpha") = alpha,
    Rcpp::Named("eta") = eta,
    Rcpp::Named("beta") = beta,
    Rcpp::Named("gamma") = gamma,
    Rcpp::Named("delta") = delta,
    Rcpp::Named("nu_vec") = nu_vec
  );
  
  if(return_tree == true)
  {
    my_tree.compute_varphi_post();  
    my_tree.representative_tree();     
    vector<unsigned short> levels = my_tree.get_level_nodes();
    vector<double> alt_probs = my_tree.get_alt_prob_nodes();
    vector<vec> effect_sizes = my_tree.get_effect_size_nodes(); 
    vector<int> directions = my_tree.get_direction_nodes();
    vector< vector<double> > sides = my_tree.get_sides_nodes(a, b);   
    vector< Col< unsigned > > data_points =  my_tree.get_data_points_nodes();
    vector<unsigned short> node_idx = my_tree.get_idx_nodes();
    
    my_tree.clear();
        
    Rcpp::List rep_tree = Rcpp::List::create( 
      Rcpp::Named( "Levels") = levels,
      Rcpp::Named( "AltProbs") = alt_probs,
      Rcpp::Named( "EffectSizes") = effect_sizes,
      Rcpp::Named( "Directions") = directions,
      Rcpp::Named( "Regions") = sides,
      Rcpp::Named( "DataPoints") = data_points,
      Rcpp::Named( "Ids") = node_idx 
    );
    
    return Rcpp::List::create(  
      Rcpp::Named( "PostGlobNull") = post_glob_null,
      Rcpp::Named( "PriorGlobNull") = prior_glob_null,
      Rcpp::Named( "LogLikelihood") = loglike,
      Rcpp::Named( "RepresentativeTree") = rep_tree,
      Rcpp::Named( "Data" ) = data,
      Rcpp::Named( "OtherInfo" ) = other_info
    ) ;   
  }
  else
  {
    my_tree.clear();
    Rcpp::List rep_tree = Rcpp::List::create();
        
    return Rcpp::List::create(  
      Rcpp::Named( "PostGlobNull") = post_glob_null,
      Rcpp::Named( "PriorGlobNull") = prior_glob_null,
      Rcpp::Named( "LogLikelihood") = loglike,
      Rcpp::Named( "RepresentativeTree") = rep_tree,
      Rcpp::Named( "Data" ) = data,
      Rcpp::Named( "OtherInfo" ) = other_info
    ) ;   

  }

    
}

