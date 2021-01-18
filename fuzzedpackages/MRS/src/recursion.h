#ifndef RECURSION_H
#define RECURSION_H

using namespace Rcpp;
using namespace arma;
using namespace std;

class class_tree
{
  public:
  Mat<unsigned int> X; // the data
  arma::vec G; // the group labels of each observation
  arma::vec H; // the sub-group labels of each observation (used for dANOVA)
  int n_tot; // total number of observations
  int p;    // dimension of the sample space
  int n_states; // number of hidden states
  vec init_state; // initial state of the hidden process
  int n_groups;  // number of groups
  Col<int> n_subgroups; // number of subgrous within each group (used for dANOVA)
  Col<int> cum_subgroups; // cumulative sum of n_subgroups; (used for dANOVA)
  int K;    // maximum depth of the tree
  arma::vec nu_vec;  // discrete prior on the parameter nu
  double alpha;   // pseudo-counts
  double beta, gamma, delta, eta;  // parameters of the transition probability matrix
  bool return_global_null, return_tree;
  int min_n_node;   // Node in the tree is returned if there are more than min_n_node data-points in it.
  //constructor
  class_tree( Mat<unsigned int> X, 
              arma::vec G,
              arma::vec H,
              vec init_state, 
              int n_groups, 
              Col<int> n_subgroups,
              int K, 
              arma::vec nu_vec,
              double alpha = 0.5,
              double beta = 1.0,
              double gamma = 0.3,
              double delta = 0.3,
              double eta = 0.3,
              bool return_global_null = true,
              bool return_tree = true,
              int min_n_node = 0
            );
            
  // compute posterior recursively          
  void update();      
  // compute PMAPs recursively
  void compute_varphi_post();
  // extract the MAP tree
  void representative_tree(); 


  
  // getters
  vector< Col< unsigned > > get_data_points_nodes();
  vector<unsigned short> get_level_nodes();
  vector<double> get_alt_prob_nodes();
  vector<vec> get_effect_size_nodes();
  vector<int> get_direction_nodes();
  vector< vector<double> > get_sides_nodes(vec a, vec b);
  vector<unsigned short> get_idx_nodes();
  double get_posterior_global_null();
  double get_prior_global_null(); 
  double get_marginal_loglikelihood(); 



  
  // clear memory 
  void clear();
            
  private:
  // observations in a tree structure
  int **data;
  // node info in a tree structure
  double **xi_post;
  double **chi;
  double **upsilon;
  double **psi_post;
  double **lambda_post;
  double **varphi_post;
  double **nu_post;
  double **theta_post;
  int **map;
  unsigned long  *modelscount;
  
  // save the representative tree and nodes information
  result_cubes_type result_cubes;
  
  // for each state compute the most likely children's states and cutting direction
  void compute_map(INDEX_TYPE& I, int level, arma::mat lambda_post);
  
  // compute posterior cutting direction probabilities
  arma::mat compute_lambda_post(  INDEX_TYPE& I, 
                                  int level, 
                                  arma::vec log_lambda, 
                                  arma::mat kappa, 
                                  arma::vec chi  );

  // compute posterior psi 
  double compute_post_psi(  INDEX_TYPE& I, 
                            int level, 
                            double null_prob, 
                            double prune_prob, 
                            arma::vec lambda_post );

  // compute posterior transition probabilities at individual node given chi
  arma::mat compute_post_xi(arma::vec chi, int level);

  // compute chi for each node given kappa 
  arma::vec compute_chi(arma::mat kappa, arma::vec log_lambda);

  // compute kappa at individual node for each state and cutting direction
  arma::mat compute_kappa(INDEX_TYPE& I, int level);

  // compute marginal likelihood at individual node for each state
  arma::vec compute_m(INDEX_TYPE& I, int level, int d);
  
  // compute marginal likelihood for each value of nu, and the theta hat for each group under the alternative
  Rcpp::List compute_m_anova(INDEX_TYPE& I, int level, int d);


  // return prior transition probability matrix
  arma::mat prior_transition_matrix(int level);
  // return prior transition probability
  double prior_transition(int s, int t, int level);
  
  
  void representative_subtree(  INDEX_TYPE& I, 
                                int level, 
                                unsigned short node_index, 
                                Mat<unsigned int> X_binary,
                                Col<unsigned int> data_indices,
                                Col<unsigned int> cut_counts,
                                uword state_star ); 
                                
  //compute effect size for mrs model
  Rcpp::List mrs_effect_size(INDEX_TYPE& I, int level, int top_direction);
  // compute effect size for anova model
  Rcpp::List anova_effect_size(INDEX_TYPE& I, int level, int top_direction);  

                                
  void save_index(  INDEX_TYPE& I, 
                    int level,                
                    double alt_prob, 
                    vec effect_size, 
                    int direction,
                    Col<unsigned int> data_indices,
                    unsigned short node_index );                               

  
  // organize observations in the tree
  void init();
  // subfunction to organize observations in the tree
  void add_data_to_subtree( INDEX_TYPE I, 
                            int level, 
                            int x_curr, 
                            int part_count, 
                            Col<unsigned int> obs,
                            unsigned int group);
                            
                            
  // functions to retreive information about a node (or a child of a node)
  double * get_node_nu_post(INDEX_TYPE& I, int level);
  double * get_node_theta_post(INDEX_TYPE& I, int level);
  double * get_node_xi_post(INDEX_TYPE& I, int level);
  double * get_node_varphi_post(INDEX_TYPE& I, int level);
  double * get_node_psi_post(INDEX_TYPE& I, int level);
  double * get_node_lambda_post(INDEX_TYPE& I, int level);
  double * get_node_chi(INDEX_TYPE& I, int level);
  double * get_node_upsilon(INDEX_TYPE& I, int level);
  int * get_node_data(INDEX_TYPE& I, int level);
  int * get_node_map(INDEX_TYPE& I, int level);
  
  double * get_child_xi_post(INDEX_TYPE& I, int i, int level, unsigned short which);
  double * get_child_varphi_post(INDEX_TYPE& I, int i, int level, unsigned short which);
  double * get_child_psi_post(INDEX_TYPE& I, int i, int level, unsigned short which);
  double * get_child_lambda_post(INDEX_TYPE& I, int i, int level, unsigned short which);
  double * get_child_chi(INDEX_TYPE& I, int i, int level, unsigned short which);
  double * get_child_upsilon(INDEX_TYPE& I, int i, int level, unsigned short which);
  int * get_child_data(INDEX_TYPE& I, int i, int level, unsigned short which);               
  int * get_child_map(INDEX_TYPE& I, int i, int level, unsigned short which);               
  

  
};

#endif
