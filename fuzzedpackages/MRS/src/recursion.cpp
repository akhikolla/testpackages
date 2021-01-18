#include <RcppArmadillo.h>
#include "helpers.h"   
#include "recursion.h"  

using namespace Rcpp;
using namespace arma;
using namespace std;

class_tree::class_tree( Mat<unsigned int> X,
                        vec G,
                        vec H,
                        vec init_state,
                        int n_groups,
                        Col<int> n_subgroups,
                        int K,
                        vec nu_vec,
                        double alpha,
                        double beta,
                        double gamma,
                        double delta,
                        double eta,
                        bool return_global_null,
                        bool return_tree,
                        int min_n_node):
                        X(X),
                        G(G),
                        H(H),
                        init_state(init_state),
                        n_groups(n_groups),
                        n_subgroups(n_subgroups),
                        K(K),
                        nu_vec(nu_vec),
                        alpha(alpha),
                        beta(beta),
                        gamma(gamma),
                        delta(delta),
                        eta(eta),
                        return_global_null(return_global_null),
                        return_tree(return_tree),
                        min_n_node(min_n_node)
{
  n_tot = X.n_rows;
  p = X.n_cols;
  n_states = init_state.n_elem;
  cum_subgroups.set_size(n_subgroups.n_elem + 1);
  cum_subgroups(0) = 0;
  cum_subgroups.subvec(1, n_subgroups.n_elem) = cumsum(n_subgroups);
  init();
}



void class_tree::update()
{
  double *CHI_CURR, *LAMBDA_CURR;
  double *XI_CURR, *PSI_CURR, *UPSILON_CURR;
  double *NU_CURR, *THETA_CURR;
  int *DATA_CURR;
  int num_data_points_node;
  INDEX_TYPE I;
  vec log_lambda(p); log_lambda.fill( (-1.0) * log((double)p) );  
  mat lambda_post(n_states, p);
  mat xi_post(n_states, n_states);
  mat kappa(n_states, p);
  vec chi(n_states);
  int it;
  
  for(int level = (K+1); level >= 0; level--)
  {
    unsigned count = 0;
    I = init_index(level);
    while(count < modelscount[level])
    {
      CHI_CURR = get_node_chi(I, level);  
      DATA_CURR = get_node_data(I, level);
      if( (level < K + 1) && (return_tree == true) )
      {
        XI_CURR = get_node_xi_post(I, level);
        LAMBDA_CURR = get_node_lambda_post(I, level);  
        
        if (sum(n_subgroups) != n_groups) {
          NU_CURR = get_node_nu_post(I, level);
          THETA_CURR = get_node_theta_post(I, level);
        }
      }
      if(return_global_null == true)
        PSI_CURR = get_node_psi_post(I, level);  

      if(return_tree == true)
        UPSILON_CURR = get_node_upsilon(I, level);
      
      for(int j=0; j<(int)pow2(level); j++)
      {
        I.var[MAXVAR] = j;
        num_data_points_node = sum_elem(DATA_CURR, sum(n_subgroups));
        
        // XI_CURR[0] : log P(null -> null | data )
        // XI_CURR[1] : log P(null -> alternative | data )
        // XI_CURR[2] : log P(null -> prune | data )
        // XI_CURR[3] : log P(alternative -> null | data )
        // XI_CURR[4] : log P(alternative -> alternative | data )
        // XI_CURR[5] : log P(alternative -> prune | data )
        // XI_CURR[6] : log P(prune -> null | data )
        // XI_CURR[7] : log P(prune -> alternative | data )
        // XI_CURR[8] : log P(prune -> prune | data )
        
        // CHI_CURR[0] : log marginal likelihood when state is null
        // CHI_CURR[1] : log marginal likelihood when state  is alternative
        // CHI_CURR[2] : log marginal likelihood when state is prune
        
        // s = hidden state, d = cutting direction
        // LAMBDA_CURR[s=0 d=0], LAMBDA_CURR[s=1 d=0], ..., LAMBDA_CURR[s=1 d=1], ... 
        
        // PSI_CURR[0] : global null probability when the parent's state is null
        
        // UPSILON_CURR[0] : MAP likelihood when state is null
        // UPSILON_CURR[1] : MAP likelihood when state is alternative
        // UPSILON_CURR[2] : MAP likelihood when state is prune
        
        // MAP_CURR[0] : most likely left-child's state if parent's state is null
        // MAP_CURR[1] : most likely right-child's state if parent's state is null
        // MAP_CURR[2] : most likely cutting direction if the state is null        
        // MAP_CURR[3] : most likely left-child's state if parent's state is alternative
        // ...
        // MAP_CURR[8] : most likely cutting direction if the state is prune
        
        // NU_CURR[0] : P( nu_1 | S = 1, d = 0 ), P( nu_2 | S = 1, d = 0 ), ... P( nu_1 | S = 1, d = 1 ), ...
        
        // THETA_CURR: \hat \theta for each group and for each nu.
        
        if( level == (K+1) )
        {
          for(int s = 0; s < n_states; s++)
          {
            CHI_CURR[s] = num_data_points_node*level*log(2.0);
            if(return_tree == true)
              UPSILON_CURR[s] = 0.0;
          }         
          
          if(return_global_null == true)
            PSI_CURR[0] = 0;  
        } else if( (num_data_points_node == 0) || (num_data_points_node == 1) )
        {
          for(int s = 0; s < n_states; s++)
            CHI_CURR[s] = num_data_points_node*level*log(2.0);  

          if(return_global_null == true)
          {
            PSI_CURR[0] = compute_post_psi( I, 
                                            level, 
                                            prior_transition(0, 0, level),
                                            prior_transition(0, 2, level),
                                            log_lambda );            
          }

          if(return_tree == true)
          {
            // save xi
            it = 0;
            for(int s = 0; s < n_states; s++)
            {
              for(int t = 0; t < n_states; t++)
              {
                  XI_CURR[it] = prior_transition(s, t, level);
                  it++;   
              }
            }
            // save lambda
            it = 0;
            for(int s = 0; s < n_states; s++)
            {
              for(int d = 0; d < p; d++)
              {
                LAMBDA_CURR[it] = log_lambda(d);
                it++;
              }              
            }  
            // compute the MAP quantities
            // NOTE: MAP_CURR and UPSILON_CURR are updated within the function          
            for(int s=0; s < n_states; s++)
              lambda_post.row(s) = log_lambda.t();  
            compute_map(I, level, lambda_post);   
            
            // save thetas and nus if anova model
            if (sum(n_subgroups) != n_groups) {
              int n_nu_grid = nu_vec.n_elem;
              it = 0;
              for (int d = 0; d < p; d++) {
                for (int ii = 0; ii < n_nu_grid; ii++) {
                  NU_CURR[it] = -log(n_nu_grid);
                  it++;
                }
              }
              
              it = 0;
              for (int d = 0; d < p; d++) {
                for (int ii = 0; ii < n_nu_grid; ii++) {
                  for (int jj = 0; jj < n_groups; jj++) {
                    THETA_CURR[it] = 0.5;
                    it++;
                  }
                }
              }
            }
          }
                            
        } else
        {
          // NOTE: NU_CURR and THETA_CURR are update within the function (only for ANDOVA)
          kappa = compute_kappa(I, level);
          chi = compute_chi(kappa, log_lambda);
          for(int s = 0; s < n_states; s++)
            CHI_CURR[s] = chi(s);
          
          if(return_global_null == true || return_tree == true)
          {
            xi_post = compute_post_xi(chi, level);  
              
            lambda_post  = compute_lambda_post( I, 
                                                level, 
                                                log_lambda, 
                                                kappa, 
                                                chi  );
          }


          if(return_global_null == true)
          {
            // save global null
            PSI_CURR[0] = compute_post_psi( I, 
                                            level, 
                                            xi_post(0, 0),
                                            xi_post(0, 2),
                                            lambda_post.row(0).t() );            
          }

          if(return_tree == true)
          {
            // save xi
            it = 0;
            for(int s = 0; s < n_states; s++)
            {
              for(int t = 0; t < n_states; t++)
              {
                  XI_CURR[it] = xi_post(s, t);
                  it++;                
              }
            }   
            // save lambda            
            it = 0;
            for(int s = 0; s < n_states; s++)
            {
              for(int d = 0; d < p; d++)
              {
                LAMBDA_CURR[it] = lambda_post(s, d);
                it++;
              }              
            }                 
            // compute the MAP quantities
            // NOTE: MAP_CURR and UPSILON_CURR are updated within the function
            compute_map(I, level, lambda_post);
          }
        
        }
        
        CHI_CURR += n_states;
        DATA_CURR += sum(n_subgroups);
        if( return_global_null == true)
          PSI_CURR++;
        if( return_tree == true)
        {
          UPSILON_CURR += n_states;
          if(level < K + 1)
          {
            XI_CURR += n_states*n_states;
            LAMBDA_CURR += n_states*p;    
            
            if (sum(n_subgroups) != n_groups) {
              int n_nu_grid = nu_vec.n_elem;
              NU_CURR += p * n_nu_grid;
              THETA_CURR += p * n_nu_grid * n_groups;
            }

          }
        }                  
      }      
      I = get_next_node(I, p, level);
      count++;
    }    
  }  
}


void class_tree::compute_map(INDEX_TYPE& I, int level, arma::mat lambda_post)
{
  double *UPSILON_CURR = get_node_upsilon(I, level);
  int *MAP_CURR = get_node_map(I, level);
  int it;
  uword  row, col, slice;
  double max_val;
  if(level == K)
  {
    it = 0;
    vec like_vec(p);
    for(int s = 0; s < n_states; s++)  
    {           
      // for each cutting direction...
      like_vec = lambda_post.row(s).t();
      // ... and pick the "best" cutting direction for fixed s   
      max_val = like_vec.max(slice);
      UPSILON_CURR[s] = max_val;
      MAP_CURR[it + 2] = slice;
      it += 3;      
    }
  }
  else
  {
    double *XI_CHILD_0, *XI_CHILD_1;
    double *UPSILON_CHILD_0, *UPSILON_CHILD_1;
    // organize chi (the transition probability matrices) of the children in a cube  
    cube cube_xi_child_0(n_states, n_states, p);
    cube cube_xi_child_1(n_states, n_states, p);    
    for(int d = 0; d < p; d++)
    {
      it = 0;
      XI_CHILD_0 = get_child_xi_post(I,d,level,0);
      XI_CHILD_1 = get_child_xi_post(I,d,level,1);     
      for(int s = 0; s < n_states; s++)
      {
        for(int t = 0; t < n_states; t++)
        {
          cube_xi_child_0(s,t,d) = XI_CHILD_0[it];
          cube_xi_child_1(s,t,d) = XI_CHILD_1[it];
          it++;
        }
      }
    }
    
    // organize Upsilon ( the MAP likelihood vector) of the children in a matrix
    mat mat_upsilon_child_0(n_states, p);
    mat mat_upsilon_child_1(n_states, p);
    
    for(int d = 0; d < p; d++)
    {
      it = 0;
      UPSILON_CHILD_0 = get_child_upsilon(I,d,level,0);
      UPSILON_CHILD_1 = get_child_upsilon(I,d,level,1);     
      for(int s = 0; s < n_states; s++)
      {
        mat_upsilon_child_0(s,d) = UPSILON_CHILD_0[it];
        mat_upsilon_child_1(s,d) = UPSILON_CHILD_1[it];
        it++;
      }
    }
    
    it = 0;
    // for each current state
    cube like_cube(n_states, n_states, p);
    uword  row, col, slice;
    double max_val;
    for(int s = 0; s < n_states; s++)  
    {           
      // for each cutting direction...
      for(int d = 0; d < p; d++)
      {
        // for each left child state (child_0)
        for(int ss = 0; ss < n_states; ss++)
        {
          // for each right child state (child_1)
          for(int tt = 0; tt < n_states; tt++)
          {
            // ... compute the likelihood                  
            like_cube(ss, tt, d) = lambda_post(s, d) + cube_xi_child_0(s,ss,d) + cube_xi_child_1(s,tt,d) 
              + mat_upsilon_child_0(ss,d) + mat_upsilon_child_1(tt,d);
          }
        } 
      }
      // ... and pick the "best" (d, ss, tt) combination for fixed s   
      max_val = like_cube.max(row, col, slice);
      UPSILON_CURR[s] = max_val;
      MAP_CURR[it] = row;
      MAP_CURR[it + 1] = col;
      MAP_CURR[it + 2] = slice;
      it += 3;
    }    
  }
}

void class_tree::compute_varphi_post()
{
  double *VARPHI_CURR, *XI_CURR;
  double *VARPHI_PARENT, *LAMBDA_PARENT;
  INDEX_TYPE I;  
  bool first_time;  
  vec temp(n_states);
  int level = 0;
  I = init_index(level);
  mat xi_post(n_states, n_states);
  vec parent_varphi_post(n_states);
  mat parent_lambda_post(n_states, p);
  int it=0;
  VARPHI_CURR = get_node_varphi_post(I, level);
  XI_CURR = get_node_xi_post(I, level);
  
  for(int s = 0; s < n_states; s++)
    VARPHI_CURR[s] = XI_CURR[s];
    
  for(level=1; level<=K; level++)
  {
    unsigned count = 0;
    I = init_index(level);    
    while(count < modelscount[level])
    {            
      XI_CURR = get_node_xi_post(I, level);
      VARPHI_CURR = get_node_varphi_post(I, level);
      for(int j=0; j<(int)pow2(level); j++)
      {
        I.var[MAXVAR] = j;
        // retrieve transition probability matrix
        it = 0;
        for(int s = 0; s < n_states; s++)
        {
          for(int t = 0; t < n_states; t++)
          {
            xi_post(s,t) = XI_CURR[it];
            it++;
          }  
        }
    
        first_time = true;
        for(int i = 0; i < p; i++)
        {
          for(int which = 0; which <= 1; which++)
          {
            std::pair<bool, INDEX_TYPE>  Ip = make_parent_index(I, i, level, which);
            if( Ip.first == true )
            {
              VARPHI_PARENT = get_node_varphi_post(Ip.second, level - 1);             
              LAMBDA_PARENT = get_node_lambda_post(Ip.second, level - 1);           
                            
              it = 0;
              for(int s = 0; s < n_states; s++)
              {
                parent_varphi_post(s) = VARPHI_PARENT[s]; 
                for(int d = 0; d < p; d++)
                {
                  parent_lambda_post(s,d) = LAMBDA_PARENT[it];
                  it++;
                }
              }                
              
              temp.fill(log(0.0));
              for(int s = 0; s < n_states; s++)
              {
                for(int t = 0; t < n_states; t++)
                  temp(s) = log_exp_x_plus_exp_y( temp(s), 
                    parent_varphi_post(t) + parent_lambda_post(t, i) + xi_post(t, s) );
              }
                
              if(first_time == true)
              {
                first_time = false;
                for(int s = 0; s < n_states; s++)
                  VARPHI_CURR[s] = temp(s);                                   
              }
              else
              {
                for(int s = 0; s < n_states; s++)
                  VARPHI_CURR[s] = log_exp_x_plus_exp_y( VARPHI_CURR[s], temp(s)); 
              }                
            }              
          }
        }
        XI_CURR += (n_states*n_states);
        VARPHI_CURR += n_states; 
      }      
      I = get_next_node(I,p,level);
      count++;
    }    
  }    
}

arma::mat class_tree::compute_lambda_post(  INDEX_TYPE& I, 
                                            int level, 
                                            arma::vec log_lambda, 
                                            arma::mat kappa, 
                                            arma::vec chi  )
{
  mat output(n_states, p);
  for(int s = 0; s < n_states; s++)
  {
    for(int d = 0; d < p; d++)
      output(s, d) = log_lambda(d) + kappa(s,d) - chi(s);
  }
  return output;
}


double class_tree::compute_post_psi(  INDEX_TYPE& I, 
                                      int level, 
                                      double null_prob, 
                                      double prune_prob, 
                                      arma::vec lambda_post )
{
  double *PSI_CHILD_0, *PSI_CHILD_1;
  double temp = log(0.0);
  for(int d = 0; d < p; d++)
  {
    PSI_CHILD_0 = get_child_psi_post(I, d, level, 0);
    PSI_CHILD_1 = get_child_psi_post(I, d, level, 1);
    temp = log_exp_x_plus_exp_y( temp, lambda_post(d) + PSI_CHILD_0[0] + PSI_CHILD_1[0] );
  }              
  return log_exp_x_plus_exp_y( prune_prob, null_prob + temp ) ;
}

arma::mat class_tree::compute_post_xi(arma::vec chi, int level)
{
  mat output(n_states, n_states);
  double temp;
  mat xi_prior = prior_transition_matrix(level);
  for(int s = 0; s < n_states; s++)
  {
    temp = xi_prior(s, 0) + chi(0);
    for(int t = 1; t < n_states; t++)
      temp = log_exp_x_plus_exp_y(temp, xi_prior(s, t) + chi(t) );      
    output.row(s) = xi_prior.row(s) + chi.t() -  temp;
  }
  return output;
}

arma::vec class_tree::compute_chi(arma::mat kappa, arma::vec log_lambda)
{
  vec output(n_states);  
  for(int s = 0; s < n_states; s++)
  {
    output(s) = log_lambda(0) + kappa(s, 0);
    for(int d = 1; d < p; d++)
      output(s) = log_exp_x_plus_exp_y( output(s), log_lambda(d) + kappa(s, d) );
  }
  return output;
}

arma::mat class_tree::compute_kappa(INDEX_TYPE& I, int level)
{
  mat output(n_states, p);
  double *CHI_CHILD_0, *CHI_CHILD_1;
  mat xi_prior = prior_transition_matrix(level+1);
  vec chi_vec_child_0(n_states), chi_vec_child_1(n_states), m(n_states); 
  int n_nu_grid = nu_vec.n_elem;
  vec temp_nu_vec(n_nu_grid + 1);
  int it;
  double temp_0, temp_1;
  
  for(int d = 0; d < p; d++)
  {
    it = 0;
    CHI_CHILD_0 = get_child_chi(I,d,level,0);
    CHI_CHILD_1 = get_child_chi(I,d,level,1);
    if( sum(n_subgroups) == n_groups)
      m = compute_m(I, level, d);
    else {
      List list_m_thetas = compute_m_anova(I, level, d);
      temp_nu_vec = Rcpp::as<vec>(list_m_thetas["m_nus"]);
      mat thetas = Rcpp::as<mat>(list_m_thetas["thetas"]);
      m(0) = temp_nu_vec(0);
      m(2) = temp_nu_vec(0);
      m(1) = log(0.0);
      for (int ii = 0; ii < n_nu_grid; ii++) {
        m(1) = log_exp_x_plus_exp_y(m(1), temp_nu_vec(ii + 1));
      }
      
      if (return_tree == true) {
        double *NU_POST = get_node_nu_post(I, level);
        double *THETA_POST = get_node_theta_post(I, level);
        for (int ii = 0; ii < n_nu_grid; ii++) {
          NU_POST[d * n_nu_grid + ii] = temp_nu_vec(ii + 1) - m(1);
          for (int jj = 0; jj < n_groups; jj++) {
            THETA_POST[d * n_nu_grid * n_groups + ii * n_groups + jj] = thetas(jj, ii);
          }
        }
      }

    }
    for(int s = 0; s < n_states; s++)
    {
      chi_vec_child_0(s) = CHI_CHILD_0[it];
      chi_vec_child_1(s) = CHI_CHILD_1[it];
      it++;
    }
    
    for(int s = 0; s < n_states; s++)
    {
      temp_0 = xi_prior(s, 0) + chi_vec_child_0(0);
      temp_1 = xi_prior(s, 0) + chi_vec_child_1(0);
      for(int t = 1; t < n_states; t++)
      {
        temp_0 = log_exp_x_plus_exp_y( temp_0, xi_prior(s, t) + chi_vec_child_0(t) );
        temp_1 = log_exp_x_plus_exp_y( temp_1, xi_prior(s, t) + chi_vec_child_1(t) );
      }
      output(s, d) = m(s) + temp_0 + temp_1; 
    }    
  }
  
  return output;
}

arma::vec class_tree::compute_m(INDEX_TYPE& I, int level, int d)
{
  vec output(n_states);
  int *DATA_CHILD_0, *DATA_CHILD_1;
  int all_data_0, all_data_1;

  DATA_CHILD_0 = get_child_data(I,d,level,0);
  DATA_CHILD_1 = get_child_data(I,d,level,1);
    
  all_data_0 = sum_elem(DATA_CHILD_0, n_groups);
  all_data_1 = sum_elem(DATA_CHILD_1, n_groups);
  
  output(0) = lgamma(all_data_0 + alpha) + lgamma(all_data_1 + alpha) 
        - lgamma(all_data_0 + all_data_1  + 2 * alpha) 
        - ( 2 * lgamma(alpha) - lgamma(2*alpha) );
        
  output(1) = - n_groups*( 2 * lgamma(alpha) - lgamma(2*alpha) );
  for(int i = 0; i < n_groups; i++)
    output(1) += lgamma(DATA_CHILD_0[i]+alpha) + lgamma(DATA_CHILD_1[i] + alpha) 
        - lgamma(DATA_CHILD_0[i] + DATA_CHILD_1[i] + 2.0*alpha);
 
  output(2) = output(0);    
  return output;

}


Rcpp::List class_tree::compute_m_anova(INDEX_TYPE& I, int level, int d)
{
  int *DATA_CHILD_0, *DATA_CHILD_1;
  DATA_CHILD_0 = get_child_data(I,d,level,0);
  DATA_CHILD_1 = get_child_data(I,d,level,1);
    
  int n_grid = nu_vec.n_elem;
  vec m_nus(n_grid + 1); m_nus.fill(log(0.0));
  mat thetas(n_groups, n_grid); thetas.fill(0.5);
  int n_grid_theta = 4;
  vec theta0(n_grid_theta);
  theta0 << 0.125 << 0.375 << 0.625 << 0.875 ;

  vec data_0( sum(n_subgroups) ); 
  vec data_1( sum(n_subgroups) );
  for(int i = 0; i < sum(n_subgroups); i++)
  {
    data_0(i) = DATA_CHILD_0[i]; 
    data_1(i) = DATA_CHILD_1[i];
  }  
  
  //under the null
  if(sum(data_0) == 0 || sum(data_1) == 0)
  {
    for(int g = 0; g < n_grid; g++)
    {
      for(int h = 0; h < n_grid_theta; h++)
        m_nus(0) = log_exp_x_plus_exp_y( m_nus(0), 
          eval_h(theta0(h), data_0, data_1, nu_vec(g), alpha ) - log(n_grid) - log(n_grid_theta) ); 
    }
  }
  else
  {
    for(int g = 0; g < n_grid; g++)
    {
      vec tt = newtonMethod(data_0, data_1, nu_vec(g), alpha);
      if(std::isnan(tt(1)))
      {
        for(int h = 0; h < n_grid_theta; h++)
          m_nus(0) = log_exp_x_plus_exp_y( m_nus(0), 
            eval_h(theta0(h), data_0, data_1, nu_vec(g), alpha ) - log(n_grid) - log(n_grid_theta) ); 
      }
      else
        m_nus(0) = log_exp_x_plus_exp_y( m_nus(0), tt(1)  - log(n_grid) );    
    }
      
    
  }
  //under the alternative
  mat temp_mat(n_groups, n_grid); temp_mat.fill(log(0.0));
  for(int j = 0; j < n_groups; j++)
  {
    vec temp_data_0 = data_0.subvec( cum_subgroups(j), cum_subgroups(j+1) - 1  );
    vec temp_data_1 = data_1.subvec( cum_subgroups(j), cum_subgroups(j+1) - 1  );
    if(sum(temp_data_0) == 0 &&  sum(temp_data_1) == 0)
    {
      for(int g = 0; g < n_grid; g++)
        temp_mat(j,g) = 0;
    }      
    else if(sum(temp_data_0) == 0 || sum(temp_data_1) == 0)
    {
      for(int g = 0; g < n_grid; g++)
      {
        for(int h = 0; h < n_grid_theta; h++)
          temp_mat(j,g) = log_exp_x_plus_exp_y( temp_mat(j,g), 
            eval_h(theta0(h), temp_data_0, temp_data_1, nu_vec(g), alpha ) - log(n_grid_theta) ); 
      }
    }
    else
    {
      for(int g = 0; g < n_grid; g++)
      {
        vec tt = newtonMethod(temp_data_0, temp_data_1, nu_vec(g), alpha);       
        if(std::isnan(tt(1)))
        {
          for(int h = 0; h < n_grid_theta; h++)
            temp_mat(j,g) = log_exp_x_plus_exp_y( temp_mat(j,g), 
              eval_h(theta0(h), temp_data_0, temp_data_1, nu_vec(g), alpha ) - log(n_grid_theta) );      
        }
        else {
          temp_mat(j,g) = tt(1);
          thetas(j,g) = tt(0);
        }
      }
         
      
    }
  }
  for(int g = 0; g < n_grid; g++)
    m_nus(g + 1) = sum(temp_mat.col(g)) - log(n_grid);
    
  return Rcpp::List::create(  
  Rcpp::Named( "m_nus" ) = m_nus,
  Rcpp::Named( "thetas" ) = thetas);  
}




arma::mat class_tree::prior_transition_matrix(int level)
{
  mat xi_prior(n_states, n_states);
  for(int s = 0; s < n_states; s++)
  {
    for(int t = 0; t < n_states; t++)
      xi_prior(s,t) = prior_transition(s, t, level);
  }
  return xi_prior;
}

double class_tree::prior_transition(int s, int t, int level)
{
  if(level == 0)
  {
    return log(init_state(t));
  }
  else
  {
    if( s == 1 )   // alternative
    {
        if( t == 1 )  // alternative
            return( log( 1.0 - eta) + log(delta) -level*log(beta)  );
            
        else if( t == 0 )   // null
            return( log(1.0 - eta) + log( 1.0 -  delta*pow( beta , -level ) ) ); 
            
        else
          return log(eta);        
   
    }
    else if( s == 0 )   // null
    {
        if( t == 1 )    // alternative
            return( log( 1.0 - eta) + log(gamma) -level*log(2.0)  );      
        else if( t == 0 )  // null
            return( log(1.0 - eta) + log( 1.0 -  gamma*pow( 2.0 , -level ) )  ); 
        else    
          return log(eta);
    }
    else
    {
        if( t == 2 )    // prune
          return log(1.0);
        else
          return log(0.0);
    }
  }
}


void class_tree::representative_tree() 
{   
  // function to find nested sequence of regions with probability 

  INDEX_TYPE I_root = init_index(0); 
  Col< unsigned int > data_indices(n_tot);
  Col< unsigned int > cut_counts(p); cut_counts.fill(0);
  for(int i=0; i < n_tot; i++)  
    data_indices(i) = i+1;
  representative_subtree(I_root, 0, 0, X, data_indices, cut_counts, 0);

}

void class_tree::representative_subtree(  INDEX_TYPE& I, 
                                          int level, 
                                          unsigned short node_index, 
                                          Mat<unsigned int> X_binary,
                                          Col<unsigned int> data_indices,
                                          Col<unsigned int> cut_counts,
                                          uword state_star
                                        ) 
{
    if ( level <= K ) 
    {        
        int state_0, state_1, top_direction;        
        int *MAP_CURR = get_node_map(I, level);    
        INDEX_TYPE J_0, J_1;
        unsigned short new_node_index_0, new_node_index_1;      
    
        if(level == 0)
        {       
          double *VARPHI_CURR = get_node_varphi_post(I, level);
          vec loglike(n_states);
          for(int s = 0; s < n_states; s++)
            loglike(s) = VARPHI_CURR[s];
          double max_val = loglike.max(state_star);          
        }
        
        int *DATA_CURR = get_node_data(I, level);
        int  num_data_points_node = sum_elem(DATA_CURR, sum(n_subgroups));
        
        // find the most likely state for the children states
        state_0 = MAP_CURR[state_star*n_states];
        state_1 = MAP_CURR[state_star*n_states+1];        
        // find most likely cutting direction
        top_direction = MAP_CURR[state_star*n_states+2];
        
        vec effect_size(n_groups);
        double alt_state_prob;
        
        if (sum(n_subgroups) != n_groups) {
          List temp_list = anova_effect_size(I, level, top_direction);  
          effect_size =  Rcpp::as<vec>(temp_list["effect_size"]);
          alt_state_prob = Rcpp::as<double>(temp_list["alt_state_prob"]);
        } else {
          List temp_list = mrs_effect_size(I, level, top_direction);  
          effect_size =  Rcpp::as<vec>(temp_list["effect_size"]);
          alt_state_prob = Rcpp::as<double>(temp_list["alt_state_prob"]);          
        }

        bool flag = false;
        if( num_data_points_node > min_n_node )
        {
          flag = true;
        
          vector<unsigned int> left_indices, right_indices;
          for(int i=0; i<(int)X_binary.n_rows; i++)
          {
            if( ( X_binary(i,top_direction) >> cut_counts(top_direction))  & 1 )
              right_indices.push_back(i);
            else
              left_indices.push_back(i);
          }
          cut_counts(top_direction)++;
          
          Col<unsigned int> l_indices = conv_to< Col<unsigned int> >::from(left_indices);
          Col<unsigned int> r_indices = conv_to< Col<unsigned int> >::from(right_indices);
                    
          // move to the 2 children of the set    
          J_0 = make_child_index(I, top_direction, level, 0 );
    	    new_node_index_0 =  node_index << 1 ;
    	    representative_subtree(J_0, level+1, new_node_index_0, 
            X_binary.rows(l_indices), data_indices(l_indices), cut_counts, state_0); 
	
	        J_1 = make_child_index(I, top_direction, level, 1 );
    	    new_node_index_1 =  ( node_index << 1 ) | 1 ;	
		      representative_subtree(J_1, level+1, new_node_index_1, 
            X_binary.rows(r_indices), data_indices(r_indices), cut_counts, state_1); 
        }
        
        if(flag)
        {
          save_index(I, level, alt_state_prob, effect_size, 
            (top_direction+1), data_indices, node_index);  
        } 
    }   
}

Rcpp::List class_tree::mrs_effect_size(INDEX_TYPE& I, int level, int top_direction) {  
  vec effect_size(n_groups); 
  
  int *DATA_CHILD_0, *DATA_CHILD_1;
  DATA_CHILD_0 = get_child_data(I,top_direction,level,0);
  DATA_CHILD_1 = get_child_data(I,top_direction,level,1);
  double *VARPHI_CURR = get_node_varphi_post(I, level);
  
  int n_0, n_1;
  int n_sample = 1000;
  mat theta(n_sample, n_groups);
  int it = 0;
  for(int j = 0; j < n_groups; j++)
  {
    n_0 = 0;
    n_1 = 0;
    for(int i = 0; i < n_subgroups(j); i++)
    {
      n_0 += DATA_CHILD_0[it];
      n_1 += DATA_CHILD_1[it];
      it++;
    }
    
    vec temp = rbeta(n_sample, (double)n_0 + alpha, (double)n_1 + alpha );
    theta.col(j) = temp;
  }
        
  for(int j = 0; j < n_groups; j++)
    effect_size(j) = mean(log(theta.col(j)) - log(1.0 - theta.col(j)) 
      - log( (sum(theta,1) - theta.col(j) )/(n_groups - 1.0)  )
      + log( 1.0 - (sum(theta,1) - theta.col(j) )/(n_groups - 1.0)  ) ); 

  double den = log(0.0);
  for(int s = 0; s < n_states; s++)
    den = log_exp_x_plus_exp_y(den, VARPHI_CURR[s]);
  effect_size = effect_size * exp(VARPHI_CURR[1] - den);
  
  return Rcpp::List::create(  
    Rcpp::Named( "effect_size" ) = effect_size,
    Rcpp::Named( "alt_state_prob" ) = exp(VARPHI_CURR[1] - den)) ;    
}


Rcpp::List class_tree::anova_effect_size(INDEX_TYPE& I, int level, int top_direction) {  
  
  vec effect_size(n_groups); effect_size.fill(0);

  double *VARPHI_CURR = get_node_varphi_post(I, level);
  double *NU_CURR = get_node_nu_post(I, level);
  double *THETA_CURR = get_node_theta_post(I, level);
  int n_nu_grid = nu_vec.n_elem;

  mat thetas(n_groups, n_nu_grid);
  vec nus(n_nu_grid);


  int it = 0;
  for (int i = 0; i < n_nu_grid; i++) {
    nus(i) = exp(NU_CURR[i]);
    for(int j = 0; j < n_groups; j++) {
      thetas(j,i) = THETA_CURR[it];
      it++;
    }
  }    
  
  for (int i = 0; i < n_nu_grid; i++) {
    for (int j = 0; j < n_groups; j++) {
      effect_size(j) += nus(i) * (log(thetas(j, i)) - log(1.0 - thetas(j, i)) 
        - log((sum(thetas.col(i)) - thetas(j, i))/(n_groups - 1.0))
        + log(1.0 - (sum(thetas.col(i)) - thetas(j, i))/(n_groups - 1.0)));
    }
  }
  
  double den = log(0.0);
  for(int s = 0; s < n_states; s++)
    den = log_exp_x_plus_exp_y(den, VARPHI_CURR[s]);
  effect_size = effect_size * exp(VARPHI_CURR[1] - den);
  
  return Rcpp::List::create(  
    Rcpp::Named( "effect_size" ) = effect_size,
    Rcpp::Named( "alt_state_prob" ) = exp(VARPHI_CURR[1] - den)) ;    
}


void class_tree::save_index(  INDEX_TYPE& I, 
                              int level,                
                              double alt_prob, 
                              vec effect_size, 
                              int direction,
                              Col<unsigned int> data_indices,
                              unsigned short node_index) 
{
  unsigned short x_curr = 0;
  unsigned short index_prev_var = 0;
  unsigned short lower = 0;
  short x_curr_count = -1;
  cube_type new_cube;
  vector<side_type> new_sides;
  int i;  
  for ( i = 0; i < level; i++) 
  {
    if ( I.var[i] - index_prev_var - 1 > 0 ) 
    { // next variable
      side_type new_side;
      new_side.var = x_curr; 
      new_side.extremes[0] = lower; 
      new_side.extremes[1] = lower + ((unsigned int) 1 << (K-x_curr_count - 1)) - 1; 
      new_sides.push_back(new_side);
      lower = 0;
      x_curr_count = 0;
    }    
    else 
      x_curr_count++;
        
    x_curr += I.var[i] - index_prev_var - 1;
    lower |= (((I.var[MAXVAR] >> i) & (unsigned int) 1)) << (K-x_curr_count - 1);  
    index_prev_var = I.var[i];
  }
  
  if (level > 0) 
  {
    side_type new_side;
    new_side.var = x_curr; 
    new_side.extremes[0] = lower; 
    new_side.extremes[1] =lower + ((unsigned int) 1 << (K-x_curr_count - 1)) - 1; 
    new_sides.push_back(new_side);
  }
  
  new_cube.sides = new_sides;
  new_cube.alt_prob = alt_prob;
  new_cube.level = level;
  new_cube.effect_size = effect_size;
  new_cube.direction = direction;
  new_cube.node_idx = node_index;
  new_cube.data_points = data_indices;  
  result_cubes.push_back(new_cube);  
}




vector< Col< unsigned > > class_tree::get_data_points_nodes()
{
  result_cubes_type::iterator it;
  vector< Col< unsigned > > v ;
  for(it=result_cubes.begin(); it<result_cubes.end(); it++)
    v.push_back(it->data_points);
  return v;
}

vector<unsigned short> class_tree::get_level_nodes()
{
  result_cubes_type::iterator it;
  vector<unsigned short> v ;
  for(it=result_cubes.begin(); it<result_cubes.end(); it++)
    v.push_back(it->level);
  return v;
}

vector<unsigned short> class_tree::get_idx_nodes()
{
  result_cubes_type::iterator it;
  vector<unsigned short> v ;
  for(it=result_cubes.begin(); it<result_cubes.end(); it++)
    v.push_back(it->node_idx);
  return v;
}

vector<double> class_tree::get_alt_prob_nodes()
{
  result_cubes_type::iterator it;
  vector<double> v ;
	for(it=result_cubes.begin(); it<result_cubes.end(); it++)
    {
		v.push_back(it->alt_prob);
	}
	return v;
}


vector< vec > class_tree::get_effect_size_nodes()
{
  result_cubes_type::iterator it;
	vector< vec > v ;
	for(it=result_cubes.begin(); it<result_cubes.end(); it++)
    {
		v.push_back(it->effect_size);
	}
	return v;
}

vector<int> class_tree::get_direction_nodes()
{
	result_cubes_type::iterator it;
	vector<int> v ;
	for(it=result_cubes.begin(); it<result_cubes.end(); it++)
    {
		v.push_back(it->direction);
	}
	return v;
}



vector< vector<double> > class_tree::get_sides_nodes(vec a, vec b)
  {
  result_cubes_type::iterator it;
  vector<side_type>::iterator km;
  vector<vector<double> > v;
  unsigned short actual_var = 0;
  
  for(it=result_cubes.begin(); it < result_cubes.end(); it++)
  {
    vector<double> w;
    actual_var = 0;	
  
    for(km = it->sides.begin(); km < it->sides.end()  ; km++)
    {
      while(km->var > actual_var)
      {      // add uncutted variables
        w.push_back( -b(actual_var) / a(actual_var) );
        w.push_back( (1.0 - b(actual_var)) / a(actual_var) );
        actual_var ++;
      }			       
      w.push_back( ( km->extremes[0]/((double)pow2(K)) - b(actual_var))/(a(actual_var)) );
      w.push_back( ( (1.0+km->extremes[1])/((double)pow2(K)) - b(actual_var))/(a(actual_var)) );	
      actual_var++;	
    }
    while( actual_var < p )
    {
      w.push_back( -b(actual_var) / a(actual_var) );
      w.push_back( (1.0 - b(actual_var)) / a(actual_var) );
      actual_var ++;
    }
  // add final uncutted variables
    v.push_back(w);
  }
  return v;
}


double class_tree::get_posterior_global_null() 
{ 
  return exp( psi_post[0][0] );
} 

// Compute the marginal prior probability of the null hypothesis 
// that the distributions are equal
double class_tree::get_prior_global_null() 
{ 
  double log_Psi = 0;
  for (int level = K; level>0; level--) 
    log_Psi = log_exp_x_plus_exp_y ( prior_transition(0,2,level), prior_transition(0,0,level)  + 2.0*log_Psi ) ;

  log_Psi = log_exp_x_plus_exp_y ( log(init_state(2)), log(init_state(0))  + 2.0*log_Psi ) ;

  return exp(log_Psi);    
}


double class_tree::get_marginal_loglikelihood() 
{ 
  double temp = log(0.0);
  for(int s = 0; s < n_states; s++)
  {
    temp = log_exp_x_plus_exp_y(  temp, init_state(s) + chi[0][s] );
  }
  return temp;
}   



void class_tree::init()
{
  unsigned long  j, l;
  data = new int*[K+2];
  modelscount = new unsigned long int[K+2];
  chi = new double*[K+2];
  
  if(return_global_null == true)
    psi_post = new double*[K+2];
  
  if(return_tree == true)
  {
    xi_post = new double*[K+1];
    varphi_post = new double*[K+1];
    lambda_post = new double*[K+1];  
    upsilon = new double*[K+2];
    if (sum(n_subgroups) != n_groups) {
      nu_post = new double*[K+1];
      theta_post = new double*[K+1];
    }
    map = new int*[K+1];
  }
  
  for(int i = 0; i <= K + 1; i++ )
  {
    modelscount[i] = Choose(p + i - 1, i);    
    data[i] = new int[(unsigned long )(modelscount[i]*sum(n_subgroups)) << i ];
    if( (i <= K) && (return_tree == true) )
    {
      xi_post[i] = new double[(unsigned long )(modelscount[i]*n_states*n_states) << i ];
      varphi_post[i] = new double[(unsigned long )(modelscount[i]*n_states) << i ];
      lambda_post[i] = new double[(unsigned long )(modelscount[i]*n_states*p) << i ];
      if (sum(n_subgroups) != n_groups) {
        nu_post[i] = new double[(unsigned long )(modelscount[i]*nu_vec.n_elem*p) << i ];
        theta_post[i] = new double[(unsigned long )(modelscount[i]*nu_vec.n_elem*p*n_groups) << i ];
      }
      map[i] = new int[(unsigned long )(modelscount[i]*n_states*3) << i ];
      
    } 
    
    chi[i] = new double[(unsigned long )(modelscount[i]*n_states) << i ];
    if(return_global_null == true)
      psi_post[i] = new double[(unsigned long )(modelscount[i]) << i ];
    if(return_tree == true)
      upsilon[i] = new double[(unsigned long )(modelscount[i]*n_states) << i ];
    
    for(j = 0; j < modelscount[i]; j++)
    {
      for(l = 0; l < pow2(i); l++)
      {
        for(int km = 0; km < sum(n_subgroups); km++ )
        {
          // counts of data-points in node for each group and sub-group
          data[i][(j*pow2(i)+l)*sum(n_subgroups)+km] = 0;
        }          
      }
    }
  }
  
  INDEX_TYPE I_root = init_index(0);  
  for(int i=0; i < n_tot; i++)
    add_data_to_subtree(I_root, 0, 1, 0, X.row(i).t(), cum_subgroups(G(i)-1) + H(i) - 1 );
      
}

void class_tree::add_data_to_subtree( INDEX_TYPE I, 
                                      int level, 
                                      int x_curr, 
                                      int part_count, 
                                      Col<unsigned int> obs,
                                      unsigned int group)
{
  int *NODE_CURR;
  INDEX_TYPE I_child;
  NODE_CURR = get_node_data(I, level);
  NODE_CURR[group] += 1;  // add one observation to the node for that group 
  int i = 0;
  
  if(level < K + 1)
  {
    i = x_curr - 1;
    I_child = make_child_index(I, i, level, (obs(i) >> part_count) & 1);
    add_data_to_subtree(I_child, level+1, x_curr, part_count+1, obs, group);
    for(i=x_curr; i<p; i++)
    {
      I_child = make_child_index(I, i, level, obs(i) & 1);
      add_data_to_subtree(I_child, level+1, i+1, 1, obs, group);
    }
    
  }
  
}

double * class_tree::get_node_nu_post(INDEX_TYPE& I, int level)
{
    return &nu_post[level][(get_node_index(I, level, nu_vec.n_elem * p))];
}

double * class_tree::get_node_theta_post(INDEX_TYPE& I, int level)
{
    return &theta_post[level][(get_node_index(I, level, nu_vec.n_elem * p * n_groups))];
}

double * class_tree::get_node_xi_post(INDEX_TYPE& I, int level)
{
    return &xi_post[level][(get_node_index(I, level, n_states*n_states))];
}

double * class_tree::get_node_varphi_post(INDEX_TYPE& I, int level)
{
    return &varphi_post[level][(get_node_index(I, level, n_states))];
}

double * class_tree::get_node_lambda_post(INDEX_TYPE& I, int level)
{
    return &lambda_post[level][(get_node_index(I, level, n_states*p))];
}

double * class_tree::get_node_psi_post(INDEX_TYPE& I, int level)
{
    return &psi_post[level][(get_node_index(I, level, 1))];
}

double * class_tree::get_node_chi(INDEX_TYPE& I, int level)
{
    return &chi[level][(get_node_index(I, level, n_states))];
}

double * class_tree::get_node_upsilon(INDEX_TYPE& I, int level)
{
    return &upsilon[level][(get_node_index(I, level, n_states))];
}

int * class_tree::get_node_data(INDEX_TYPE& I, int level)
{
    return &data[level][(get_node_index(I, level, sum(n_subgroups) ))];
}

int * class_tree::get_node_map(INDEX_TYPE& I, int level)
{
    return &map[level][(get_node_index(I, level, 3*n_states))];
}



double * class_tree::get_child_xi_post(INDEX_TYPE& I, int i, int level, unsigned short which)
{
    INDEX_TYPE child_index = make_child_index(I,i,level,which);
    return &xi_post[level+1][(get_node_index(child_index,level+1, n_states*n_states))];
}

double * class_tree::get_child_varphi_post(INDEX_TYPE& I, int i, int level, unsigned short which)
{
    INDEX_TYPE child_index = make_child_index(I,i,level,which);
    return &varphi_post[level+1][(get_node_index(child_index,level+1, n_states))];
}

double * class_tree::get_child_lambda_post(INDEX_TYPE& I, int i, int level, unsigned short which)
{
    INDEX_TYPE child_index = make_child_index(I,i,level,which);
    return &lambda_post[level+1][(get_node_index(child_index,level+1, n_states*p))];
}

double * class_tree::get_child_psi_post(INDEX_TYPE& I, int i, int level, unsigned short which)
{
    INDEX_TYPE child_index = make_child_index(I,i,level,which);
    return &psi_post[level+1][(get_node_index(child_index,level+1, 1))];
}

double * class_tree::get_child_chi(INDEX_TYPE& I, int i, int level, unsigned short which)
{
    INDEX_TYPE child_index = make_child_index(I,i,level,which);
    return &chi[level+1][(get_node_index(child_index,level+1, n_states))];
}

double * class_tree::get_child_upsilon(INDEX_TYPE& I, int i, int level, unsigned short which)
{
    INDEX_TYPE child_index = make_child_index(I,i,level,which);
    return &upsilon[level+1][(get_node_index(child_index,level+1, n_states))];
}

int * class_tree::get_child_data(INDEX_TYPE& I, int i, int level, unsigned short which)
{
    INDEX_TYPE child_index = make_child_index(I,i,level,which);
    return &data[level+1][(get_node_index(child_index,level+1, sum(n_subgroups) ))];
}

int * class_tree::get_child_map(INDEX_TYPE& I, int i, int level, unsigned short which)
{
    INDEX_TYPE child_index = make_child_index(I,i,level,which);
    return &map[level+1][(get_node_index(child_index,level+1, 3*n_states))];
}

void class_tree::clear() 
{
  for (int i = 0; i <= (K+1); i++) 
  {
      delete [] chi[i];         
      delete [] data[i];
      
      if(return_global_null == true)
        delete [] psi_post[i];
        
      if(return_tree == true)
      {
        delete [] upsilon[i]; 
        if( i <= K )
        {
          delete [] lambda_post[i];          
          delete [] map[i];
          delete [] varphi_post[i];
          delete [] xi_post[i];
          if (sum(n_subgroups) != n_groups) {
            delete [] nu_post[i];
            delete [] theta_post[i];
          }

        }               
      }
  }
  delete [] chi; chi = NULL;
  delete [] data; data = NULL;
  delete [] modelscount; modelscount= NULL;
  
  if( return_global_null == true)
    delete [] psi_post; 
  psi_post = NULL;
  
  if(return_tree == true)
  {
    delete [] lambda_post; 
    delete [] map; 
    delete [] varphi_post; 
    delete [] xi_post; 
    delete [] upsilon; 
    if (sum(n_subgroups) != n_groups) {
      delete [] nu_post;
      delete [] theta_post;
    }
  }
  lambda_post = NULL;
  map = NULL;
  varphi_post = NULL;
  xi_post = NULL;  
  upsilon = NULL;
  nu_post = NULL;
  theta_post = NULL;
  
}


 
 
 
