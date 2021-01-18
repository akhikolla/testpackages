# include <RcppParallel.h>
# include <RcppArmadillo.h>
# include <sys/time.h>
# include <vector>
# include <algorithm>
# include <cstdlib>
# include <iostream>
#include <RcppArmadilloExtensions/sample.h>
#include <iomanip>

#include <Rcpp.h>

using namespace std;
using namespace Rcpp;
// [[Rcpp::depends(RcppParallel, RcppArmadillo)]]

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
arma::mat boosting_stratify_core(arma::vec& sample, arma::vec& delta, arma::vec& facility, int& num_facility,
                                 arma::mat& X , int& M_stop, double& rate, int adj_variables) {
  int p = X.n_cols, n_sample = sample.size();
  int p_adj = adj_variables;
  arma::vec L2(p_adj);
  arma::vec L1(p_adj);
  arma::rowvec beta = arma::zeros<arma::rowvec>(p);
  arma::vec E = arma::ones(n_sample);
  std::vector<std::vector<int> > facility_idx(num_facility);
  for (int i = 0; i < n_sample; i++)
    facility_idx[facility(sample(i)) - 1].push_back(i);
  arma::mat selection_df(M_stop,p+p_adj);
  
  for (int m = 0; m < M_stop; m++) {
    arma::vec S0 = arma::zeros(n_sample);
    arma::vec S1 = arma::zeros(n_sample);
    // Find j_star, compute S0 and S1, update S1 dynamically
    int j_star = 0;
    double L1_max = 0.0;
    for (int j = 0; j < p; j++) {
      if(p_adj!=0){
        if(j==0){
          j = p_adj;
        }
        L1.zeros();
        arma::mat S1_temp = arma::zeros(n_sample,p_adj);
        arma::uvec adj_var_loc = arma::linspace<arma::uvec>(0,p_adj-1,p_adj);
        for (int g = 0; g < num_facility; g++) {
          for (unsigned int i = facility_idx[g].size() - 1; i--;) {
            int index = facility_idx[g][i];
            if (S0(index) == 0)
              S0(index) = E(index) + (i == facility_idx[g].size() - 1 ? 0.0 : S0(facility_idx[g][i + 1]));
            arma::rowvec X_samp_idx = X.row(sample(index));
            S1_temp.row(index) = X_samp_idx(adj_var_loc) * E(index) + (i == facility_idx[g].size() - 1 ? 0.0 : S1_temp(facility_idx[g][i + 1]));
            if (delta(sample(index)) == 0) continue;
            L1 += (X_samp_idx(adj_var_loc) - S1_temp.row(index) / S0(index));
          }
        }
        
      }
      
      double L1_cur = 0.0;
      arma::vec S1_temp = arma::zeros(n_sample);
      for (int g = 0; g < num_facility; g++) {
        for (unsigned int i = facility_idx[g].size() - 1; i--;) {
          int index = facility_idx[g][i];
          if (S0(index) == 0)
            S0(index) = E(index) + (i == facility_idx[g].size() - 1 ? 0.0 : S0(facility_idx[g][i + 1]));
          S1_temp(index) = X(sample(index), j) * E(index) + (i == facility_idx[g].size() - 1 ? 0.0 : S1_temp(facility_idx[g][i + 1]));
          if (delta(sample(index)) == 0) continue;
          L1_cur += (X(sample(index), j) - S1_temp(index) / S0(index));
        }
      }
      if (std::abs(L1_cur) > std::abs(L1_max)) {
        j_star = j;
        L1_max = L1_cur;
        S1 = S1_temp;
      }
    }
    // Compute S2, L2_star with j_star
    arma::vec S2 = arma::zeros(n_sample);
    double L2_star = 0.0;
    for (int g = 0; g < num_facility; g++) {
      for (unsigned int i = facility_idx[g].size() - 1; i--;) {
        int index = facility_idx[g][i];
        S2(index) = X(sample(index), j_star) * X(sample(index), j_star) * E(index) + (i == facility_idx[g].size() - 1 ? 0.0 : S2(facility_idx[g][i + 1]));
        if (delta(sample(index)) == 0) continue;
        L2_star += (S2(index) / S0(index) - S1(index) * S1(index) / (S0(index) * S0(index)));
      }
    }
    // find L2 to estimate beta for adjustment variables
    if(p_adj!=0){
      arma::mat S2_adj(n_sample,p_adj);
      L2.zeros();
      arma::uvec adj_var_loc = arma::linspace<arma::uvec>(0,p_adj-1,p_adj);
      for (int g = 0; g < num_facility; g++) {
        for (int i = facility_idx[g].size() - 1; i >= 0; i--) {
          int index = facility_idx[g][i];
          arma::rowvec X_row = X.row(sample(index));
          
          S2.row(index) = trans(X_row(adj_var_loc)) * E(index) * X_row(adj_var_loc)  + (i == facility_idx[g].size() - 1 ? 0.0 : S2(facility_idx[g][i + 1])) ;
          if (delta(sample(index)) == 0) continue;
          L2 += (S2(index) / S0(index) - S1(index) * S1(index) / (S0(index) * S0(index)));
        }
      }
    }
    
    
    // Update beta and E
    double temp = rate * L1_max / L2_star;
    beta(j_star) += temp;
    for (int i = 0; i < n_sample; i++)
      E(i) *= exp(X(sample(i), j_star) * temp);
    
    // Update beta for adjustment variables
    if(p_adj!=0){
      arma::uvec adj_var_loc = arma::linspace<arma::uvec>(0,p_adj-1,p_adj);
      beta(adj_var_loc) = beta(adj_var_loc) + rate * L1 / L2;
    }
    selection_df.row(m) = beta;
  }
  return selection_df;
}

arma::vec reverse_vec(arma::vec x) {
  std::reverse(x.begin(), x.end());
  return x;
}

// [[Rcpp::export]]
arma::mat boosting_stratify_path(arma::vec& sample, arma::vec& delta, arma::vec& facility, int& num_facility, arma::mat& X, int& M_stop , double& rate, int adj_variables) {
  int p = X.n_cols-adj_variables, n_sample = sample.size();
  int p_adj = adj_variables;
  arma::rowvec L2(p_adj);
  arma::rowvec L1(p_adj);
  arma::rowvec beta = arma::zeros<arma::rowvec>(p);
  arma::vec E = arma::ones(n_sample);
  std::vector<std::vector<int> > facility_idx(num_facility);
  for (int i = 0; i < n_sample; i++)
    facility_idx[facility(sample(i)) - 1].push_back(i);
  arma::mat selection_df(M_stop,p+p_adj);
  
  for (int m = 0; m < M_stop; m++)  {
    arma::vec S0 = arma::zeros(n_sample);
    arma::vec S1 = arma::zeros(n_sample);
    arma::mat S1_temp_adj = arma::zeros(n_sample,p_adj);
    // Find j_star, compute S0 and S1, update S1 dynamically
    int j_star = 0;
    double L1_max = 0.0;
    for (int j = 0; j < p; j++) {
      if(p_adj!=0){
        if(j==0){
          j = p_adj;
        }
        L1.zeros();
        arma::uvec adj_var_loc = arma::linspace<arma::uvec>(0,p_adj-1,p_adj);
        for (int g = 0; g < num_facility; g++) {
          for (unsigned int i = facility_idx[g].size() - 1; i--;) {
            int index = facility_idx[g][i];
            if (S0(index) == 0)
              S0(index) = E(index) + (i == facility_idx[g].size() - 1 ? 0.0 : S0(facility_idx[g][i + 1]));
            arma::rowvec X_samp_idx = X.row(sample(index));
            S1_temp_adj.row(index) = trans(X_samp_idx(adj_var_loc) * E(index) + (i == facility_idx[g].size() - 1 ? 0.0 : S1_temp_adj(facility_idx[g][i + 1])));
            if (delta(sample(index)) == 0) continue;
            L1 += (trans(X_samp_idx(adj_var_loc)) - S1_temp_adj.row(index) / S0(index));
          }
        }
      }
      double L1_cur = 0.0;
      arma::vec S1_temp = arma::zeros(n_sample);
      for (int g = 0; g < num_facility; g++) {
        for (unsigned int i = facility_idx[g].size() - 1; i--;) {
          int index = facility_idx[g][i];
          if (S0(index) == 0)
            S0(index) = E(index) + (i == facility_idx[g].size() - 1 ? 0.0 : S0(facility_idx[g][i + 1]));
          S1_temp(index) = X(sample(index), j) * E(index) + (i == facility_idx[g].size() - 1 ? 0.0 : S1_temp(facility_idx[g][i + 1]));
          if (delta(sample(index)) == 0) continue;
          L1_cur += (X(sample(index), j) - S1_temp(index) / S0(index));
        }
      }
      if (std::abs(L1_cur) > std::abs(L1_max)) {
        j_star = j;
        L1_max = L1_cur;
        S1 = S1_temp;
      }
    }
    // Compute S2, L2_star with j_star
    arma::vec S2 = arma::zeros(n_sample);
    double L2_star = 0.0;
    for (int g = 0; g < num_facility; g++) {
      for (int i = facility_idx[g].size() - 1; i >= 0; i--) {
        int index = facility_idx[g][i];
        S2(index) = X(sample(index), j_star) * X(sample(index), j_star) * E(index) + (i == facility_idx[g].size() - 1 ? 0.0 : S2(facility_idx[g][i + 1]));
        if (delta(sample(index)) == 0) continue;
        L2_star += (S2(index) / S0(index) - S1(index) * S1(index) / (S0(index) * S0(index)));
      }
    }
    // find L2 to estimate beta for adjustment variables
    if(p_adj!=0){
      arma::vec S2_adj = arma::zeros(n_sample);
      L2.zeros();
      arma::uvec adj_var_loc = arma::linspace<arma::uvec>(0,p_adj-1,p_adj);
      for (int g = 0; g < num_facility; g++) {
        for (unsigned int i = facility_idx[g].size() - 1; i--;) {
          int index = facility_idx[g][i];
          arma::rowvec X_row = X.row(sample(index));
          S2_adj.row(index) = trans(X_row(adj_var_loc)) * E(index) * X_row(adj_var_loc)  + (i == facility_idx[g].size() - 1 ? 0.0 : S2_adj(facility_idx[g][i + 1])) ;
          if (delta(sample(index)) == 0) continue;
          L2 += (S2_adj(index) / S0(index) - S1(index) * S1(index) / (S0(index) * S0(index)));
        }
      }
    }
    // Update beta and E
    double temp = rate * L1_max / L2_star;
    beta(j_star) += temp;
    for (int i = 0; i < n_sample; i++)
      E(i) *= exp(X(sample(i), j_star) * temp);
    
    // Update beta for adjustment variables
    if(p_adj!=0){
      arma::uvec adj_var_loc = arma::linspace<arma::uvec>(0,p_adj-1,p_adj);
      beta(adj_var_loc) = beta(adj_var_loc) + trans(rate * L1 / L2);
    }
    selection_df.row(m) = beta;
  }
  return selection_df;
}

// [[Rcpp::export]]
Rcpp::List boosting_stratify_numselected1(arma::vec& sample, arma::vec& delta, arma::vec& facility, int& num_facility, arma::mat& X, double& num_selected , double& rate, int adj_variables) {
  int p = X.n_cols-adj_variables, n_sample = sample.size();
  int p_adj = adj_variables;
  arma::rowvec L2(p_adj);
  arma::rowvec L1(p_adj);
  arma::rowvec beta = arma::zeros<arma::rowvec>(p);
  arma::vec E = arma::ones(n_sample);
  std::vector<std::vector<int> > facility_idx(num_facility);
  arma::mat selection_df;
  
  for (int i = 0; i < n_sample; i++)
    facility_idx[facility(sample(i)) - 1].push_back(i);
  int selected_key = 0;
  int num_iterations = 0;
  while(selected_key < num_selected+1) {
    num_iterations++;
    arma::vec S0 = arma::zeros(n_sample);
    arma::vec S1 = arma::zeros(n_sample);
    arma::mat S1_temp_adj = arma::zeros(n_sample,p_adj);
    // Find j_star, compute S0 and S1, update S1 dynamically
    int j_star = 0;
    double L1_max = 0.0;
    for (int j = 0; j < p; j++) {
      if(p_adj!=0){
        if(j==0){
          j = p_adj;
        }
        L1.zeros();
        arma::uvec adj_var_loc = arma::linspace<arma::uvec>(0,p_adj-1,p_adj);
        for (int g = 0; g < num_facility; g++) {
          for (unsigned int i = facility_idx[g].size() - 1; i--;) {
            int index = facility_idx[g][i];
            if (S0(index) == 0)
              S0(index) = E(index) + (i == facility_idx[g].size() - 1 ? 0.0 : S0(facility_idx[g][i + 1]));
            arma::rowvec X_samp_idx = X.row(sample(index));
            S1_temp_adj.row(index) = trans(X_samp_idx(adj_var_loc) * E(index) + (i == facility_idx[g].size() - 1 ? 0.0 : S1_temp_adj(facility_idx[g][i + 1])));
            if (delta(sample(index)) == 0) continue;
            L1 += (trans(X_samp_idx(adj_var_loc)) - S1_temp_adj.row(index) / S0(index));
          }
        }
        
      }
      
      double L1_cur = 0.0;
      arma::vec S1_temp = arma::zeros(n_sample);
      for (int g = 0; g < num_facility; g++) {
        for (unsigned int i = facility_idx[g].size() - 1; i--;) {
          int index = facility_idx[g][i];
          if (S0(index) == 0)
            S0(index) = E(index) + (i == facility_idx[g].size() - 1 ? 0.0 : S0(facility_idx[g][i + 1]));
          S1_temp(index) = X(sample(index), j) * E(index) + (i == facility_idx[g].size() - 1 ? 0.0 : S1_temp(facility_idx[g][i + 1]));
          if (delta(sample(index)) == 0) continue;
          L1_cur += (X(sample(index), j) - S1_temp(index) / S0(index));
        }
      }
      if (std::abs(L1_cur) > std::abs(L1_max)) {
        j_star = j;
        L1_max = L1_cur;
        S1 = S1_temp;
      }
    }
    // Compute S2, L2_star with j_star
    arma::vec S2 = arma::zeros(n_sample);
    double L2_star = 0.0;
    for (int g = 0; g < num_facility; g++) {
      for (unsigned int i = facility_idx[g].size() - 1; i--;) {
        int index = facility_idx[g][i];
        S2(index) = X(sample(index), j_star) * X(sample(index), j_star) * E(index) + (i == facility_idx[g].size() - 1 ? 0.0 : S2(facility_idx[g][i + 1]));
        if (delta(sample(index)) == 0) continue;
        L2_star += (S2(index) / S0(index) - S1(index) * S1(index) / (S0(index) * S0(index)));
      }
    }
    // find L2 to estimate beta for adjustment variables
    if(p_adj!=0){
      arma::vec S2_adj = arma::zeros(n_sample);
      L2.zeros();
      arma::uvec adj_var_loc = arma::linspace<arma::uvec>(0,p_adj-1,p_adj);
      for (int g = 0; g < num_facility; g++) {
        for (unsigned int i = facility_idx[g].size() - 1; i--;) {
          int index = facility_idx[g][i];
          arma::rowvec X_row = X.row(sample(index));
          S2_adj.row(index) = trans(X_row(adj_var_loc)) * E(index) * X_row(adj_var_loc)  + (i == facility_idx[g].size() - 1 ? 0.0 : S2_adj(facility_idx[g][i + 1])) ;
          if (delta(sample(index)) == 0) continue;
          L2 += (S2_adj(index) / S0(index) - S1(index) * S1(index) / (S0(index) * S0(index)));
        }
      }
    }
    // Update beta and E
    double temp = rate * L1_max / L2_star;
    beta(j_star) += temp;
    for (int i = 0; i < n_sample; i++)
      E(i) *= exp(X(sample(i), j_star) * temp);
    
    // Update beta for adjustment variables
    if(p_adj!=0){
      arma::uvec adj_var_loc = arma::linspace<arma::uvec>(0,p_adj-1,p_adj);
      beta(adj_var_loc) = beta(adj_var_loc) + trans(rate * L1 / L2);
    }
    selection_df.insert_rows(selection_df.n_rows, beta);
    arma::vec num_selected_temp = (nonzeros(beta));
    selected_key = (num_selected_temp).n_elem - p_adj;
  }
  num_iterations = num_iterations - 1;
  return List::create(Named("beta")=beta, Named("num_iterations")=num_iterations, Named("selection_df")=selection_df);
}

arma::vec loglik(int n, arma::vec delta, arma::mat z, arma::vec beta){
  arma::vec exp_z_beta = exp(z*beta);
  arma::vec rev1 = reverse_vec(exp_z_beta);
  arma::vec cumulative_rev1 = cumsum(rev1);
  arma::vec S0 = reverse_vec(cumulative_rev1);
  arma::vec partial_likelihood = delta%((z*beta)-log(S0));
  
  return(partial_likelihood);
}
// [[Rcpp::export]]
Rcpp::List boosting_stratify_likelihood1(arma::vec& sample, arma::vec& delta, arma::vec& facility, int& num_facility, arma::mat& X, double& rate, double delta_likelihood, int adj_variables) {
  int p = X.n_cols-adj_variables, n_sample = sample.size();
  int p_adj = adj_variables;
  arma::rowvec L2(p_adj);
  arma::rowvec L1(p_adj);
  arma::rowvec beta = arma::zeros<arma::rowvec>(p);
  arma::vec E = arma::ones(n_sample);
  std::vector<std::vector<int> > facility_idx(num_facility);
  for (int i = 0; i < n_sample; i++)
    facility_idx[facility(sample(i)) - 1].push_back(i);
  int num_iterations = 0;
  int key = -1;
  int delay = 3;
  bool mstop = FALSE;
  int maxit = 1*X.n_cols;
  
  NumericVector log_lik;
  NumericVector log_dist;
  arma::mat selection_df;
  
  while(key< maxit) {
    key++;
    num_iterations++;
    arma::vec S0 = arma::zeros(n_sample);
    arma::vec S1 = arma::zeros(n_sample);
    arma::mat S1_temp_adj = arma::zeros(n_sample,p_adj);
    // Find j_star, compute S0 and S1, update S1 dynamically
    int j_star = 0;
    double L1_max = 0.0;
    for (int j = 0; j < p; j++) {
      if(p_adj!=0){
        if(j==0){
          j = p_adj;
        }
        L1.zeros();
        arma::uvec adj_var_loc = arma::linspace<arma::uvec>(0,p_adj-1,p_adj);
        for (int g = 0; g < num_facility; g++) {
          for (unsigned int i = facility_idx[g].size() - 1; i--;) {
            int index = facility_idx[g][i];
            if (S0(index) == 0)
              S0(index) = E(index) + (i == facility_idx[g].size() - 1 ? 0.0 : S0(facility_idx[g][i + 1]));
            arma::rowvec X_samp_idx = X.row(sample(index));
            S1_temp_adj.row(index) = trans(X_samp_idx(adj_var_loc) * E(index) + (i == facility_idx[g].size() - 1 ? 0.0 : S1_temp_adj(facility_idx[g][i + 1])));
            if (delta(sample(index)) == 0) continue;
            L1 += (trans(X_samp_idx(adj_var_loc)) - S1_temp_adj.row(index) / S0(index));
          }
        }
      }
      double L1_cur = 0.0;
      arma::vec S1_temp = arma::zeros(n_sample);
      for (int g = 0; g < num_facility; g++) {
        for (unsigned int i = facility_idx[g].size() - 1; i--;) {
          int index = facility_idx[g][i];
          if (S0(index) == 0)
            S0(index) = E(index) + (i == facility_idx[g].size() - 1 ? 0.0 : S0(facility_idx[g][i + 1]));
          S1_temp(index) = X(sample(index), j) * E(index) + (i == facility_idx[g].size() - 1 ? 0.0 : S1_temp(facility_idx[g][i + 1]));
          if (delta(sample(index)) == 0) continue;
          L1_cur += (X(sample(index), j) - S1_temp(index) / S0(index));
        }
      }
      if (std::abs(L1_cur) > std::abs(L1_max)) {
        j_star = j;
        L1_max = L1_cur;
        S1 = S1_temp;
      }
    }
    // Compute S2, L2_star with j_star
    arma::vec S2 = arma::zeros(n_sample);
    double L2_star = 0.0;
    for (int g = 0; g < num_facility; g++) {
      for (unsigned int i = facility_idx[g].size() - 1; i--;) {
        int index = facility_idx[g][i];
        S2(index) = X(sample(index), j_star) * X(sample(index), j_star) * E(index) + (i == facility_idx[g].size() - 1 ? 0.0 : S2(facility_idx[g][i + 1]));
        if (delta(sample(index)) == 0) continue;
        L2_star += (S2(index) / S0(index) - S1(index) * S1(index) / (S0(index) * S0(index)));
      }
    }
    // find L2 to estimate beta for adjustment variables
    if(p_adj!=0){
      arma::vec S2_adj = arma::zeros(n_sample);
      L2.zeros();
      arma::uvec adj_var_loc = arma::linspace<arma::uvec>(0,p_adj-1,p_adj);
      for (int g = 0; g < num_facility; g++) {
        for (unsigned int i = facility_idx[g].size() - 1; i--;) {
          int index = facility_idx[g][i];
          arma::rowvec X_row = X.row(sample(index));
          S2_adj.row(index) = trans(X_row(adj_var_loc)) * E(index) * X_row(adj_var_loc)  + (i == facility_idx[g].size() - 1 ? 0.0 : S2_adj(facility_idx[g][i + 1])) ;
          if (delta(sample(index)) == 0) continue;
          L2 += (S2_adj(index) / S0(index) - S1(index) * S1(index) / (S0(index) * S0(index)));
        }
      }
    }
    // Update beta and E
    double temp = rate * L1_max / L2_star;
    beta(j_star) += temp;
    for (int i = 0; i < n_sample; i++)
      E(i) *= exp(X(sample(i), j_star) * temp);
    selection_df.insert_rows(selection_df.n_rows, beta);
    
    // Update beta for adjustment variables
    if(p_adj!=0){
      arma::uvec adj_var_loc = arma::linspace<arma::uvec>(0,p_adj-1,p_adj);
      beta(adj_var_loc) = beta(adj_var_loc) + trans(rate * L1 / L2);
    }
    log_lik.push_back(sum(loglik(delta.size(),delta,X,trans(beta))));
    if(key < delay)
      log_dist.push_back(1);
    else{
      log_dist.push_back((log_lik(key)-log_lik(key-1))/(log_lik(key)-log_lik(delay-1)));
      arma::vec idx_arma = arma::linspace<arma::vec>(key-delay,key,key-delay+1);
      IntegerVector idx = as<IntegerVector>(wrap(idx_arma));
      NumericVector subset = log_dist[idx];
      double max_value = max(subset);
      mstop = max_value<delta_likelihood;
    }
    
    if(mstop==TRUE)
      break;
  }
  
  return List::create(Named("beta")=beta, Named("num_iterations")=num_iterations, Named("selection_df")=selection_df);
}
// [[Rcpp::export]]
Rcpp::List boosting_stratify_BIC1(arma::vec& sample, arma::vec& delta, arma::vec& facility, int& num_facility, arma::mat& X, double& rate, bool& early_stop,  int adj_variables, double gamma, bool aic=false) {
  
  int p = X.n_cols-adj_variables, n_sample = sample.size();
  int p_adj = adj_variables;
  arma::rowvec L2(p_adj);
  arma::rowvec L1(p_adj);
  arma::rowvec beta = arma::zeros<arma::rowvec>(p);
  arma::vec E = arma::ones(n_sample);
  std::vector<std::vector<int> > facility_idx(num_facility);
  for (int i = 0; i < n_sample; i++)
    facility_idx[facility(sample(i)) - 1].push_back(i);
  int num_iterations = 0;
  int key = -1;
  int delay = 10;
  bool mstop = FALSE;
  int maxit = 10000;
  arma::mat selection_df;
  
  NumericVector log_lik;
  NumericVector BIC;
  // count number of uncensored events in the data
  arma::uvec uncensored_events = find(delta == 1);
  int num_uncensored = uncensored_events.size();
  
  if(aic==true){
    num_uncensored = exp(2);
    gamma = 0;
  }
  
  
  while(key< maxit) {
    key++;
    num_iterations++;
    arma::vec S0 = arma::zeros(n_sample);
    arma::vec S1 = arma::zeros(n_sample);
    arma::mat S1_temp_adj = arma::zeros(n_sample,p_adj);
    // Find j_star, compute S0 and S1, update S1 dynamically
    int j_star = 0;
    double L1_max = 0.0;
    for (int j = 0; j < p; j++) {
      if(p_adj!=0){
        if(j==0){
          j = p_adj;
        }
        L1.zeros();
        arma::uvec adj_var_loc = arma::linspace<arma::uvec>(0,p_adj-1,p_adj);
        for (int g = 0; g < num_facility; g++) {
          for (unsigned int i = facility_idx[g].size() - 1; i--;) {
            int index = facility_idx[g][i];
            if (S0(index) == 0)
              S0(index) = E(index) + (i == facility_idx[g].size() - 1 ? 0.0 : S0(facility_idx[g][i + 1]));
            arma::rowvec X_samp_idx = X.row(sample(index));
            S1_temp_adj.row(index) = trans(X_samp_idx(adj_var_loc) * E(index) + (i == facility_idx[g].size() - 1 ? 0.0 : S1_temp_adj(facility_idx[g][i + 1])));
            if (delta(sample(index)) == 0) continue;
            L1 += (trans(X_samp_idx(adj_var_loc)) - S1_temp_adj.row(index) / S0(index));
          }
        }
      }
      double L1_cur = 0.0;
      arma::vec S1_temp = arma::zeros(n_sample);
      for (int g = 0; g < num_facility; g++) {
        for (unsigned int i = facility_idx[g].size() - 1; i--;) {
          int index = facility_idx[g][i];
          if (S0(index) == 0)
            S0(index) = E(index) + (i == facility_idx[g].size() - 1 ? 0.0 : S0(facility_idx[g][i + 1]));
          S1_temp(index) = X(sample(index), j) * E(index) + (i == facility_idx[g].size() - 1 ? 0.0 : S1_temp(facility_idx[g][i + 1]));
          if (delta(sample(index)) == 0) continue;
          L1_cur += (X(sample(index), j) - S1_temp(index) / S0(index));
        }
      }
      if (std::abs(L1_cur) > std::abs(L1_max)) {
        j_star = j;
        L1_max = L1_cur;
        S1 = S1_temp;
      }
    }
    // Compute S2, L2_star with j_star
    arma::vec S2 = arma::zeros(n_sample);
    double L2_star = 0.0;
    for (int g = 0; g < num_facility; g++) {
      for (int i = facility_idx[g].size() - 1; i >= 0; i--) {
        int index = facility_idx[g][i];
        S2(index) = X(sample(index), j_star) * X(sample(index), j_star) * E(index) + (i == facility_idx[g].size() - 1 ? 0.0 : S2(facility_idx[g][i + 1]));
        if (delta(sample(index)) == 0) continue;
        L2_star += (S2(index) / S0(index) - S1(index) * S1(index) / (S0(index) * S0(index)));
      }
    }
    // find L2 to estimate beta for adjustment variables
    if(p_adj!=0){
      arma::vec S2_adj = arma::zeros(n_sample);
      L2.zeros();
      arma::uvec adj_var_loc = arma::linspace<arma::uvec>(0,p_adj-1,p_adj);
      for (int g = 0; g < num_facility; g++) {
        for (unsigned int i = facility_idx[g].size() - 1; i--;) {
          int index = facility_idx[g][i];
          arma::rowvec X_row = X.row(sample(index));
          S2_adj.row(index) = trans(X_row(adj_var_loc)) * E(index) * X_row(adj_var_loc)  + (i == facility_idx[g].size() - 1 ? 0.0 : S2_adj(facility_idx[g][i + 1])) ;
          if (delta(sample(index)) == 0) continue;
          L2 += (S2_adj(index) / S0(index) - S1(index) * S1(index) / (S0(index) * S0(index)));
        }
      }
    }
    // Update beta and E
    double temp = rate * L1_max / L2_star;
    beta(j_star) += temp;
    for (int i = 0; i < n_sample; i++)
      E(i) *= exp(X(sample(i), j_star) * temp);
    selection_df.insert_rows(selection_df.n_rows, beta);
    
    // Update beta for adjustment variables
    if(p_adj!=0){
      arma::uvec adj_var_loc = arma::linspace<arma::uvec>(0,p_adj-1,p_adj);
      beta(adj_var_loc) = beta(adj_var_loc) + trans(rate * L1 / L2);
    }
    log_lik.push_back(sum(loglik(delta.size(),delta,X,trans(beta))));
    // count number of parameters in current model
    arma::uvec nonzero_parameters = find(beta != 0);
    int num_parameters = nonzero_parameters.size();
    if(early_stop == TRUE){
      if(key < delay)
        BIC.push_back(-2*log_lik(key) + num_parameters*log(num_uncensored) - 2*gamma*lgamma(num_parameters+1) - 2*gamma*lgamma(p-num_parameters+1));
      else{
        BIC.push_back(-2*log_lik(key) + num_parameters*log(num_uncensored) - 2*gamma*lgamma(num_parameters+1) - 2*gamma*lgamma(p-num_parameters+1));
        arma::vec idx_arma = arma::linspace<arma::vec>(key-delay+1,key,key-delay);
        IntegerVector idx = as<IntegerVector>(wrap(idx_arma));
        NumericVector subset = BIC[idx];
        double min_value = min(subset);
        mstop = min_value>BIC[key-delay];
      }
      if(mstop==TRUE)
        break;
    }
    else{
      if(key < 1000)
        BIC.push_back(-2*log_lik(key) + num_parameters*log(num_uncensored) + 2*gamma*lgamma(p+1) - 2*gamma*lgamma(num_parameters+1) - 2*gamma*lgamma(p-num_parameters+1));
      else{
        num_iterations = which_min(BIC);
        break;
      }
    }
  }
  return List::create(Named("beta")=beta, Named("num_iterations")=num_iterations, Named("Information Criteria")=BIC, Named("selection_df")=selection_df);
}

arma::rowvec dloglik_stratify2_Cpp(int n, arma::vec delta, arma::mat z, arma::vec beta, arma::vec facility, int number_facility=1){
  arma::vec S0(delta.size());
  S0.zeros();
  arma::mat S1(n,beta.size());
  S1.zeros();
  arma::vec temp = exp(z*beta);
  arma::mat S1_pre(z.n_rows, z.n_cols);
  for(unsigned int i=0;i<z.n_rows;i++){
    S1_pre.row(i) = z.row(i)*temp(i);
  }
  int F_key = 0;
  while(F_key < number_facility){
    F_key++;
    arma::uvec loc = find(facility==F_key);
    arma::vec temp2 = reverse_vec(cumsum(reverse_vec(temp(loc))));
    S0(loc) = temp2;
    arma::Col<double> temp3;
    arma::mat S1_row_loc = S1_pre.rows(loc);
    arma::mat S1_loc(loc.size(), beta.size());
    
    for(unsigned int i=0;i<S1_row_loc.n_cols;i++){
      arma::vec S1_col = reverse_vec(cumsum(reverse_vec(S1_row_loc.col(i))));
      S1_loc.col(i) = S1_col;
    }
    S1.rows(loc) = S1_loc;
  }
  
  arma::mat S1_S0(n,beta.size());
  for(int i=0; i< S1.n_cols;i++){
    S1_S0.col(i) = (S1.col(i))/S0;
  }
  
  arma::mat delta_z_s1_s0 = trans(delta)*(z-S1_S0);
  arma::rowvec L1 = sum(delta_z_s1_s0,0);
  
  return(L1);
}

arma::mat cvkfold(int n, int k){
  int fl = floor(n/k);
  arma::mat folds(0,k);
  int index=-1;
  for(int i=0;i<fl;i++){
    for(int j=0;j<k;j++){
      index++;
      arma::vec newrow(k);
      newrow.ones();
      newrow(j)=0;
      folds.insert_rows(index,newrow.t());
    }
  }
  arma::mat folds_shuffle = shuffle(folds,0);
  return(folds_shuffle);
}

arma::vec normalize(arma::vec x){
  arma::vec y = sqrt(x.size())*(x-mean(x))/sqrt(sum(square(x-mean(x))));
  return(y);
}

List boosting_stratify_core_update(arma::vec& sample,
                                   arma::vec& delta,
                                   arma::vec& facility,
                                   int& num_facility,
                                   arma::mat& X,
                                   int& M_start, arma::rowvec& beta_start, arma::vec& L1_start, arma::vec& L2_start, arma::vec& E_start,
                                   int& M_stop, double& rate, int adj_variables) {
  int p = X.n_cols, n_sample = sample.size();
  int p_adj = adj_variables;
  std::vector<std::vector<int> > facility_idx(num_facility);
  for (int i = 0; i < n_sample; i++)
    facility_idx[facility(sample(i)) - 1].push_back(i);
  
  arma::vec L2 = L2_start;
  arma::vec L1 = L1_start;
  arma::rowvec beta = beta_start;
  arma::vec E = E_start;
  
  for (int m = M_start; m < M_stop; m++) {
    arma::vec S0 = arma::zeros(n_sample);
    arma::vec S1 = arma::zeros(n_sample);
    // Find j_star, compute S0 and S1, update S1 dynamically
    int j_star = 0;
    double L1_max = 0.0;
    for (int j = 0; j < p; j++) {
      if(p_adj!=0){
        j=p_adj;
        arma::mat S1_temp = arma::zeros(n_sample,p_adj);
        arma::uvec adj_var_loc = arma::linspace<arma::uvec>(0,p_adj-1,p_adj);
        for (int g = 0; g < num_facility; g++) {
          for (unsigned int i = facility_idx[g].size() - 1; i--;) {
            int index = facility_idx[g][i];
            if (S0(index) == 0)
              S0(index) = E(index) + (i == facility_idx[g].size() - 1 ? 0.0 : S0(facility_idx[g][i + 1]));
            arma::rowvec X_samp_idx = X.row(sample(index));
            S1_temp.row(index) = X_samp_idx(adj_var_loc) * E(index) + (i == facility_idx[g].size() - 1 ? 0.0 : S1_temp(facility_idx[g][i + 1]));
            if (delta(sample(index)) == 0) continue;
            L1 += (X_samp_idx(adj_var_loc) - S1_temp.row(index) / S0(index));
          }
        }
      }
      double L1_cur = 0.0;
      arma::vec S1_temp = arma::zeros(n_sample);
      for (int g = 0; g < num_facility; g++) {
        for (unsigned int i = facility_idx[g].size() - 1; i--;) {
          int index = facility_idx[g][i];
          if (S0(index) == 0)
            S0(index) = E(index) + (i == facility_idx[g].size() - 1 ? 0.0 : S0(facility_idx[g][i + 1]));
          S1_temp(index) = X(sample(index), j) * E(index) + (i == facility_idx[g].size() - 1 ? 0.0 : S1_temp(facility_idx[g][i + 1]));
          if (delta(sample(index)) == 0) continue;
          L1_cur += (X(sample(index), j) - S1_temp(index) / S0(index));
        }
      }
      if (std::abs(L1_cur) > std::abs(L1_max)) {
        j_star = j;
        L1_max = L1_cur;
        S1 = S1_temp;
      }
    }
    // Compute S2, L2_star with j_star
    arma::vec S2 = arma::zeros(n_sample);
    double L2_star = 0.0;
    for (int g = 0; g < num_facility; g++) {
      for (int i = facility_idx[g].size() - 1; i >= 0; i--) {
        int index = facility_idx[g][i];
        S2(index) = X(sample(index), j_star) * X(sample(index), j_star) * E(index) + (i == facility_idx[g].size() - 1 ? 0.0 : S2(facility_idx[g][i + 1]));
        if (delta(sample(index)) == 0) continue;
        L2_star += (S2(index) / S0(index) - S1(index) * S1(index) / (S0(index) * S0(index)));
      }
    }
    // find L2 to estimate beta for adjustment variables
    if(p_adj!=0){
      arma::mat S2_adj(n_sample,p_adj);
      L2.zeros();
      arma::uvec adj_var_loc = arma::linspace<arma::uvec>(0,p_adj-1,p_adj);
      for (int g = 0; g < num_facility; g++) {
        for (unsigned int i = facility_idx[g].size() - 1; i--;) {
          int index = facility_idx[g][i];
          arma::rowvec X_row = X.row(sample(index));
          S2.row(index) = trans(X_row(adj_var_loc)) * E(index) * X_row(adj_var_loc)  + (i == facility_idx[g].size() - 1 ? 0.0 : S2(facility_idx[g][i + 1])) ;
          if (delta(sample(index)) == 0) continue;
          L2 += (S2(index) / S0(index) - S1(index) * S1(index) / (S0(index) * S0(index)));
        }
      }
    }
    // Update beta and E
    double temp = rate * L1_max / L2_star;
    beta(j_star) += temp;
    for (int i = 0; i < n_sample; i++)
      E(i) *= exp(X(sample(i), j_star) * temp);
    // Update beta for adjustment variables
    if(p_adj!=0){
      arma::uvec adj_var_loc = arma::linspace<arma::uvec>(0,p_adj-1,p_adj);
      beta(adj_var_loc) = beta(adj_var_loc) + rate * L1 / L2;
    }
  }
  return List::create(Named("beta") = beta, Named("L1") = L1, Named("L2") = L2, Named("E") = E);
}

double boosting_stratify_core_vec(arma::vec& sample, arma::vec& delta, arma::vec& facility, int& num_facility, arma::vec& X, int& M_stop, double& rate) {
  int n_sample = sample.size();
  double beta = 0;
  arma::vec E = arma::ones(n_sample);
  std::vector<std::vector<int> > facility_idx(num_facility);
  for (int i = 0; i < n_sample; i++)
    facility_idx[facility(sample(i)) - 1].push_back(i);
  
  for (int m = 0; m < M_stop; m++) {
    arma::vec S0 = arma::zeros(n_sample);
    arma::vec S1 = arma::zeros(n_sample);
    // Find j_star, compute S0 and S1, update S1 dynamically
    int j_star = 0;
    double L1_max = 0.0;
    double L1_cur = 0.0;
    arma::vec S1_temp = arma::zeros(n_sample);
    for (int g = 0; g < num_facility; g++) {
      for (unsigned int i = facility_idx[g].size() - 1; i--;) {
        int index = facility_idx[g][i];
        if (S0(index) == 0)
          S0(index) = E(index) + (i == facility_idx[g].size() - 1 ? 0.0 : S0(facility_idx[g][i + 1]));
        S1_temp(index) = X(sample(index)) * E(index) + (i == facility_idx[g].size() - 1 ? 0.0 : S1_temp(facility_idx[g][i + 1]));
        if (delta(sample(index)) == 0) continue;
        L1_cur += (X(sample(index)) - S1_temp(index) / S0(index));
      }
    }
    if (std::abs(L1_cur) > std::abs(L1_max)) {
      L1_max = L1_cur;
      S1 = S1_temp;
    }
    // Compute S2, L2_star with j_star
    arma::vec S2 = arma::zeros(n_sample);
    double L2_star = 0.0;
    for (int g = 0; g < num_facility; g++) {
      for (unsigned int i = facility_idx[g].size() - 1; i--;) {
        int index = facility_idx[g][i];
        S2(index) = X(sample(index)) * X(sample(index)) * E(index) + (i == facility_idx[g].size() - 1 ? 0.0 : S2(facility_idx[g][i + 1]));
        if (delta(sample(index)) == 0) continue;
        L2_star += (S2(index) / S0(index) - S1(index) * S1(index) / (S0(index) * S0(index)));
      }
    }
    // Update beta and E
    double temp = rate * L1_max / L2_star;
    beta += temp;
    for (int i = 0; i < n_sample; i++)
      E(i) *= exp(X(sample(i), j_star) * temp);
  }
  return beta;
}

arma::rowvec boosting_stratify_core_beta(arma::vec& sample, arma::vec& delta, arma::vec& facility, int& num_facility,
                                         arma::mat& X , int& M_stop, double& rate, int adj_variables) {
  int p = X.n_cols, n_sample = sample.size();
  int p_adj = adj_variables;
  arma::vec L2(p_adj);
  arma::vec L1(p_adj);
  arma::rowvec beta = arma::zeros<arma::rowvec>(p);
  arma::vec E = arma::ones(n_sample);
  std::vector<std::vector<int> > facility_idx(num_facility);
  for (int i = 0; i < n_sample; i++)
    facility_idx[facility(sample(i)) - 1].push_back(i);
  for (int m = 0; m < M_stop; m++) {
    arma::vec S0 = arma::zeros(n_sample);
    arma::vec S1 = arma::zeros(n_sample);
    // Find j_star, compute S0 and S1, update S1 dynamically
    int j_star = 0;
    double L1_max = 0.0;
    for (int j = 0; j < p; j++) {
      if(p_adj!=0){
        j=p_adj;
        arma::mat S1_temp = arma::zeros(n_sample,p_adj);
        arma::uvec adj_var_loc = arma::linspace<arma::uvec>(0,p_adj-1,p_adj);
        for (int g = 0; g < num_facility; g++) {
          for (unsigned int i = facility_idx[g].size() - 1; i--;) {
            int index = facility_idx[g][i];
            if (S0(index) == 0)
              S0(index) = E(index) + (i == facility_idx[g].size() - 1 ? 0.0 : S0(facility_idx[g][i + 1]));
            arma::rowvec X_samp_idx = X.row(sample(index));
            S1_temp.row(index) = X_samp_idx(adj_var_loc) * E(index) + (i == facility_idx[g].size() - 1 ? 0.0 : S1_temp(facility_idx[g][i + 1]));
            if (delta(sample(index)) == 0) continue;
            L1 += (X_samp_idx(adj_var_loc) - S1_temp.row(index) / S0(index));
          }
        }
        
      }
      
      double L1_cur = 0.0;
      arma::vec S1_temp = arma::zeros(n_sample);
      for (int g = 0; g < num_facility; g++) {
        for (unsigned int i = facility_idx[g].size() - 1;  i--;) {
          int index = facility_idx[g][i];
          if (S0(index) == 0)
            S0(index) = E(index) + (i == facility_idx[g].size() - 1 ? 0.0 : S0(facility_idx[g][i + 1]));
          S1_temp(index) = X(sample(index), j) * E(index) + (i == facility_idx[g].size() - 1 ? 0.0 : S1_temp(facility_idx[g][i + 1]));
          if (delta(sample(index)) == 0) continue;
          L1_cur += (X(sample(index), j) - S1_temp(index) / S0(index));
        }
      }
      if (std::abs(L1_cur) > std::abs(L1_max)) {
        j_star = j;
        L1_max = L1_cur;
        S1 = S1_temp;
      }
    }
    // Compute S2, L2_star with j_star
    arma::vec S2 = arma::zeros(n_sample);
    double L2_star = 0.0;
    for (int g = 0; g < num_facility; g++) {
      for (unsigned int i = facility_idx[g].size() - 1;  i--;) {
        int index = facility_idx[g][i];
        S2(index) = X(sample(index), j_star) * X(sample(index), j_star) * E(index) + (i == facility_idx[g].size() - 1 ? 0.0 : S2(facility_idx[g][i + 1]));
        if (delta(sample(index)) == 0) continue;
        L2_star += (S2(index) / S0(index) - S1(index) * S1(index) / (S0(index) * S0(index)));
      }
    }
    // find L2 to estimate beta for adjustment variables
    if(p_adj!=0){
      arma::mat S2_adj(n_sample,p_adj);
      L2.zeros();
      arma::uvec adj_var_loc = arma::linspace<arma::uvec>(0,p_adj-1,p_adj);
      for (int g = 0; g < num_facility; g++) {
        for (unsigned int i = facility_idx[g].size() - 1; i--;) {
          int index = facility_idx[g][i];
          arma::rowvec X_row = X.row(sample(index));
          
          S2.row(index) = trans(X_row(adj_var_loc)) * E(index) * X_row(adj_var_loc)  + (i == facility_idx[g].size() - 1 ? 0.0 : S2(facility_idx[g][i + 1])) ;
          if (delta(sample(index)) == 0) continue;
          L2 += (S2(index) / S0(index) - S1(index) * S1(index) / (S0(index) * S0(index)));
        }
      }
    }
    
    // Update beta and E
    double temp = rate * L1_max / L2_star;
    beta(j_star) += temp;
    for (int i = 0; i < n_sample; i++)
      E(i) *= exp(X(sample(i), j_star) * temp);
    
    // Update beta for adjustment variables
    if(p_adj!=0){
      arma::uvec adj_var_loc = arma::linspace<arma::uvec>(0,p_adj-1,p_adj);
      beta(adj_var_loc) = beta(adj_var_loc) + rate * L1 / L2;
    }
  }
  return beta;
}
// [[Rcpp::export]]
Rcpp::List cross_validation_func_update(int K, arma::vec time, arma::vec delta, arma::mat z, arma::vec facility,
                                        double rate=0.01, int track=10, int M_stop=100, int adj_variables=0){
  arma::mat selection_df;
  arma::vec fac_unique = arma::unique(facility);
  int num_facility = fac_unique.size();
  arma::vec IDvector(z.n_rows);
  for(unsigned int i=0; i< z.n_rows; i++) IDvector(i) = i;
  arma::rowvec boosting_coef = boosting_stratify_core_beta(IDvector, delta, facility, num_facility, z, M_stop, rate, adj_variables);
  arma::vec coef_nonzero = boosting_coef(find(abs(boosting_coef)>(rate/2)));
  
  int p1 = coef_nonzero.size();
  int p=z.n_cols;
  int p_adj = adj_variables;
  arma::vec beta(p);
  beta.zeros();
  int N = time.size();
  arma::Col<double> cvrisk_all;
  arma::mat beta_cv(K,p);
  beta_cv.zeros();
  arma::mat weight = cvkfold(N,K);
  arma::vec unique_facility = arma::unique(facility);
  
  List N_k;
  List delta_k;
  List time_k;
  List z_k;
  List F_k;
  List facility_k;
  List IDvector_k;
  
  for (int k=0;k<K;k++) {
    arma::uvec weight_k = arma::find(weight.col(k)>0);
    N_k.push_back(sum(weight.col(k)>0));
    arma::vec delta_weight_k = delta(weight_k);
    delta_k.push_back(delta_weight_k);
    arma::vec time_weight_k = time(weight_k);
    time_k.push_back(time_weight_k);
    arma::mat z_k_temp = z.rows(weight_k);
    arma::vec  id = arma::conv_to<arma::vec>::from(weight_k);
    IDvector_k.push_back(id);
    
    int j=p1;
    arma::mat k_temp(weight_k.size(), p-p1);
    while (j<p){
      j++;
      k_temp.col(j-p1-1) = normalize(z_k_temp.col(j-1));
    }
    if(p1>0){
      arma::uvec zero_p1_list = arma::linspace<arma::uvec>(0,p1-1,p1);
      arma::mat z_k_temp_modified = z_k_temp.cols(zero_p1_list);
      z_k.push_back(join_rows(z_k_temp_modified,k_temp));
    }
    else{
      z_k.push_back(k_temp);
    }
    arma::vec facility_k_temp=facility(find(weight.col(k)>0));
    arma::vec unique_provfs = arma::unique(facility_k_temp);
    F_k.push_back(unique_provfs.size());
    
    j=0;
    while (j<(int)unique_provfs.size()){
      arma::uvec location = find(facility_k_temp==unique_provfs(j));
      for(unsigned int i=0;i<location.size();i++){
        facility_k_temp(location(i))=j+1;
      }
      j++;
    }
    facility_k.push_back(facility_k_temp);
    
  }
  bool stage2_key=FALSE;
  
  List beta_start;
  List L1_start;
  List L2_start;
  List E_start;
  for(int k =0; k<K; k++){
    arma::rowvec beta_temp = arma::zeros<arma::rowvec>(p);
    arma::vec L2_temp(p_adj);
    arma::vec L1_temp(p_adj);
    int N_kk = N_k[k];
    arma::vec E_temp = arma::ones(N_kk);
    
    beta_start.push_back(beta_temp);
    L2_start.push_back(L2_temp);
    L1_start.push_back(L1_temp);
    E_start.push_back(E_temp);
  }
  
  int m=0;
  int m1;
  while (stage2_key==FALSE){
    double cvrisk=0;
    for (int k=0; k<K;k++) {
      arma::mat zk_k = z_k[k];
      arma::vec facility_kk = facility_k[k];
      arma::vec delta_kk = delta_k[k];
      arma::vec time_kk = time_k[k];
      arma::vec unique_facility_kk = arma::unique(facility_kk);
      arma::vec id_kk = IDvector_k[k];
      arma::rowvec beta_temp = beta_start[k];
      arma::vec L2_temp = L2_start[k];
      arma::vec L1_temp = L1_start[k];
      arma::vec E_temp = E_start[k];
      m1 = m+10L;
      
      List cv_res_k = boosting_stratify_core_update(id_kk, delta, facility, num_facility, z,
                                                    m, beta_temp,  L1_temp, L2_temp, E_temp,
                                                    m1, rate, adj_variables);
      
      beta_start[k] = cv_res_k["beta"];
      L2_start[k] = cv_res_k["L2"];
      L1_start[k] = cv_res_k["L1"];
      E_start[k] = cv_res_k["E"];
      
      int N_kk = N_k[k];
      int F_kk = F_k[k];
      arma::rowvec beta_cv_rowk = beta_start[k];
      arma::rowvec diff = dloglik_stratify2_Cpp(N_kk,delta_kk,zk_k,trans(beta_cv_rowk),facility_kk,F_kk);
      
      arma::mat diff_abs = abs(diff);
      double j_star_value= diff_abs.max();
      arma::uvec j_star = find(diff_abs==j_star_value);
      arma::vec z_j_star = z.cols(j_star);
      
      beta_cv_rowk(j_star) = beta_cv_rowk(j_star) + rate*boosting_stratify_core_vec(id_kk, delta, facility, num_facility, z_j_star, m1, rate);
      beta_cv = join_cols(beta_cv,beta_cv_rowk);
      double cv_temp=-(sum(loglik(N,delta,z,trans(beta_cv_rowk))))+sum(loglik(N_kk,delta_kk,zk_k,trans(beta_cv_rowk)));
      cvrisk=cvrisk+cv_temp;
    }
    cvrisk_all.insert_rows(m/10,1);
    cvrisk_all(m/10) = cvrisk;
    
    if(m>=(100+(20*track))){
      arma::uvec compare_max1 = arma::linspace<arma::uvec>(track,(m/10)-10,(m/10)-9-track); // change m when iterating by more than 1
      arma::uvec compare_max2 = arma::linspace<arma::uvec>((m/10)-9,m/10,10);
      arma::vec cvrisk_compare1 = cvrisk_all(compare_max1);
      arma::vec cvrisk_compare2 = cvrisk_all(compare_max2);
      
      if ( min(cvrisk_compare2) >= min(cvrisk_compare1) ) {
        stage2_key=TRUE;
      }
    }
    else{
      stage2_key=FALSE;
    }
    m=m+10;
    if(m>5000){
      stage2_key=TRUE;
    }
  }
  return List::create(Named("mstop")=m1, Named("cvrisk") = cvrisk_all);
}
