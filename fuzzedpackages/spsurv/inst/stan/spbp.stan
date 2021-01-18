// include Functions block.
#include /include/loglikbp.stan

// Data block (important).
data{
  // setting observed data:
  int<lower=1> n;
  int<lower=1> m;
  int<lower=1> q;
  real tau;
  int<lower=0, upper=1> approach;
  int<lower=0, upper=3> dist;
  int<lower=0, upper=1> null;
  int<lower=0, upper=2> M;
  vector<lower=0, upper = 1>[n] status;
  int<lower=0> id [n];
  vector<lower=0>[n] z; // for coding purposes
  vector<lower=0>[n] time;

  // observed arrays:
  matrix[n,q] X;
  matrix[n,m] b;
  matrix[n,m] B;

  // setting hyperparameters:
  int<lower=0> priordist [2];
  vector [4] priorpars;

  // beta (vector)
  int<lower=0> priordist_beta [q];
  vector [q] location_beta;
  vector<lower=0> [q] scale_beta;

  // Standard quantities
  vector<lower=0>[q] std; // feature standard deviatons
  vector[q] means;        // feature means
}

// Parametes block (important).
parameters{
  vector[q] beta_scaled;           // feature effect
  vector<lower=0>[m] gamma_scaled; //BP basis effect
}

transformed parameters{
    //// Declare statement
      vector[n] log_lik;              // declare log likelihood
    // scaled coefficients
      vector[q] beta;                 // declare scaled beta
      vector<lower=0>[m] gamma;       // declare scaled BP gamma

    // reparametrized beta (due to HM dynamics)
      vector[q] beta_std = (beta_scaled - to_vector(location_beta)) ./ to_vector(scale_beta);

    //// Define declared variables
      beta = beta_scaled ./ std;      // define beta to original scale

       // definition if model is AFT
      if(M == 2){
          gamma = gamma_scaled * exp(sum(beta_scaled .* means ./ std));
      }// if model is PO or PH
      else{
          gamma = gamma_scaled * exp(-sum(beta_scaled .* means ./ std));
      }
        // definition if model is null
      if(null == 1){
          log_lik = loglik_null(beta_scaled, gamma_scaled, status, X, b, B, M, dist, id, z);
      } // definition if model is P0
      else{
        if(M == 0){
            log_lik = loglik_po(beta_scaled, gamma_scaled, status, X, b, B, dist, id, z);
        } // definition if model PH
        else if( M == 1){
            log_lik = loglik_ph(beta_scaled, gamma_scaled, status, X, b, B, dist, id, z);
        } // definition if model is AFT
        else{
            log_lik = loglik_aft(time, beta_scaled, gamma_scaled, status, X, b, B, dist, id, z);
        }
      }
  }

// Model block (important).
model{

  if(approach == 1){ // priors
  ////// Beta Prior
      for(i in 1:q){
        if(priordist_beta[i] == 0){
          beta_std ~ normal(0, 1);
        }
        else{
          beta_std ~ cauchy(0, 1);
        }
      }
  ////// Gamma Prior
      if(priordist[1] == 1){
        gamma_scaled ~ gamma(priorpars[1], priorpars[2]);
      }
      else if(priordist[1] == 2){
        gamma_scaled ~ inv_gamma(priorpars[1], priorpars[2]);
      }
      else{
        gamma_scaled ~ lognormal(priorpars[1], priorpars[2]);
      }
   }
     target += sum(log_lik);
 }

// Final line empty to avoid warnings.
