// include Functions block.
#include /include/loglikbp.stan

// Data block (important).
data{
  // setting observed data:
  int<lower=1> n;
  int<lower=1> m;
  int<lower=1> q;
  real tau;
  int<lower=0, upper=3> dist;
  int<lower=0, upper=1> null;
  int<lower=0, upper=2> M;
  vector<lower=0, upper = 1>[n] status;
  int<lower=0> id [n];
  vector<lower=0>[n] time;

  // observed arrays:
  matrix[n,q] X;
  matrix[n,m] b;
  matrix[n,m] B;

  // setting hyperparameters:
  int<lower=0> priordist [2];
  real priorpars [4];

  // beta (vector)
  int<lower=0> priordist_beta [q];
  real location_beta [q];
  real<lower=0> scale_beta [q];

  // Standard quantities
  vector<lower=0>[q] std; // feature standard deviatons
  vector[q] means; // feature means
}

// Parametes block (important).
parameters{
  vector[q] beta; // feature effect
  vector[n] z;    // random frailty
  vector<lower=0>[m] gamma; //BP basis effect
  //// hyperparameters
    real <lower=0> kappa; // multiplicative frailty precision
}

transformed parameters{

    vector[n] log_lik; //log likelihood
    vector[m] nu; // exp of the BP basis effect

    real<lower=0> sigma; // stdev of the additive frailty

    vector[q] beta_std; // standardized feature effect
    vector<lower=0>[m] gamma_std; // standardized BP basis effect

    beta_std = beta ./ std;
    if(M == 2){
        gamma_std = gamma * exp(sum(beta .* means ./ std));
    }
    else{
        gamma_std = gamma * exp(-sum(beta .* means ./ std));
    }

    nu = log(gamma);
    sigma = inv_sqrt(kappa);

    if(null == 1){
        log_lik = loglik_null(beta, gamma, status, X, b, B, M, dist, id, z);
    }
    else{
      if(M == 0){
          log_lik = loglik_po(beta, gamma, status, X, b, B, dist, id, z);
      }
      else if( M == 1){
          log_lik = loglik_ph(beta, gamma, status, X, b, B, dist, id, z);
      }
      else{
          log_lik = loglik_aft(time, beta, gamma, status, X, b, B, dist, id, z);
      }
    }
}

// Model block (important).
model{

////// Beta Prior
    for(i in 1:q){
      print("location = ", location_beta[i], " scale = ", scale_beta[i]);

      if(priordist_beta[i] == 0){
        beta ~ normal(location_beta[i], scale_beta[i]);
      }
      else{
        beta ~ cauchy(location_beta[i], scale_beta[i]);
      }
    }
////// Gamma Prior
    print("h1 = ", priorpars[1], " h2 = ", priorpars[2]);
    if(priordist[1] == 1){
      gamma ~ gamma(priorpars[1], priorpars[2]);
    }
    else if(priordist[1] == 2){
      gamma ~ inv_gamma(priorpars[1], priorpars[2]);
    }
    else{
      gamma ~ lognormal(priorpars[1], priorpars[2]);
    }

////// Frailty Priors
	  //frailty
	  if(dist == 1){
	    z ~ gamma(kappa, kappa); // gamma frailty
	  }
	  else if(dist == 2){
	    z ~ normal(0, sigma); // gaussian frailty
	  }
	  else if(dist == 3){
	    z ~ cauchy(0, sigma); // t-student frailty
	  }
	  else{
	    // print("dist ", dist);
	    z ~ normal(0, 0.000001);
	  }
	  // kappa
	  if(dist > 0 ){
	        if(priordist[2] == 0){
      kappa ~ normal(priorpars[3], priorpars[4]);
    }
    else if(priordist[2] == 1){
      kappa ~ gamma(priorpars[3], priorpars[4]);
    }
    else if(priordist[2] == 2){
      kappa ~ inv_gamma(priorpars[3], priorpars[4]);
    }
    else{
      kappa ~ lognormal(priorpars[3], priorpars[4]);
    }
  }
	 target += sum(log_lik);
}

// Final line empty to avoid warnings.
