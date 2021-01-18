
#include chunks/mylib.stan

data {
  int<lower=0> n;
  int<lower=0> p;
  int<lower=0> q;
  int<lower=0> y[n];
  matrix[n, p] X;
  matrix[n, q] Z;
  row_vector[p] x_mean;
  vector<lower=0>[p] x_sd;
  row_vector[q] z_mean;
  vector<lower=0>[q] z_sd;  
  int<lower=0, upper=1> approach;
  real mu_beta;
  real<lower=0> sigma_beta;
  real mu_psi;
  real<lower=0> sigma_psi;
}

parameters {
  vector[q] psi_std;
  vector[p] beta_std;
}

transformed parameters{
  vector[q] psi;
  vector[p] beta;
  
  if(p==1){
    beta[1] = beta_std[1]/x_sd[1];
  }else{
    beta[2:p] = beta_std[2:p] ./ x_sd[2:p];
    beta[1] = beta_std[1]/x_sd[1] - x_mean[2:p]*beta[2:p];
  }
  
  if(q==1){
    psi[1] = psi_std[1]/z_sd[1];
  }else{
    psi[2:q] = psi_std[2:q] ./ z_sd[2:q];
    psi[1] = psi_std[1]/z_sd[1] - z_mean[2:q]*psi[2:q];
  }
}

model{
    // likelihood:
    vector[n] eta1;
    vector[n] eta2;
    real mu[n];
    real theta[n];
    vector[n] omega;
    eta1 = Z*psi_std;
    eta2 = X*beta_std;
    for(t in 1:n){
      omega[t] = inv_logit(eta1[t]);
      mu[t] = exp(eta2[t]);
      theta[t] = lambertW(mu[t]);
      if(y[t] == 0)
        target += log_sum_exp(bernoulli_lpmf(1 | omega[t]),
                              bernoulli_lpmf(0 | omega[t]) + bell_lpmf(y[t] | theta[t]));
        else
          target += bernoulli_lpmf(0 | omega[t]) + bell_lpmf(y[t] | theta[t]);
    }
    if(approach==1){
      // prior distributions:
      beta_std ~ normal(mu_beta, sigma_beta);
      psi_std ~ normal(mu_psi, sigma_psi);
    }
}


