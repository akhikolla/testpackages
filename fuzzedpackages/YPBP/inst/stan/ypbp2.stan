
#include /chunks/logliksbp.stan

data{
  int<lower=1> n;
  int<lower=1> m;
  int<lower=1> q;
  int<lower=1> p;
  vector[n] time;
  vector[n] status;
  matrix[n, q] Z;
  matrix[n, p] X;
  matrix[n, m] g;
  matrix[n, m] G;
  real<lower=0> tau;
  real h1_gamma;
  real h2_gamma;
  real mu_psi;
  real mu_phi;
  real mu_beta;
  real<lower=0> sigma_psi;
  real<lower=0> sigma_phi;
  real<lower=0> sigma_beta;
  int<lower=3, upper=4> M;
  int<lower=0, upper=1> approach;
}

parameters{
  vector[q] psi;
  vector[q] phi;
  vector[p] beta;
  vector<lower=0>[m] gamma;
}


transformed parameters{
  vector[n] loglik;
  if(M==3){
    loglik = loglik3_bp(status, Z, X, g, G, tau, gamma, psi, phi, beta);
  }else{
    loglik = loglik4_bp(status, Z, X, g, G, tau, gamma, psi, phi, beta);
  }
}


model{
  target += sum(loglik);
  if(approach==1){
    gamma ~ lognormal(h1_gamma, h2_gamma);
    psi ~ normal(mu_psi, sigma_psi);
    phi ~ normal(mu_phi, sigma_phi);
    beta ~ normal(mu_beta, sigma_beta);
  }
}





