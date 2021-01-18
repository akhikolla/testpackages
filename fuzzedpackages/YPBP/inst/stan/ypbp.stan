
#include /chunks/logliksbp.stan

data{
  int<lower=1> n;
  int<lower=1> m;
  int<lower=1> q;
  vector[n] status;
  matrix[n,q] Z;
  matrix[n,m] g;
  matrix[n,m] G;
  real<lower=0> tau;
  real h1_gamma;
  real h2_gamma;
  real mu_psi;
  real mu_phi;
  real<lower=0> sigma_psi;
  real<lower=0> sigma_phi;
  int<lower=1, upper=2> M;
  int<lower=0, upper=1> approach;
}


parameters{
  vector[q] psi;
  vector[q] phi;
  vector<lower=0>[m] gamma;
}


transformed parameters{
  vector[n] loglik;
  if(M==1){
    loglik = loglik1_bp(status, Z, g, G, tau, gamma, psi, phi);
  }else{
    loglik = loglik2_bp(status, Z, g, G, tau, gamma, psi, phi);
  }
}


model{
  target += sum(loglik);
  if(approach==1){
    gamma ~ lognormal(h1_gamma, h2_gamma);
    psi ~ normal(0, sigma_psi);
    phi ~ normal(0, sigma_phi);
  }
}





