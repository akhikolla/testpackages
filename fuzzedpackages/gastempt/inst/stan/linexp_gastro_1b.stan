// A minimal program to fit linexp gastric emptying curves with Stan
// Conventional method with loop indexing
// The prior mean of the initial volume can be set from the calling function.
// All other priors are fixed.

data{
  real prior_v0;
  int<lower=0> n; // Number of data
  int<lower=0> n_record; // Number of records
  int record[n];
  vector[n] minute;
  vector[n] volume;
}

parameters{
  vector<lower=0>[n_record] v0;
  vector<lower=0>[n_record] kappa;
  vector<lower=0>[n_record] tempt;
  real <lower=0> sigma;
  real <lower=0> mu_kappa;
  real <lower=0> sigma_kappa;
}


model{
  int  reci;
  real v0r;
  real kappar;
  real temptr;
  real vol[n];
  mu_kappa ~ normal(1.5,0.5);
  sigma_kappa ~ normal(1,0.5);

  v0    ~ normal(prior_v0, 100);
  kappa ~ lognormal(mu_kappa, sigma_kappa);
  tempt ~ normal(60, 20);
  sigma ~ gamma(20, 0.5);

for (i in 1:n){
   reci = record[i];
   v0r = v0[reci];
   kappar = kappa[reci];
   temptr = tempt[reci];
   vol[i] = v0r*(1+kappar*minute[i]/temptr)*exp(-minute[i]/temptr);
  }
  volume ~ normal(vol, sigma);
}


