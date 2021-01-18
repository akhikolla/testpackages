
#include chunks/mylib.stan

data{
  int<lower=1> n;
  int<lower=1> p;
  int y[n];
  matrix[n, p] X;
  row_vector[p] x_mean;
  vector<lower=0>[p] x_sd;
  int<lower=0, upper=1> approach;
  real mu_beta;
  real<lower=0> sigma_beta;
}

parameters{
  vector[p] beta_std;
}


transformed parameters{
  vector[p] beta;
  if(p==1){
    beta[1] = beta_std[1]/x_sd[1];
  }else{
    beta[2:p] = beta_std[2:p] ./ x_sd[2:p];
    beta[1] = beta_std[1]/x_sd[1] - x_mean[2:p]*beta[2:p];
  }
}

model{
  real theta[n];
  vector[n] eta = X*beta_std;
  real mu[n];
  for(i in 1:n){
    mu[i] = exp(eta[i]);
    theta[i] = lambertW(mu[i]);
  }

  target += loglik_bell(y, theta);
  if(approach==1){
    beta_std ~ normal(mu_beta, sigma_beta);
  }

}
