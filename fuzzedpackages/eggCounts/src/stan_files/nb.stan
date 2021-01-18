data {
  int J; // number of animals
  int ystarraw[J]; // McMaster count
  real CF[J];
}
parameters{
  real<lower=0> kappa;
  real<lower=0> mu;
  real<lower=0> mui[J];
}
transformed parameters{
  real lambda[J];
  for (i in 1:J){
    lambda[i] = mui[i]/CF[i];
  }
}
model{
  mu ~ gamma(1, 0.001);    // priors
  kappa ~ gamma(1, 0.7);
  mui ~ gamma(kappa, kappa/mu);       // likelihoods, gamma(shape, rate)
  ystarraw ~ poisson(lambda);
}

