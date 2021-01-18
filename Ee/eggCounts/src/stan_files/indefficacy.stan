data{
  int J; // number of animals
  int ystararaw[J]; // after treatment McMaster count
  int ystarbraw[J]; // before treatment McMaster count
  real fpre[J];
  real fpost[J];
}
parameters{
  real<lower=0> kappa;
  real<lower=0> mu;
  real<lower=0> delta[J];
  real<lower=0> delta_shape;
  real<lower=0, upper=1> delta_mu;
  real<lower=0> mub[J];
}
transformed parameters{
  real lambdaa[J];
  real lambdab[J];
  for (i in 1:J){
    lambdab[i] = mub[i]/fpre[i];
    lambdaa[i] = delta[i]*mub[i]/fpost[i];
  }
}
model{
  mu ~ gamma(1,0.001);    // prior
  kappa ~ gamma(1,0.7);
  delta ~ gamma(delta_shape, delta_shape/delta_mu); // shape, rate
  delta_shape ~ normal(2, 1);
  delta_mu ~ beta(1,1);
  mub ~ gamma(kappa, kappa/mu);   // likelihoods 
  ystararaw ~ poisson(lambdaa);
  ystarbraw ~ poisson(lambdab);
}
