data{
  int J; // number of animals
  int ystararaw[J]; // after treatment McMaster count
  int ystarbraw[J]; // before treatment McMaster count
  real fpre[J];
  real fpost[J];
}

parameters{
  real<lower=0,upper=1> delta;
  real<lower=0> mu;
}

transformed parameters{
  real lambdaa[J];
  real lambdab[J];
  for (i in 1:J){
    lambdab[i] = mu/fpre[i];
    lambdaa[i] = delta*mu/fpost[i];
    }
}

model{
  mu ~ gamma(1, 0.001);    // prior
  delta ~ beta(1,1);
  ystararaw ~ poisson(lambdaa);         // likelihoods
  ystarbraw ~ poisson(lambdab);
}
