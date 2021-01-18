data{
  int Ja; // number of animals
  int Jb;
  int ystararaw[Ja]; // after treatment McMaster count
  int ystarbraw[Jb]; // before treatment McMaster count
  real fpost[Ja];
  real fpre[Jb];
}
parameters{
  real<lower=0> kappa;
  real<lower=0> mu;
  real<lower=0,upper=1> delta;
  real<lower=0> mub[Jb]; // true epg before treatment
  real<lower=0> mua[Ja]; // true epg after treatment
}
transformed parameters{
  real lambdaa[Ja];
  real lambdab[Jb];
  for (i in 1:Jb){
    lambdab[i] = mub[i]/fpre[i];
  }
  for (i in 1:Ja){
    lambdaa[i] = delta*mua[i]/fpost[i];
  }
}
model{
  mu ~ gamma(1,0.001);    // priors
  kappa ~ gamma(1,0.7);
  delta ~ beta(1,1);
  mub ~ gamma(kappa, kappa/mu); // likelihoods
  mua ~ gamma(kappa, kappa/mu);   
  ystarbraw ~ poisson(lambdab);
  ystararaw ~ poisson(lambdaa);
}
