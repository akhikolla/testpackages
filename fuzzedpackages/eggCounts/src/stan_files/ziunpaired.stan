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
  real<lower=0> mua[Ja];
  real<lower=0> mub[Jb];
  real<lower=0,upper=1> phi;
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
  mu ~ gamma(1,0.001);    // prior
  kappa ~ gamma(1,0.7);
  delta ~ beta(1,1);
  phi ~ beta(1,1);
  // likelihoods ZIPoGa
    mub ~ gamma(kappa, kappa/mu); 
    mua ~ gamma(kappa, kappa/mu); 
  for (n in 1:Jb){
    if (ystarbraw[n] == 0)
      target += log_sum_exp(bernoulli_lpmf(1 | phi), bernoulli_lpmf(0 | phi)+poisson_lpmf(ystarbraw[n] | lambdab[n]));
    else
      target += bernoulli_lpmf(0 | phi) + poisson_lpmf(ystarbraw[n] | lambdab[n]);
  }
  for (n in 1:Ja){
    if (ystararaw[n] == 0)
      target += log_sum_exp(bernoulli_lpmf(1 | phi), bernoulli_lpmf(0 | phi)+poisson_lpmf(ystararaw[n] | lambdaa[n]));
    else 
      target += bernoulli_lpmf(0 | phi) + poisson_lpmf(ystararaw[n] | lambdaa[n]);
  }
}
