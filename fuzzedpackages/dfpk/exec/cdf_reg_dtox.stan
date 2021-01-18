data {
  int<lower=0> N;         // number of previous patients
  int y[N];               // binary response
  real dose[N];           // log dose 
  real beta0mean[2];
  real beta1mean[2];
}
parameters {
  real<lower=beta0mean[1], upper=beta0mean[2]> beta0; 
  real<lower=beta1mean[1], upper=beta1mean[2]> beta1;
}
model {
  real p[N];          // probabilities
  vector[N] z;        // logistic transformation
  for (n in 1:N){
  z[n] = -beta0 + dose[n]*beta1;
  p[n] = normal_cdf(z[n],0,1);
  }
  y ~ bernoulli(p);
  beta0 ~ uniform(beta0mean[1], beta0mean[2]);
  beta1 ~ uniform(beta1mean[1], beta1mean[2]);
}
