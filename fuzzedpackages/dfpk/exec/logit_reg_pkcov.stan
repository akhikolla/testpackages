data {
  int<lower=0> N;         // number of previous patients
  int y[N];               // binary response
  real dose[N];           // log dose 
  vector[N] dauc;
  real beta0mean;
  real beta1mean[2];
}
parameters {
  real<lower=0, upper=25> beta0;
  real<lower=beta1mean[1], upper=beta1mean[2]> beta1;
  real <lower=0, upper = 5> beta2;
}
model {
  real p[N];          // probabilities
  vector[N] z;        // logistic transformation
  for (n in 1:N){
  z[n] = beta0 - dose[n]*beta1 - beta2*dauc[n];
  p[n] = 1 / (1 + exp(z[n]));
  }
  y ~ bernoulli(p);
  beta0 ~ uniform(0.0,25);
  beta1 ~ uniform(beta1mean[1], beta1mean[2]);
  beta2 ~ uniform(0.0,5);
}
