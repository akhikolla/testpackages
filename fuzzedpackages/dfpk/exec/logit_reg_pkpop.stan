data {
  int<lower=0> N;         // number of previous patients
  int y[N];               // binary response
  real dose[N];           // log dose 
  real beta3mean;
  real beta4mean;
}
parameters {
  real<lower=0, upper=beta3mean> beta3; 
  real<lower=0, upper=beta4mean> beta4;
}
model {
  real p[N];          // probabilities
  vector[N] z;        // logistic transformation
  for (n in 1:N){
  z[n] = beta3 - dose[n]*beta4;
  p[n] = 1 / (1 + exp(z[n]));
  }
  y ~ bernoulli(p);
  beta3 ~ uniform(0.0, beta3mean);
  beta4 ~ uniform(0.0, beta4mean);
}
