data {
  int<lower=0> N;         // number of previous patients
  int y[N];               // binary response
  real dose[N];           // log dose
  real beta2mean;
  real beta3mean;
}
parameters {
  real<lower=0, upper=beta2mean> beta2; 
  real<lower=0, upper=beta3mean> beta3;
}
model {
  real p[N];          // probabilities
  vector[N] z;        // logistic transformation
  for (n in 1:N){
  z[n] = beta2 - dose[n]*beta3;
  p[n] = 1 / (1 + exp(z[n]));
  }
  y ~ bernoulli(p);
  beta2 ~ uniform(0.0, beta2mean);
  beta3 ~ uniform(0.0, beta3mean);
}
