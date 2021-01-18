data {
  int<lower=0> N; 			// number of patient
  vector[N] auc; 			  // log auc
  matrix[N,2] dose; 		// log dose + intercept
  real beta0;
  real mu;
}
parameters {
  vector[2] b;
  real<lower=0,upper=1> sigma;
}
model{
  auc ~ normal(dose*b, sigma);
  sigma ~ beta(1, 1);
  b[1] ~ normal(mu, beta0);
  b[2] ~ normal(1, beta0);
}
