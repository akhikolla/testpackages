//
// model 2:
//       full stratified model
//
//      theta_g ~ N(mu_g, sigma^2)
//      mu_g    ~ N(0, 10^3)
//

data {
	int<lower=0>  SIZE;
	vector[SIZE]  Y;
	vector[SIZE]  SIGY;
  real<lower=0> B;
	real<lower=0> DELTA;
  int<lower=0, upper=1> PRIORSIG;
}

parameters {
	vector[SIZE] mu;
	vector<lower=0, upper=1>[SIZE] uvs;
	vector[SIZE] nvs;
}

transformed parameters {
	vector<lower=0>[SIZE] vs;

  if (0 == PRIORSIG) {
    vs = exp(log(SIGY) + (uvs * 2 - 1) * DELTA);
  } else {
    vs = exp(log(SIGY) + nvs * sqrt(DELTA));
  }
}

model {
	mu  ~ normal(0, sqrt(B));
  uvs ~ uniform(0,1);
  nvs ~ normal(0,1);
  Y   ~ normal(mu, vs);
}

generated quantities {
  vector[SIZE] log_lik;
  for (i in 1:SIZE) {
    log_lik[i] = normal_lpdf(Y[i] | mu[i], vs[i]);    
  }
}
