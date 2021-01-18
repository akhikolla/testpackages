//
// model 1:
//    no subgroup effect
//
//      thetag ~ N(mu, sigma^2)
//      mu     ~ N(0, B)
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
	real tau;
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
	tau ~ normal(0, sqrt(B));
  uvs ~ uniform(0,1);
  nvs ~ normal(0,1);
  Y   ~ normal(tau, vs);
}

generated quantities {
  vector[SIZE] mu;
  vector[SIZE] log_lik;

  for (i in 1:SIZE) {
    mu[i] = tau;
  }

  for (i in 1:SIZE) {
    log_lik[i] = normal_lpdf(Y[i] | mu[i], vs[i]);    
  }
}

