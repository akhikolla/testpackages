//
// model 4:
//    simple shrinkage
//
//     theta_g ~ N(b+phi_g, sigma^2)
//     phi_g ~ N(0, omega^2)
//
//     omega^2 ~ N(0,1)
//     
//

data {
	int<lower=0>  SIZE;
	vector[SIZE]  Y;
	vector[SIZE]  SIGY;
  real<lower=0> B;
	real<lower=0> D;
	real<lower=0> DELTA;
  int<lower=0, upper=1> PRIORSIG;
}

parameters {
  real b0;
  real<lower=0> omega;
 	vector<lower=0, upper=1>[SIZE] uvs;
	vector[SIZE] nvs;
  vector[SIZE] nphi;
}

transformed parameters{
	vector<lower=0>[SIZE] vs;
  vector[SIZE] phi;
	vector[SIZE] mu;

  if (0 == PRIORSIG) {
    vs = exp(log(SIGY) + (uvs * 2 - 1) * DELTA);
  } else {
    vs = exp(log(SIGY) + nvs * sqrt(DELTA));
  }

  // non-centralization
  phi = nphi * omega;
  mu  = b0 + phi;
}

model {
  b0    ~ normal(0, sqrt(B));
  nphi  ~ normal(0,1);
  uvs   ~ uniform(0,1);
  nvs   ~ normal(0,1);

  if (0 == D) {
    //jeffreys
    target += -log(omega);
  } else {
    //half normal
    omega ~ normal(0, sqrt(D));
  }

	Y ~ normal(mu, vs);
}

generated quantities {
  vector[SIZE] log_lik;
  for (i in 1:SIZE) {
    log_lik[i] = normal_lpdf(Y[i] | mu[i], vs[i]);    
  }
}
