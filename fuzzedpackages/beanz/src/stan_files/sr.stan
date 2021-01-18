//
// model 3:
//    simple regression
//
//     theta_g ~ N(tau, sigma^2)
//         tau = b0+b1*X1+..+bp*Xp
//          b  ~ N(0, 10^3)
//

data {
	int<lower=0>  SIZE;
	int<lower=0>  NX;
	vector[SIZE]  Y;
	vector[SIZE]  SIGY;
	matrix[SIZE, NX] X;

  real<lower=0> B;
	real<lower=0> C;
	real<lower=0> DELTA;
  int<lower=0, upper=1> PRIORSIG;
}


parameters {
	real       b0;
	vector[NX] bgamma;
 	vector<lower=0, upper=1>[SIZE] uvs;
	vector[SIZE] nvs;
}

transformed parameters{
	vector<lower=0>[SIZE] vs;
	vector[SIZE] mu;

  if (0 == PRIORSIG) {
    vs = exp(log(SIGY) + (uvs * 2 - 1) * DELTA);
  } else {
    vs = exp(log(SIGY) + nvs * sqrt(DELTA));
  }

  mu = b0+X*bgamma;
}

model {
	b0     ~ normal(0, sqrt(B));
	bgamma ~ normal(0, sqrt(C));
	uvs    ~ uniform(0,1);
  nvs    ~ normal(0,1);
	Y      ~ normal(mu, vs);
}

generated quantities {
  vector[SIZE] log_lik;
  for (i in 1:SIZE) {
    log_lik[i] = normal_lpdf(Y[i] | mu[i], vs[i]);    
  }
}
