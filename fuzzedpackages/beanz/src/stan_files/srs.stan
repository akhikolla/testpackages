//
// model 5:
//    simple shrinkage+regression
//
//     theta_g  ~ N(tau, sigma^2)
//          tau = b0+b1X1+..+bpXp+phi_g
//     phi_g ~ N(0, omega^2)
//       beta~ N(0,1000)
//     tau^2 ~ HalfN(0,1)

data {
  int<lower=0>     SIZE;
	int<lower=0>     NX;
	vector[SIZE]     Y;
	vector[SIZE]     SIGY;
	matrix[SIZE, NX] X;

  real<lower=0> B;
  real<lower=0> C;
	real<lower=0> D;
	real<lower=0> DELTA;
  int<lower=0, upper=1> PRIORSIG;
}

parameters {
  real<lower=0> omega;
  real b0;
  vector[NX] bgamma;
  vector<lower=0, upper=1>[SIZE] uvs;
	vector[SIZE] nvs;
  vector[SIZE] nphi;
}

transformed parameters{
	vector<lower=0>[SIZE] vs;
	vector[SIZE] mu;
  vector[SIZE] phi;

  if (0 == PRIORSIG) {
    vs = exp(log(SIGY) + (uvs * 2 - 1) * DELTA);
  } else {
    vs = exp(log(SIGY) + nvs * sqrt(DELTA));
  }

  phi = nphi * omega;
  mu  = b0 + X * bgamma + phi;
}


model {
  b0      ~ normal(0, sqrt(B));
  bgamma  ~ normal(0, sqrt(C));
  nphi    ~ normal(0,1);
  uvs     ~ uniform(0,1);
  nvs     ~ normal(0,1);
  //phi     ~ normal(0, omega);

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
