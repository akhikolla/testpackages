//
// model 7:
//     Extended dixon and Simon
//
//     theta_g ~ N(theta, sigma^2)
//     theta = b0+b1X1+..+bpXp+c1X1X2+..+cX1Xp+..+zX1X2...Xp
//        a ~ N(0,1000)
//        b ~ N(0, omega[1]^2)
//        c ~ N(0, omega[2]^2)
//             ...
//     omega^2 ~ HalfN(0,1)

data {
  int<lower=0> SIZE;
  int<lower=0> NX;
  int<lower=0> NTAU;

  vector[SIZE]     Y;
	vector[SIZE]     SIGY;
  matrix[SIZE, NX] X;
  int<lower=0>     TAUINX[NX];

  real<lower=0> B;
	real<lower=0> D;
	real<lower=0> DELTA;
  int<lower=0, upper=1> PRIORSIG;
}

parameters {
  real b0;
  real<lower=0> omega[NTAU];
  vector<lower=0, upper=1>[SIZE] uvs;
	vector[SIZE] nvs;
  vector[NX] nomega;

}

transformed parameters{
	vector<lower=0>[SIZE] vs;
  vector[NX] bgamma;
	vector[SIZE] mu;

  if (0 == PRIORSIG) {
    vs = exp(log(SIGY) + (uvs * 2 - 1) * DELTA);
  } else {
    vs = exp(log(SIGY) + nvs * sqrt(DELTA));
  }

  for (i in 1:NX) {
    bgamma[i] = omega[TAUINX[i]] * nomega[i];
  }

  mu = b0+X*bgamma;
}

model {
  b0     ~ normal(0, sqrt(B));
  nomega ~ normal(0,1);
  uvs    ~ uniform(0,1);
  nvs    ~ normal(0,1);

  if (0 == D) {
    //jeffreys
    for (i in 1:NTAU) {
      target += -log(omega[i]);
    }
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
