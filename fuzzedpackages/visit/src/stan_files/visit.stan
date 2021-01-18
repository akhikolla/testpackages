functions {
  vector ftheta(real ptox, real pres, real rho) {
    vector[4] theta;
    if (1 == rho) {
      theta[1] = (1-ptox) * (1-pres); //00
      theta[2] = (1-ptox) * pres;     //01
      theta[3] = ptox     * (1-pres); //10
      theta[4] = ptox     * pres;     //11
    } else {
      theta[4] = -(sqrt((ptox+pres-ptox*rho-pres*rho-1)^2-4*(rho-1)*ptox*pres*rho)
                    + (ptox+pres-ptox*rho-pres*rho-1))/2/(rho-1);

      if (theta[4] < 0)
        theta[4] = 0;

      theta[3] = ptox - theta[4];
      if (theta[3] < 0)
        theta[3] = 0;

      theta[2] = pres - theta[4];
      if (theta[2] < 0)
        theta[2] = 0;

      theta[1] = rho*theta[2]*theta[3]/theta[4];
    }

    theta = theta*(1/sum(theta));
    return(theta);
  }
}

data {
  //total dose levels
  int<lower=0>  NDOSE;
  //benchmark tox rate for each level
  real<lower=0> TAU[NDOSE];
  //parameters for prior dist
  real          PAR[NDOSE,4]; //a,b,c,d
  //std for alpha prior
  real<lower=0> SDALPHA;

  //number of dose levels with observations
  int<lower=0> NOBSLEVEL;
  //observed data: 00 01 10 11
  int<lower=0> OBSY[NOBSLEVEL+1,4]; 

  //model selection
  int<lower=0>  SINGLEALPHA;
  int<lower=0>  SINGLERHO;
  int<lower=0>  FIXRHOAT1;
}

parameters {
  real          alpha[NDOSE];
  real<lower=0> rho[NDOSE];
  real<lower=0.0000001, upper=0.99999999> pres[NDOSE];
}

transformed parameters{
  real<lower=0.0000001, upper=0.99999999> ptox[NDOSE];
  vector<lower=0, upper=1>[4] theta[NDOSE];

  for (i in 1:NDOSE) {
    real r;
    real a;

    if (1 == FIXRHOAT1) {
      r = 1;
    } else if (1 == SINGLERHO) {
      r = rho[1];
    } else {
      r = rho[i];
    }

    if (1 == SINGLEALPHA) {
      a = alpha[1];
    } else {
      a = alpha[i];
    }

    ptox[i]  = pow(TAU[i], exp(a));
    theta[i] = ftheta(ptox[i], pres[i], r);
  }

}

model {
  //prior
  alpha[1] ~ normal(0, SDALPHA);
  for (i in 2:NDOSE) {
    alpha[i] ~ normal(alpha[i-1], SDALPHA);
  }

  for (i in 1:NDOSE) {
    pres[i] ~ beta(PAR[i,1], PAR[i,2]);
    rho[i]  ~ lognormal(PAR[i,3], PAR[i,4]);
  }

  //likelihood
  if (0 < NOBSLEVEL) {
    for (i in 1:NOBSLEVEL) {
      OBSY[i] ~ multinomial(theta[i]);
    }
  }
}

