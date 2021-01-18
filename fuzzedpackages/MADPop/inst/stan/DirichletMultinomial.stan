// STAN model for the Dirichlet-Multinomial Distribution
// {mlysy,pwjkim}@uwaterloo.ca, january 2015
//
// model is:
// X_l | rho_l ~ Multinomial(N_l, rho_l)
// rho_l ~ Dirichlet(eta * alpha),
// where sum(alpha = 1).
// thus the prior expectation is E[rho_l] = alpha.
//
// for now, let's assume that eta is given, though ideally it will be estimated from the
// data.

// [[devtools.stan::export("hUM.mod")]]

functions {
  // unnormalized dirichlet-multinomial distribution
  // (for likelihood only)
  real dirichlet_multinomial_lpmf(int[] x, vector eta) {
    real ans;
    ans = 0.0;
    for(ii in 1:num_elements(x)) {
      ans += lgamma(x[ii] + eta[ii]) - lgamma(eta[ii]);
    }
    return ans + lgamma(sum(eta)) - lgamma(sum(x)+sum(eta));
  }
  // same thing but vectorized, i.e. accepts matrix X
  // current bug prevents overloading...
  real Dirichlet_Multinomial_lpmf(int[,] X, vector eta) {
    int D[2];
    real ans;
    real seta;
    real slgeta;
    D = dims(X);
    // eta constants
    seta = sum(eta);
    slgeta = 0.0;
    for(jj in 1:D[2]) {
      slgeta += lgamma(eta[jj]);
    }
    // all of X
    ans = 0.0;
    for(ii in 1:D[1]) {
      for(jj in 1:D[2]) {
	ans += lgamma(X[ii,jj] + eta[jj]);
      }
      ans -= lgamma(sum(X[ii])+seta);
    }
    ans += D[1] * (lgamma(seta) - slgeta);
    return ans;
  }
}

data {
  int<lower=1> nG; // number of observed genotype categories
  int<lower=1> nL; // number of lakes
  int<lower=0> X[nL,nG]; // vector of counts for each lake
  //real<lower=0> eta; // hyper prior precision parameter
  int<lower=0,upper=nL> nLrho; // number of lakes for which to generate samples from rho
  int<lower=1,upper=nL> iLrho[nLrho]; // indices of these lakes
}

parameters {
  simplex[nG] alpha; // probability parameters for dirichlet-multinomial
  real<lower=0> eta; // precision parameter
}

model {
  //vector<lower=0> kappa[nG]; // dirichlet parameters
  //kappa = eta * alpha;
  X ~ Dirichlet_Multinomial(eta * alpha);
  // avoid ridiculously small dirichlet variance
  // using p(1/sum(eta)) \propto 1
  //increment_log_prob(-2*log(1+eta));
  target += -2*log(1+eta);
}

generated quantities {
  // individual lake multinomial parameters
  simplex[nG] rho[nLrho];
  if(nLrho > 0) {
    for(ii in 1:nLrho) {
      rho[ii] = dirichlet_rng(to_vector(X[iLrho[ii]]) + eta * alpha);
    }
  }
}
