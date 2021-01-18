// copied with minor modifications from https://github.com/luisdamiano/stancon18
// CC-BY 4.0

functions {
  vector normalize(vector x) {
    return x / sum(x);
  }
}

data {
  int<lower=1> T;                   // number of observations (length)
  int<lower=1> K;                   // number of hidden states
  real x_t[T];                      // observations
  int<lower=0> est_sigma;           // flag, whether to estimate sigma (1) or use values passed in (0)
  real sigma_t[T];               // estimated sigma for each observation
}

parameters {
  // Discrete state model
  simplex[K] p_1k;                  // initial state probabilities
  simplex[K] A_ij[K];               // transition probabilities
                                    // A_ij[i][j] = p(z_t = j | z_{t-1} = i)
  // Continuous observation model
  ordered[K] mu_k;                  // observation means
  real<lower=0> sigma_k[K];         // observation standard deviations, optionally estimated if est_sigma == 1. Can the quantity K * est_sigma be used to dimension sigma_k?
}

transformed parameters {
  vector[K] unalpha_tk[T];

  { // Forward algorithm log p(z_t = j | x_{1:t})
    real accumulator[K];

    if(est_sigma == 1) {
      // use estimated sigma values
      unalpha_tk[1] = log(p_1k) + normal_lpdf(x_t[1] | mu_k, sigma_k);
    } else {
      // otherwise use values passed in by user, fixed
      unalpha_tk[1] = log(p_1k) + normal_lpdf(x_t[1] | mu_k, sigma_t[1]);
    }

    for (t in 2:T) {
      for (j in 1:K) { // j = current (t)
        for (i in 1:K) { // i = previous (t-1)
                         // Murphy (2012) Eq. 17.48
                         // belief state      + transition prob + local evidence at t
            if(est_sigma == 1) {
              // use estimated sigma values
              accumulator[i] = unalpha_tk[t-1, i] + log(A_ij[i, j]) + normal_lpdf(x_t[t] | mu_k[j], sigma_k[j]);
            } else {
              // otherwise use values passed in by user, fixed
              accumulator[i] = unalpha_tk[t-1, i] + log(A_ij[i, j]) + normal_lpdf(x_t[t] | mu_k[j], sigma_t[t]);
            }

        }
        unalpha_tk[t, j] = log_sum_exp(accumulator);
      }
    }
  } // Forward
}

model {
  sigma_k ~ student_t(3, 0, 1);
  mu_k ~ student_t(3, 0, 3);
  target += log_sum_exp(unalpha_tk[T]); // Note: update based only on last unalpha_tk
}

generated quantities {
  vector[K] unbeta_tk[T];
  vector[K] ungamma_tk[T];
  vector[K] alpha_tk[T];
  vector[K] beta_tk[T];
  vector[K] gamma_tk[T];
  vector[T] log_lik; // added to store log-likelihood for calculation of LOOIC
  int<lower=1, upper=K> zstar_t[T];
  real logp_zstar_t;

  { // Forward algortihm
    for (t in 1:T)
      alpha_tk[t] = softmax(unalpha_tk[t]);
  } // Forward

  { // Backward algorithm log p(x_{t+1:T} | z_t = j)
    real accumulator[K];

    for (j in 1:K)
      unbeta_tk[T, j] = 1;

    for (tforward in 0:(T-2)) {
      int t;
      t = T - tforward;

      for (j in 1:K) { // j = previous (t-1)
        for (i in 1:K) { // i = next (t)
                         // Murphy (2012) Eq. 17.58
                         // backwards t    + transition prob + local evidence at t
            if(est_sigma == 1) {
              accumulator[i] = unbeta_tk[t, i] + log(A_ij[j, i]) + normal_lpdf(x_t[t] | mu_k[i], sigma_k[i]);
            } else {
              accumulator[i] = unbeta_tk[t, i] + log(A_ij[j, i]) + normal_lpdf(x_t[t] | mu_k[i], sigma_t[t]);
            }

          }
        unbeta_tk[t-1, j] = log_sum_exp(accumulator);
      }
    }

    for (t in 1:T)
      beta_tk[t] = softmax(unbeta_tk[t]);
  } // Backward

  { // Forwards-backwards algorithm log p(z_t = j | x_{1:T})
    for(t in 1:T) {
        ungamma_tk[t] = alpha_tk[t] .* beta_tk[t];
        gamma_tk[t] = normalize(ungamma_tk[t]);
    }

    for(t in 1:T) {
      // gamma_tk is vector of normalized probability of state given all data, p(z_t = j | x_{1:T})

      log_lik[t] = 0; // initialize
      // log_lik accumulator. need to sum to integrate over states,
      // p(x_t) = p(x_t | z_t = 1) * p(z_t = 1)...
      // gamma_tk is p(x_t | z_t = k), alpha_tk is p(z_t = k | x[1:T])
      //if(est_sigma == 1) {
      for (j in 1:K) {
          log_lik[t] = log_lik[t] + gamma_tk[t,j]*alpha_tk[t,j];
      }
      //} else {
      //  for (j in 1:K)
      //    log_lik[t] = log_lik[t] + gamma_tk[t,j]*alpha_tk[t,j];
      //}
      log_lik[t] = log(log_lik[t]);
    }

  } // Forwards-backwards

  { // Viterbi algorithm
    int a_tk[T, K];                 // backpointer to the most likely previous state on the most probable path
    real delta_tk[T, K];            // max prob for the seq up to t
                                    // with final output from state k for time t
    if(est_sigma == 1) {
    for (j in 1:K)
      delta_tk[1, K] = normal_lpdf(x_t[1] | mu_k[j], sigma_k[j]);
    } else {
    for (j in 1:K)
      delta_tk[1, K] = normal_lpdf(x_t[1] | mu_k[j], sigma_t[1]);
    }

    for (t in 2:T) {
      for (j in 1:K) { // j = current (t)
        delta_tk[t, j] = negative_infinity();
        for (i in 1:K) { // i = previous (t-1)
          real logp;
          if(est_sigma == 1) {
          logp = delta_tk[t-1, i] + log(A_ij[i, j]) + normal_lpdf(x_t[t] | mu_k[j], sigma_k[j]);
          } else {
            logp = delta_tk[t-1, i] + log(A_ij[i, j]) + normal_lpdf(x_t[t] | mu_k[j], sigma_t[t]);
          }
          if (logp > delta_tk[t, j]) {
            a_tk[t, j] = i;
            delta_tk[t, j] = logp;
          }
        }
      }
    }

    logp_zstar_t = max(delta_tk[T]);

    for (j in 1:K)
      if (delta_tk[T, j] == logp_zstar_t)
        zstar_t[T] = j;

    for (t in 1:(T - 1)) {
      zstar_t[T - t] = a_tk[T - t + 1, zstar_t[T - t + 1]];
    }
  }
}
