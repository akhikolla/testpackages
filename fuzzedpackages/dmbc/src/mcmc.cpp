// mcmc.cpp

#include "dmbc.h"

void dmbc_mcmc_binom(
  double* z_chain,
  double* alpha_chain,
  double* eta_chain,
  double* sigma2_chain,
  double* lambda_chain,
  double* prob_chain,
  double* x_chain,
  double* x_ind_chain,
  double* accept,
  double* loglik,
  double* logprior,
  double* logpost,
  int* Dm,
  double* z,
  int* x,
  int* ng,
  double* alpha,
  double* eta,
  double* sigma2,
  double* lambda,
  const double* hyper_eta_a,
  const double* hyper_eta_b,
  const double* hyper_lambda,
  const double gamma_z,
  const double gamma_alpha,
  const double hyper_sigma2_a,
  const double hyper_sigma2_b,
  int totiter,
  int n,
  int p,
  int S,
  int G,
  int verbose){
  int niter = 0, t = 0, ifix = 0;
  int m = n*(n - 1)/2;

  double* z_prop = new double[n*p];
  double* z_g = new double[n*p];
  double* z_old = new double[p];
  double* z_new = new double[p];
  int* ng_tmp = new int[G];
  for(int g = 0; g < G; g++){
    ng_tmp[g] = 0;
  } 
  int* Dm_g = new int[m*S];
  int* Dm_s = new int[m];
  double* delta = new double[m];
  double* pi_g = new double[m];
  double* alpha_plus_delta = new double[m];
  double* d_sigma2 = new double[G];

  double* accept_z_i = new double[G];
  double* accept_alpha = new double[G];
  double sum_z_i = 0, sum_alpha = 0;
  for(int g = 0; g < G; g++){
    accept[g] = 0;
    accept[g + G] = 0;
    accept_z_i[g] = 0;
    accept_alpha[g] = 0;
    d_sigma2[g] = NA_REAL;
  }
  int* iperm = new int[n];
  double* sigma = new double[p*p];
  double* mean_g = new double[p];
  double* sigma_g = new double[p*p];
  for(int j = 0; j < p; j++){
    mean_g[j] = 0;
    for(int i = 0; i < p; i++){
      sigma[i + j*p] = 0;
      sigma_g[i + j*p] = 0;
      if(i == j){
        sigma[i + p*j] = gamma_z*gamma_z;
      }
    }
  }
  double sum_z2 = 0;
  double* prob = new double[G];
  double maxprob = 0, prob_sum = 0;

  double* hyper_lambda_plus_ng = new double[G];

  double z_A = 0, z_B = 0, ran_unif = 0;
  double* dmn_z_new = new double;
  double* dmn_z_old = new double;
  double alpha_old = 0, alpha_new = 0;
  double alpha_A = 0, alpha_B = 0;
  double* loglik_z_prop = new double;
  double* loglik_z = new double;
  double* loglik_alpha_prop = new double;
  double* loglik_alpha = new double;
  double* lprior_z = new double[n];
  double lprior_alpha = 0;
  double* lprior_eta = new double;
  double* lprior_sigma2 = new double;
  double* lprior_lambda = new double;
  double sum_lprior_z = 0;

  GetRNGstate();

  while (niter < totiter){
    niter++;
    for(int g = 0; g < G; g++){
      // generate Z_g by random walk Metropolis-Hastings
      t = 0;
      tableC(ng_tmp, x, S, G);
      for(int s = 0; s < S; s++){
        if(x[s] == (g + 1)){
          for(int h = 0; h < m; h++){
            Dm_g[t] = Dm[h + m*s];
            t++;
          }
        }
      }
      for(int i = 0; i < n; i++){
        iperm[i] = i;  // or reshuffle randomly
      }
      for(int i = 0; i < n; i++){
        ifix = iperm[i];
        for(int j = 0; j < p; j++){
          z_old[j] = z[ifix + n*j + n*p*g];
          for(int k = 0; k < n; k++){
            z_prop[k + n*j] = z[k + n*j + n*p*g];
            z_g[k + n*j] = z[k + n*j + n*p*g];
          }
        }
        rmultinorm(z_new, 1, z_old, sigma, p);
        for(int j = 0; j < p; j++){
          z_prop[ifix + n*j] = z_new[j];
        }

        loglik_rbmds_binom(loglik_z_prop, Dm_g, z_prop, alpha[g], n, p, ng_tmp[g]);
        loglik_rbmds_binom(loglik_z, Dm_g, z_g, alpha[g], n, p, ng_tmp[g]);
        z_A = *loglik_z_prop - *loglik_z;
        for(int j = 0; j < p; j++){
          sigma_g[j + p*j] = eta[g];
        }
        dmultinorm(dmn_z_new, z_new, mean_g, sigma_g, 1, p, 1);
        dmultinorm(dmn_z_old, z_old, mean_g, sigma_g, 1, p, 1);
        z_B = *dmn_z_new - *dmn_z_old;

        ran_unif = R::runif(0, 1);
        if(ran_unif < exp(z_A + z_B)){
          for(int j = 0; j < p; j++){
            z_old[j] = z_new[j];
          }
          sum_z_i++;
          accept_z_i[g]++;
        }
        for(int j = 0; j < p; j++){
          z[ifix + n*j + n*p*g] = z_old[j];
        }
      }
      for(int i = 0; i < n; i++){
        for(int j = 0; j < p; j++){
          for(int h = 0; h < G; h++){
            z_chain[(niter - 1) + totiter*i + totiter*n*j + totiter*n*p*h] = z[i + n*j + n*p*h];
          }
          z_g[i + n*j] = z[i + n*j + n*p*g];
        }
      }

      // generate alpha_g by using random walk Metropolis-Hastings
      alpha_old = alpha[g];
      alpha_new = R::rnorm(alpha_old, gamma_alpha);

      loglik_rbmds_binom(loglik_alpha_prop, Dm_g, z_g, alpha_new, n, p, ng_tmp[g]);
      loglik_rbmds_binom(loglik_alpha, Dm_g, z_g, alpha_old, n, p, ng_tmp[g]);
      alpha_A = *loglik_alpha_prop - *loglik_alpha;
      alpha_B = R::dnorm(alpha_new, 0, sqrt(sigma2[g]), 1) - R::dnorm(alpha_old, 0, sqrt(sigma2[g]), 1);

      ran_unif = R::runif(0, 1);
      if(ran_unif < exp(alpha_A + alpha_B)){
        alpha_old = alpha_new;
        accept_alpha[g]++;
        sum_alpha++;
      }

      alpha[g] = alpha_old;
      alpha_chain[(niter - 1) + totiter*g] = alpha[g];

      // generate eta_g using its full conditional posterior distribution
      sum_z2 = 0;
      for(int j = 0; j < p; j++){
        for(int k = 0; k < n; k++){
          sum_z2 += z_g[k + n*j]*z_g[k + n*j];
        }
      }
      rinvgamma(&eta[g], 1, (hyper_eta_a[g] + n*p/2.0), (hyper_eta_b[g] + sum_z2/2.0));

      eta_chain[(niter - 1) + totiter*g] = eta[g];

      // generate sigma2_g using its full conditional posterior distribution
      rinvgamma(&sigma2[g], 1, (hyper_sigma2_a + 1), (hyper_sigma2_b + alpha[g]*alpha[g]/2.0));

      sigma2_chain[(niter - 1) + totiter*g] = sigma2[g];
    }

    // generate lambda using its full conditional posterior distribution
    for(int g = 0; g < G; g++){
      hyper_lambda_plus_ng[g] = hyper_lambda[g] + ng[g];
    }
    rdirichlet(lambda, 1, hyper_lambda_plus_ng, G);

    for(int g = 0; g < G; g++){
      lambda_chain[(niter - 1) + totiter*g] = lambda[g];
    }

    // generate x using its full conditional posterior distribution
    for(int s = 0; s < S; s++){
      for(int k = 0; k < m; k++){
        Dm_s[k] = Dm[k + m*s];
      }
      for(int g = 0; g < G; g++){
        for(int j = 0; j < p; j++){
          for(int i = 0; i < n; i++){
            z_g[i + n*j] = z[i + n*j + n*p*g];
          }
        }
        dist(delta, z_g, n, p);
        for(int i = 0; i < m; i++){
          alpha_plus_delta[i] = alpha[g] + delta[i];
        }
        expit(pi_g, alpha_plus_delta, m);
        dprodber(&prob[g], Dm_s, pi_g, m, 1);
        // if(any_na_nan(prob, G)){
        //  REprintf("### niter = %d  - Some NaN in prob[g]! ###\n", niter);
        //  REprintf("### niter = %d ###\n", niter);
        //  REprintf("prob[g] = %1.9f\n", prob[g]);
        //  for(int k = 0; k < m; k++){
        //    if((pi_g[k] == 0) || (pi_g[k] == 1)){
        //      REprintf("   k = %d - pi_g[k] = %1.9f\n", k, pi_g[k]);
        //      REprintf("   apd = %1.9f - a = %1.9f - d = %1.9f\n", alpha_plus_delta[k], alpha[g], delta[k]);
        //    }
        //  }
        // }
      }
      maxprob = *std::max_element(prob, prob + G);
      for(int g = 0; g < G; g++){
        prob[g] -= maxprob;
        prob[g] = lambda[g]*exp(prob[g]);
        prob_sum += prob[g];
      }
      for(int g = 0; g < G; g++){
        prob[g] /= prob_sum;
      }
      prob_sum = 0;
      for(int g = 0; g < G; g++){
        prob_chain[(niter - 1) + totiter*s + totiter*S*g] = prob[g];
      }

      sample_no_rep(G, prob, iperm, 1, &x[s]);
      x_chain[(niter - 1) + totiter*s] = x[s];
    }
    tableC(ng, x, S, G);
    for(int g = 0; g < G; g++){
      for(int s = 0; s < S; s++){
        x_ind_chain[(niter - 1) + totiter*s + totiter*S*g] = (x[s] == (g + 1)) ? 1 : 0;
      }
    }

    // print the information
    if(((niter % 500) == 0) && verbose){
      REprintf("   iteration %d/%d ==> acceptance z_i (avg): %1.4f - acceptance alpha (avg): %1.4f\n", niter, totiter, sum_z_i/(G*n*niter),
        sum_alpha/(G*niter));
    }

    // calculate the loglikelihood, logprior and logposterior for the sampled parameter values
    logprior[niter - 1] = 0;
    for(int g = 0; g < G; g++){
      for(int j = 0; j < p; j++){
        for(int i = 0; i < n; i++){
          z_g[i + n*j] = z[i + n*j + n*p*g];
        }
      }
      for(int j = 0; j < p; j++){
        sigma_g[j + p*j] = eta[g];
      }
      dmultinorm(lprior_z, z_g, mean_g, sigma_g, n, p, 1);
      sum_lprior_z = 0;
      for(int i = 0; i < n; i++){
        sum_lprior_z += lprior_z[i];
      }
      lprior_alpha = R::dnorm(alpha[g], 0, sqrt(sigma2[g]), 1);
      sum_z2 = 0;
      for(int j = 0; j < p; j++){
        for(int i = 0; i < n; i++){
          sum_z2 += z_g[i + n*j]*z_g[i + n*j];
        }
      }
      dinvgamma(lprior_eta, &eta[g], hyper_eta_a[g], hyper_eta_b[g], 1, 1);
      dinvgamma(lprior_sigma2, &sigma2[g], hyper_sigma2_a, hyper_sigma2_b, 1, 1);
      logprior[niter - 1] += sum_lprior_z + lprior_alpha + *lprior_eta + *lprior_sigma2;
    }
    ddirichlet(lprior_lambda, lambda, hyper_lambda, 1, G, 1);
    logprior[niter - 1] += *lprior_lambda;
    loglik_dmbc(&loglik[niter - 1], Dm, z, alpha, d_sigma2, lambda, x, n, p, S, G, "binomial");
    logpost[niter - 1] = loglik[niter - 1] + logprior[niter - 1];
    R_CheckUserInterrupt();
  }

  PutRNGstate();

  for(int g = 0; g < G; g++){
    accept[g] = accept_z_i[g]/(n*totiter);
    accept[g + G] = accept_alpha[g]/totiter;
  }

  delete[] z_prop;
  delete[] z_g;
  delete[] z_old;
  delete[] z_new;

  delete[] accept_z_i;
  delete[] accept_alpha;
  delete[] iperm;
  delete[] sigma;
  delete[] mean_g;
  delete[] sigma_g;
  delete[] prob;
  delete[] ng_tmp;
  delete[] Dm_g;
  delete[] Dm_s;
  delete[] delta;
  delete[] alpha_plus_delta;
  delete[] pi_g;
  delete[] hyper_lambda_plus_ng;
  delete[] d_sigma2;

  delete   dmn_z_new;
  delete   dmn_z_old;
  delete   loglik_z_prop;
  delete   loglik_z;
  delete   loglik_alpha_prop;
  delete   loglik_alpha;
  delete[] lprior_z;
  delete   lprior_eta;
  delete   lprior_sigma2;
  delete   lprior_lambda;
}
