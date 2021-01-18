#include <cmath>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <rgen.h>
// [[Rcpp::depends(rgen)]]
#include "kalman.h"
#include "sample.h"
#include "distributions.h"
#include "ide_helpers.h"
#include "misc_helpers.h"

using namespace Rcpp;

// [[Rcpp::export]]
List eof(arma::mat Y, arma::mat F, arma::mat G_0, arma::mat Sigma_G_inv,
         arma::colvec m_0, arma::mat C_0, arma::mat scale_W,
         NumericVector params, CharacterVector proc_model,
         const int n_samples, const bool verbose) {
  
  // Extract scalar parameters
  const bool AR = proc_model(0) == "AR", DENSE = proc_model(0) == "Dense";
  const int P = G_0.n_rows, T = Y.n_cols, df_W = params["df_W"];
  const double alpha_sigma2 = params["alpha_sigma2"];
  const double beta_sigma2  = params["beta_sigma2"];
  const double alpha_lambda = params["alpha_lambda"];
  const double beta_lambda  = params["beta_lambda"];
  double sigma2_i = params["sigma2"];
  const bool sample_sigma2  = sigma2_i == NA, Discount = df_W == NA;
  
  // Create matrices and cubes for FFBS
  Y.insert_cols(0, 1); // make Y true-indexed; i.e. index 1 is t_1
  arma::cube theta(P, T+1, n_samples), R_inv(P, P, T+1), C(P, P, T+1), G, C_T;
  arma::mat a(P, T+1), m(P, T+1);
  m.col(0) = m_0;
  C.slice(0) = C_0;
  
  // Process matrix: G
  if (AR) {
    G_0 = G_0.diag();
    G = arma::zeros(P, P, n_samples+1);
    G.slice(0).diag() = G_0;
  } else if (DENSE) {
    G = arma::zeros(P, P, n_samples+1);
    G.slice(0) = G_0;
    G_0.reshape(P*P, 1);
  } else {
    G.set_size(P, P, 1);
    G.slice(0).eye();
  }
  
  // Find starting values for thetas
  theta.slice(0) = F.t() * Y;
  
  // Observation error: sigma2
  arma::vec sigma2;
  if (sample_sigma2) {
    sigma2.set_size(n_samples+1);
    sampleSigma2(sigma2.at(0), alpha_sigma2, beta_sigma2,
                 Y, F, theta.slice(0));
  }
  
  // Process error: W or lambda
  arma::vec lambda;
  arma::cube W;
  if (Discount) {
    lambda.set_size(n_samples+1);
    lambda.at(0) = rigamma(alpha_lambda, beta_lambda);
    C_T.set_size(P, P, n_samples);
  } 
  else { // Sample W from inverse-Wishart distribution
    W.set_size(P, P, n_samples+1);
    sampleW(W.slice(0), theta.slice(0), G.slice(0), scale_W, df_W);
  }
  
  // Sampling loop
  int G_idx = 0; // Incremented each iteration for AR and Dense models
  for (int i=0; i<n_samples; ++i) {
    checkUserInterrupt();
    
    // Set sigma2_i for FFBS
    if (sample_sigma2) sigma2_i = sigma2.at(i);
    
    // FFBS
    if (verbose) Rcout << "Filtering sample number " << i+1 << std::endl;
    if (Discount) {
      kalman(m, C, a, R_inv, Y, F, G.slice(G_idx), sigma2_i, lambda.at(i));
      C_T.slice(i) = C.slice(T); // Save for predictions
    } else {
      kalman(m, C, a, R_inv, Y, F, G.slice(G_idx), sigma2_i, -1, W.slice(i));
    }
    
    if (verbose) Rcout << "Drawing sample number " << i+1 << std::endl;
    backwardSample(theta.slice(i), m, a, C, G.slice(G_idx), R_inv);
    
    // Sample Sigma2
    if (sample_sigma2) {
      sampleSigma2(sigma2.at(i+1), alpha_sigma2, beta_sigma2,
                   Y, F, theta.slice(i));
    }
    
    // Sample W
    if (Discount) {
      sampleLambda(lambda.at(i+1), alpha_lambda, beta_lambda,
                   G.slice(G_idx), C, theta.slice(i));
    } else {
      sampleW(W.slice(i+1), theta.slice(i), G.slice(G_idx), scale_W, df_W);
    }
    
    // Sample G
    if (AR) {
      if (Discount) { // Update depends on W
        sampleAR(G.slice(G_idx+1), R_inv, theta.slice(i),
                 Sigma_G_inv, G_0, true, lambda.at(i+1));
      } else {
        sampleAR(G.slice(G_idx+1), W.slices(i+1, i+1), 
                 theta.slice(i), Sigma_G_inv, G_0);
      }
    } else if (DENSE) {
      if (Discount) { // Update depends on W
        sampleG(G.slice(G_idx+1), R_inv, theta.slice(i),
                Sigma_G_inv, G_0, true, lambda.at(i+1));
      } else {
        sampleG(G.slice(G_idx+1), W.slices(i+1, i+1), 
                theta.slice(i), Sigma_G_inv, G_0);
      }
    } // Nothing to update for RW case
    
    // Increment G_idx for AR and Dense models
    if (AR || DENSE) ++G_idx;
  }
  
  // Process results
  G.shed_slice(0);
  if (sample_sigma2) sigma2 = sigma2.subvec(1, n_samples);
  if (Discount) lambda = lambda.subvec(1, n_samples);
  
  List results;
  results["F"] = F; // Observation matrix
  results["theta"]  = theta; // State vector
  if (AR || DENSE) results["G"] = G; // Process matrix
  
  // Observation error
  if (sample_sigma2) results["sigma2"] = sigma2;
  else results["sigma2"] = sigma2_i;
  
  // Process error
  if (Discount) {
    results["C_T"] = C_T;
    results["lambda"] = lambda; 
  } else {
    W.shed_slice(0);
    results["W"] = W;
  }
  
  return results;
}

// [[Rcpp::export]]
List ide(arma::mat Y, arma::mat locs, arma::colvec m_0, arma::mat C_0,
         arma::mat mean_mu_kernel, arma::mat var_mu_kernel, arma::cube K,
         arma::cube scale_Sigma_kernel, arma::mat scale_W, NumericVector params, 
         const int n_samples, const bool verbose) {
  
  // Extract/set scalar parameters
  const double J = params["J"], L = params["L"], df_W = params["df_W"];
  const double alpha_sigma2 = params["alpha_sigma2"];
  const double beta_sigma2 = params["beta_sigma2"];
  const double alpha_lambda = params["alpha_lambda"];
  const double beta_lambda = params["beta_lambda"];
  const double proposal_factor_mu = params["proposal_factor_mu"];
  const double proposal_factor_Sigma = params["proposal_factor_Sigma"];
  const double df_Sigma_kernel = params["df_Sigma_kernel"];
  double sigma2_i = params["sigma2"], mh_ratio = 0.0;
  const int P = (2*J + 1) * (2*J + 1) , T = Y.n_cols, S = Y.n_rows;
  const int locs_dim = locs.n_cols, n_knots = K.n_cols;
  const int kernel_samples_per_iter = params["kernel_samples_per_iter"];
  int mu_acceptances = 0, Sigma_acceptances = 0;
  const bool sample_sigma2 = sigma2_i == NA, Discount = df_W == NA;
  const bool SV = params["SV"] > 0;
  const double Sigma_kernel_proposal_df = locs_dim + df_Sigma_kernel/proposal_factor_Sigma;
  const double Sigma_kernel_adjustment = Sigma_kernel_proposal_df - locs_dim - 1;
  
  // Create matrices and cubes for FFBS
  Y.insert_cols(0, 1); // make Y true-indexed; i.e. index 1 is t_1
  arma::cube theta(P, T+1, n_samples), G(P, P, n_samples+1);
  arma::cube R_inv(P, P, T+1), C(P, P, T+1), C_T;
  arma::mat a(P, T+1), m(P, T+1);
  m.col(0) = m_0;
  C.slice(0) = C_0;
  
  // Create objects for storing sampled mu_kernel and Sigma_kernel
  double u;
  
  // G
  arma::mat G_proposal, G_current;
  
  // mu_kernel
  arma::cube mu_kernel, mu_kernel_knots;
  arma::mat mu_kernel_proposal, mu_kernel_knots_proposal;
  arma::mat mu_kernel_current, mu_kernel_knots_current;
  arma::mat mu_kernel_proposal_var = proposal_factor_mu * proposal_factor_mu
                                                        * var_mu_kernel;
  
  // Sigma_kernel
  arma::cube Sigma_kernel_proposal, Sigma_kernel_knots_proposal;
  arma::cube Sigma_kernel_current, Sigma_kernel_knots_current;
  arma::field<arma::cube> Sigma_kernel(n_samples+1);
  arma::field<arma::cube> Sigma_kernel_knots(n_samples+1);
  
  if (SV) {
    // Set size of kernel parameter objects
    mu_kernel.set_size(S, locs_dim, n_samples+1);
    mu_kernel_knots.set_size(n_knots, locs_dim, n_samples+1);
    Sigma_kernel_proposal.set_size(locs_dim, locs_dim, S);
    Sigma_kernel_knots_proposal.set_size(locs_dim, locs_dim, n_knots);
    for (int i=0; i<=n_samples; ++i) {
      Sigma_kernel.at(i).set_size(locs_dim, locs_dim, S);
      Sigma_kernel_knots.at(i).set_size(locs_dim, locs_dim, n_knots);
    }
    
    // Set initial values
    mu_kernel_knots.slice(0) = mean_mu_kernel;
    mu_kernel.slice(0) = K.slice(0) * mu_kernel_knots.slice(0);
    Sigma_kernel_knots.at(0) = scale_Sigma_kernel / (df_Sigma_kernel-locs_dim-1);
    mapSigma(Sigma_kernel.at(0), Sigma_kernel_knots.at(0), K);
  }
  else {
    mu_kernel.set_size(locs_dim, 1, n_samples+1);
    mu_kernel.slice(0) = mean_mu_kernel;
    Sigma_kernel_proposal.set_size(locs_dim, locs_dim, 1);
    for (int i=0; i<=n_samples; ++i) {
      Sigma_kernel.at(i).set_size(locs_dim, locs_dim, 1);
    }
    Sigma_kernel.at(0) = scale_Sigma_kernel / (df_Sigma_kernel-locs_dim-1);
  }
  
  // Create observation matrix (F) and initial process matrix (G)
  arma::mat w_for_B = makeW(J, L);
  arma::mat F = makeF(locs, w_for_B, L);
  const arma::mat FtFiFt = arma::solve(F.t() * F, F.t());
  arma::mat B(S, P);
  makeB(B, mu_kernel.slice(0), Sigma_kernel.at(0), locs, w_for_B, L);
  G.slice(0) = FtFiFt * B;
  
  // Observation error
  arma::vec sigma2;
  if (sample_sigma2) {
    sigma2.set_size(n_samples+1);
    sigma2.at(0) = rigamma(alpha_sigma2, beta_sigma2);
  }
  
  // Process error
  arma::vec lambda;
  arma::cube W;
  if (Discount) {
    lambda.set_size(n_samples+1);
    lambda.at(0) = rigamma(alpha_lambda, beta_lambda);
    C_T.set_size(P, P, n_samples);
  } else { // Sample W from inverse-Wishart
    W.set_size(P, P, n_samples+1);
    W.slice(0) = rgen::riwishart(df_W, scale_W);
  }
  
  // Begin MCMC
  for (int i=0; i<n_samples; ++i) {
    checkUserInterrupt();
    
    // FFBS
    if (verbose) Rcout << "Filtering sample number " << i+1 << std::endl;
    if (sample_sigma2) sigma2_i = sigma2.at(i);
    
    if (Discount) {
      kalman(m, C, a, R_inv, Y, F, G.slice(i), sigma2_i, lambda.at(i));
      C_T.slice(i) = C.slice(T); // Save for predictions
    } else {
      kalman(m, C, a, R_inv, Y, F, G.slice(i), sigma2_i, -1, W.slice(i));
    }
    
    if (verbose) Rcout << "Drawing sample number " << i+1 << std::endl;
    backwardSample(theta.slice(i), m, a, C, G.slice(i), R_inv);
    
    // Sigma2
    if (sample_sigma2) {
      sampleSigma2(sigma2.at(i+1), alpha_sigma2, beta_sigma2, Y, F, theta.slice(i));
    }
    
    // Process error
    if (Discount) {
      sampleLambda(lambda.at(i+1), alpha_lambda, beta_lambda,
                   G.slice(i), C, theta.slice(i));
    } else {
      sampleW(W.slice(i+1), theta.slice(i), G.slice(i), scale_W, df_W);
    }
    
    // Begin MH step for kernel parameters
    
    // Set previous values for kernel parameters
    G_current = G.slice(i);
    mu_kernel_current = mu_kernel.slice(i);
    Sigma_kernel_current = Sigma_kernel.at(i);
    if (SV) {
      mu_kernel_knots_current = mu_kernel_knots.slice(i);
      Sigma_kernel_knots_current = Sigma_kernel_knots.at(i);
    }
    
    // Kernel parameter sampling loop (only saves last values)
    for (int j=0; j<kernel_samples_per_iter; ++j) {
      
      // MH step for mu
      // Propose mu and calculate probabilities under prior
      if (SV) {
        mu_kernel_knots_proposal = proposeMu(mu_kernel_knots_current, mu_kernel_proposal_var);
        mh_ratio  = ldmvnorm(mu_kernel_knots_proposal, mu_kernel_knots.slice(0),
                             var_mu_kernel);
        mh_ratio -= ldmvnorm(mu_kernel_knots_current, mu_kernel_knots.slice(0),
                             var_mu_kernel);
        mu_kernel_proposal = K.slice(0) * mu_kernel_knots_proposal; // map to spatial locations
      } else {
        mu_kernel_proposal = proposeMu(mu_kernel_current, mu_kernel_proposal_var);
        mh_ratio  = ldmvnorm(mu_kernel_proposal, mu_kernel.slice(0), var_mu_kernel);
        mh_ratio -= ldmvnorm(mu_kernel_current, mu_kernel.slice(0), var_mu_kernel);
      }
      
      // Calculate implied proposal for G and likelihood
      makeB(B, mu_kernel_proposal, Sigma_kernel_current, locs, w_for_B, L);
      G_proposal = FtFiFt * B; 
      
      if (Discount) {
        mh_ratio += kernelLikelihoodDiscount(G_proposal, theta.slice(i), C, lambda.at(i+1));
        mh_ratio -= kernelLikelihoodDiscount(G_current, theta.slice(i), C, lambda.at(i+1));
      } else {
        mh_ratio += kernelLikelihood(G_proposal, theta.slice(i), W.slice(i+1));
        mh_ratio -= kernelLikelihood(G_current, theta.slice(i), W.slice(i+1));
      }
      
      // Accept according to mh-ratio
      u = R::runif(0, 1);
      if (log(u) < mh_ratio) {
        ++mu_acceptances;
        G_current = G_proposal;
        mu_kernel_current = mu_kernel_proposal;
        if (SV) mu_kernel_knots_current = mu_kernel_knots_proposal;
      } // Otherwise, no need to update current values
      
      // MH step for Sigma
      if (SV) {
        // Propose Sigma at all knot locations
        for (int k=0; k<n_knots; ++k) {
          Sigma_kernel_knots_proposal.slice(k) = rgen::riwishart(
                                                   Sigma_kernel_proposal_df,
                                                   Sigma_kernel_knots_current.slice(k) *
                                                   Sigma_kernel_adjustment);
        }
        // Probabilities under prior
        mh_ratio  = ldiwishart(Sigma_kernel_knots_proposal, df_Sigma_kernel, 
                               scale_Sigma_kernel);
        mh_ratio -= ldiwishart(Sigma_kernel_knots_current, df_Sigma_kernel,
                               scale_Sigma_kernel);
        
        // Transition probabilities
        mh_ratio -= ldiwishart(Sigma_kernel_knots_proposal, Sigma_kernel_proposal_df,
                               Sigma_kernel_knots_current * Sigma_kernel_adjustment);
        mh_ratio += ldiwishart(Sigma_kernel_knots_current, Sigma_kernel_proposal_df,
                               Sigma_kernel_knots_proposal * Sigma_kernel_adjustment);
        
        // Map to all spatial locations
        mapSigma(Sigma_kernel_proposal, Sigma_kernel_knots_proposal, K.slice(0));
      } else {
        // Propose new Sigma
        Sigma_kernel_proposal.slice(0) = rgen::riwishart(Sigma_kernel_proposal_df,
                                                         Sigma_kernel_current.slice(0) * 
                                                         Sigma_kernel_adjustment);
        // Probabilities under prior
        mh_ratio  = ldiwishart(Sigma_kernel_proposal, df_Sigma_kernel, 
                               scale_Sigma_kernel);
        mh_ratio -= ldiwishart(Sigma_kernel_current, df_Sigma_kernel,
                               scale_Sigma_kernel);
        
        // Transition probabilities
        mh_ratio -= ldiwishart(Sigma_kernel_proposal, Sigma_kernel_proposal_df,
                               Sigma_kernel_current * Sigma_kernel_adjustment);
        mh_ratio += ldiwishart(Sigma_kernel_current, Sigma_kernel_proposal_df,
                               Sigma_kernel_proposal * Sigma_kernel_adjustment);
      }
      
      // Calculate likelihood
      makeB(B, mu_kernel_current, Sigma_kernel_proposal, locs, w_for_B, L);
      G_proposal = FtFiFt * B; 
      
      if (Discount) {
        mh_ratio += kernelLikelihoodDiscount(G_proposal, theta.slice(i), C, lambda.at(i+1));
        mh_ratio -= kernelLikelihoodDiscount(G_current, theta.slice(i), C, lambda.at(i+1));
      } else {
        mh_ratio += kernelLikelihood(G_proposal, theta.slice(i), W.slice(i+1));
        mh_ratio -= kernelLikelihood(G_current, theta.slice(i), W.slice(i+1));
      }
      
      // Accept according to mh-ratio
      u = R::runif(0, 1);
      if (log(u) < mh_ratio) {
        ++Sigma_acceptances;
        G_current = G_proposal;
        Sigma_kernel_current = Sigma_kernel_proposal;
        if (SV) Sigma_kernel_knots_current = Sigma_kernel_knots_proposal;
      } // Otherwise, no need to update current values
      // end iteration of kernel parameter sampling loop
    } // end kernel parameter sampling loop
    
    G.slice(i+1) = G_current;
    mu_kernel.slice(i+1) = mu_kernel_current;
    Sigma_kernel.at(i+1) = Sigma_kernel_current;
    if (SV) {
      mu_kernel_knots.slice(i+1) = mu_kernel_knots_current;
      Sigma_kernel_knots.at(i+1) = Sigma_kernel_knots_current;
    }
  } // end sampling loop
  
  // Drop starting values
  G.shed_slice(0);
  mu_kernel.shed_slice(0);
  if (sample_sigma2) sigma2 = sigma2.subvec(1, n_samples);
  if (Discount) lambda = lambda.subvec(1, n_samples);
  
  // Last element of Sigma_kernel is removed in dstm_ide in dstm.R
  
  // Save results to list
  List results;
  results["theta"]  = theta;
  results["F"] = F;
  results["G"] = G;
  results["mu_kernel"] = mu_kernel;
  results["Sigma_kernel"] = Sigma_kernel;
  if (SV) {
    results["mu_kernel_knots"] = mu_kernel_knots;
    results["Sigma_kernel_knots"] = Sigma_kernel_knots;
    results["process_convolution_map"] = K;
  }
  
  if (Discount) {
    results["lambda"] = lambda;
    results["C_T"] = C_T;
  }
  else {
    W.shed_slice(0);
    results["W"] = W;
  }
  
  if (sample_sigma2) results["sigma2"] = sigma2;
  else results["sigma2"] = sigma2_i;
  
  results["mu_acceptances"] = mu_acceptances;
  results["Sigma_acceptances"] = Sigma_acceptances;
  
  return results;
}
