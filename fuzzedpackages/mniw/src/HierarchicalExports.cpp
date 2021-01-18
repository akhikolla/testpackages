/// @file HierarchicalExports.cpp
///
/// @brief Exported Rcpp functions for Hierarchical Normal Unequal-Variance model. 

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Eigen;
using namespace Rcpp;
//#include <iostream>
#include "mniw/TriUtils.h"
#include "mniw/Hierarchical.h"
using namespace mniw;

/// Bayesian inference for the Hierarchical Normal Unequal-Variance model.
///
/// Generate `nSamples` draws from a Gibbs sampler on the posterior random-effect and hyperparameter distribution of the Hierarchical Normal Unequal-Variance model.
///
/// @param [in] nSamples Integer number of posterior draws to return.
/// @param [in] nBurn Integer number of samples to initially discard as burn-in.
/// @param [in] Y Response matrix of size `N x q`, each row is an observation.
/// @param [in] X Covariate matrix of size `N x q`, each row is an observation.
/// @param [in] V Matrix of size `q x (Nq)` containing the variance of each response about its mean.
/// @param [in] Lambda Prior mean matrix of size `p x q`.
/// @param [in] Omega Prior row-wise precision matrix of size `p x p`.
/// @param [in] Psi Prior column-wise scale matrix of size `q x q`.
/// @param [in] nu Prior shape parameter.
/// @param [in] Beta0 Initial random-effects coefficient matrix of size `p x q`.
/// @param [in] iSigma0 Initial random-effects precision matrix of size `q x q`.
/// @param [in] Mu0 Initial random-effects matrix of size `N x q`.
/// @param [in] updateBetaSigma Boolean; whether or not to update the random-effects mean and variance matrices `Beta` and `Sigma`.
/// @param [in] updateMu Boolean; whether or not to update the random-effects `Mu`.
/// @param [in] storeBetaSigma Boolean; whether or not to return `Beta` and `Sigma`.
/// @param [in] storeMu Boolean; whether or not to return `Mu`.
///
/// @return `Rcpp::List` with elements `Beta`, `Sigma`, and `Mu`, being matrices of size `p x (q*nSamples)`, `q x (q*nSamples)` and `q x (N*nSamples)`.
// [[Rcpp::export]]
List HierUneqVModelGibbs(int nSamples, int nBurn,
			 Eigen::MatrixXd Y, Eigen::MatrixXd X,
			 Eigen::MatrixXd V,
			 Eigen::MatrixXd Lambda, Eigen::MatrixXd Omega,
			 Eigen::MatrixXd Psi, double nu,
			 Eigen::MatrixXd Beta0, Eigen::MatrixXd iSigma0,
			 Eigen::MatrixXd Mu0,
			 bool updateBetaSigma, bool updateMu,
			 bool storeBetaSigma, bool storeMu) {
  // dimensions of the problem
  int N = Y.rows();
  int p = X.cols();
  int q = Y.cols();
  // output variables
  MatrixXd BetaOut(p, (storeBetaSigma ? nSamples*q : q));
  MatrixXd SigmaOut(q, (storeBetaSigma ? nSamples*q: q));
  MatrixXd MuOut(q, (storeMu ? nSamples*N : N));
  if(!storeBetaSigma) {
    BetaOut.setZero();
    SigmaOut.setZero();
  }
  if(!storeMu) {
    MuOut.setZero();
  }
  // internal variables
  HierNormal hiernorm(Y, X, V, Lambda, Omega, Psi, nu);
  hiernorm.setMu(Mu0);
  hiernorm.setBeta(Beta0);
  hiernorm.setISigma(iSigma0);
  // Gibbs sampler
  for(int ii=-nBurn; ii<nSamples; ii++) {
    // Rprintf("iteration %i\n", ii);
    if(updateBetaSigma) {
      hiernorm.UpdateBetaSigma();
    }
    if(updateMu) {
      hiernorm.UpdateMu();
    }
    // storage
    if(ii >= 0) {
      if(storeBetaSigma) {
	hiernorm.getBeta(BetaOut.block(0,ii*q,p,q));
	hiernorm.getSigma(SigmaOut.block(0,ii*q,q,q));
      }
      if(storeMu) {
	hiernorm.getMut(MuOut.block(0,ii*N,q,N));
      }
    }
  }
  return List::create(_["Beta"] = BetaOut,
		      _["Sigma"] = SigmaOut,
		      _["Mu"] = MuOut);
}
