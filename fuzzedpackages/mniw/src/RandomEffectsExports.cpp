/// @file RandomEffectsExports.cpp
///
/// @brief Exported Rcpp functions for the Multivariate Random-Effects Normal distribution.

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Eigen;
#include "mniw/TriUtils.h"
#include "mniw/RandomEffects.h"
using namespace mniw;

//////////////////////////////////////////////////////////////////

/// Generate a random sample from the Random-Effects Normal distribution.
///
/// Generate `N` independent draws from \f$p(\boldsymbol{\mu} \mid \boldsymbol{x})\f$, where
/// \f[
/// \begin{aligned}
/// \boldsymbol{\mu} & \sim \mathcal N(\boldsymbol{\lambda}, \boldsymbol{\Sigma}) \\
/// \boldsymbol{x} \mid \boldsymbol{\mu} & \sim \mathcal N(\boldsymbol{\mu}, \boldsymbol{V}).
/// \end{aligned}
/// \f]
/// Each argument can be vectorized, meaning that it can have length either `N` or `1`, denoted here as `n`.
///
/// @param [in] N Integer number of random draws.
/// @param [in] x Matrix of size `q x n` of observations.
/// @param [in] V Matrix of size `q x nq` of observation variances.
/// @param [in] lambda Matrix of size `q x n` of prior means.
/// @param [in] Sigma Matrix of size `q x nq` of prior variances.
///
/// @return Matrix of size `q x Nq` of random draws of observation means `mu`.
//[[Rcpp::export]]
Eigen::MatrixXd GenerateRandomEffectsNormal(int N,
					    Eigen::MatrixXd x,
					    Eigen::MatrixXd V,
					    Eigen::MatrixXd lambda,
					    Eigen::MatrixXd Sigma) {
  int q = lambda.rows();
  bool singleLambda = (lambda.cols() == 1);
  bool singleX = (x.cols() == 1);
  bool singleV = (V.cols() == q);
  bool singleSigma = (Sigma.cols() == q);
  // output variables
  MatrixXd mu(q,N);
  // internal variables
  // TempPQ *tmp = new TempPQ(1,q);
  LLT<MatrixXd> lltq(q);
  MatrixXd Iq = MatrixXd::Identity(q,q);
  MatrixXd C = MatrixXd::Zero(q,q);
  MatrixXd Omega = MatrixXd::Zero(q,q);
  VectorXd z = VectorXd::Zero(q);
  RanfxNormal renorm(q);
  if(singleV) {
    // tmp->lltq.compute(V);
    // C = tmp->lltq.solve(tmp->Iq);
    lltq.compute(V);
    C = lltq.solve(Iq);
  }
  if(singleSigma) {
    // tmp->lltq.compute(Sigma);
    // Omega = tmp->lltq.solve(tmp->Iq);
    lltq.compute(Sigma);
    Omega = lltq.solve(Iq);
  }
  if(singleLambda && singleX) {
    z = lambda - x;
  }
  for(int ii=0; ii<N; ii++) {
    if(!singleV) {
      // tmp->lltq.compute(V.block(0,ii*q,q,q));
      // C = tmp->lltq.solve(tmp->Iq);
      lltq.compute(V.block(0,ii*q,q,q));
      C = lltq.solve(Iq);
    }
    if(!singleSigma) {
      // tmp->lltq.compute(Sigma.block(0,ii*q,q,q));
      // Omega = tmp->lltq.solve(tmp->Iq);
      lltq.compute(Sigma.block(0,ii*q,q,q));
      Omega = lltq.solve(Iq);
    }
    if(!(singleLambda && singleX)) {
      z = lambda.col(ii*(!singleLambda)) - x.col(ii*(!singleX));
    }
    // GenerateRandomEffectsO(mu.col(ii), z, C, Omega, tmp);
    renorm.GenerateO(mu.col(ii), z, C, Omega);
    mu.col(ii) += x.col(ii*(!singleX));
  }
  // delete tmp;
  return(mu);
}


