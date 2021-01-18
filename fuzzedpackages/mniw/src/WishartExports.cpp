/// @file WishartExports.cpp
///
/// @brief Exported Rcpp functions for the %Wishart distribution.

#include <Rcpp.h>
using namespace Rcpp;
#include <RcppEigen.h>
using namespace Eigen;
//#include <iostream>
#include "mniw/Wishart.h"
using namespace mniw;

// -----------------------------------------------------------------------------

/// Log-density of the Wishart and Inverse-Wishart distributions.
///
/// Evaluate the log-density of `N` observations of a `q x q` Wishart or Inverse-Wishart distribution.  Each argument except `inverse` can be vectorized, meaning that it can have length either `N` or `1`, denoted here as `n`.
///
/// @param [in] X Matrix of size `q x nq` of random variable observations.
/// @param [in] Psi Matrix of size `q x nq` of scale matrices.
/// @param [in] nu Vector of `n` degrees-of-freedom parameters.
/// @param [in] inverse Whether to evaluate the log-density of the Wishart or Inverse-Wishart distribution.
///
/// @return A vector of `N` log-density evaluations.
//[[Rcpp::export]]
Eigen::VectorXd LogDensityWishart(Eigen::MatrixXd X, Eigen::MatrixXd Psi,
				  Eigen::VectorXd nu, bool inverse = false) {
  int q = X.rows();
  int N = X.cols()/q;
  N = std::max<int>(N, Psi.cols()/q);
  N = std::max<int>(N, nu.size());
  // output variables
  VectorXd logDens(N);
  // internal variables
  Wishart wish(q);
  LLT<MatrixXd> cholX(q);
  LLT<MatrixXd> cholPsi(q);
  double ldPsi;
  double ldX;
  bool singleX = X.cols() == q;
  bool singlePsi = Psi.cols() == q;
  bool singleNu = nu.size() == 1;
  if(singleX) {
    cholX.compute(X);
    ldX = logDetCholV(cholX);
  }
  if(singlePsi) {
    cholPsi.compute(Psi);
    ldPsi = logDetCholV(cholPsi);
  }
  for(int ii=0; ii<N; ii++) {
    if(!singleX) {
      cholX.compute(X.block(0,ii*q,q,q));
      ldX = logDetCholV(cholX);
    }
    if(!singlePsi) {
      cholPsi.compute(Psi.block(0,ii*q,q,q));
      ldPsi = logDetCholV(cholPsi);
    }
    logDens(ii) = wish.LogDens(X.block(0,ii*(!singleX)*q,q,q), cholX, ldX,
			       Psi.block(0,ii*(!singlePsi)*q,q,q), cholPsi,
			       ldPsi, nu(ii*(!singleNu)), inverse);
  }
  return logDens;
}


/// Generate a random sample from the Wishart or Inverse-Wishart distribution.
///
/// Generate `N` independent draws from a `q x q` Wishart or Inverse-Wishart.  Each argument except `inverse` can be vectorized, meaning that it can have length either `N` or `1`, denoted here as `n`.
///
/// @param [in] N Integer number of random draws to produce
/// @param [in] Psi Matrix of size `q x nq` of scale matrices.
/// @param [in] nu Vector of `n` degrees-of-freedom parameters.
/// @param [in] inverse Whether to evaluate the log-density of the Wishart or Inverse-Wishart distribution.
///
/// @return A `q x Nq` matrix of `N` random draws.
//[[Rcpp::export]]
Eigen::MatrixXd GenerateWishart(int N, Eigen::MatrixXd Psi, Eigen::VectorXd nu,
				bool inverse = false) {
  int q = Psi.rows();
  bool singlePsi = (Psi.cols() == q);
  bool singleNu = (nu.size() == 1);
  // output variables
  MatrixXd X(q,N*q);
  // internal variables
  // TempPQ *tmp = new TempPQ(1,q);
  LLT<MatrixXd> lltq(q);
  MatrixXd Lq = MatrixXd::Zero(q,q);
  MatrixXd Uq = MatrixXd::Zero(q,q);
  MatrixXd Iq = MatrixXd::Identity(q,q);
  MatrixXd PsiL = MatrixXd::Zero(q,q);
  MatrixXd VL = MatrixXd::Zero(q,q);
  Wishart wish(q);
  // wishart lower triangular factor
  if(singlePsi) {
    if(!inverse) {
      // tmp->lltq.compute(Psi);
      // PsiL = tmp->lltq.matrixL();
      lltq.compute(Psi);
      PsiL = lltq.matrixL();
    }
    else {
      // ReverseCholesky(PsiL, Psi, tmp->lltq);
      ReverseCholesky(PsiL, Psi, lltq);
    }
  }
  // for-loop
  for(int ii=0; ii<N; ii++) {
    if(!inverse) {
      if(!singlePsi) {
	// tmp->lltq.compute(Psi.block(0,ii*q,q,q));
	// PsiL = tmp->lltq.matrixL();
	lltq.compute(Psi.block(0,ii*q,q,q));
	PsiL = lltq.matrixL();
      }
      // GenerateWishartLowerTri(VL, PsiL, nu(ii*(!singleNu)), tmp->Lq);
      // CrossProdLLt(X.block(0,ii*q,q,q), VL, tmp->Uq);
      wish.GenerateLowerTri(VL, PsiL, nu(ii*(!singleNu)));
      CrossProdLLt(X.block(0,ii*q,q,q), VL, Uq);
    }
    else {
      if(!singlePsi) {
	// ReverseCholesky(PsiL, Psi.block(0,ii*q,q,q), tmp->lltq);
	ReverseCholesky(PsiL, Psi.block(0,ii*q,q,q), lltq);
      }
      // GenerateWishartLowerTriXi(VL, PsiL, nu(ii*(!singleNu)));
      wish.GenerateLowerTriXi(VL, PsiL, nu(ii*(!singleNu)));
      // InverseLLt(X.block(0,ii*q,q,q), VL, tmp->Lq, tmp->Uq, tmp->Iq);
      // wish.GenerateLowerTriXi(VL, PsiL, nu(ii*(!singleNu)));
      InverseLLt(X.block(0,ii*q,q,q), VL, Lq, Uq, Iq);
    }
  }
  // delete tmp;
  return X;
}
