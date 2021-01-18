/// @file MatrixTExports.cpp
///
/// @brief Exported Rcpp functions for the Matrix-T distribution.

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Eigen;
#include "mniw/MatrixT.h"
using namespace mniw;

/// Log-density of the Matrix-T distribution.
///
/// Evaluate the log-density of `N` observations of a `p x q` dimensional Matrix-T distribution.  Each argument can be vectorized, meaning that it can have length either `N` or `1`, denoted here as `n`.
///
/// @param [in] X Matrix of `p x nq` random matrix observations.
/// @param [in] Lambda Matrix of `p x nq` mean matrices.
/// @param [in] SigmaR Matrix of `p x np` row-wise variance matrices.
/// @param [in] SigmaC Matrix of `q x nq` column-wise variance matrices.
/// @param [in] nu Vector of `n` degrees-of-freedom parameters.
///
/// @return Vector of `N` log-density evaluations.
//[[Rcpp::export]]
Eigen::VectorXd LogDensityMatrixT(Eigen::MatrixXd X,
				  Eigen::MatrixXd Lambda,
				  Eigen::MatrixXd SigmaR,
				  Eigen::MatrixXd SigmaC,
				  Eigen::VectorXd nu) {
  // dimensions of the problem
  int p = SigmaR.rows();
  int q = SigmaC.rows();
  int N = X.cols()/q;
  N = std::max<int>(N, Lambda.cols()/q);
  N = std::max<int>(N, SigmaR.cols()/p);
  N = std::max<int>(N, SigmaC.cols()/q);
  N = std::max<int>(N, nu.size());
  // output variables
  VectorXd logDens(N);
  // internal variables
  LLT<MatrixXd> cholSigmaR(p);
  LLT<MatrixXd> cholSigmaC(q);
  double ldSigmaR, ldSigmaC;
  bool singleX = (X.cols() == q);
  bool singleLambda = (Lambda.cols() == q);
  bool singleSigmaR = (SigmaR.cols() == p);
  bool singleSigmaC = (SigmaC.cols() == q);
  bool singleNu = (nu.size() == 1);
  MatrixT matT(p,q);
  // precompute Cholesky factors and log-determinants
  if(singleSigmaR) {
    cholSigmaR.compute(SigmaR);
    ldSigmaR = logDetCholV(cholSigmaR);
  }
  if(singleSigmaC) {
    cholSigmaC.compute(SigmaC);
    ldSigmaC = logDetCholV(cholSigmaC);
  }
  // main loop
  for(int ii=0; ii<N; ii++) {
    // compute Chol/log-det if not precomputed
    if(!singleSigmaR) {
      cholSigmaR.compute(SigmaR.block(0,p*ii,p,p));
      ldSigmaR = logDetCholV(cholSigmaR);
    }
    if(!singleSigmaC) {
      cholSigmaC.compute(SigmaC.block(0,q*ii,q,q));
      ldSigmaC = logDetCholV(cholSigmaC);
    }
    // density evaluation
    logDens(ii) = matT.LogDens(X.block(0,q*ii*(!singleX),p,q),
			       Lambda.block(0,q*ii*(!singleLambda),p,q),
			       SigmaR.block(0,p*ii*(!singleSigmaR),p,p),
			       cholSigmaR, ldSigmaR,
			       SigmaC.block(0,q*ii*(!singleSigmaC),q,q),
			       cholSigmaC, ldSigmaC,
			       nu(ii*(!singleNu)));
  }
  return logDens;
}

/// Generate a random sample from the Matrix-T distribution.
///
/// Generate `N` independent draws from a `p x q` dimensional Matrix-T distribution.  Each argument can be vectorized, meaning that it can have length either `N` or `1`, denoted here as `n`.
///
/// @param [in] N Integer number of random draws
/// @param [in] Lambda Matrix of `p x nq` mean matrices.
/// @param [in] SigmaR Matrix of `p x np` row-wise variance or precision matrices.
/// @param [in] SigmaC Matrix of `q x nq` column-wise variance matrices.
/// @param [in] nu Vector of `n` degrees-of-freedom parameters.
/// @param [in] inverse Boolean; `true/false` indicates that `SigmaR` is on the precision/variance scale.
/// @return Matrix of `p x Nq` random draws.
// [[Rcpp::export]]
Eigen::MatrixXd GenerateMatrixT(int N, Eigen::MatrixXd Lambda,
				Eigen::MatrixXd SigmaR,
				Eigen::MatrixXd SigmaC,
				Eigen::VectorXd nu,
				bool inverse = false) {
  // problem dimensions
  int p = SigmaR.rows();
  int q = SigmaC.rows();
  bool singleLambda = (Lambda.cols() == q);
  bool singleSigmaR = (SigmaR.cols() == p);
  bool singleSigmaC = (SigmaC.cols() == q);
  bool singleNu = (nu.size() == 1);
  // internal variables
  LLT<MatrixXd> lltq(q);
  LLT<MatrixXd> lltp(p);
  MatrixXd SigmaRL = MatrixXd::Zero(p,p);
  MatrixXd OmegaRU = MatrixXd::Zero(p,p);
  // MatrixXd CL = MatrixXd::Zero(q,q);
  MatrixXd OmegaCL = MatrixXd::Zero(q,q);
  MatrixT matT(p,q);
  // output variables
  MatrixXd X(p,N*q);
  // precomputations
  if(singleSigmaC) {
    ReverseCholesky(OmegaCL, SigmaC, lltq);
  }
  if(singleSigmaR) {
    lltp.compute(SigmaR);
    if(!inverse) {
      SigmaRL = lltp.matrixL();
    }
    else {
      OmegaRU = lltp.matrixU();
    }
  }
  // main loop
  for(int ii=0; ii<N; ii++) {
    if(!singleSigmaC) {
      ReverseCholesky(OmegaCL, SigmaC.block(0,ii*q,q,q), lltq);
    }
    if(!singleSigmaR) {
      lltp.compute(SigmaR.block(0,ii*p,p,p));
    }
    if(!inverse) {
      if(!singleSigmaR) {
	SigmaRL = lltp.matrixL();
      }
      matT.GenerateRowSColO(X.block(0,ii*q,p,q),
			       Lambda.block(0,ii*q*(!singleLambda),p,q),
			    SigmaRL, OmegaCL, nu(ii*(!singleNu)));
    }
    else {
      if(!singleSigmaR) {
	OmegaRU = lltp.matrixU();
      }
      matT.GenerateRowOColO(X.block(0,ii*q,p,q),
			    Lambda.block(0,ii*q*(!singleLambda),p,q),
			    OmegaRU, OmegaCL, nu(ii*(!singleNu)));
    }
  }
  return X;
}
