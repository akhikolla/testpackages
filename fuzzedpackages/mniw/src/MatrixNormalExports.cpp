/// @file MatrixNormalExports.cpp
///
/// @brief Exported Rcpp functions for the Matrix-Normal distribution.

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Eigen;
//#include <iostream>
#include "mniw/MatrixNormal.h"
using namespace mniw;

/// Log-density of the Matrix-Normal distribution.
///
/// Evaluate the log-density of `N` observations of a `p x q` dimensional Matrix-Normal distribution.  Each argument can be vectorized, meaning that it can have length either `N` or `1`, denoted here as `n`.
///
/// @param [in] X Matrix of `p x nq` random matrix observations.
/// @param [in] Lambda Matrix of `p x nq` mean matrices.
/// @param [in] SigmaR Matrix of `p x np` row-wise variance matrices.
/// @param [in] SigmaC Matrix of `q x nq` column-wise variance matrices.
///
/// @return Vector of `N` log-density evaluations.
//[[Rcpp::export]]
Eigen::VectorXd LogDensityMatrixNormal(Eigen::MatrixXd X, Eigen::MatrixXd Lambda,
				       Eigen::MatrixXd SigmaR,
				       Eigen::MatrixXd SigmaC) {
  // dimensions of the problem
  int q = SigmaC.rows();
  int p = SigmaR.rows();
  int N = X.cols()/q;
  N = std::max<int>(N, Lambda.cols()/q);
  N = std::max<int>(N, SigmaR.cols()/p);
  N = std::max<int>(N, SigmaC.cols()/q);
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
  MatrixNormal matnorm(p,q);
  if(singleSigmaR) {
    cholSigmaR.compute(SigmaR);
    ldSigmaR = logDetCholV(cholSigmaR);
  }
  if(singleSigmaC) {
    cholSigmaC.compute(SigmaC);
    ldSigmaC = logDetCholV(cholSigmaC);
  }
  for(int ii=0; ii<N; ii++) {
    if(!singleSigmaR) {
      cholSigmaR.compute(SigmaR.block(0,p*ii,p,p));
      ldSigmaR = logDetCholV(cholSigmaR);
    }
    if(!singleSigmaC) {
      cholSigmaC.compute(SigmaC.block(0,q*ii,q,q));
      ldSigmaC = logDetCholV(cholSigmaC);
    }
    logDens(ii) = matnorm.LogDens(X.block(0,q*ii*(!singleX),p,q),
				  Lambda.block(0,q*ii*(!singleLambda),p,q),
				  cholSigmaR, ldSigmaR, cholSigmaC, ldSigmaC);
  }
  return logDens;
}

/// Generate a random sample from the Matrix-Normal distribution.
///
/// Generate `N` independent draws from a `p x q` dimensional Matrix-Normal distribution.  Each argument can be vectorized, meaning that it can have length either `N` or `1`, denoted here as `n`.
///
/// @param [in] N Integer number of random draws
/// @param [in] Lambda Matrix of `p x nq` mean matrices.
/// @param [in] SigmaR Matrix of `p x np` row-wise variance matrices.
/// @param [in] SigmaC Matrix of `q x nq` column-wise variance matrices.
/// @return Matrix of `p x Nq` random draws.
//[[Rcpp::export]]
Eigen::MatrixXd GenerateMatrixNormal(int N, Eigen::MatrixXd Lambda,
				     Eigen::MatrixXd SigmaR,
				     Eigen::MatrixXd SigmaC) {
  int p = SigmaR.rows();
  int q = SigmaC.rows();
  bool singleLambda = (Lambda.cols() == q);
  bool singleSigmaR = (SigmaR.cols() == p);
  bool singleSigmaC = (SigmaC.cols() == q);
  int ii;
  // output variables
  MatrixXd X(p,N*q);
  // internal variables
  MatrixXd SigmaRL = MatrixXd::Zero(p,p);
  MatrixXd SigmaCU = MatrixXd::Zero(q,q);
  LLT<MatrixXd> lltq(q);
  LLT<MatrixXd> lltp(p);
  MatrixNormal matnorm(p,q);
  if(singleSigmaR) {
    lltp.compute(SigmaR);
    SigmaRL = lltp.matrixL();
  }
  if(singleSigmaC) {
    lltq.compute(SigmaC);
    SigmaCU = lltq.matrixU();
  }
  for(ii=0; ii<N; ii++) {
    if(!singleSigmaR) {
      lltp.compute(SigmaR.block(0,ii*p,p,p));
      SigmaRL = lltp.matrixL();
    }
    if(!singleSigmaC) {
      lltq.compute(SigmaC.block(0,ii*q,q,q));
      SigmaCU = lltq.matrixU();
    }
    matnorm.GenerateRowSColS(X.block(0,ii*q,p,q),
			     Lambda.block(0,ii*q*(!singleLambda),p,q),
			     SigmaRL, SigmaCU);
  }
  return X;
}
