/// @file MatrixNIWExports.cpp
///
/// @brief Exported Rcpp functions for the Matrix-Normal Inverse-Wishart distribution.

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Eigen;
using namespace Rcpp;
//#include <iostream>
#include "mniw/TriUtils.h"
#include "mniw/Wishart.h"
#include "mniw/MatrixNormal.h"
using namespace mniw;

//////////////////////////////////////////////////////////////////


/// Generate a random sample from the Matrix-Normal Inverse-Wishart distribution.
///
/// Generate `N` independent draws from a `p x q`/`q x q` dimensional Matrix-Normal Inverse-Wishart distribution.  Each argument can be vectorized, meaning that it can have length either `N` or `1`, denoted here as `n`.
///
/// @param [in] N Integer number of random draws
/// @param [in] Lambda Matrix of `p x nq` mean matrices.
/// @param [in] Sigma Matrix of `p x np` row-wise variance or precision matrices.
/// @param [in] Psi Matrix of `q x nq` scale matrices.
/// @param [in] nu Vector of `n` degrees-of-freedom parameters.
/// @param [in] inverse Boolean; `true/false` indicates that `Sigma` is on the precision/variance scale.
/// @return List with elements `X` and `V`, consisting of matrices of size `p x Nq` and `q x Nq` random draws respectively.
//[[Rcpp::export]]
List GenerateMatrixNIW(int N,
		       Eigen::MatrixXd Lambda, Eigen::MatrixXd Sigma,
		       Eigen::MatrixXd Psi, Eigen::VectorXd nu,
		       bool inverse = false) {
  // problem dimensions
  int p = Lambda.rows();
  int q = Psi.rows();
  bool singleLambda = (Lambda.cols() == q);
  bool singleSigma = (Sigma.cols() == p);
  bool singlePsi = (Psi.cols() == q);
  bool singleNu = (nu.size() == 1);
  // internal variables
  LLT<MatrixXd> lltq(q);
  MatrixXd Lq = MatrixXd::Zero(q,q);
  MatrixXd Uq = MatrixXd::Zero(q,q);
  MatrixXd Iq = MatrixXd::Identity(q,q);
  MatrixXd CL = MatrixXd::Zero(q,q);
  MatrixXd SigmaL = MatrixXd::Zero(p,p);
  MatrixXd OmegaU = MatrixXd::Zero(p,p);
  MatrixXd XiL = MatrixXd::Zero(q,q);
  LLT<MatrixXd> lltp(p);
  Wishart wish(q);
  MatrixNormal matnorm(p,q);
  // output variables
  MatrixXd X(p,N*q);
  MatrixXd V(q,N*q);
  // precomputations
  if(singlePsi) {
    ReverseCholesky(XiL, Psi, lltq);
  }
  if(singleSigma) {
    lltp.compute(Sigma);
    if(!inverse) {
      SigmaL = lltp.matrixL();
    }
    else {
      OmegaU = lltp.matrixU();
    }
  }
  // main loop
  for(int ii=0; ii<N; ii++) {
    if(!singlePsi) {
      ReverseCholesky(XiL, Psi.block(0,ii*q,q,q), lltq);
    }
    wish.GenerateLowerTriXi(CL, XiL, nu(ii*(!singleNu)));
    if(!singleSigma) {
      lltp.compute(Sigma.block(0,ii*p,p,p));
    }
    if(!inverse) {
      if(!singleSigma) {
	SigmaL = lltp.matrixL();
      }
      matnorm.GenerateRowSColO(X.block(0,ii*q,p,q),
			       Lambda.block(0,ii*q*(!singleLambda),p,q),
			       SigmaL, CL);
    }
    else {
      if(!singleSigma) {
	OmegaU = lltp.matrixU();
      }
      matnorm.GenerateRowOColO(X.block(0,ii*q,p,q),
			       Lambda.block(0,ii*q*(!singleLambda),p,q),
			       OmegaU, CL);
    }
    InverseLLt(V.block(0,ii*q,q,q), CL, Lq, Uq, Iq);
  }
  return List::create(_["X"] = wrap(X),
		      _["V"] = wrap(V));
}
