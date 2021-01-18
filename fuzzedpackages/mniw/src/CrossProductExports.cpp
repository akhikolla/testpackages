/// @file CrossProductExports.cpp
///
/// @brief Exported Rcpp functions for matrix-weighted inner products.
///
/// A matrix-weighted inner product is of the form
/// \f[
/// \boldsymbol{W}_{q \times r} = \boldsymbol{X}_{q \times p}'\boldsymbol{V}_{p \times p} \boldsymbol{Y}_{p \times r},
/// \f]
/// where \f$\boldsymbol{V}\f$ is a symmetric positive-definite matrix.

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Eigen;
//#include <iostream>

//////////////////////////////////////////////////////////////////

/// Matrix-weighted self inner-product.
///
/// Calculate inner products of the form `t(X) V X` or `t(X) V^{-1} X`.  Both `X` and `V` can be vectorized, i.e., of size `N` or `1`, referred to here as `n`.
///
/// @param [in] X Matrix of size `p x qn`.
/// @param [in] V Matrix of size `p x pn`, where each `p x p` block is symmetric positive-definite.
/// @param [in] p Integer specifying the dimension of `V`.
/// @param [in] q Integer specifying the column dimension of `X`.
/// @param [in] inverse Boolean; whether to use `V` or `V^{-1}` for the matrix wweight.
///
/// @return Matrix of size `q x qn`, where each block of size `q x q` is an inner product.
//[[Rcpp::export]]
Eigen::MatrixXd CrossProdVXX(Eigen::MatrixXd X, Eigen::MatrixXd V,
			     int p, int q, bool inverse = false) {
  bool singleX = (X.cols() == q);
  bool singleV = (V.cols() == p);
  int N = singleX ? V.cols()/p : X.cols()/q;
  // output variables
  MatrixXd Y(q,q*N);
  // internal variables
  LLT<MatrixXd> lltp(p);
  MatrixXd VX(p,q);
  MatrixXd Xt(q,p);
  if(singleX) {
    Xt = X.adjoint();
  }
  if(inverse & singleV) {
    lltp.compute(V);
  }
  for(int ii=0; ii<N; ii++) {
    if(!inverse) {
      VX.noalias() = V.block(0,ii*p*(!singleV),p,p) * X.block(0,ii*q*(!singleX),p,q);
    }
    else {
      if(!singleV) {
	lltp.compute(V.block(0,ii*p,p,p));
      }
      VX = lltp.solve(X.block(0,ii*q*(!singleX),p,q));
    }
    if(singleX) {
      Y.block(0,ii*q,q,q).noalias() = Xt * VX;
    }
    else {
      Y.block(0,ii*q,q,q).noalias() = X.block(0,ii*q,p,q).adjoint() * VX;
    }
  }
  return(Y);
}

/// Matrix-weighted inner-product.
///
/// Calculate inner products of the form `t(X) V Y` or `t(X) V^{-1} Y`.  Each of `X`, `V`, and `Y` can be vectorized, i.e., of size `N` or `1`, referred to here as `n`.
///
/// @param [in] X Matrix of size `p x qn`.
/// @param [in] Y Matrix of size `p x rn`.
/// @param [in] V Matrix of size `p x pn`, where each `p x p` block is symmetric positive-definite.
/// @param [in] p Integer specifying the dimension of `V`.
/// @param [in] q Integer specifying the column dimension of `X`.
/// @param [in] r Integer specifying the column dimension of `Y`.
/// @param [in] inverse Boolean; whether to use `V` or `V^{-1}` for the matrix wweight.
///
/// @return Matrix of size `q x rn`, where each block of size `q x r` is an inner product.
//[[Rcpp::export]]
Eigen::MatrixXd CrossProdVXY(Eigen::MatrixXd X, Eigen::MatrixXd Y,
			     Eigen::MatrixXd V,
			     int p, int q, int r, bool inverse = false) {
  bool singleX = (X.cols() == q);
  bool singleY = (Y.cols() == r);
  bool singleV = (V.cols() == p);
  int N = singleX ? V.cols()/p : X.cols()/q;
  N = singleY ? N : Y.cols()/r;
  // output variables
  MatrixXd W(q,r*N);
  // internal variables
  LLT<MatrixXd> lltp(p);
  MatrixXd VY(p,r);
  MatrixXd Xt(q,p);
  // if both X,V or both Y,V are single, compute that product once
  if(inverse & singleV) {
    lltp.compute(V);
  }
  if(singleX) {
    if(singleV) {
      // X'V or X'V^{-1}
      if(!inverse) {
	Xt.noalias() = X.adjoint() * V;
      }
      else {
	Xt = lltp.solve(X).adjoint();
      }
    } else {
      // just transpose X
      Xt = X.adjoint();
    }
  } else if (singleY & singleV) {
    // VY or V^{-1}Y
    if(!inverse) {
      VY.noalias() = V * Y;
    }
    else {
      VY = lltp.solve(Y);
    }
  }
  //std::cout << "t(X) = " << Xt << std::endl;
  //std::cout << "V = " << V << std::endl;
  //std::cout << "Y = " << Y << std::endl;
  for(int ii=0; ii<N; ii++) {
    if(singleV & singleX) {
      // single X'V
      W.block(0,ii*r,q,r).noalias() = Xt * Y.block(0,ii*r*(!singleY),p,r);
    } else if(singleV & singleY) {
      // single VY
      W.block(0,ii*r,q,r).noalias() = X.block(0,ii*q*(!singleX),p,q).adjoint() * VY;
    } else {
      // V is not single, or neither X nor Y are single
      if(!inverse) {
	VY.noalias() = V.block(0,ii*p*(!singleV),p,p) * Y.block(0,ii*r*(!singleY),p,r);
      } else {
	if(!singleV) {
	  lltp.compute(V.block(0,ii*p,p,p));
	}
	VY = lltp.solve(Y.block(0,ii*r*(!singleY),p,r));
      }
      if(singleX) {
	W.block(0,ii*r,q,r).noalias() = Xt * VY;
      } else {
	W.block(0,ii*r,q,r).noalias() = X.block(0,ii*q,p,q).adjoint() * VY;
      }
    }
  }
  return(W);
}
