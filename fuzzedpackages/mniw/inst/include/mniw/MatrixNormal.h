/// @file MatrixNormal.h
///
/// @brief Density evaluation and random number generation for the Matrix Normal distribution.

#ifndef MatrixNormal_h
#define MatrixNormal_h 1

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "mniw/TriUtils.h"

namespace mniw {

  using namespace Eigen;

  /// The Matrix Normal distribution.
  class MatrixNormal {
  private:
    // storage
    int p_;
    int q_;
    MatrixXd Z_;
    LLT<MatrixXd> cholSigmaR_;
    // double ldSigmaR;
    LLT<MatrixXd> cholSigmaC_;
    // double ldSigmaC;
  public:
    /// Constructor.
    MatrixNormal(int p, int q);
    /// Log-density evaluation.
    double LogDens(const Ref<const MatrixXd>& X, const Ref<const MatrixXd>& Lambda,
		   const Ref<const MatrixXd>& SigmaR, const Ref<const MatrixXd>& SigmaC);
    /// Log-density evaluation with precomputations.
    double LogDens(const Ref<const MatrixXd>& X,
		   const Ref<const MatrixXd>& Lambda,
		   LLT<MatrixXd>& cholSigmaR, double ldSigmaR, 
		   LLT<MatrixXd>& cholSigmaC, double ldSigmaC);
    /// Random draw with RowV/ColV on the variance/variance scale.
    void GenerateRowSColS(Ref<MatrixXd> X,
			  const Ref<const MatrixXd>& Lambda,
			  const Ref<const MatrixXd>& SigmaRL,
			  const Ref<const MatrixXd>& SigmaCU);
    /// Random number generation with RowV/ColV on the variance/precision scale.
    void GenerateRowSColO(Ref<MatrixXd> X,
			  const Ref<const MatrixXd>& Lambda,
			  const Ref<const MatrixXd>& SigmaRL,
			  const Ref<const MatrixXd>& OmegaCL);
    /// Random number generation with RowV/ColV on the precision/precision scale.
    void GenerateRowOColO(Ref<MatrixXd> X,
			  const Ref<const MatrixXd>& Lambda,
			  const Ref<const MatrixXd>& OmegaRU,
			  const Ref<const MatrixXd>& OmegaCL);
  };

  /// @param [in] p Number of rows of Matrix Normal distribution.
  /// @param [in] q Number of columns of Matrix Normal distribution.
  inline MatrixNormal::MatrixNormal(int p, int q) {
    p_ = p;
    q_ = q;
    Z_ = MatrixXd::Zero(p,q);
    cholSigmaR_.compute(MatrixXd::Identity(q_,q_));
    cholSigmaC_.compute(MatrixXd::Identity(p_,p_));
  }

  // --- log-density -------------------------------------------------------------

  /// @param [in] X Observation matrix of size `p x q`.
  /// @param [in] Lambda Mean matrix of size `p x q`.
  /// @param [in] SigmaR Row-variance matrix of size `p x p`.
  /// @param [in] SigmaC Column-variance matrix of size `q x q`.
  ///
  /// @return The log-density evaluated at `X`.
  inline double MatrixNormal::LogDens(const Ref<const MatrixXd>& X,
				      const Ref<const MatrixXd>& Lambda,
				      const Ref<const MatrixXd>& SigmaR,
				      const Ref<const MatrixXd>& SigmaC) {
    cholSigmaR_.compute(SigmaR);
    cholSigmaC_.compute(SigmaC);
    return LogDens(X, Lambda,
		   cholSigmaR_, logDetCholV(cholSigmaR_),
		   cholSigmaC_, logDetCholV(cholSigmaC_));
  }

  /// Identical to the shorter `LogDens`, except with pre-computed Cholesky factors and log-determinants for `SigmaR` and `SigmaC` (faster in a for-loop where one of these is held fixed).
  ///
  /// @param [in] X Observation matrix of size `p x q`.
  /// @param [in] Lambda Mean matrix of size `p x q`.
  /// @param [in] cholSigmaR Lower Cholesky factor of row-variance matrix, as an `Eigen::LLT` object.
  /// @param [in] ldSigmaR Log-determinant of `cholSigmaR`.
  /// @param [in] cholSigmaC Lower Cholesky factor of column-variance matrix, as an `Eigen::LLT` object.
  ///
  /// @param [in] ldSigmaC Log-determinant of `cholSigmaC`.
  ///
  /// @return The log-density evaluated as `X`.
  inline double MatrixNormal::LogDens(const Ref<const MatrixXd>& X,
				      const Ref<const MatrixXd>& Lambda,
				      LLT<MatrixXd>& cholSigmaR, double ldSigmaR, 
				      LLT<MatrixXd>& cholSigmaC, double ldSigmaC) {
    //double ldens;
    // double p = X.rows();
    // double q = X.cols();
    Z_ = X-Lambda;
    // if(CalcSigmaR) {
    //   cholSigmaR.compute(SigmaR);
    //   ldSigmaR = 0.0;
    //   for(int ii=0; ii<p; ii++) {
    //     ldSigmaR += log(cholSigmaR.matrixL()(ii,ii));
    //   }
    // }
    // if(CalcSigmaC) {
    //   cholSigmaC.compute(SigmaC);
    //   ldSigmaC = 0.0;
    //   for(int ii=0; ii<q; ii++) {
    //     ldSigmaC += log(cholSigmaC.matrixL()(ii,ii));
    //   }
    // }
    cholSigmaR.matrixL().solveInPlace(Z_);
    cholSigmaC.matrixU().solveInPlace<OnTheRight>(Z_);
    return -(.5 * Z_.squaredNorm() + q_ * ldSigmaR + p_ * ldSigmaC + p_*q_ * M_LN_SQRT_2PI);
  }

  // --- rng ---------------------------------------------------------------------

  /// @param [out] X Matrix of size `p x q` in which to store the random draw.
  /// @param [in] Lambda Mean matrix of size `p x q`.
  /// @param [in] SigmaRL Lower Cholesky factor of the row-variance matrix `SigmaR` (a matrix of size `p x p`).
  /// @param [in] SigmaCU Upper Cholesky factor of the column-variance matrix `SigmaC` (a matrix of size `q x q`).
  inline void MatrixNormal::GenerateRowSColS(Ref<MatrixXd> X,
					     const Ref<const MatrixXd>& Lambda,
					     const Ref<const MatrixXd>& SigmaRL,
					     const Ref<const MatrixXd>& SigmaCU) {
    int ii, jj;
    // int p = Lambda.rows();
    // int q = Lambda.cols();
    // populate Z_ with iid normals
    for(ii=0; ii<p_; ii++) {
      for(jj=0; jj<q_; jj++) {
	Z_(ii,jj) = norm_rand();
      }
    }
    // scaling
    triMultXU(X, Z_, SigmaCU);
    triMultLX(Z_, SigmaRL, X);
    X = Z_ + Lambda;
    return;
  }

  /// @param [out] X Matrix of size `p x q` in which to store the random draw.
  /// @param [in] Lambda Mean matrix of size `p x q`.
  /// @param [in] SigmaRL Lower Cholesky factor of the row-variance matrix `SigmaR` (a matrix of size `p x p`).
  /// @param [in] OmegaCL Lower Cholesky factor of the column-precision matrix `SigmaC^{-1}` (a matrix of size `q x q`).
  inline void MatrixNormal::GenerateRowSColO(Ref<MatrixXd> X,
					     const Ref<const MatrixXd>& Lambda,
					     const Ref<const MatrixXd>& SigmaRL,
					     const Ref<const MatrixXd>& OmegaCL) {
    int ii, jj;
    // int p = Lambda.rows();
    // int q = Lambda.cols();
    // populate Z_ with iid normals
    for(ii=0; ii<p_; ii++) {
      for(jj=0; jj<q_; jj++) {
	Z_(ii,jj) = norm_rand();
      }
    }
    // scaling
    triMultXLi(Z_, OmegaCL);
    triMultLX(X, SigmaRL, Z_);
    X += Lambda;
    return;
  }

  /// @param [out] X Matrix of size `p x q` in which to store the random draw.
  /// @param [in] Lambda Mean matrix of size `p x q`.
  /// @param [in] OmegaRU Upper Cholesky factor of the row-precision matrix `SigmaR^{-1}` (a matrix of size `p x p`).
  /// @param [in] OmegaCL Lower Cholesky factor of the column-precision matrix `SigmaC^{-1}` (a matrix of size `q x q`).
  inline void MatrixNormal::GenerateRowOColO(Ref<MatrixXd> X,
					     const Ref<const MatrixXd>& Lambda,
					     const Ref<const MatrixXd>& OmegaRU,
					     const Ref<const MatrixXd>& OmegaCL) {
    int ii, jj;
    // int p = Lambda.rows();
    // int q = Lambda.cols();
    // populate Z with iid normals
    for(ii=0; ii<p_; ii++) {
      for(jj=0; jj<q_; jj++) {
	Z_(ii,jj) = norm_rand();
      }
    }
    // scaling
    triMultXLi(Z_, OmegaCL);
    triMultUiX(OmegaRU, Z_);
    X = Z_ + Lambda;
    return;
  }

} // namespace mniw
#endif
