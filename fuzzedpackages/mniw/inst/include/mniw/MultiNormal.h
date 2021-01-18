/// @file MultiNormal.h
///
/// @brief Density evaluation and random number generation for the Multivariate Normal distribution.

#ifndef MultiNormal_h
#define MultiNormal_h 1

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "mniw/TriUtils.h"

namespace mniw {

  using namespace Eigen;
  
  /// The Multivariate Normal distribution.
  class MultiNormal {
  private:
    int q_;
    LLT<MatrixXd> cholSigma_;
    MatrixXd SigmaL_;
    VectorXd z_;
    MatrixXd Z_;
  public:
    /// Constructor.
    MultiNormal(int q);
    /// Log-density with precomputations.
    double LogDens(const Ref<const VectorXd>& x,
		   const Ref<const VectorXd>& mu,
		   LLT<MatrixXd>& cholSigma, double ldSigma);
    /// Log-density.
    double LogDens(const Ref<const VectorXd>& x,
		   const Ref<const VectorXd>& mu,
		   const Ref<const MatrixXd>& Sigma);
    /// Random number generation with pre-computations.
    void Generate(Ref<VectorXd> x, const Ref<const VectorXd>& mu,
		  LLT<MatrixXd>& cholSigma, Ref<MatrixXd> SigmaL);
    /// Random number generation.
    void Generate(Ref<VectorXd> x, const Ref<const VectorXd>& mu,
		  const Ref<const MatrixXd>& Sigma);
  };

  /// @param [in] q Number of dimensions of the Multivariate Normal.
  inline MultiNormal::MultiNormal(int q) {
    q_ = q;
    z_ = VectorXd::Zero(q_);
    SigmaL_ = MatrixXd::Zero(q_,q_);
    cholSigma_.compute(MatrixXd::Identity(q_,q_));
  }

  /// Identical to the shorter version of `LogDens`, but with the Cholesky decomposition and its log-determinant pre-computed.
  ///
  /// @param [in] x Vector of observations.
  /// @param [in] mu Mean vector.
  /// @param [in] cholSigma Cholesky decomposition of the variance matrix, supplied as an `Eigen::LLT` object.
  /// @param [in] ldSigma Log-determinant of the Cholesky factor of the variance.
  ///
  /// @return The log-density evaluated at `x`.
  inline double MultiNormal::LogDens(const Ref<const VectorXd>& x,
				     const Ref<const VectorXd>& mu,
				     LLT<MatrixXd>& cholSigma, double ldSigma) {
    //double ldens;
    // double q = x.size();
    z_ = x-mu;
    // if(CalcSigma) {
    //   cholSigma.compute(Sigma);
    //   ldSigma = 0.0;
    //   for(int ii=0; ii<q; ii++) {
    //     ldSigma += log(cholSigma.matrixL()(ii,ii));
    //   }
    // }
    cholSigma.matrixL().solveInPlace(z_);
    return -(.5 * z_.squaredNorm() + ldSigma + q_ * M_LN_SQRT_2PI);
  }

  /// @param [in] x Vector of observations.
  /// @param [in] mu Mean vector.
  /// @param [in] Sigma Variance matrix.
  ///
  /// @return The log-density evaluated at `x`.
  inline double MultiNormal::LogDens(const Ref<const VectorXd>& x,
				     const Ref<const VectorXd>& mu,
				     const Ref<const MatrixXd>& Sigma) {
    cholSigma_.compute(Sigma);
    return LogDens(x, mu, cholSigma_, logDetCholV(cholSigma_));
  }

  /// @param [out] x Vector to which random draw is assigned.
  /// @param [in] mu Mean vector.
  /// @param [in] cholSigma Cholesky decomposition of the variance, supplied as an `Eigen::LLT` object.
  /// @param [in] SigmaL Dense representation of this object.  This is because `TriMultLX` does not accept `Eigen::TriangularView` objects and hopefully can be removed in the future...
  inline void MultiNormal::Generate(Ref<VectorXd> x,
				    const Ref<const VectorXd>& mu,
				    LLT<MatrixXd>& cholSigma,
				    Ref<MatrixXd> SigmaL) {
    for(int ii=0; ii<q_; ii++) {
      z_(ii) = norm_rand(); // generate iid standard normals
    }
    triMultLX(x, SigmaL, z_); // multiply by Sigma^{1/2}
    x += mu; // add mu
    return;
  }

  /// @param [out] x Vector to which random draw is assigned.
  /// @param [in] mu Mean vector.
  /// @param [in] Sigma Variance matrix.
  inline void MultiNormal::Generate(Ref<VectorXd> x,
				    const Ref<const VectorXd>& mu,
				    const Ref<const MatrixXd>& Sigma) {
    cholSigma_.compute(Sigma);
    SigmaL_ = cholSigma_.matrixL();
    Generate(x, mu, cholSigma_, SigmaL_);
  }

} // namespace mniw

#endif
