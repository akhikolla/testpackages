/// @file RandomEffects.h
///
/// @brief Multivariate Random-Effects Normal distribution.

#ifndef RandomEffects_h
#define RandomEffects_h 1

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "mniw/TriUtils.h"
// #include <iostream>

namespace mniw {

  using namespace Eigen;  
  
  /// The Multivariate Random-Effects Normal distribution.
  class RanfxNormal {
  private: // storage
    int q_;
    MatrixXd Mq_;
    MatrixXd Uq_;
    LLT<MatrixXd> lltq_;
  public:
    /// Constructor.
    RanfxNormal(int q);
    /// Random number generation with prior precision matrix on the Cholesky scale.
    void GenerateOL(Ref<VectorXd> mu, const Ref<const VectorXd>& x,
		    const Ref<const MatrixXd>& C,
		    const Ref<const MatrixXd>& OmegaL);
    /// Random number generation with prior precision matrix on the regular scale.
    void GenerateO(Ref<VectorXd> mu, const Ref<const VectorXd>& x,
		   const Ref<const MatrixXd>& C,
		   const Ref<const MatrixXd>& Omega);
  };

  /// @param [in] q Number of dimensions of the Multivariate Normal.
  inline RanfxNormal::RanfxNormal(int q) {
    q_ = q;
    Mq_ = MatrixXd::Zero(q_,q_);
    Uq_ = MatrixXd::Zero(q_,q_);
    lltq_.compute(MatrixXd::Identity(q_,q_));
  }

  /// @param [out] mu Vector in which to assign the random draw.
  /// @param [in] x Vector of observed data.
  /// @param [in] C Precision matrix of the observed data, `C = V^{-1}`.
  /// @param [in] OmegaL Cholesky factor of the prior precision matrix, `Sigma^{-1} = OmegaL * t(OmegaL)`.
  inline void RanfxNormal::GenerateOL(Ref<VectorXd> mu,
				      const Ref<const VectorXd>& x,
				      const Ref<const MatrixXd>& C,
				      const Ref<const MatrixXd>& OmegaL) {
    // std::cout << "mu = \n" << mu << std::endl;
    // std::cout << "x = \n" << x << std::endl;
    // std::cout << "C = \n" << C << std::endl;
    // std::cout << "OmegaL = \n" << OmegaL << std::endl;
    // int q = mu.size();
    // mu = Q * x
    CrossProdLLt(Mq_, OmegaL, Uq_);
    // std::cout << Mq_ << std::endl;
    mu.noalias() = Mq_ * x;
    // mu = G_L^{-1} z
    Mq_ += C;
    lltq_.compute(Mq_);
    lltq_.matrixL().solveInPlace(mu);
    // mu = mu + N(0,I)
    for(int ii=0; ii<q_; ii++) {
      mu[ii] += norm_rand();
    }
    // mu = G_U^{-1} z + y
    lltq_.matrixU().solveInPlace(mu);
    return;
  }

  /// @param [out] mu Vector in which to assign the random draw.
  /// @param [in] x Vector of observed data.
  /// @param [in] C Precision matrix of the observed data, `C = V^{-1}`.
  /// @param [in] Omega Prior precision matrix, `Omega = Sigma^{-1}`.
  inline void RanfxNormal::GenerateO(Ref<VectorXd> mu,
				     const Ref<const VectorXd>& x,
				     const Ref<const MatrixXd>& C,
				     const Ref<const MatrixXd>& Omega) {
    // int q = mu.size();
    // mu = Q * x
    mu.noalias() = Omega * x;
    // mu = G_L^{-1} z
    Mq_ = Omega + C;
    lltq_.compute(Mq_);
    lltq_.matrixL().solveInPlace(mu);
    // mu = mu + N(0,I)
    for(int ii=0; ii<q_; ii++) {
      mu[ii] += norm_rand();
    }
    // mu = G_U^{-1} z + y
    lltq_.matrixU().solveInPlace(mu);
    return;
  }

} // namespace mniw
  
#endif
