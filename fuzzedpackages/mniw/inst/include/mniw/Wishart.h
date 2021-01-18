/// @file Wishart.h
///
/// @brief Density evaluation and random number generation for the (Inverse) Wishart distribution.

#ifndef Wishart_h
#define Wishart_h 1

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "mniw/TriUtils.h"

namespace mniw {

#define gammaln R::lgammafn
#define chisq_rand Rf_rchisq
  using namespace Eigen;

  /// Multivariate log-gamma function.

  /// Computes the multivariate log-gamma function, which for real \f$\alpha\f$ and integer \f$q\f$ is defined as the logarithm of
  ///
  /// \f[
  /// \Gamma_q(\alpha) = \pi^{q(q-1)/4} \prod_{j=1}^q \Gamma(\alpha + (1-j)/2).
  /// \f]
  ///
  /// @param [in] alpha Real scalar.
  /// @param [in] q Integer scalar.
  /// @return The multivariate log-gamma function \f$\log \Gamma_q(\alpha)\f$.
  ///
  inline double logMultiGamma(double alpha, int q) {
    double lmg = 0.0;
    for(int ii=0; ii<q; ii++) {
      lmg += gammaln(alpha - .5*ii);
    }
    lmg += .5*q*(q-1) * M_LN_SQRT_PI;
    return lmg;
  }

  /// The %Wishart and Inverse %Wishart distributions.

  /// The %Wishart distribution is on a random positive-definite matrix \f$\bm{X}_{q\times q}\f$ is is denoted \f$X \sim \mathrm{Wish}(\Psi, \nu)\f$, and defined as \f$X = (L Z)(L Z)'\f$, where 
  /// 
  /// - \f$\Psi_{q\times q} = LL'\f$ is the positive-definite matrix scale parameter,
  /// - \f$\nu > q\f$ is the shape parameter, 
  /// 
  /// and \f$Z_{q\times q}\f$ is a random lower-triangular matrix with elements
  /// \f[
  ///   Z_{ij}
  ///   \begin{cases} 
  ///     \stackrel{\mathrm{iid}}{\sim} \mathcal N(0,1) & i < j \\
  ///     \stackrel{\mathrm{ind}}{\sim} \chi^2_{(\nu-i+1)} & i = j \\
  ///     = 0 & i > j.
  ///   \end{cases}
  /// \f]
  /// 
  /// The log-density of the %Wishart distribution is
  /// \f[
  /// \log p(X \mid \Psi, \nu) = -\tfrac{1}{2} \left[\mathrm{tr}(\Psi^{-1} X) + (q+1-\nu)\log |X| + \nu \log |\Psi|  + \nu q \log(2) + 2 \log \Gamma_q(\tfrac \nu 2)\right],
  /// \f]
  /// where  \f$\Gamma_n(x)\f$ is the multivariate Gamma function defined as
  /// \f[
  /// \Gamma_n(x) = \pi^{n(n-1)/4} \prod_{j=1}^n \Gamma\big(x + \tfrac 1 2 (1-j)\big).
  /// \f]
  /// 
  /// The Inverse-%Wishart distribution \f$X \sim \mathrm{InvWish}(\Psi, \nu)\f$ is defined as \f$X^{-1} \sim \mathrm{Wish}(\Psi^{-1}, \nu)\f$.  Its log-density is given by
  /// \f[
  /// \log p(X \mid \Psi, \nu) = -\tfrac 1 2 \left[\mathrm{tr}(\Psi X^{-1}) + (\nu+q+1) \log |X| - \nu \log |\Psi| + \nu q \log(2) + 2 \log \Gamma_q(\tfrac \nu 2)\right].
  /// \f]
  class Wishart {
  private:
    // storage
    int q_;
    LLT<MatrixXd> cholX_;
    LLT<MatrixXd> cholPsi_;
    MatrixXd Z_;
    MatrixXd XL_;
  public:
    /// Constructor
    Wishart(int q);
    /// Log-density
    double LogDens(const Ref<const MatrixXd>& X,
		   const Ref<const MatrixXd>& Psi,
		   double nu, bool inv);
    /// Log-density with pre-computations
    double LogDens(const Ref<const MatrixXd>& X,
		   LLT<MatrixXd>& cholX, double ldX,
		   const Ref<const MatrixXd>& Psi,
		   LLT<MatrixXd>& cholPsi, double ldPsi,
		   double nu, bool inv);
    /// Random draw with scale matrix input
    void GenerateLowerTri(Ref<MatrixXd> XL,
			  const Ref<const MatrixXd>& PsiL, double nu);
    /// Random draw with precision matrix input
    void GenerateLowerTriXi(Ref<MatrixXd> XL,
			    const Ref<const MatrixXd>& XiL, double nu);
  };

  /// @param [in] q Size of the Wishart/Inverse-Wishart distribution.
  inline Wishart::Wishart(int q) {
    q_ = q;
    cholX_.compute(MatrixXd::Identity(q_,q_));
    cholPsi_.compute(MatrixXd::Identity(q_,q_));
    Z_ = MatrixXd::Zero(q_,q_);
    XL_ = MatrixXd::Zero(q_,q_);
  }

  /// @param [in] X Observation matrix of size `q x q`.
  /// @param [in] Psi Scale matrix parameter of size `q x q`.
  /// @param [in] nu Shape parameter.
  /// @param [in] inv Whether or not the Inverse-Wishart distribution is desired.
  ///
  /// @return The log-density of the distribution.
  inline double Wishart::LogDens(const Ref<const MatrixXd>& X,
				 const Ref<const MatrixXd>& Psi,
				 double nu, bool inv) {
    cholX_.compute(X);
    cholPsi_.compute(Psi);
    return LogDens(X, cholX_, logDetCholV(cholX_),
		   Psi, cholPsi_, logDetCholV(cholPsi_), nu, inv);
  }

  /// Identical to Wishart::LogDens(const Ref<const MatrixXd>&,const Ref<const MatrixXd>&,double,bool), except with pre-computed Cholesky factors and log-determinants for `X` and `Psi` (faster in a for-loop where one of these is held fixed).
  ///
  /// @param [in] X Observation matrix of size `q x q`.
  /// @param [in] cholX Cholesky decomposition of a variance matrix `X`.  This is represented by an `Eigen::LLT` object, for which `cholX.compute()` has already been applied.
  /// @param [in] ldX Log-determinant of `cholX`.
  /// @param [in] Psi Scale matrix parameter of size `q x q`.
  /// @param [in] cholPsi Cholesky decomposition of the scale matrix parameter `Psi`.
  /// @param [in] ldPsi Log-determinant of `cholPsi`.
  /// @param [in] nu Shape parameter.
  /// @param [in] inv Whether or not the Inverse-Wishart distribution is desired.
  ///
  /// @return The log-density of the distribution.
  inline double Wishart::LogDens(const Ref<const MatrixXd>& X,
				 LLT<MatrixXd>& cholX, double ldX,
				 const Ref<const MatrixXd>& Psi,
				 LLT<MatrixXd>& cholPsi, double ldPsi,
				 double nu, bool inv) {
    double ldens;
    // trace
    if(!inv) {
      Z_.noalias() = cholPsi.solve(X);
    }
    else {
      Z_.noalias() = cholX.solve(Psi);
    }
    // std::cout << cholX.matrixLLT() << std::endl;
    // std::cout << Z_ << std::endl;
    // rest of calculation
    if(!inv) {
      ldens = (nu-q_-1) * ldX - nu * ldPsi;
    }
    else {
      ldens = -(nu+q_+1) * ldX + nu * ldPsi;
    }
    ldens -= .5*(Z_.trace() + q_*nu * M_LN2) + logMultiGamma(.5*nu, q_);
    return ldens;
  }

  /// Simulate a random draw from the lower Cholesky factor of a Wishart distribution with given scale matrix and shape parameters.
  ///
  /// @param [out] XL Matrix of size `q x q` containing the random draw, which lives only on the lower triangular half of the matrix.
  /// @param [in] PsiL Matrix of size `q x q` containing the lower Cholesky factor of the scale matrix parameter `Psi`.
  /// @param [in] nu Shape parameter.
  inline void Wishart::GenerateLowerTri(Ref<MatrixXd> XL,
					const Ref<const MatrixXd>& PsiL,
					double nu) {
    int ii, jj;
    for(ii=0; ii<q_; ii++) {
      // diagonal
      XL_(ii,ii) = sqrt(chisq_rand(nu-ii));
      // off-diagonals
      for(jj=0; jj<ii; jj++) {
	XL_(ii,jj) = norm_rand();
      }
    }
    // multiply by precision matrix XL = PsiL * XL
    triMultLL(XL, PsiL, XL_);
    return;
  }

  /// Simulate a random draw from the lower Cholesky factor of a Wishart distribution with given precision matrix and shape parameters. 
  ///
  /// @param [out] XL Matrix of size `q x q` containing the random draw, which lives only on the lower triangular half of the matrix.
  /// @param [in] XiL Matrix containing the inverse of the lower Cholesky factor of the scale matrix parameter `Psi`, namely `Psi^{-1} = XiL' * XiL`.
  /// @param [in] nu Shape parameter.
  inline void Wishart::GenerateLowerTriXi(Ref<MatrixXd> XL,
					  const Ref<const MatrixXd>& XiL,
					  double nu) {
    int ii, jj;
    for(ii=0; ii<q_; ii++) {
      // diagonal
      XL(ii,ii) = sqrt(chisq_rand(nu-ii));
      // off-diagonals
      for(jj=0; jj<ii; jj++) {
	XL(ii,jj) = norm_rand();
      }
    }
    // multiply by precision matrix
    triMultLiX(XL, XiL);
    return;
  }

} // namespace mniw
  
#endif
