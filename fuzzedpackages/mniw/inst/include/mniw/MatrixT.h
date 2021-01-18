/// @file MatrixT.h
///
/// @brief Density evaluation and random number generation for the Matrix-t distribution.

#ifndef MatrixT_h
#define MatrixT_h 1

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "mniw/TriUtils.h"
#include "mniw/Wishart.h"
#include "mniw/MatrixNormal.h"

namespace mniw {

  using namespace Eigen;

  /// The Matrix-t distribution.
  class MatrixT {
  private:
    // storage
    int p_;
    int q_;
    int pq_; // p * q
    bool pMin_; // isTRUE(p < q)
    MatrixXd Z_;
    LLT<MatrixXd> cholSigmaR_;
    LLT<MatrixXd> cholSigmaC_;
    LLT<MatrixXd> lltd_;
    MatrixXd A_;
    MatrixXd B_;
    MatrixXd CL_;
    Wishart *wish_;
    MatrixNormal *matnorm_;
  public:
    /// Constructor.
    MatrixT(int p, int q);
    /// Destructor.
    ~MatrixT();
    /// Log-density evaluation.
    double LogDens(const Ref<const MatrixXd>&  X,
		   const Ref<const MatrixXd>& Mu,
		   const Ref<const MatrixXd>& SigmaR,
		   const Ref<const MatrixXd>& SigmaC,
		   double nu);
    /// Log-density evaluation with precomputations.
    double LogDens(const Ref<const MatrixXd>& X,
		   const Ref<const MatrixXd>& Mu,
		   const Ref<const MatrixXd>& SigmaR,
		   LLT<MatrixXd>& cholSigmaR, double ldSigmaR,
		   const Ref<const MatrixXd>& SigmaC,
		   LLT<MatrixXd>& cholSigmaC, double ldSigmaC, double nu);
    /// Random draw with RowV/ColV on the variance/precision scale.
    void GenerateRowSColO(Ref<MatrixXd> X,
			  const Ref<const MatrixXd>& Mu,
			  const Ref<const MatrixXd>& SigmaRL,
			  const Ref<const MatrixXd>& OmegaCL, double nu);
    /// Random draw with RowV/ColV on the precision/precision scale.
    void GenerateRowOColO(Ref<MatrixXd> X,
			  const Ref<const MatrixXd>& Mu,
			  const Ref<const MatrixXd>& OmegaRU,
			  const Ref<const MatrixXd>& OmegaCL, double nu);
  };

  /// @param [in] p Number of rows of Matrix-t distribution.
  /// @param [in] q Number of columns of Matrix-t distribution.
  inline MatrixT::MatrixT(int p, int q) {
    // problem dimensions
    p_ = p;
    q_ = q;
    // precomputations
    pMin_ = p_ < q_;
    pq_ = p_*q_; 
    // storage
    Z_ = MatrixXd::Zero(p_,q_);
    cholSigmaR_.compute(MatrixXd::Identity(p_,p_));
    cholSigmaC_.compute(MatrixXd::Identity(q_,q_));
    if(pMin_) {
      A_ = MatrixXd::Zero(q_,p_);
      B_ = MatrixXd::Zero(p_,p_);
      lltd_.compute(MatrixXd::Identity(p_,p_));
    }
    else {
      A_ = MatrixXd::Zero(p_,q_);
      B_ = MatrixXd::Zero(q_,q_);
      lltd_.compute(MatrixXd::Identity(q_,q_));
    }
    CL_ = MatrixXd::Zero(q_,q_);
    wish_ = new Wishart(q_);
    matnorm_ = new MatrixNormal(p_,q_);
  }

  inline MatrixT::~MatrixT() {
    delete wish_;
    delete matnorm_;
  }

  // --- log-density --------------------------------------------------------

  /// @param [in] X Observation matrix of size `p x q`.
  /// @param [in] Mu Mean matrix of size `p x q`.
  /// @param [in] SigmaR Row-variance matrix of size `p x p`.
  /// @param [in] SigmaC Column-variance matrix of size `q x q`.
  /// @param [in] nu Shape parameter.
  ///
  /// @return The log-density evaluated at `X`.
  inline double MatrixT::LogDens(const Ref<const MatrixXd>&  X,
				 const Ref<const MatrixXd>& Mu,
				 const Ref<const MatrixXd>& SigmaR,
				 const Ref<const MatrixXd>& SigmaC,
				 double nu) {
    cholSigmaR_.compute(SigmaR);
    cholSigmaC_.compute(SigmaC);
    return LogDens(X, Mu,
		   SigmaR, cholSigmaR_, logDetCholV(cholSigmaR_),
		   SigmaC, cholSigmaC_, logDetCholV(cholSigmaC_), nu);
  }

  /// Identical to the shorter `LogDens`, except with pre-computed Cholesky factor and log-determinants for `SigmaR` and `SigmaC` (faster in a for-loop where one of these is held fixed).
  ///
  /// @param [in] X Observation matrix of size `p x q`.
  /// @param [in] Mu Mean matrix of size `p x q`.
  /// @param [in] SigmaR Row-variance matrix of size `p x p`.
  /// @param [in] cholSigmaR Lower Cholesky factor of row-variance matrix, as an `Eigen::LLT` object.
  /// @param [in] ldSigmaR Log-determinant of `cholSigmaR`.
  /// @param [in] SigmaC Column-variance matrix of size `q x q`.
  /// @param [in] cholSigmaC Lower Cholesky factor of column-variance matrix, as an `Eigen::LLT` object.
  /// @param [in] ldSigmaC Log-determinant of `cholSigmaC`.
  /// @param [in] nu Shape parameter.
  ///
  /// @return The log-density evaluated at `X`.
  inline double MatrixT::LogDens(const Ref<const MatrixXd>& X,
				 const Ref<const MatrixXd>& Mu,
				 const Ref<const MatrixXd>& SigmaR,
				 LLT<MatrixXd>& cholSigmaR, double ldSigmaR,
				 const Ref<const MatrixXd>& SigmaC,
				 LLT<MatrixXd>& cholSigmaC, double ldSigmaC,
				 double nu) {
    double ldens;
    double nuq = nu + q_ - 1.0;
    double nupq = nuq + p_;
    // calculate log-determinant by Sylvester formula
    Z_.noalias() = X - Mu;
    if(pMin_) {
      A_ = cholSigmaC.solve(Z_.adjoint());
      B_.noalias() = SigmaR + Z_ * A_;
      lltd_.compute(B_);
      ldens = logDetCholV(lltd_) - ldSigmaR;
    }
    else {
      A_ = cholSigmaR.solve(Z_);
      B_.noalias() = SigmaC + Z_.adjoint() * A_;
      lltd_.compute(B_);
      ldens = logDetCholV(lltd_) - ldSigmaC;
    }
    // add components
    ldens = nupq * ldens + q_ * ldSigmaR + p_ * ldSigmaC + pq_ * M_LN_SQRT_PI;
    return -ldens + logMultiGamma(.5 * nupq, q_) - logMultiGamma(.5 * nuq, q_);
  }

  // --- rng ----------------------------------------------------------------

  /// @param [out] X Matrix of size `p x q` in which to store the random draw.
  /// @param [in] Mu Mean matrix of size `p x q`.
  /// @param [in] SigmaRL Lower Cholesky factor of the row-variance matrix `SigmaR` (a matrix of size `p x p`).
  /// @param [in] OmegaCL Lower Cholesky factor of the column-precision matrix: `SigmaC^{-1} = OmegaCL * t(OmegaCL)` (a matrix of size `q x q`).
  /// @param [in] nu Shape parameter.
  inline void MatrixT::GenerateRowSColO(Ref<MatrixXd> X,
					const Ref<const MatrixXd>& Mu,
					const Ref<const MatrixXd>& SigmaRL,
					const Ref<const MatrixXd>& OmegaCL,
					double nu) {
    wish_->GenerateLowerTriXi(CL_, OmegaCL, nu);
    matnorm_->GenerateRowSColO(X, Mu, SigmaRL, CL_);
    return;
  }

  /// @param [out] X Matrix of size `p x q` in which to store the random draw.
  /// @param [in] Mu Mean matrix of size `p x q`.
  /// @param [in] OmegaRU Upper Cholesky factor of the row-precision matrix: `SigmaR^{-1} = t(OmegaRU) * OmegaRU` (a matrix of size `p x p`).
  /// @param [in] OmegaCL Lower Cholesky factor of the column-precision matrix: `SigmaC^{-1} = OmegaCL * t(OmegaCL)` (a matrix of size `q x q`).
  /// @param [in] nu Shape parameter.
  inline void MatrixT::GenerateRowOColO(Ref<MatrixXd> X,
					const Ref<const MatrixXd>& Mu,
					const Ref<const MatrixXd>& OmegaRU,
					const Ref<const MatrixXd>& OmegaCL,
					double nu) {
    wish_->GenerateLowerTriXi(CL_, OmegaCL, nu);
    matnorm_->GenerateRowOColO(X, Mu, OmegaRU, CL_);
    return;
  }

} // namespace mniw

#endif
