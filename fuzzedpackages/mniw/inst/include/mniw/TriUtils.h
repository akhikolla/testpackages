/// @file TriUtils.h
///
/// @brief Utility functions for triangular matrices.
///
/// Various types of products, solvers, and decompositions involving lower and upper triangular matrices.  Almost all arguments to these functions can be subsets of larger matrices; see "Writing Functions Taking Eigen Types as Parameters" in the Eigen documentation.

#ifndef TriUtils_h
#define TriUtils_h 1

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
//#include <iostream>

namespace mniw {
  
  using namespace Eigen;

  // --- matrix multiplication ---------------------------------------------------

  /// @brief Multiplication of two lower triangular matrices
  ///
  /// Performs the matrix multiplication `X = L1 * L2`, where `L1` and `L2` are lower triangular matrices.
  ///
  /// @param [in] L1 First `n x n` lower triangular matrix.
  /// @param [in] L2 Second `n x n` lower triangular matrix.
  /// @param [out] X Matrix of size `n x n` containing the product of `L1` and `L2`.
  ///
  /// @note Does not work properly if any of the inputs arguments are also outputs.
  /// @note For safest results, all "off-triangular" elements of triangular arguments should be set to zero.
  template <typename T1, typename T2>
  void triMultLL(const Eigen::MatrixBase<T1>& X,
		 const Ref<const MatrixXd>& L1,
		 const Eigen::MatrixBase<T2>& L2) {
    // hack to use triangularView: see Eigen documentation
    //"Writing Functions Taking Eigen Types as Parameters"
    // #define _VL const_cast<Eigen::MatrixBase<T1>& >(VL)
    // #define _XL const_cast<Eigen::MatrixBase<T2>& >(XL)
#define _X const_cast<Eigen::MatrixBase<T1>& >(X)
#define _L2 const_cast<Eigen::MatrixBase<T2>& >(L2)
    _X.template triangularView<Eigen::Lower>() = L1 * _L2.template triangularView<Eigen::Lower>();
#undef _X
#undef _L2
    return;
  }

  /// @brief Right-multiplication by upper triangular matrix
  ///
  /// Performs the matrix multiplication `Y = X * U`, where `U` is an upper triangular matrix.
  ///
  /// @param [in] X Matrix of size `n x p` to be multiplied.
  /// @param [in] U Lower triangular matrix of size `n x n` right-multiplying `X`.
  /// @param [out] Y Matrix of size `n x p` containing the product of `X` and `U`.
  ///
  /// @note Does not work properly if any of the inputs arguments are also outputs.
  /// @note For safest results, all "off-triangular" elements of triangular arguments should be set to zero.
  template <typename T1>
  void triMultXU(Ref<MatrixXd> Y, const Ref<const MatrixXd>& X,
		 const Eigen::MatrixBase<T1>& U) {
    Y.noalias() = X * U.template triangularView<Eigen::Upper>();
    return;
  }

  /// @brief Left-multiplication by a lower triangular matrix
  ///
  /// Performs the matrix multiplication `Y = L * X`, where `L` is a lower triangular matrix.
  ///
  /// @param [in] X Matrix of size `n x p` to be multiplied.
  /// @param [in] L Lower triangular matrix of size `n x n` left-multiplying `X`.
  /// @param [out] Y Matrix of size `n x p` containing the product of `L` and `X`.
  ///
  /// @note Does not work properly if any of the inputs arguments are also outputs.
  /// @note For safest results, all "off-triangular" elements of triangular arguments should be set to zero.
  template <typename T1>
  void triMultLX(Ref<MatrixXd> Y, const Eigen::MatrixBase<T1>& L,
		 const Ref<const MatrixXd>& X) {
    Y.noalias() = L.template triangularView<Eigen::Lower>() * X;
    return;
  }

  // --- solution of linear system -----------------------------------------------

  /// @brief In-place solution of a lower triangular system
  ///
  /// Performs the matrix multiplication `X = L^{-1} * X`, where `L` is a lower triangular matrix.
  ///
  /// @param [in,out] X Matrix of size `n x p` on RHS of linear system.  This is also where solution will be returned, i.e., in-place.
  /// @param [in] L Lower triangular matrix of size `n x n` on LHS of linear system.
  ///
  /// @note For safest results, all "off-triangular" elements of triangular arguments should be set to zero.
  template <typename T1>
  void triMultLiX(Ref<MatrixXd> X, const Eigen::MatrixBase<T1>& L) {
#ifdef _MSC_VER
    L.template triangularView<Eigen::Lower>().solveInPlace(X);
#else
    L.template triangularView<Eigen::Lower>().template solveInPlace(X);
#endif
    return;
  }

  /// @brief In-place solution of a reverse lower triangular system
  ///
  /// Performs the multiplication `X = X * L^{-1}`, where `L` is a lower triangular matrix.
  ///
  /// @param [in,out] X Matrix of size `n x p` on RHS of linear system.  This is also where solution will be returned, i.e., in-place.
  /// @param [in] L Lower triangular matrix of size `n x n` on LHS of linear system.
  ///
  /// @note Does not work properly if any of the inputs arguments are also outputs.
  /// @note For safest results, all "off-triangular" elements of triangular arguments should be set to zero.
  template <typename T1>
  void triMultXLi(Ref<MatrixXd> X, const Eigen::MatrixBase<T1>& L) {
#ifdef _MSC_VER
    L.template triangularView<Eigen::Lower>().solveInPlace<OnTheRight>(X);
#else
    L.template triangularView<Eigen::Lower>().template solveInPlace<OnTheRight>(X);
#endif
    return;
  }

  /// @brief In-place solution of an upper triangular system
  ///
  /// Performs the multiplication `X = U^{-1} * X`, where `U` is an upper triangular matrix.
  ///
  /// @param [in,out] X Matrix of size `n x p` on RHS of linear system.  This is also where solution will be returned, i.e., in-place.
  /// @param [in] U Upper triangular matrix of size `n x n` on LHS of linear system.
  ///
  /// @note Does not work properly if any of the inputs arguments are also outputs.
  /// @note For safest results, all "off-triangular" elements of triangular arguments should be set to zero.
  template <typename T1>
  void triMultUiX(const Eigen::MatrixBase<T1>& U, Ref<MatrixXd> X) {
#ifdef _MSC_VER
    U.template triangularView<Eigen::Upper>().solveInPlace(X);
#else
    U.template triangularView<Eigen::Upper>().template solveInPlace(X);
#endif
    return;
  }


  // --- transpose-products ------------------------------------------------------

  /// @brief Transpose-product of upper triangular matrices
  ///
  /// Performs the multiplication `X = U&apos; * U`, where `U` is an upper triangular matrix.
  ///
  /// @param [out] X Matrix of size `n x n` containing the transpose-product of `U`.
  /// @param [in] U Upper triangular matrix of size `n x n`.
  /// @param [in] L Lower triangular matrix of size `n x n` used for intermediate calculations.
  ///
  /// @note Does not work properly if any of the inputs arguments are also outputs.
  /// @note For safest results, all "off-triangular" elements of triangular arguments should be set to zero.
  template <typename T1, typename T2>
  void crossProdUtU(const Eigen::MatrixBase<T1>& X,
		    const Ref<const MatrixXd>& U,
		    const Eigen::MatrixBase<T2>& L) {
    // hack
#define _X const_cast<Eigen::MatrixBase<T1>& >(X)
#define _L const_cast<Eigen::MatrixBase<T2>& >(L)
    _L.template triangularView<Eigen::Lower>() = U.adjoint();
    _X.template triangularView<Eigen::Upper>() = L.template triangularView<Eigen::Lower>() * U;
    _X.template triangularView<Eigen::Lower>() = X.adjoint();
#undef _X
#undef _L
    return;
  }

  /// @brief Reverse transpose-product of lower triangular matrices
  ///
  /// Performs the multiplication `X = L * L&apos;`, where `L` is a lower triangular matrix.
  ///
  /// @param [out] X Matrix of size `n x n` containing the reverse transpose-product of `L`.
  /// @param [in] L Lower triangular matrix of size `n x n`.
  /// @param [in] U Upper triangular matrix of size `n x n` used for intermediate calculations.
  ///
  /// @note Does not work properly if any of the inputs arguments are also outputs.
  /// @note For safest results, all "off-triangular" elements of triangular arguments should be set to zero.
  template <typename T1, typename T2>
  void CrossProdLLt(const Eigen::MatrixBase<T1>& X,
		    const Ref<const MatrixXd>& L,
		    const Eigen::MatrixBase<T2>& U) {
    // hack
#define _X const_cast<Eigen::MatrixBase<T1>& >(X)
#define _U const_cast<Eigen::MatrixBase<T2>& >(U)
    _U.template triangularView<Eigen::Upper>() = L.adjoint();
    _X.template triangularView<Eigen::Upper>() = L * U.template triangularView<Eigen::Upper>();
    _X.template triangularView<Eigen::Lower>() = X.adjoint();
#undef _X
#undef _U
    return;
  }

  /// @brief Transpose-product of lower triangular matrices
  ///
  /// Performs the multiplication `X = L&apos; * L`, where `L` is a lower triangular matrix.
  ///
  /// @param [out] X Matrix of size `n x n` containing the transpose-product of `L`.
  /// @param [in] L Lower triangular matrix of size `n x n`.
  /// @param [in] U Upper triangular matrix of size `n x n` used for intermediate calculations.
  ///
  /// @note Does not work properly if any of the inputs arguments are also outputs.
  /// @note For safest results, all "off-triangular" elements of triangular arguments should be set to zero.
  template <typename T1, typename T2>
  void CrossProdLtL(const Eigen::MatrixBase<T1>& X,
		    const Ref<const MatrixXd>& L,
		    const Eigen::MatrixBase<T2>& U) {
    // hack
#define _X const_cast<Eigen::MatrixBase<T1>& >(X)
#define _U const_cast<Eigen::MatrixBase<T2>& >(U)
    _U.template triangularView<Eigen::Upper>() = L.adjoint();
    _X.template triangularView<Eigen::Upper>() = U.template triangularView<Eigen::Upper>() * L;
    _X.template triangularView<Eigen::Lower>() = X.adjoint();
#undef _X
#undef _U
    return;
  }

  /// @brief Inverse of reverse transpose-product
  ///
  /// Performs the calculation `X = (L&apos; * L)^{-1}`, where `L` is a lower triangular matrix.
  ///
  /// @param [out] X Matrix of size `n x n` containing the inverse of the reverse transpose-product of `L`.
  /// @param [in] L Lower triangular matrix of size `n x n`.
  /// @param [in] U Upper triangular matrix of size `n x n` used for intermediate calculations.
  ///
  /// @param [in] L2 Lower triangular matrix of size `n x n` used for intermediate calculations.
  /// @param [in] I Identity matrix of size `n x n` used for intermediate calculations.
  ///
  /// @note Does not work properly if any of the inputs arguments are also outputs.
  /// @note For safest results, all "off-triangular" elements of triangular arguments should be set to zero.
  template <typename T1>
  void InverseLLt(Ref<MatrixXd> X,
		  const Eigen::MatrixBase<T1>& L,
		  Ref<MatrixXd> L2,
		  Ref<MatrixXd> U,
		  const Ref<const MatrixXd>& I) {
    // L2 = L^{-1}
    L2 = L.template triangularView<Eigen::Lower>().solve(I);
    // X = L2' * L2 = L'{-1} * L^{-1} = (L * L')^{-1}
    CrossProdLtL(X, L2, U);
    return;
  }

  // --- reverse-cholesky --------------------------------------------------------

  // /// Anti-transpose of a matrix
  // ///
  // /// Transpose the lower triangular elements of a square matrix across the anti-diagonal (the remaining elements of the output are left untouched).  So for example, we would have
  // ///
  // /// \f[
  // /// \begin{bmatrix} 1 & 2 & 3 \\ 4 & 5 & 6 \\ 7 & 8 & 9 \end{bmatrix}
  // /// \qquad \longrightarrow \qquad
  // /// \begin{bmatrix} 9 & & \\ 8 & 5 & \\ 7 & 4 & 1 \end{bmatrix}.
  // /// \f]
  // ///
  // /// @param [out] Y Matrix of size `n x n` to which anti-transpose is output.
  // /// @param [in] X Matrix of size `n x n` of which the lower triangular elements will be anti-transposed.
  // ///
  // /// @note Does not work properly if any of the inputs arguments are also outputs.
  // inline void AntiTransposeLowerTri(Ref<MatrixXd> Y, const Ref<const MatrixXd>& X) {
  //   int q = X.cols();
  //   int ii, jj;
  //   // anti-transpose lower part of X into Y
  //   for(ii=0; ii<q; ii++) {
  //     for(jj=0; jj<=ii; jj++) {
  //       Y(q-1-jj,q-1-ii) = X(ii,jj);
  //     }
  //   }
  //   return;
  // }

  /// @brief Reverse-Cholesky decomposition of a positive-definite matrix
  ///
  /// Calculates the lower triangular matrix `L` satisfying `V = L&apos;L`, where `V` is a positive-definite matrix.
  ///
  /// @param [out] L Lower triangular matrix of size `n x n`.
  /// @param [in] V Positive-definite matrix of size `n x n`.
  /// @param [in] llt Cholesky solver via reference to `Eigen::LLT` object used for intermediate calculations.
  ///
  /// @note These calculations leave the upper triangular half of `L` unchanged, such that the output is truly triangular only if the upper triangular half of `L` was initialized to zero.
  ///
  inline void ReverseCholesky(Ref<MatrixXd> L,
			      const Ref<const MatrixXd>& V,
			      LLT<MatrixXd>& llt) {
    int q = L.cols();
    int ii, jj;
    // L = anti_transpose(lower_tri(V))
    for(ii=0; ii<q; ii++) {
      for(jj=0; jj<=ii; jj++) {
	L(q-1-jj,q-1-ii) = V(ii,jj);
      }
    }
    // llt "=" chol(L)
    llt.compute(L.triangularView<Eigen::Lower>());
    // L "=" anti_transpose(llq)
    for(ii=0; ii<q; ii++) {
      for(jj=0; jj<=ii; jj++) {
	L(q-1-jj,q-1-ii) = llt.matrixL()(ii,jj);
      }
    }
    // // upper_tri(L) = 0
    // if(upperZero) {
    L = L.triangularView<Eigen::Lower>();
    // }
    return;
  }


  // --- log-determinant ---------------------------------------------------------

  /// @brief Logarithm of the determinant of a Cholesky decomposition
  ///
  /// Calculates `log|L|`, where `L` is the lower triangular factor of a Cholesky decomposition.
  ///
  /// @param [in] cholV Cholesky factor of size `n x n`.  This is represented by an object of class `Eigen::LLT`, such that the triangular factor is given by `L = cholV.matrixL()`.
  /// @return The logarithm of the determinant of `L`.
  inline double logDetCholV(LLT<MatrixXd>& cholV) {
    double ldC = 0.0;
    for(int ii=0; ii<cholV.cols(); ii++) {
      ldC += log(cholV.matrixL()(ii,ii));
    }
    return ldC;
  }

  /// @brief Logarithm of the determinant of a positive-definite matrix
  ///
  /// Calculates `log|V|`, where `V` is a positive-definite matrix.
  ///
  /// @param [in] V Positive-definite matrix of size `n x n`.
  /// @param [in] cholV Cholesky solver of size `n x n` required for intermediate calculations.
  /// @return The logarithm of the determinant of `V`.
  inline double logDetV(MatrixXd V, LLT<MatrixXd> cholV) {
    double ldV = 0.0;
    cholV.compute(V);
    return 2.0 * logDetCholV(cholV);
    // for(int ii=0; ii<V.cols(); ii++) {
    //   ldV += log(cholV.matrixL()(ii,ii));
    // }
    // return 2.0*ldV;
  }

} // namespace mniw
  
#endif
