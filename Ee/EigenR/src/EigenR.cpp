// -*- mode: C++; c-indent-level: 2; c-basic-offset: 2; indent-tabs-mode: nil;
// -*-

#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

/* -------------------------------------------------------------------------- */
Eigen::MatrixXcd matricesToMatrixXcd(const Eigen::MatrixXd& Re,
                                     const Eigen::MatrixXd& Im) {
  const std::complex<double> I_ {0.0, 1.0};
  return Re.cast<std::complex<double>>() + I_ * Im.cast<std::complex<double>>();
}

Eigen::VectorXcd vectorsToVectorXcd(const Eigen::VectorXd& Re,
                                    const Eigen::VectorXd& Im) {
  const std::complex<double> I_ {0.0, 1.0};
  return Re.cast<std::complex<double>>() + I_ * Im.cast<std::complex<double>>();
}

Rcpp::List cplxMatrixToList(const Eigen::MatrixXcd& M) {
  return Rcpp::List::create(Rcpp::Named("real") = M.real(),
                            Rcpp::Named("imag") = M.imag());
}

Rcpp::List cplxVectorToList(const Eigen::VectorXcd& V) {
  return Rcpp::List::create(Rcpp::Named("real") = V.real(),
                            Rcpp::Named("imag") = V.imag());
}

/* -------------------------------------------------------------------------- */
/*
Rcpp::List cplxRcppMatrixToList(const Rcpp::ComplexMatrix M) {
  Rcpp::NumericMatrix realPart(M.nrow(), M.ncol());
  Rcpp::NumericMatrix imagPart(M.nrow(), M.ncol());
  for(auto i = 0; i < M.nrow(); i++) {
    for(auto j = 0; j < M.ncol(); j++) {
      const std::complex<double> z = M(i, j);
      realPart(i, j) = real(z);
      imagPart(i, j) = imag(z);
    }
  }
  return Rcpp::List::create(Rcpp::Named("real") = realPart,
                            Rcpp::Named("imag") = imagPart);
}
*/

Rcpp::ComplexVector cplxMatrixToRcpp(const Eigen::MatrixXcd& M) {
  Eigen::MatrixXd Mreal = M.real();
  Eigen::MatrixXd Mimag = M.imag();
  SEXP MrealS = Rcpp::wrap(Mreal);
  SEXP MimagS = Rcpp::wrap(Mimag);
  Rcpp::NumericMatrix outReal(MrealS);
  Rcpp::NumericMatrix outImag(MimagS);
  Rcpp::ComplexMatrix outRealCplx(outReal);
  Rcpp::ComplexMatrix outImagCplx(outImag);
  Rcomplex I; I.r = 0.0; I.i = 1.0;
  Rcpp::ComplexVector out = outRealCplx + I * outImagCplx;
  out.attr("dim") = Rcpp::Dimension(M.rows(), M.cols());
  return out;
}

Rcpp::NumericMatrix dblMatrixToRcpp(const Eigen::MatrixXd& M) {
  SEXP Ms = Rcpp::wrap(M);
  Rcpp::NumericMatrix out(Ms);
  return out;
}

/* Sparse stuff ------------------------------------------------------------- */
Eigen::SparseMatrix<double> realSparseMatrix(const std::vector<size_t>& i,
                                             const std::vector<size_t>& j,
                                             const std::vector<double>& Mij,
                                             const size_t nrows,
                                             const size_t ncols) {
  Eigen::SparseMatrix<double> out(nrows, ncols);
  out.reserve(Mij.size());
  for(auto k = 0; k < i.size(); k++) {
    out.insert(i[k], j[k]) = Mij[k];
  }
  return out;
}

Eigen::SparseMatrix<std::complex<double>> cplxSparseMatrix(
    const std::vector<size_t>& i,
    const std::vector<size_t>& j,
    const std::vector<std::complex<double>>& Mij,
    const size_t nrows,
    const size_t ncols) {
  Eigen::SparseMatrix<std::complex<double>> out(nrows, ncols);
  out.reserve(Mij.size());
  for(auto k = 0; k < i.size(); k++) {
    out.insert(i[k], j[k]) = Mij[k];
  }
  return out;
}

/* determinant -------------------------------------------------------------- */
template <typename Number>
Number determinant(
    const Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>& M) {
  return M.determinant();
}

template <typename Number>
Number determinant_sparse(Eigen::SparseMatrix<Number>& M) {
  Eigen::SparseLU<Eigen::SparseMatrix<Number>> solver;
  M.makeCompressed();
  solver.analyzePattern(M);
  solver.factorize(M);
  if(solver.info() != Eigen::Success) {
    throw Rcpp::exception("LU factorization has failed.");
  }
  return solver.determinant();
}

// [[Rcpp::export]]
double EigenR_det_real(const Eigen::MatrixXd& M) {
  return determinant<double>(M);
}

// [[Rcpp::export]]
std::complex<double> EigenR_det_cplx(const Eigen::MatrixXd& Re,
                                     const Eigen::MatrixXd& Im) {
  const Eigen::MatrixXcd M = matricesToMatrixXcd(Re, Im);
  return determinant<std::complex<double>>(M);
}

// [[Rcpp::export]]
double EigenR_det_sparse_real(const std::vector<size_t>& i,
                              const std::vector<size_t>& j,
                              const std::vector<double>& Mij,
                              const size_t nrows,
                              const size_t ncols) {
  Eigen::SparseMatrix<double> M = realSparseMatrix(i, j, Mij, nrows, ncols);
  return determinant_sparse<double>(M);
}

// [[Rcpp::export]]
std::complex<double> EigenR_det_sparse_cplx(
    const std::vector<size_t>& i,
    const std::vector<size_t>& j,
    const std::vector<std::complex<double>>& Mij,
    const size_t nrows,
    const size_t ncols) {
  Eigen::SparseMatrix<std::complex<double>> M =
      cplxSparseMatrix(i, j, Mij, nrows, ncols);
  return determinant_sparse<std::complex<double>>(M);
}

/* rank --------------------------------------------------------------------- */
template <typename Number>
unsigned rank(const Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>& M) {
  return M.colPivHouseholderQr().rank();
}

// [[Rcpp::export]]
unsigned EigenR_rank_real(const Eigen::MatrixXd& M) {
  return rank<double>(M);
}

// [[Rcpp::export]]
unsigned EigenR_rank_cplx(const Eigen::MatrixXd& Re,
                          const Eigen::MatrixXd& Im) {
  const Eigen::MatrixXcd M = matricesToMatrixXcd(Re, Im);
  return rank<std::complex<double>>(M);
}

/* inverse ------------------------------------------------------------------ */
template <typename Number>
Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic> inverse(
    const Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>& M) {
  return M.inverse();
}

// [[Rcpp::export]]
Eigen::MatrixXd EigenR_inverse_real(const Eigen::MatrixXd& M) {
  return inverse<double>(M);
}

// [[Rcpp::export]]
Rcpp::List EigenR_inverse_cplx(const Eigen::MatrixXd& Re,
                               const Eigen::MatrixXd& Im) {
  const Eigen::MatrixXcd M = matricesToMatrixXcd(Re, Im);
  const Eigen::MatrixXcd Minv = inverse<std::complex<double>>(M);
  return cplxMatrixToList(Minv);
}

/* kernel COD --------------------------------------------------------------- */
template <typename Number>
Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic> kernel_COD(
    const Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>& M) {
  // https://stackoverflow.com/a/53598471/1100107
  Eigen::CompleteOrthogonalDecomposition<
      Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>>
      cod(M);
  const Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic> P =
      cod.colsPermutation();
  const Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic> V =
      cod.matrixZ().transpose();
  const Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic> Kernel =
      P * V.rightCols(V.cols() - cod.rank());
  return Kernel;
}

// [[Rcpp::export]]
Eigen::MatrixXd EigenR_kernel_COD_real(const Eigen::MatrixXd& M) {
  return kernel_COD<double>(M);
}

// [[Rcpp::export]]
Rcpp::List EigenR_kernel_COD_cplx(const Eigen::MatrixXd& Re,
                                  const Eigen::MatrixXd& Im) {
  const Eigen::MatrixXcd M = matricesToMatrixXcd(Re, Im);
  const Eigen::MatrixXcd Kernel = kernel_COD<std::complex<double>>(M);
  return cplxMatrixToList(Kernel);
}

/* kernel LU ---------------------------------------------------------------- */
template <typename Number>
Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic> kernel_LU(
    const Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>& M) {
  const Eigen::FullPivLU<Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>>
      lu(M);
  return lu.kernel();
}

// [[Rcpp::export]]
Eigen::MatrixXd EigenR_kernel_LU_real(const Eigen::MatrixXd& M) {
  return kernel_LU<double>(M);
}

// [[Rcpp::export]]
Rcpp::List EigenR_kernel_LU_cplx(const Eigen::MatrixXd& Re,
                                 const Eigen::MatrixXd& Im) {
  const Eigen::MatrixXcd M = matricesToMatrixXcd(Re, Im);
  const Eigen::MatrixXcd Kernel = kernel_LU<std::complex<double>>(M);
  return cplxMatrixToList(Kernel);
}

/* image LU ----------------------------------------------------------------- */
template <typename Number>
Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic> image_LU(
    const Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>& M) {
  const Eigen::FullPivLU<Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>>
      lu(M);
  return lu.image(M);
}

// [[Rcpp::export]]
Eigen::MatrixXd EigenR_image_LU_real(const Eigen::MatrixXd& M) {
  return image_LU<double>(M);
}

// [[Rcpp::export]]
Rcpp::List EigenR_image_LU_cplx(const Eigen::MatrixXd& Re,
                                const Eigen::MatrixXd& Im) {
  const Eigen::MatrixXcd M = matricesToMatrixXcd(Re, Im);
  const Eigen::MatrixXcd Image = image_LU<std::complex<double>>(M);
  return cplxMatrixToList(Image);
}

/* image QR ----------------------------------------------------------------- */
template <typename Number>
Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic> image_QR(
    const Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>& M) {
  const Eigen::ColPivHouseholderQR<
      Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>>
      qr = M.colPivHouseholderQr();
  const Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic> Q =
      qr.householderQ().setLength(qr.nonzeroPivots());
  return Q.leftCols(qr.rank());
}

// [[Rcpp::export]]
Eigen::MatrixXd EigenR_image_QR_real(const Eigen::MatrixXd& M) {
  return image_QR<double>(M);
}

// [[Rcpp::export]]
Rcpp::List EigenR_image_QR_cplx(const Eigen::MatrixXd& Re,
                                const Eigen::MatrixXd& Im) {
  const Eigen::MatrixXcd M = matricesToMatrixXcd(Re, Im);
  const Eigen::MatrixXcd Image = image_QR<std::complex<double>>(M);
  return cplxMatrixToList(Image);
}

/* image COD ---------------------------------------------------------------- */
template <typename Number>
Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic> image_COD(
    const Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>& M) {
  Eigen::CompleteOrthogonalDecomposition<
      Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>>
      cod(M);
  const Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic> Q =
      cod.householderQ();
  return Q.leftCols(cod.rank());
}

// [[Rcpp::export]]
Eigen::MatrixXd EigenR_image_COD_real(const Eigen::MatrixXd& M) {
  return image_COD<double>(M);
}

// [[Rcpp::export]]
Rcpp::List EigenR_image_COD_cplx(const Eigen::MatrixXd& Re,
                                 const Eigen::MatrixXd& Im) {
  const Eigen::MatrixXcd M = matricesToMatrixXcd(Re, Im);
  const Eigen::MatrixXcd Image = image_COD<std::complex<double>>(M);
  return cplxMatrixToList(Image);
}

/* QR ----------------------------------------------------------------------- */
template <typename Number>
std::vector<Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>> QR(
    const Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>& M) {
  const Eigen::HouseholderQR<
      Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>>
      qr = M.householderQr();
  const Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic> R =
      qr.matrixQR().template triangularView<Eigen::Upper>();
  const Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic> Q =
      qr.householderQ();
  return {Q, R};
}

// [[Rcpp::export]]
Rcpp::List EigenR_QR_real(const Eigen::MatrixXd& M) {
  const std::vector<Eigen::MatrixXd> QRdecomp = QR<double>(M);
  return Rcpp::List::create(Rcpp::Named("Q") = QRdecomp[0],
                            Rcpp::Named("R") = QRdecomp[1]);
}

// [[Rcpp::export]]
Rcpp::List EigenR_QR_cplx(const Eigen::MatrixXd& Re,
                          const Eigen::MatrixXd& Im) {
  const Eigen::MatrixXcd M = matricesToMatrixXcd(Re, Im);
  const std::vector<Eigen::MatrixXcd> QRdecomp = QR<std::complex<double>>(M);
  return Rcpp::List::create(Rcpp::Named("Q") = cplxMatrixToList(QRdecomp[0]),
                            Rcpp::Named("R") = cplxMatrixToList(QRdecomp[1]));
}

/* Cholesky ----------------------------------------------------------------- */
template <typename Number>
struct Cholesky {
  Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic> U;
  Number determinant;
};

template <typename Number>
Cholesky<Number> chol(
    const Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>& M) {
  const Eigen::LLT<Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>>
      lltOfM(M);
  if(lltOfM.info() != Eigen::Success) {
    throw Rcpp::exception("The matrix is not positive definite.");
  }
  const Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic> 
    U = lltOfM.matrixU();
  Cholesky<Number> out;
  out.U = U;
  out.determinant = pow(U.diagonal().prod(), 2);
  return out;
}

/*
template <typename Matrix, typename Number>
Matrix chol_sparse(
    Eigen::SparseMatrix<Number>& M) {
  Eigen::SimplicialLLT<Eigen::SparseMatrix<Number>> solver;
  M.makeCompressed();
  solver.analyzePattern(M);
  solver.factorize(M);
  if(solver.info() != Eigen::Success) {
    throw Rcpp::exception("LU factorization has failed.");
  }
  //Matrix U;

  if(std::is_same<Number, std::complex<double>>::value) {
    Rcomplex I;
    I.r = 0.0; I.i = 1.0;
    Rcpp::NumericMatrix Ureal0 = Rcpp::wrap(solver.matrixU().real());
    Rcpp::NumericMatrix Uimag0 = Rcpp::wrap(solver.matrixU().imag());
    Matrix Ureal(Ureal0); Matrix Uimag(Uimag0);
    Matrix U(Ureal + I * Uimag);
    return U;
  }else {
    Matrix U = Rcpp::wrap(solver.matrixU());
    return U;
  }
  U.attr("determinant") = solver.determinant();
  return U;
}
*/

template <typename Number>
Cholesky<Number> chol_sparse(Eigen::SparseMatrix<Number>& M) {
  Eigen::SimplicialLLT<Eigen::SparseMatrix<Number>> solver;
  M.makeCompressed();
  solver.analyzePattern(M);
  solver.factorize(M);
  if(solver.info() != Eigen::Success) {
    throw Rcpp::exception("LU factorization has failed.");
  }
  Cholesky<Number> out;
  out.U = solver.matrixU();
  out.determinant = solver.determinant();
  return out;
}

/*
Rcpp::NumericMatrix chol_sparse_real(Eigen::SparseMatrix<double>& M) {
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
  M.makeCompressed();
  solver.analyzePattern(M);
  solver.factorize(M);
  if(solver.info() != Eigen::Success) {
    throw Rcpp::exception("LU factorization has failed.");
  }
  Eigen::MatrixXd U = solver.matrixU();
  SEXP s = Rcpp::wrap(U);
  Rcpp::NumericMatrix out(s);
  out.attr("determinant") = solver.determinant();
  return out;
}

Rcpp::ComplexVector chol_sparse_cplx(
    Eigen::SparseMatrix<std::complex<double>>& M) {
  Eigen::SimplicialLLT<Eigen::SparseMatrix<std::complex<double>>> solver;
  M.makeCompressed();
  solver.analyzePattern(M);
  solver.factorize(M);
  if(solver.info() != Eigen::Success) {
    throw Rcpp::exception("LU factorization has failed.");
  }
  Eigen::MatrixXcd U = solver.matrixU();
  Eigen::MatrixXd Ureal = U.real();
  Eigen::MatrixXd Uimag = U.imag();
  SEXP UrealS = Rcpp::wrap(Ureal);
  SEXP UimagS = Rcpp::wrap(Uimag);
  Rcpp::NumericMatrix outReal(UrealS);
  Rcpp::NumericMatrix outImag(UimagS);
  Rcpp::ComplexMatrix outRealCplx(outReal);
  Rcpp::ComplexMatrix outImagCplx(outImag);
  Rcomplex I;
  I.r = 0.0; I.i = 1.0;
  Rcpp::ComplexVector out = outRealCplx + I * outImagCplx;
  out.attr("dim") = Rcpp::Dimension(U.rows(), U.cols());
  out.attr("determinant") = solver.determinant();
  return out;
}
*/

// [[Rcpp::export]]
Rcpp::NumericMatrix EigenR_chol_real(const Eigen::MatrixXd& M) {
  Cholesky<double> cholesky = chol<double>(M);
  Rcpp::NumericMatrix U = dblMatrixToRcpp(cholesky.U);
  U.attr("determinant") = cholesky.determinant;
  return U;
}

// [[Rcpp::export]]
Rcpp::ComplexVector EigenR_chol_cplx(const Eigen::MatrixXd& Re,
                                     const Eigen::MatrixXd& Im) {
  const Eigen::MatrixXcd M = matricesToMatrixXcd(Re, Im);
  Cholesky<std::complex<double>> cholesky = chol<std::complex<double>>(M);
  Rcpp::ComplexVector U = cplxMatrixToRcpp(cholesky.U);
  U.attr("determinant") = cholesky.determinant;
  return U;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix EigenR_chol_sparse_real(const std::vector<size_t>& i,
                                            const std::vector<size_t>& j,
                                            const std::vector<double>& Mij,
                                            const size_t nrows,
                                            const size_t ncols) {
  Eigen::SparseMatrix<double> M = realSparseMatrix(i, j, Mij, nrows, ncols);
  Cholesky<double> cholesky = chol_sparse<double>(M);
  Rcpp::NumericMatrix U = dblMatrixToRcpp(cholesky.U);
  U.attr("determinant") = cholesky.determinant;
  return U;
}

// [[Rcpp::export]]
Rcpp::ComplexVector EigenR_chol_sparse_cplx(
    const std::vector<size_t>& i,
    const std::vector<size_t>& j,
    const std::vector<std::complex<double>>& Mij,
    const size_t nrows,
    const size_t ncols) {
  Eigen::SparseMatrix<std::complex<double>> M =
      cplxSparseMatrix(i, j, Mij, nrows, ncols);
  Cholesky<std::complex<double>> cholesky =
      chol_sparse<std::complex<double>>(M);
  Rcpp::ComplexVector U = cplxMatrixToRcpp(cholesky.U);
  U.attr("determinant") = cholesky.determinant;
  return U;
}

/* UtDU --------------------------------------------------------------------- */
template <typename Number>
Rcpp::List UtDU(
    const Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>& M) {
  const Eigen::LDLT<Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>>
      ldltOfM(M);
  if(ldltOfM.info() != Eigen::Success) {
    throw Rcpp::exception("Factorization has failed.");
  }
  const Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic> U =
      ldltOfM.matrixU();
  const Eigen::Matrix<Number, Eigen::Dynamic, 1> D = ldltOfM.vectorD();
  const Eigen::Transpositions<Eigen::Dynamic> T = ldltOfM.transpositionsP();
  Eigen::VectorXi perm(T.size());
  for(auto i = 0; i < T.size(); i++) {
    perm(i) = i;
  }
  Rcpp::List out =
      Rcpp::List::create(Rcpp::Named("U") = U, Rcpp::Named("D") = D,
                         Rcpp::Named("perm") = T * perm);
  bool positive = ldltOfM.isPositive();
  double rcond = ldltOfM.rcond();
  out.attr("positive") = positive;
  out.attr("rcond") = rcond;
  return out;
}

/*
template <typename Number>
Rcpp::List UtDU_sparse(
    Eigen::SparseMatrix<Number>& M) {
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<Number>> ldltOfM;
  M.makeCompressed();
  ldltOfM.analyzePattern(M);
  ldltOfM.factorize(M);
  if(ldltOfM.info() != Eigen::Success) {
    throw Rcpp::exception("Factorization has failed.");
  }
  const Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic> U =
    ldltOfM.matrixU();
  const Eigen::Matrix<Number, Eigen::Dynamic, 1> D = ldltOfM.vectorD();
  const Eigen::Transpositions<Eigen::Dynamic> T = ldltOfM.transpositionsP();
  Eigen::VectorXi perm(T.size());
  for(auto i = 0; i < T.size(); i++) {
    perm(i) = i;
  }
  const Rcpp::List out =
    Rcpp::List::create(Rcpp::Named("U") = U, Rcpp::Named("D") = D,
                       Rcpp::Named("perm") = T * perm);
  return out;
}
*/

// [[Rcpp::export]]
Rcpp::List EigenR_UtDU_real(const Eigen::MatrixXd& M) {
  return UtDU<double>(M);
}

// [[Rcpp::export]]
Rcpp::List EigenR_UtDU_cplx(const Eigen::MatrixXd& Re,
                            const Eigen::MatrixXd& Im) {
  const Eigen::MatrixXcd M = matricesToMatrixXcd(Re, Im);
  const Rcpp::List utdu = UtDU<std::complex<double>>(M);
  Rcpp::List out =
      Rcpp::List::create(Rcpp::Named("U") = cplxMatrixToList(utdu["U"]),
                         Rcpp::Named("D") = cplxVectorToList(utdu["D"]),
                         Rcpp::Named("perm") = utdu["perm"]);
  out.attr("positive") = utdu.attr("positive");
  out.attr("rcond") = utdu.attr("rcond");
  return out;
}

/*
// [[Rcpp::export]]
Rcpp::List EigenR_UtDU_sparse_real(const std::vector<size_t>& i,
                                   const std::vector<size_t>& j,
                                   const std::vector<double>& Mij,
                                   const size_t nrows,
                                   const size_t ncols) {
  Eigen::SparseMatrix<double> M = realSparseMatrix(i, j, Mij, nrows, ncols);
  return UtDU_sparse<double>(M);
}

// [[Rcpp::export]]
Rcpp::List EigenR_UtDU_sparse_cplx(const std::vector<size_t>& i,
                                   const std::vector<size_t>& j,
                                   const std::vector<std::complex<double>>& Mij,
                                   const size_t nrows,
                                   const size_t ncols) {
  Eigen::SparseMatrix<std::complex<double>> M =
    cplxSparseMatrix(i, j, Mij, nrows, ncols);
  const Rcpp::List utdu = UtDU_sparse<std::complex<double>>(M);
  Rcpp::List out =
    Rcpp::List::create(Rcpp::Named("U") = cplxMatrixToList(utdu["U"]),
                       Rcpp::Named("D") = cplxVectorToList(utdu["D"]),
                       Rcpp::Named("perm") = utdu["perm"]);
  return out;
}
*/

/* Least-squares ------------------------------------------------------------ */
template <typename Number>
Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic> lsSolve(
    const Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>& A,
    const Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>& b) {
  return A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
}

// [[Rcpp::export]]
Eigen::MatrixXd EigenR_lsSolve_real(const Eigen::MatrixXd& A,
                                    const Eigen::MatrixXd& b) {
  return lsSolve<double>(A, b);
}

// [[Rcpp::export]]
Rcpp::List EigenR_lsSolve_cplx(const Eigen::MatrixXd& ReA,
                               const Eigen::MatrixXd& ImA,
                               const Eigen::MatrixXd& Reb,
                               const Eigen::MatrixXd& Imb) {
  const Eigen::MatrixXcd A = matricesToMatrixXcd(ReA, ImA);
  const Eigen::MatrixXcd b = matricesToMatrixXcd(Reb, Imb);
  const Eigen::MatrixXcd X = lsSolve<std::complex<double>>(A, b);
  return cplxMatrixToList(X);
}
