// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::MappedSparseMatrix<double> MSpMat;
typedef Eigen::Map<Eigen::MatrixXd> MapMatd;

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//
// [[Rcpp::export]]
Eigen::MatrixXd pcg_dense(const Eigen::Map<Eigen::MatrixXd> A, const Eigen::Map<Eigen::MatrixXd> B, const double tol) {

  
  Eigen::ConjugateGradient<Eigen::MatrixXd, Eigen::Lower|Eigen::Upper> cg;
  cg.setTolerance(tol);
  cg.compute(A);
  Eigen::MatrixXd x = cg.solve(B);
  
  return x;
}

// [[Rcpp::export]]
Eigen::MatrixXd pcg_sparse(const Eigen::MappedSparseMatrix<double> A, const Eigen::Map<Eigen::MatrixXd> B, const double tol) {
  
  Eigen::ConjugateGradient<SpMat, Eigen::Lower|Eigen::Upper > solver_cg;
  solver_cg.setTolerance(tol);
  solver_cg.compute(A);
  Eigen::MatrixXd x = solver_cg.solve(B);
  
  return x;
}
