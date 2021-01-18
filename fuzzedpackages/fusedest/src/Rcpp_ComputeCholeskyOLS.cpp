#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)
// For more on using Rcpp click the Help button on the editor toolbar

Eigen::VectorXd ComputeCholeskyOLS(Eigen::SparseMatrix<double> LT_s, Eigen::VectorXd XTy_s, Eigen::VectorXd b){

  using Eigen::Map;
  using Eigen::MatrixXd;
  using Eigen::SparseMatrix;
  using Eigen::MappedSparseMatrix;

  int p_s(XTy_s.size());

  VectorXd theta = VectorXd(p_s).setZero();
  VectorXd beta = VectorXd(p_s).setZero();

  theta = LT_s.transpose().triangularView<Lower>().solve(XTy_s + b);

  beta = LT_s.triangularView<Upper>().solve(theta);

  return beta;
  }
