#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)


SEXP ComputeBlockOLS(SEXP m, SEXP p, Eigen::SparseMatrix<double> inv_XTX_s, Eigen::VectorXd XTy_s, Eigen::MatrixXd b_s){

  using Eigen::Map;
  using Eigen::MatrixXd;

  int p_s = Rcpp::as<int>(p);
  int m_s = Rcpp::as<int>(m);

  MatrixXd beta = MatrixXd(m_s, p_s).setZero();

  int ind_j = 0;

  for(int j = 0; j < m_s; j++){

    ind_j = j*p_s;

    beta.row(j) = inv_XTX_s.block(ind_j, ind_j, p_s, p_s)*VectorXd(XTy_s.segment(ind_j, p_s) + VectorXd(b_s.row(j)));

  }

  return Rcpp::wrap(beta);
  }
