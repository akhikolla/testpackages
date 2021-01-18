#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)

Eigen::SparseMatrix<double> ComputeXTXrhoGTG(Eigen::SparseMatrix<double> XTX_s, Eigen::SparseMatrix<double> G_s, SEXP rho){

  using Eigen::Map;
  using Eigen::MatrixXd;
  using Eigen::SparseMatrix;
  using Eigen::MappedSparseMatrix;


  double rho_s = Rcpp::as<double>(rho);

  Eigen::SparseMatrix<double> rhoGTG_s = rho_s*G_s.adjoint()*G_s;
  Eigen::SparseMatrix<double> XTX_rhoGTG = XTX_s + rhoGTG_s;

  return XTX_rhoGTG;
  }
