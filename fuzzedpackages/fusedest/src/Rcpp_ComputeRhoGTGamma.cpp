#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)

Eigen::VectorXd ComputeRhoGTGamma(Eigen::SparseMatrix<double> G_s, Eigen::VectorXd gamma_s, SEXP rho){

  using Eigen::Map;
  using Eigen::MatrixXd;
  using Eigen::SparseMatrix;
  using Eigen::MappedSparseMatrix;


  double rho_s = Rcpp::as<double>(rho);

  Eigen::VectorXd rhoGTGamma = rho_s*G_s.adjoint()*gamma_s;

  /*inv_XTX_rhoGTG = XTX_rhoGTG.llt().solve(Id);*/

  return rhoGTGamma;
  }
