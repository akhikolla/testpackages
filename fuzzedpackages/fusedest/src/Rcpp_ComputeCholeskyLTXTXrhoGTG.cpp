#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)
// For more on using Rcpp click the Help button on the editor toolbar

Eigen::SparseMatrix<double> ComputeCholeskyLTXTXrhoGTG(Eigen::SparseMatrix<double> XTX_s,  Eigen::SparseMatrix<double> G_s, SEXP rho){

  using Eigen::Map;
  using Eigen::MatrixXd;
  using Eigen::SparseMatrix;
  using Eigen::MappedSparseMatrix;


  double rho_s = Rcpp::as<double>(rho);

  Eigen::SparseMatrix<double> XTX_rhoGTG = XTX_s+ rho_s*G_s.adjoint()*G_s;

  SimplicialLLT<SparseMatrix<double>, Lower, NaturalOrdering<int> > llt(XTX_rhoGTG);

  return llt.matrixL().adjoint();
  }
