#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

SEXP ComputeInvXTXrhoGTG(Eigen::SparseMatrix<double> XTX_s, Eigen::SparseMatrix<double> G_s, SEXP rho){

  using Eigen::Map;
  using Eigen::MatrixXd;
  using Eigen::SparseMatrix;
  using Eigen::MappedSparseMatrix;


  double rho_s = Rcpp::as<double>(rho);

  const int m_p(XTX_s.cols());

  Eigen::SparseMatrix<double> XTX_rhoGTG = XTX_s+rho_s*G_s.adjoint()*G_s;

  /*inv_XTX_rhoGTG = XTX_rhoGTG.llt().solve(Id);*/

  SimplicialLLT<SparseMatrix<double> > solver;
  solver.compute(XTX_rhoGTG);
  SparseMatrix<double> Id(m_p,m_p);
  Id.setIdentity();

  auto inv_XTX_rhoGTG = solver.solve(Id);

  return Rcpp::wrap(inv_XTX_rhoGTG);
  }
