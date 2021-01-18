#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)

SEXP ComputeGTGamma_Eta_ADMM(Eigen::SparseMatrix<double> G_s, Eigen::VectorXd gamma_s, Eigen::VectorXd eta_s, SEXP rho){

  using Eigen::Map;
  using Eigen::VectorXd;

  double rho_s = Rcpp::as<double>(rho);

  const int m_p(G_s.cols());

  VectorXd b = VectorXd(m_p).setZero();

  b = rho_s*G_s.adjoint()*(gamma_s-eta_s);


  return Rcpp::List::create(Rcpp::Named("b") = Rcpp::wrap(b));

                           }
