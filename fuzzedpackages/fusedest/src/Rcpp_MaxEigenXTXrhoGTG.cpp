#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)

SEXP MaxEigenXTXrhoGTG(SEXP X, SEXP b, SEXP max_iter, SEXP tol_err){

  using Eigen::Map;
  using Eigen::MatrixXd;
  using Eigen::VectorXd;

  const Map<MatrixXd> X_s(as<Map<MatrixXd> >(X));
  const Map<VectorXd> b_s(as<Map<VectorXd> >(b)); /* random vector */
  int max_iter_s = Rcpp::as<int>(max_iter); /* Maximum number of iterations */
  double tol_err_s = Rcpp::as<double>(tol_err);  /* Tolerance error */


  /*Perform eigenvalue decomposition!!*/

  MatrixXd XTX = X_s.adjoint()*X_s;
  const int p(XTX.cols());
  VectorXd b_old = b_s;
  VectorXd b_new = VectorXd(p).setZero(); /* Initial values for betahat */
  VectorXd XTXb = VectorXd(p).setZero();


  int l = 0;
  double err = 1e6;

  while(err > tol_err_s && l < max_iter_s){

    XTXb = XTX*b_old;
    double Z = VectorXd(XTXb).norm();
    b_new =  XTXb/Z;
    err = VectorXd(b_old - b_new).norm();
    l+= 1;
    b_old = b_new;

  }

  double mu_power_d = double(b_old.squaredNorm());
  double mu_power_n = double(b_old.adjoint()*XTX*b_old);
  double mu_power = mu_power_n/mu_power_d;


  return  Rcpp::List::create(Rcpp::Named("max_eigenvalue") = Rcpp::wrap(mu_power));
                           }
