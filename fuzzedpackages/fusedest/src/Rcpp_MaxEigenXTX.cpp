#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)


double MaxEigenXTX(Eigen::MatrixXd X, Eigen::VectorXd z, SEXP max_iter, SEXP tol_err){

  using Eigen::Map;
  using Eigen::MatrixXd;
  using Eigen::VectorXd;

  int max_iter_s = Rcpp::as<int>(max_iter); /* Maximum number of iterations */
  double tol_err_s = Rcpp::as<double>(tol_err);  /* Tolerance error */


  /*Perform eigenvalue decomposition!!*/

  MatrixXd XTX = X.adjoint()*X;
  const int p(XTX.cols());
  VectorXd z_old = z;
  VectorXd z_new = VectorXd(p).setZero(); /* Initial values for betahat */
  VectorXd XTXz = VectorXd(p).setZero();


  int l = 0;
  double err = 1e6;
  double n_const  = 0;

  while(err > tol_err_s && l < max_iter_s){

    XTXz = XTX*z_old;
    n_const = VectorXd(XTXz).norm();
    z_new =  XTXz/n_const;
    err = VectorXd(z_old - z_new).norm();
    l+= 1;
    z_old = z_new;

  }

  double mu_power_d = double(z_old.squaredNorm());
  double mu_power_n = double(z_old.adjoint()*XTX*z_old);
  double mu_power = mu_power_n/mu_power_d;


  return mu_power;
                           }
