#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)


SEXP IRLSLogisticReg(Eigen::MatrixXd X,
                                Eigen::VectorXd y,
                                SEXP a,
                                Eigen::VectorXd b,
                                Eigen::VectorXd beta_ini,
                                SEXP max_iter,
                                SEXP tol_err){

  using Eigen::Map;
  using Eigen::MatrixXd;
  using Eigen::VectorXd;

       const double a_s = Rcpp::as<double>(a);
   const int max_iter_s = Rcpp::as<int>(max_iter);
 const double tol_err_s = Rcpp::as<double>(tol_err);

  const int n(y.size());
  const int p(X.cols());

  VectorXd beta_old = beta_ini;
  VectorXd beta = VectorXd(p).setZero();
  VectorXd theta = VectorXd(p).setZero();

  VectorXd eta = VectorXd(n).setZero();
  VectorXd exp_eta = VectorXd(n).setZero();
  VectorXd mu = VectorXd(n).setZero();
  VectorXd Psi = VectorXd(n).setZero();
  VectorXd W = VectorXd(n).setZero();
  VectorXd h = VectorXd(n).setZero();

  MatrixXd sqrt_WX = MatrixXd(n, p).setZero();
  MatrixXd XTWX_plus_a = MatrixXd(p, p).setZero();
  MatrixXd WX = MatrixXd(n, p).setZero();
  MatrixXd Id = MatrixXd(p, p).setIdentity();
  MatrixXd LT = MatrixXd(p, p).setZero();
  VectorXd err_list=VectorXd(max_iter_s).setZero();

  /* Run iteration */

  double err = 1e6;
  int l = 0;

  while(err > tol_err_s && l < max_iter_s){

    sqrt_WX.setZero();
    WX.setZero();

      eta = X*beta_old;
      exp_eta = VectorXd(eta.array().exp());
      mu = exp_eta.cwiseQuotient(VectorXd(1+ exp_eta.array()));
      Psi = VectorXd(1+exp_eta.array()).cwiseQuotient(VectorXd(mu.array()+0.000001)); /*VectorXd(1/((1-mu.array())*mu.array()));*/
        W = Psi.cwiseInverse(); /*VectorXd(mu.array()*(1-mu.array()));*/
        h = VectorXd(Psi.array()*VectorXd(y-mu).array()) + eta;

  /* Computing XTWX and WX */

  for(int i = 0; i < p; i++){
    sqrt_WX.col(i)= VectorXd(X.col(i).array()*W.array().sqrt());
    WX.col(i) = VectorXd(X.col(i).array()*W.array());
  }

  XTWX_plus_a = sqrt_WX.adjoint()*sqrt_WX + a_s*Id; /*sqrt_WX.adjoint()*sqrt_WX;*/


  /* Choleskey forward backward substitution */

  /*Eigen::LLT<MatrixXd, Lower> Chol_LLT(XTWX_plus_a);*/

  LT = XTWX_plus_a.llt().matrixL().adjoint();
  theta = LT.adjoint().triangularView<Lower>().solve(WX.adjoint()*h + b);
  beta = LT.triangularView<Upper>().solve(theta);

  LT.setZero();

     /* Computing error */

    err = VectorXd(beta_old-beta).norm()/p;
    err_list(l) = err;
    beta_old = beta;

    l += 1;
  }

  return  Rcpp::List::create(Rcpp::Named("beta") = Rcpp::wrap(beta),
                             Rcpp::Named("iter_count") = Rcpp::wrap(l),
                             Rcpp::Named("err_list") = Rcpp::wrap(err_list.segment(0,l)));
}
