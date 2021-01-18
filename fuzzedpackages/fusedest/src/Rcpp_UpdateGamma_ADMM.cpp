#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)


SEXP UpdateGamma_ADMM(SEXP q_H, SEXP p, Eigen::VectorXd Gbeta_s, Eigen::VectorXd eta_s, SEXP lambda){

  using Eigen::Map;
  using Eigen::VectorXd;

  int q_H_s = Rcpp::as<int>(q_H);
  int p_s = Rcpp::as<int>(p);
  double lambda_s = Rcpp::as<double>(lambda);

  VectorXd gamma = VectorXd(q_H_s*p_s).setZero();
  VectorXd Gbeta_plus_eta_l = VectorXd(p_s).setZero();

  double l2_norm_l = 0;

  for(int l = 0; l < q_H_s; l++){

    Gbeta_plus_eta_l = VectorXd(Gbeta_s.segment(l*p_s, p_s) + eta_s.segment(l*p_s, p_s));
    l2_norm_l = Gbeta_plus_eta_l.norm();

    if(l2_norm_l > lambda_s){

       gamma.segment(l*p_s, p_s) = (1 - lambda_s/l2_norm_l)*Gbeta_plus_eta_l;

                                  }

    l2_norm_l = 0;
    Gbeta_plus_eta_l.setZero();

     }


  return Rcpp::List::create(Rcpp::Named("gamma") = Rcpp::wrap(gamma));
                           }
