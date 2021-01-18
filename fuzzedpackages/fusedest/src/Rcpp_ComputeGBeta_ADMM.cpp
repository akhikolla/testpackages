#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

SEXP ComputeGBeta_ADMM(SEXP q_H, SEXP p, Eigen::MatrixXd beta_s, Eigen::MatrixXd H_s){

  using Eigen::Map;
  using Eigen::MatrixXd;
  using Eigen::VectorXd;

  int q_H_s = Rcpp::as<int>(q_H);
  int p_s = Rcpp::as<int>(p);

  VectorXd GB = VectorXd(q_H_s*p_s).setZero();

  int ind_i_s = 0;
  int ind_j_s = 0;

  for(int l = 0; l < q_H_s; l++){

        ind_i_s = int(H_s(l,0))-1;
        ind_j_s = int(H_s(l,1))-1;

    GB.segment(l*p_s, p_s) = VectorXd(beta_s.row(ind_i_s) - beta_s.row(ind_j_s));

        ind_i_s = 0;
        ind_j_s = 0;

     }


  return Rcpp::List::create(Rcpp::Named("Gbeta") = Rcpp::wrap(GB));
                           }
