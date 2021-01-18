#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar


Eigen::VectorXd Blockl2Norm(Eigen::VectorXd beta_i, Eigen::VectorXd beta_j, SEXP p, SEXP q_H){

using Eigen::Map;
using Eigen::VectorXd;

int p_s = Rcpp::as<int>(p);
int q_H_s = Rcpp::as<int>(q_H);

int q_Hp(beta_i.size());


VectorXd alpha = VectorXd(q_Hp).setZero();
alpha = beta_i - beta_j;

VectorXd block_l2_norm = VectorXd(q_H_s).setZero();

for(int l = 0; l < q_H_s; l++){

block_l2_norm(l) = alpha.segment(l*p_s, p_s).norm();

  }

return block_l2_norm;
}
