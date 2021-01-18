#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)


SEXP UpdateAlphal2Norm(Eigen::VectorXd theta01,
                 Eigen::VectorXd theta02,
                 Eigen::VectorXd tau,
                 SEXP p,
                 SEXP q_H,
                 SEXP lambda){

using Eigen::Map;
using Eigen::VectorXd;

int p_s = Rcpp::as<int>(p);
int q_H_s = Rcpp::as<int>(q_H);
double lambda_s = Rcpp::as<double>(lambda);

int q_Hp(theta01.size());


VectorXd alpha_old = VectorXd(q_Hp).setZero();
alpha_old = theta01 - theta02 - tau;

VectorXd alpha = VectorXd(q_Hp).setZero();

int l = 0;
double l2_norm_j = 0;

for(int j = 0; j < q_H_s; j++){

l2_norm_j = alpha_old.segment(l, p_s).norm();

  if(l2_norm_j > lambda_s){

    alpha.segment(l, p_s) = (1 - lambda_s/l2_norm_j)*alpha_old.segment(l, p_s);
  }
l += p_s;
  }

return  Rcpp::List::create(Rcpp::Named("alpha") = Rcpp::wrap(alpha));
}
