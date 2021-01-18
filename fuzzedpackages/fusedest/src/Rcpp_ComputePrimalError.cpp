#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)

double ComputePrimalError(SEXP p,
                        SEXP q_H,
                        Eigen::VectorXd theta01,
                        Eigen::VectorXd theta02,
                        Eigen::VectorXd theta03,
                        Eigen::VectorXd theta04,
                        Eigen::VectorXd beta01,
                        Eigen::VectorXd beta02,
                        Eigen::VectorXd alpha){

using Eigen::Map;
using Eigen::VectorXd;

int p_s = Rcpp::as<int>(p);
int q_H_s = Rcpp::as<int>(q_H);

double primal_err = 0;
double primal_err_j = 0;

for(int j = 0; j < q_H_s; j++){

  primal_err_j = VectorXd(theta01.segment(j*p_s, p_s)-theta03.segment(j*p_s, p_s)).norm() +
    VectorXd(theta02.segment(j*p_s, p_s)-theta04.segment(j*p_s, p_s)).norm() +
    VectorXd(beta01.segment(j*p_s, p_s)-theta03.segment(j*p_s, p_s)).norm() +
    VectorXd(beta02.segment(j*p_s, p_s)-theta04.segment(j*p_s, p_s)).norm()+
    VectorXd(alpha.segment(j*p_s, p_s)-(theta03.segment(j*p_s, p_s)-theta04.segment(j*p_s, p_s))).norm();
  primal_err += primal_err_j/(5*sqrt((double)p_s)*q_H_s);

  primal_err_j = 0;
}

return primal_err;
                           }
