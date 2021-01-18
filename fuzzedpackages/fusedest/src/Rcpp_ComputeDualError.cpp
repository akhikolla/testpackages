#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)


double ComputeDualError(SEXP p,
                        SEXP q_H,
                        Eigen::VectorXd xi01,
                        Eigen::VectorXd xi02,
                        Eigen::VectorXd tau){

using Eigen::Map;
using Eigen::VectorXd;

int p_s = Rcpp::as<int>(p);
int q_H_s = Rcpp::as<int>(q_H);

double dual_err = 0;
double dual_err_j = 0;

for(int j = 0; j < q_H_s; j++){

  dual_err_j = VectorXd(xi01.segment(j*p_s, p_s) + tau.segment(j*p_s, p_s)).norm()+
    VectorXd(xi02.segment(j*p_s, p_s) - tau.segment(j*p_s, p_s)).norm();
  dual_err += dual_err_j/(2*sqrt((double)p_s)*q_H_s);

  dual_err_j = 0;
}

return dual_err;
                           }
