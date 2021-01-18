#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)

SEXP UpdateEta_ADMM(Eigen::VectorXd eta01,
               Eigen::VectorXd Gbeta,
               Eigen::VectorXd gamma){

using Eigen::Map;
using Eigen::VectorXd;

int q_H_p(eta01.size());
VectorXd eta = VectorXd(q_H_p).setZero();

eta = eta01 + Gbeta - gamma;

return  Rcpp::List::create(Rcpp::Named("eta") = Rcpp::wrap(eta));
                           }
