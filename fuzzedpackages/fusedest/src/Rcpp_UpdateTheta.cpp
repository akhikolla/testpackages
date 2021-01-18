#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)


Eigen::VectorXd UpdateTheta(Eigen::VectorXd theta01,
                 Eigen::VectorXd alpha,
                 Eigen::VectorXd tau,
                 Eigen::VectorXd beta,
                 Eigen::VectorXd xi,
                 Eigen::VectorXd theta02){

using Eigen::Map;
using Eigen::VectorXd;

int q_Hp(theta01.size());
VectorXd theta = VectorXd(q_Hp).setZero();
theta = (theta01 + alpha + tau + beta + xi + theta02)/3;

return  theta;
                           }
