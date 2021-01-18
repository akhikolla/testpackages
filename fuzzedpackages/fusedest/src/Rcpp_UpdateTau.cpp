#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)


Eigen::VectorXd UpdateTau(Eigen::VectorXd tau01,
               Eigen::VectorXd theta01,
               Eigen::VectorXd theta02,
               Eigen::VectorXd alpha){

using Eigen::Map;
using Eigen::VectorXd;

int q_Hp(tau01.size());
VectorXd tau = VectorXd(q_Hp).setZero();

tau = tau01 - (theta01 - theta02) + alpha;

return tau;
                           }
