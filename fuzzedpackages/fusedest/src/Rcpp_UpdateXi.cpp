#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)


Eigen::VectorXd UpdateXi(Eigen::VectorXd xi01,
               Eigen::VectorXd beta,
               Eigen::VectorXd theta){

using Eigen::Map;
using Eigen::VectorXd;

int q_Hp(xi01.size());
VectorXd xi = VectorXd(q_Hp).setZero();

xi = xi01 - theta + beta;

return xi;
                           }
