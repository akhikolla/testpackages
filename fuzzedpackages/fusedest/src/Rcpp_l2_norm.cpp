#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)

SEXP Rcppl2norm(SEXP x){

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;

const Map<VectorXd> x_s(as<Map<VectorXd> >(x));

double l2_norm_x;
l2_norm_x = x_s.norm();

return Rcpp:: wrap(l2_norm_x);
}
