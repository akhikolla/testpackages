#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)

SEXP RcppWolsSolver03(SEXP invXtwX, SEXP Xtwy, SEXP b){

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::Lower;
using Eigen::Upper;
using Eigen::LLT;
using Eigen::VectorXd;
using Eigen::ArrayXXd;

const Map<MatrixXd> invXtwX_s(as<Map<MatrixXd> >(invXtwX));
const Map<VectorXd> Xtwy_s(as<Map<VectorXd> >(Xtwy));
const Map<VectorXd> b_s(as<Map<VectorXd> >(b));

const int d(invXtwX_s.cols());
VectorXd betahat = VectorXd(d);

betahat = invXtwX_s*(Xtwy_s+b_s);
return Rcpp::wrap(betahat);
}
