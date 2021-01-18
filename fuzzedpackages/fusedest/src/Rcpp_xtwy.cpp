#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)

SEXP RcppXtwy(SEXP X, SEXP y, SEXP w){

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::Lower;
using Eigen::Upper;
using Eigen::LLT;
using Eigen::VectorXd;
using Eigen::ArrayXXd;

const Map<MatrixXd> X_s(as<Map<MatrixXd> >(X));
const Map<VectorXd> y_s(as<Map<VectorXd> >(y));
const Map<VectorXd> w_s(as<Map<VectorXd> >(w));

const int d(X_s.cols());
const int n(X_s.rows());

MatrixXd wX_s2 = MatrixXd(n,d);

for(int i=0; i< d; i++){
wX_s2.col(i) = VectorXd(X_s.col(i).array()*w_s.array());
}

return Rcpp::wrap(wX_s2.adjoint()*y_s);
}
