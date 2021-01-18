#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)

SEXP Rcppl2normProx(SEXP x, SEXP l2normx, SEXP lambda){
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::ArrayXXd;
using Eigen::RowVectorXi;

const Map<VectorXd> x_s(as<Map<VectorXd> >(x));

double lambda_s = Rcpp::as<double>(lambda);
double l2normx_s = Rcpp::as<double>(l2normx);

const int p(x_s.size());

VectorXd l2_norm_x(p);
VectorXd l2_norm_prox(p);

l2_norm_x.setZero();
l2_norm_prox.setZero();


for(int j=0; j<p; j++){
l2_norm_x(j) = l2normx_s;
}

l2_norm_prox = VectorXd(VectorXd(1-lambda_s/l2_norm_x.array()).array().max(0)*x_s.array());


return Rcpp:: wrap(l2_norm_prox);
}
