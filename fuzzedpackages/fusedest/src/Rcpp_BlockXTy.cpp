#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

SEXP ComputeBlockXTy(Eigen::MatrixXd X_s,
                Eigen::VectorXd y_s,
                SEXP ind_strt,
                SEXP n_dc){

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;

const Map<VectorXd> ind_strt_s(as<Map<VectorXd> >(ind_strt));
const Map<VectorXd> n_dc_s(as<Map<VectorXd> >(n_dc));

const int m(n_dc_s.size());
const int p(X_s.cols());

VectorXd XTy = VectorXd(m*p).setZero();

int n_dc_j = 0;
int u = 0;


for(int j = 0; j < m; j++){

  u = ind_strt_s(j)-1;
  n_dc_j = n_dc_s(j);

  XTy.segment(j*p, p) = X_s.block(u, 0, n_dc_j, p).adjoint()*y_s.segment(u, n_dc_j);

}

return Rcpp::wrap(XTy);

}
