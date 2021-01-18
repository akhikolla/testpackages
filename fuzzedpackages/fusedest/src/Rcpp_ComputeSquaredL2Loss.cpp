#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)


SEXP ComputeSquaredl2Loss(Eigen::MatrixXd X_s, Eigen::VectorXd y_s, Eigen::MatrixXd beta_s, SEXP ind_strt, SEXP n_dc){

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;

const Map<VectorXd> ind_strt_s(as<Map<VectorXd> >(ind_strt));
const Map<VectorXd> n_dc_s(as<Map<VectorXd> >(n_dc));

const int m(n_dc_s.size());
const int p(X_s.cols());

double sq_l2_loss = 0;

int n_dc_j = 0;
int u = 0;


for(int j = 0; j < m; j++){

  u = ind_strt_s(j)-1;
  n_dc_j = n_dc_s(j);

  sq_l2_loss += 0.5*VectorXd(y_s.segment(u, n_dc_j)- X_s.block(u, 0, n_dc_j, p)*VectorXd(beta_s.row(j))).array().square().sum();

}

return Rcpp::wrap(sq_l2_loss);

}
