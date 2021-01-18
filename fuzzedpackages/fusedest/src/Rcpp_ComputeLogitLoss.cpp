#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)

double ComputeLogitLoss(Eigen::MatrixXd X, Eigen::VectorXd y, Eigen::MatrixXd beta, SEXP ind_strt, SEXP n_dc){

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;

const Map<VectorXd> ind_strt_s(as<Map<VectorXd> >(ind_strt));
const Map<VectorXd> n_dc_s(as<Map<VectorXd> >(n_dc));

const int m(n_dc_s.size());
const int p(X.cols());
const int max_n_j(n_dc_s.maxCoeff());
VectorXd eta = VectorXd(max_n_j).setZero();

double logit_loss = 0;
double logit_loss_1 = 0;
double logit_loss_2 = 0;

int n_dc_j = 0;
int u = 0;

for(int j = 0; j < m; j++){

  u = ind_strt_s(j)-1;
  n_dc_j = n_dc_s(j);
  eta.segment(0, n_dc_j) = X.block(u, 0, n_dc_j, p)*VectorXd(beta.row(j));
  logit_loss_1 = -y.segment(u, n_dc_j).transpose()*eta.segment(0, n_dc_j);
  logit_loss_2 = VectorXd(1 + eta.segment(0, n_dc_j).array().exp()).array().log().sum();
  logit_loss += logit_loss_1 + logit_loss_2;

  eta.setZero();
}

return logit_loss;

}
