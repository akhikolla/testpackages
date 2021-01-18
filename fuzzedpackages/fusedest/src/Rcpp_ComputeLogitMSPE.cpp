#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)

double ComputeLogitMSPE(Eigen::MatrixXd X, Eigen::VectorXd y, Eigen::MatrixXd beta, Eigen::VectorXd ind_strt, Eigen::VectorXd n_dc){

using Eigen::Map;
using Eigen::VectorXd;

const int m(n_dc.size());
const int p(X.cols());
const int max_n_j(n_dc.maxCoeff());

VectorXd eta = VectorXd(max_n_j).setZero();
VectorXd exp_eta = VectorXd(max_n_j).setZero();
VectorXd mu = VectorXd(max_n_j).setZero();
int n = int(n_dc.sum());

int n_dc_j = 0;
int u = 0;
double y_pred_k = 0;
double mspe_y = 0;

for(int j = 0; j < m; j++){

  u = ind_strt(j)-1;
  n_dc_j = n_dc(j);
  eta.segment(0, n_dc_j) = X.block(u, 0, n_dc_j, p)*VectorXd(beta.row(j));
  exp_eta.segment(0, n_dc_j) = VectorXd(eta.segment(0, n_dc_j).array().exp());
  mu.segment(0, n_dc_j) = exp_eta.segment(0, n_dc_j).cwiseQuotient(VectorXd(1+ exp_eta.segment(0, n_dc_j).array()));

  for(int k = 0; k < n_dc_j; k++){
    y_pred_k = 0;

    if(mu(k) > 0.5){
      y_pred_k = 1;
    }
   mspe_y += (y(u + k) - y_pred_k)*(y(u + k) - y_pred_k);

  }

  eta.setZero();
  exp_eta.setZero();
  mu.setZero();
}

return mspe_y/n;

}
