#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace Rcpp;

// [[Rcpp::export]]
double mu_max(NumericVector Y, List matZ){
  int n; n = Y.size();
  List resK;
  int vMax = matZ.size();
  double meany; meany = mean(Y);
  NumericVector y; y = Y-meany;
  VectorXd def; def = as<VectorXd>(y);
  NumericVector muTEMP(vMax);
  for(int v=0;v<vMax;v++){
    resK = matZ[v];
    NumericVector d; d = resK["Evalues"];
    MatrixXd Q; Q = resK["Q"];
    NumericMatrix D; D = diag(d);
    MatrixXd kv; kv = Q * as<MatrixXd>(D) * Q.transpose();
    double nrm; nrm = def.transpose()*kv*def;
    muTEMP(v) = 2*(sqrt(nrm/n));
  }
  return(max(muTEMP));
}//End mu_max
