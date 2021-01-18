#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace Rcpp;
#include "fct_kv_EvaluesQ.h"
#include "grplassoTEMP.h"
#include "penaliseTEMP.h"
// [[Rcpp::export]]
SEXP RKHSMetMod(NumericVector Y, NumericMatrix X, String kernel,
                   int Dmax, NumericVector gamma, NumericVector frc,
                   bool verbose=false){

  List Zpd; Zpd = calc_Kv(X,kernel, Dmax, true,verbose, 1e-08);
  List matZ; matZ = Zpd[0];
  List resK;
  StringVector namG; namG = Zpd[1];
  int vMax; vMax = namG.size();
  NumericVector zerosv(vMax);
  List k_v(vMax);
  int n; n = Y.size();
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
    k_v[v] = kv;
    double nrm; nrm = def.transpose()*kv*def;
    muTEMP(v) = 2*(sqrt(nrm/n));
  }
  int lm; lm = frc.size();
  double mumx = max(muTEMP);
  NumericVector mu(lm);
  List resgl(lm);
  for(int lmu=0; lmu<lm; lmu++){
    double mui; mui = mumx/frc[lmu];
    resgl[lmu] = grplasso(Y, Zpd,k_v, mui,1000,1e-4,verbose);
    mu[lmu] = mui;
  }
  List penMeMod; penMeMod = penMetaMod_cpp(Y,matZ,k_v,namG,resgl,gamma,mu,zerosv,zerosv,1000,
                                       verbose,false);
  return penMeMod;
}//End RKHSMetMod
