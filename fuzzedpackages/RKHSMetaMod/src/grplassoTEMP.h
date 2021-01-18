#ifndef grplassoTEMP_H
#define grplassoTEMP_H
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace Rcpp;
using namespace std;
struct tetafq_params
{
  NumericVector d;
  NumericVector sqd;
  MatrixXd Q;
  VectorXd R;
  double mqn;
};
double tetafq (double ro, void *params);
SEXP slvq(NumericVector d,NumericVector sqd, MatrixXd Q, VectorXd R, double mqn,
          double t0, double t1);
SEXP grplasso(NumericVector Y, List Kv, List k_v,double mu, int maxIter, double eps,
              bool verbose);
#endif
