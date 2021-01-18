#ifndef penaliseTEMP_H
#define penaliseTEMP_H
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
using Eigen::MatrixXd;
using Eigen::VectorXd;
#include <RcppGSL.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
// [[Rcpp::depends(RcppGSL)]]

using namespace Rcpp;
using namespace std;
Rcpp::NumericVector gFc(Rcpp::List liste, int l);
struct my_fct_params{
  int n;
  VectorXd Z1;
  NumericVector d;
  NumericVector sqd;
  MatrixXd Q;
  MatrixXd kv;
  double muv;
};

double my_fct (double ro, void *params);
double uniroot(int n,VectorXd Z1, NumericVector d,NumericVector sqd , MatrixXd Q,
               MatrixXd kv , double muv, double t0, double t1);
SEXP optV(int n,VectorXd Z1, NumericVector d,NumericVector sqd ,MatrixXd Q
            ,MatrixXd kv , double muv, double gamav);
struct fparams{
  NumericVector d;
  MatrixXd Q;
  MatrixXd kv;
  VectorXd R;
  double gamav;
  double muv;
};

int rvfct_f (const gsl_vector * x, void *params, gsl_vector * f);
int rvfct_df (const gsl_vector * x, void *params, gsl_matrix * J);
int rvfct_fdf (const gsl_vector * x, void *params, gsl_vector * f, gsl_matrix * J);
SEXP nleqslvgnewton(NumericVector xstart, NumericVector d, MatrixXd Q, MatrixXd kv,
                    VectorXd R, double gamav, double muv);
SEXP nleqslvhybrids(NumericVector xstart, NumericVector d, MatrixXd Q, MatrixXd kv,
                    VectorXd R, double gamav, double muv);
SEXP nleqslvbroyden(NumericVector xstart, NumericVector d, MatrixXd Q, MatrixXd kv,
                    VectorXd R, double gamav, double muv);
SEXP resiv(int n,NumericVector d, MatrixXd Q,MatrixXd kv, VectorXd R
             , double gamav, double muv, NumericVector tetav);
SEXP penMetaMod_cpp(NumericVector Y,List matZ,List k_v,StringVector namG, List resg,
                    NumericVector gamma, NumericVector mu,NumericVector gama_v,
                    NumericVector mu_v,int maxIter,
                    bool verbose, bool calcStwo);
#endif
