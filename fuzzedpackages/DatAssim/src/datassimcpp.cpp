#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
#include <algorithm>    // std::reverse
#include <string>
using namespace std;
// [[Rcpp::plugins(cpp11)]]


// This is a simple function using Rcpp that creates an R list
// containing a character vector and a numeric vector.
//
// Learn more about how to use Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//
// and browse examples of code using Rcpp at:
//
//   http://gallery.rcpp.org/
//
// [[Rcpp::export]]

List datassimcpp(arma::mat &pred, arma::mat &var, arma::mat Corr) {
  arma::mat VarDA(var);
  arma::mat PreDA(pred);
  arma::mat a(var.n_cols, var.n_rows, arma::fill::ones);
  //arma::mat Corr(var.n_cols,var.n_cols, arma::fill::ones);
      //double i = 1;
   // if (i>1) {
   // Corr = cor(PreDA);
   // } else {};
  int j;
  for (j = 0; j < var.n_rows; j++) {
    arma::mat pr = pred.row(j);
    arma::mat v = var.row(j);
    arma::mat sd = sqrt(v);
    arma::mat Cov = Corr % (sd.t()*sd);
    arma::mat a_corr(v.n_cols, v.n_cols, arma::fill::zeros);
    int i;
    for (i = 1; i < pred.n_cols; i++) {
    double k = v.n_cols-i;
      a_corr(k, k) = a((i-1),j);
  arma::mat C = Cov(i, span(0, (i-1)));
  arma::mat CC = a_corr(k, span(k, (v.n_cols-1)));
  std::reverse(CC.begin(), CC.end());
  arma::mat cov = C * CC.t();
  arma::mat z = (VarDA(j,(i-1)) - cov)/(VarDA(j,(i-1)) + v(0,i) - 2*cov);
      a(i,j) = z(0,0);
      a_corr((v.n_cols-1-i), span(k, (v.n_cols-1))) = (1-a(i,j))*a_corr(k, span(k, (v.n_cols-1)));
  arma::mat zz = (pow(a(i,j),2)*v(0,i)) + (pow((1-a(i,j)),2)*VarDA(j,(i-1))) + (2*a(i,j)*(1-a(i,j))*cov);
      VarDA(j,i) = zz(0,0);
  double tt = a(i,j)*pr(0,i) + ((1-a(i,j))*PreDA(j,(i-1)));
      PreDA(j,i) = tt;
     }
   }
  List ret;
  ret["weights"] = a.t();
  ret["PreDA"] = PreDA;
  ret["VarDA"] = VarDA;
  ret["Correlation"] = Corr;
  return ret;
}

 
