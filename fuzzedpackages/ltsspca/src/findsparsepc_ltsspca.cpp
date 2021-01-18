// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
// [[Rcpp::export]]




Rcpp::List findsparsePC(arma::mat x, arma::mat b, arma::rowvec mu, double alpha, int l,
                        int N2bis, int Npc, double tol)
{
  int it, iter;
  int n = x.n_rows;
  int p = x.n_cols;
  int h = n - floor(n*alpha);
  int k = b.n_cols;

  double s0,s,sc0,sc,delta,delta_pc;
  arma::mat xnorm;

  arma::mat xcent = x.each_row() - mu;
  arma::mat a = xcent*b;
  arma::mat resid = xcent - a*trans(b);
  arma::vec sqresid_dist = sum((resid % resid),1);
  arma::uvec x_index = sort_index(sqresid_dist,"ascend");
  arma::uvec var_index = find((b != zeros(p))==1);

  s0 = datum::inf;
  s = sum(sqresid_dist(x_index.head(h)))/h;

  arma::mat xprov;
  arma::mat xprov_cent;
  arma::mat aprov;
  arma::mat aprov_norm;
  arma::mat resid_h;
  arma::vec sqresid_h_dist;

  it = 0;  delta = 1;

  while(++it < N2bis &&  delta > tol)
  {
    xprov = x.rows(sort(x_index.head(h),"ascend"));
    xprov_cent = xprov.each_row() - mu;
    iter = 0;
    sc0 = datum::inf;
    sc = s;
    delta_pc = 1 - sc/sc0;
    while(++iter < Npc && delta_pc > tol)
    {
      aprov =  xprov_cent.cols(var_index)*b.rows(var_index);
      aprov_norm =  normalise(aprov,2);
      b =  trans(xprov_cent) * aprov_norm;
      var_index = sort_index(abs(b),"descend");
      b.rows(var_index.tail(p-l)) = zeros<mat>((p-l),k);
      mu = mean(xprov - aprov_norm*trans(b));
      xprov_cent = xprov.each_row() - mu;
      resid_h = xprov_cent;
      resid_h.cols(var_index.head(l))  = xprov_cent.cols(var_index.head(l)) - aprov_norm*trans(b.rows(var_index.head(l)));
      sqresid_h_dist = sum((resid_h % resid_h),1);
      sc = sum(sqresid_h_dist)/h;
      delta_pc = 1 - sc/sc0;
      sc0 = sc;
    }
    b = normalise(b);
    xcent = x.each_row() - mu;
    a = xcent * b;
    resid = xcent - a*trans(b);
    sqresid_dist = sum((resid % resid),1);
    x_index = sort_index(sqresid_dist,"ascend");
    s = sum(sqresid_dist(x_index.head(h)))/h;
    delta = 1-s/s0;
    s0 = s;
  }
  arma::uvec ws = sqresid_dist <= sqresid_dist(x_index(h));
  return Rcpp::List::create(Rcpp::Named("s") = s, Rcpp::Named("b") = b, Rcpp::Named("mu") = mu, Rcpp::Named("a") = a, Rcpp::Named("ws") = ws, Rcpp::Named("aprov") = aprov_norm);
}


