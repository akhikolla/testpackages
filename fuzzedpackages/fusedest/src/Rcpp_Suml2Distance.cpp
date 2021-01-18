#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)


SEXP Suml2Distance(SEXP p, SEXP q_H, Eigen::VectorXd a, Eigen::VectorXd b){

using Eigen::Map;
using Eigen::VectorXd;

int p_s = Rcpp::as<int>(p);
int q_H_s = Rcpp::as<int>(q_H);

double sum_l2_dist = 0;

for(int j = 0; j < q_H_s; j++){

  sum_l2_dist += VectorXd(a.segment(j*p_s, p_s)-b.segment(j*p_s, p_s)).norm();
}

return  Rcpp::List::create(Rcpp::Named("sum_l2_dist") = Rcpp::wrap(sum_l2_dist));
                           }
