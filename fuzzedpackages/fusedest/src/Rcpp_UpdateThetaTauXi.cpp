#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)


SEXP UpdateThetaTauXi(SEXP p,
                      SEXP q_H,
                      Eigen::VectorXd beta01,
                      Eigen::VectorXd beta02,
                      Eigen::VectorXd alpha,
                      Eigen::VectorXd tau,
                      Eigen::VectorXd theta01,
                      Eigen::VectorXd theta02,
                      Eigen::VectorXd xi01,
                      Eigen::VectorXd xi02){

using Eigen::Map;
using Eigen::VectorXd;

int p_s = Rcpp::as<int>(p);
int q_H_s = Rcpp::as<int>(q_H);

int q_Hp(theta01.size());

VectorXd tau01 = VectorXd(q_Hp).setZero();
VectorXd theta03 = VectorXd(q_Hp).setZero();
VectorXd theta04 = VectorXd(q_Hp).setZero();
VectorXd xi03 = VectorXd(q_Hp).setZero();
VectorXd xi04 = VectorXd(q_Hp).setZero();

for(int j = 0; j < q_H_s; j++){

theta03.segment(j*p_s, p_s) = (theta02.segment(j*p_s, p_s) + alpha.segment(j*p_s, p_s) + tau.segment(j*p_s, p_s) +
  beta01.segment(j*p_s, p_s) + xi01.segment(j*p_s, p_s) + theta01.segment(j*p_s, p_s))/3;

  theta04.segment(j*p_s, p_s) = (theta01.segment(j*p_s, p_s) - alpha.segment(j*p_s, p_s) - tau.segment(j*p_s, p_s) +
  beta02.segment(j*p_s, p_s) + xi02.segment(j*p_s, p_s) + theta02.segment(j*p_s, p_s))/3;

tau01.segment(j*p_s, p_s) = tau.segment(j*p_s, p_s) -(theta03.segment(j*p_s, p_s) - theta04.segment(j*p_s, p_s)) + alpha.segment(j*p_s, p_s);

xi03.segment(j*p_s, p_s) = xi01.segment(j*p_s, p_s) - theta03.segment(j*p_s, p_s) + beta01.segment(j*p_s, p_s);
xi04.segment(j*p_s, p_s) = xi02.segment(j*p_s, p_s) - theta04.segment(j*p_s, p_s) + beta02.segment(j*p_s, p_s);

}
return  Rcpp::List::create(Rcpp::Named("theta03") = Rcpp::wrap(theta03),
                             Rcpp::Named("theta04") = Rcpp::wrap(theta04),
                             Rcpp::Named("tau01") = Rcpp::wrap(tau01),
                             Rcpp::Named("xi03") = Rcpp::wrap(xi03),
                             Rcpp::Named("xi04") = Rcpp::wrap(xi04));
                           }
