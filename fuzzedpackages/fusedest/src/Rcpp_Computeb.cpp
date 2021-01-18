#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)


SEXP Computeb(SEXP p, Eigen::VectorXd theta, Eigen::VectorXd xi,  SEXP o_edge_ext, SEXP deg_dc, SEXP ind_edge_strt){

using Eigen::Map;
using Eigen::VectorXd;
using Eigen::MatrixXd;

int p_s = Rcpp::as<int>(p);
const Map<VectorXd> o_edge_ext_s(as<Map<VectorXd> >(o_edge_ext));
const Map<VectorXd> deg_dc_s(as<Map<VectorXd> >(deg_dc)); /* Vector of group sizes */
const Map<VectorXd> ind_edge_strt_s(as<Map<VectorXd> >(ind_edge_strt));

const int m(deg_dc_s.size());
const int q_H2(o_edge_ext_s.size());

MatrixXd b_list = MatrixXd(q_H2, p_s).setZero();
MatrixXd b_s = MatrixXd(m, p_s).setZero();

int u = 0;
int l = 0;

for(int i = 0; i < q_H2; i++){

  l = o_edge_ext_s(i)-1;
  b_list.row(i) = VectorXd(theta.segment(l*p_s, p_s)-xi.segment(l*p_s, p_s));
  /*l += p_s;*/
}

int deg_j = 0;
u = 0;

for(int j = 0; j < m; j++){

  u = ind_edge_strt_s(j)-1;
  deg_j = deg_dc_s(j);

  b_s.row(j) = b_list.block(u, 0, deg_j, p_s).colwise().sum();
}

return  Rcpp::List::create(Rcpp::Named("b") = Rcpp::wrap(b_s));
                           }
