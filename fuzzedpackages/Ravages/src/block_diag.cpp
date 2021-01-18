#include <RcppEigen.h>
using namespace Rcpp;

//[[Rcpp::export]]
Eigen::MatrixXd block_diag(List L){
  Eigen::MatrixXd L1 = L[0];
  if(L.size() == 1)
    return L1;
  int ngpe = L.size();
  List res(ngpe);
  //For the first group
  Eigen::MatrixXd ZeroesBottom = Eigen::MatrixXd::Zero(L1.rows()*(ngpe-1), L1.cols());
  Eigen::MatrixXd tmp(L1.rows() + ZeroesBottom.rows(), L1.cols());
  tmp << L1, ZeroesBottom;
  res[0] = tmp;
  
  for(int i = 1; i <= ngpe-1; i++){
    Eigen::MatrixXd Ltmp = L[i];
    Eigen::MatrixXd ZeroesTop = Eigen::MatrixXd::Zero(L1.rows()*i, L1.cols());
    Eigen::MatrixXd ZeroesBottom = Eigen::MatrixXd::Zero(L1.rows()*(ngpe-i-1), L1.cols());
    Eigen::MatrixXd tmp(ZeroesTop.rows() + Ltmp.rows() + ZeroesBottom.rows(), Ltmp.cols());
    tmp << ZeroesTop, Ltmp, ZeroesBottom;
    Eigen::MatrixXd tomergeLeft = res[i-1];
    Eigen::MatrixXd tmp2(tomergeLeft.rows(), tomergeLeft.cols()+tmp.cols());
    tmp2 << tomergeLeft, tmp;
    res[i] = tmp2;
  }
  return res[ngpe-1];
}

RcppExport SEXP block_diag(SEXP LSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type L(LSEXP);
    rcpp_result_gen = Rcpp::wrap(block_diag(L));
    return rcpp_result_gen;
END_RCPP
}
