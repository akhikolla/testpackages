/* ********************************************

   FIXME: sp_setXtf1 bruges i gRain (i .setRoot); skal
   re-implementeres i gRain.

 * ********************************************* */

#include <RcppEigen.h>

using namespace Rcpp;

//[[Rcpp::export(.sparse_setXtf1)]]
SEXP sparse_setXtf1 ( SEXP XX_, SEXP TF_){
  using Eigen::Map;
  using Eigen::MatrixXi;
  typedef Eigen::SparseMatrix<double> MSpMat;

  MSpMat   X(as<MSpMat>(XX_));
  const Map<MatrixXi> TF(as<Map<MatrixXi> >(TF_));
  int nrTF(TF.rows());
  int rr, ii, jj;
  for (rr=0; rr<nrTF; rr++){
    jj = (TF.coeff(rr,0))-1;
    ii = (TF.coeff(rr,1))-1;
    X.coeffRef(ii,jj) = 1;
  }
  X.makeCompressed();
  return wrap( X );
}
