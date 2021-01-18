#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)
// For more on using Rcpp click the Help button on the editor toolbar


SEXP RcppInvGram(SEXP X, SEXP w, SEXP lambda){

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::Lower;
using Eigen::Upper;
using Eigen::LLT;
using Eigen::VectorXd;
using Eigen::ArrayXXd;

const Map<MatrixXd> X_s(as<Map<MatrixXd> >(X));
const Map<VectorXd> w_s(as<Map<VectorXd> >(w));
double lambda_s = Rcpp::as<double>(lambda);

const int d(X_s.cols());
const int n(X_s.rows());

MatrixXd wX_s = MatrixXd(n,d);
MatrixXd XtwX = lambda_s*MatrixXd(d,d).setIdentity();
MatrixXd XtwX_inv = MatrixXd(d,d);
MatrixXd Id = MatrixXd(d,d).setIdentity();

for(int i=0; i< d; i++){
wX_s.col(i)= VectorXd(X_s.col(i).array()*w_s.array().sqrt());
}

XtwX+= MatrixXd(d,d).setZero().selfadjointView<Lower>().rankUpdate(wX_s.adjoint());

//const LLT<MatrixXd> llt(XtwX+lambda_s*id);

XtwX_inv = XtwX.llt().solve(Id);
return Rcpp::wrap(XtwX_inv);
}
