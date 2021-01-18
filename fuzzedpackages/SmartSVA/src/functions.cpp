#include <RcppEigen.h>

// [[Rcpp::export]]
Eigen::MatrixXd prodCpp(const Eigen::Map<Eigen::MatrixXd> & B, const Eigen::Map<Eigen::MatrixXd> & C) {

		Eigen::MatrixXd BtC(B * C);
		return BtC;
}

// [[Rcpp::export]]
Eigen::MatrixXd crossprodCpp(const Eigen::Map<Eigen::MatrixXd> & A) {
		const int n(A.cols());
    Eigen::MatrixXd AtA(Eigen::MatrixXd(n, n).setZero().selfadjointView<Eigen::Lower>().rankUpdate(A.adjoint()));
		return AtA;
}

// [[Rcpp::export]]
Eigen::MatrixXd tcrossprodCpp(const Eigen::Map<Eigen::MatrixXd> & B, const Eigen::Map<Eigen::MatrixXd> & C) {

		Eigen::MatrixXd BtC(B * C.adjoint());
		return BtC;
}
