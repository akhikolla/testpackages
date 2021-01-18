//[[Rcpp::depends(RcppEigen)]]

#include<RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

//[[Rcpp::export]]
Eigen::VectorXd lp(int n, double p0, double mu, double v0, double dt, double alpha, double c0, Eigen::VectorXd z, Eigen::VectorXd z1, Eigen::VectorXd M0, Eigen::VectorXd x0) {
	VectorXd p(n + 1), v(n + 1);
	p(0) = p0;
	v(0) = v0;
	for (int i = 1; i <= n; i++) {
		p(i) = p(i - 1) + (mu - .5*v(i - 1))*dt + sqrt(v(i - 1)*dt)*z1(i - 1) + M0(i - 1);
		v(i) = c0*(pow(z(i - 1) + sqrt(v(i - 1) / c0)*exp(-.5*alpha*dt), 2) + x0(i - 1));
	}
	return p;
}


//[[Rcpp::export]]

Eigen::VectorXd lp2(int n, double p0, double mu, double v0, double dt, double alpha, double c0, Eigen::VectorXd z, Eigen::VectorXd z1, Eigen::VectorXd x0) {
	VectorXd p(n + 1), v(n + 1);
	p(0) = p0;
	v(0) = v0;
	for (int i = 1; i <= n; i++) {
		p(i) = p(i - 1) + (mu - .5*v(i - 1))*dt + sqrt(v(i - 1)*dt)*z1(i - 1);
		v(i) = c0*(pow(z(i - 1) + sqrt(v(i - 1) / c0)*exp(-.5*alpha*dt), 2) + x0(i - 1));
	}
	return p;
}


//[[Rcpp::export]]

Eigen::VectorXd pvc(int n, double p0, double mt, double beta0, double beta1, double v0, double st, double vxs, Eigen::MatrixXd z, Eigen::VectorXd m) {
	VectorXd p(n + 1), v(n + 1);
	p(0) = p0;
	v(0) = v0;
	for (int i = 1; i <= n; i++) {
		p(i) = p(i - 1) + mt + exp(beta0 + beta1*v(i - 1))*st*z(i - 1, 0) + m(i - 1);
		v(i) = vxs*v(i - 1) + st*z(i - 1, 1);
	}
	return p;
}


//[[Rcpp::export]]

Eigen::VectorXd pvc0(int n, double p0, double mt, double beta0, double beta1, double v0, double st, double vxs, Eigen::MatrixXd z) {
	VectorXd p(n + 1), v(n + 1);
	p(0) = p0;
	v(0) = v0;
	for (int i = 1; i <= n; i++) {
		p(i) = p(i - 1) + mt + exp(beta0 + beta1*v(i - 1))*st*z(i - 1, 0);
		v(i) = vxs*v(i - 1) + st*z(i - 1, 1);
	}
	return p;
}


//[[Rcpp::export]]

Eigen::VectorXd pv2(int n, double mt, double b0, double b1, double b2, double p0, double v10, double v20, double st, Eigen::MatrixXd z, double v1xs, double v2xs, double bv2) {
	VectorXd p(n + 1), v1(n + 1), v2(n + 1);
	p(0) = p0;
	v1(0) = v10;
	v2(0) = v20;
	for (int i = 1; i <= n; i++) {
		p(i) = p(i - 1) + mt - exp(b0 + b1*v1(i - 1) + b2*v2(i - 1))*st*z(i - 1, 0);
		v1(i) = v1xs*v1(i - 1) + st*z(i - 1, 1);
		v2(i) = v2xs*v2(i - 1) + (1 + bv2*v2(i - 1))*st*z(i - 1, 2);
	}
	return p;
}
