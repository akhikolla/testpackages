#include<RcppArmadillo.h>
#include<Rmath.h>
// #include<stdio.h>
#include"BVCUtilities.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
//using namespace R;
// using Eigen::Map;                       // 'maps' rather than copies
// using Eigen::MatrixXd;                  // variable size matrix, double precision
// using Eigen::VectorXd;                  // variable size vector, double precision
// using Eigen::SelfAdjointEigenSolver;    // one of the eigenvalue solvers

arma::vec mvrnormChol(const arma::vec& mu, const arma::mat& sigma) {
   arma::vec Y = arma::randn(sigma.n_cols);
   return mu + arma::chol(sigma) * Y;
}

double rinvgaussian(double mu, double lambda){
	if(mu>1000){
		mu = 1000;
	}
	double random_sample;
    double z,y,x,u;
	z=R::rnorm(0,1);
	y=z*z;
	x=mu+0.5*mu*mu*y/lambda - 0.5*(mu/lambda)*sqrt(4*mu*lambda*y+mu*mu*y*y);
	u=R::runif(0,1);
	if(u <= mu/(mu+x)){
		random_sample = x;
	}else{
		random_sample = mu*mu/x;
	};
    return(random_sample);
}


arma::vec mvrnormCpp(const arma::vec& mu, const arma::mat& sigma, double tol){
	unsigned int p = mu.n_elem;
	arma::vec eigval;
	arma::mat eigvec;
	// double tol = std::pow(10, -6);
	arma::eig_sym(eigval, eigvec, sigma);
	if(arma::any(eigval < (-1*tol*std::abs(eigval(eigval.n_elem -1))))){
		std::string error = std::string("sigma is not positive definite");
		throw std::runtime_error(error);
	}
	arma::vec z = arma::randn(p);
	arma::vec rs = mu + eigvec * diagmat(sqrt(arma::clamp(eigval, 0, eigval.max()))) * z;
	return rs;
}

arma::vec mvrnormCpp(const arma::vec& mu, const arma::mat& sigma){
	unsigned int p = mu.n_elem;
	arma::vec eigval;
	arma::mat eigvec;
	arma::eig_sym(eigval, eigvec, sigma);
	if(arma::any(eigval <= 0)){
		std::string error = std::string("sigma is not positive definite");
		throw std::runtime_error(error);
	}
	arma::vec z = arma::randn(p);
	arma::vec rs = mu + eigvec * diagmat(sqrt(eigval)) * z;
	return rs;
}

// arma::vec mvrnormEigen(arma::vec& mu, arma::mat& sigma){ 	
	// unsigned int p = mu.n_elem;
	// Eigen::SelfAdjointEigenSolver<MatrixXd> es(Eigen::Map<Eigen::MatrixXd>(sigma.memptr(),
                                                        // sigma.n_rows,
                                                        // sigma.n_cols));
	// // Eigen::SelfAdjointEigenSolver<MatrixXd> es(sigma_eigen);
	// Eigen::MatrixXd eigenvectors = es.eigenvectors();
	// Eigen::VectorXd eigenvalues = es.eigenvalues();
	// arma::vec eigval = arma::vec(eigenvalues.data(), eigenvalues.size(), false, false);
	// arma::mat eigvec = arma::mat(eigenvectors.data(), eigenvectors.rows(), eigenvectors.cols(), false, false);
	// double tol = std::pow(10, -6);
	// // eig_sym(eigval, eigvec, sigma);
	// if(arma::any(eigval < (-1*tol*std::abs(eigval(eigval.n_elem -1))))){
		// std::string error = std::string("sigma is not positive definite");
		// throw std::runtime_error(error);
	// }
	// arma::vec z = arma::randn(p);
	// arma::vec rs = mu + eigvec * arma::diagmat(sqrt(arma::clamp(eigval, 0, eigval.max()))) * z;
	// return rs;
// }
