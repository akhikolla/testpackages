#include<RcppArmadillo.h>
#include<Rmath.h>
// #include<stdio.h>
#include"BVCUtilities.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace R;


double rtnorm(double mu, double sigma){
	double lower = -mu/sigma, ret, a, z;

	if(lower <0){  //naive
		while((ret = R::rnorm(0, 1)) < lower)
			;
	}
	else{  // accept-reject method based on an exponential distribution
		a = (lower + pow((lower*lower + 4),0.5))/2;
		do{
			ret = R::rexp(a) + lower;
			z = exp(-(ret - a)*(ret - a)/2);
		}while(z < R::runif(0,1));
	}
	return (ret * sigma + mu);
}


double rtnorm0(double mu, double sigma){
	double ret;
	while((ret = R::rnorm(mu, sigma)) <= 0);
	return (ret);
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


double rinvGauss(double mu, double lambda)
{
	double b=0.5*mu/lambda;
	double a=mu*b;
	double c=4.0*mu*lambda;
	double d=mu*mu;
	double z=R::rnorm(0,1);
	double random_sample;

	double u=R::runif(0,1);
	double v=z*z;	// Chi-square with 1 df
	double x=mu+a*v-b*sqrt(c*v+d*v*v);	// Smallest root

	if(u <= mu/(mu+x)){ // Pick x with prob mu/(mu+x), else d/x;
		random_sample = x;
	}else{
		random_sample = d/x;
	};
    return(random_sample);

}

arma::vec mvrnormCpp(const arma::vec& mu, const arma::mat& sigma, double tol){
	unsigned int p = mu.n_elem;
	arma::vec eigval;
	arma::mat eigvec;
	// double tol = std::pow(10, -6);
	arma::eig_sym(eigval, eigvec, sigma);
	// arma::eig_gen(eigval, eigvec, sigma);
	if(arma::any(eigval < (-1*tol*std::abs(eigval(eigval.n_elem -1))))){
		std::string error = std::string("covariance matrix is not positive definite");
		throw std::runtime_error(error);
	}
	arma::vec z = arma::randn(p);
	arma::vec rs = mu + eigvec * diagmat(sqrt(arma::clamp(eigval, 0, eigval.max()))) * z;
	return rs;
}

// arma::vec mvrnormCpp(const arma::vec& mu, const arma::mat& sigma){
	// unsigned int p = mu.n_elem;
	// arma::vec eigval;
	// arma::mat eigvec;
	// arma::eig_sym(eigval, eigvec, sigma);
	// double scale = std::pow(10,5);
	// if(arma::any((arma::round(eigval * scale)/scale) <= 0)){
		// std::string error = std::string("covariance matrix is not positive definite");
		// throw std::runtime_error(error);
	// }
	// arma::vec z = arma::randn(p);
	// arma::vec rs = mu + eigvec * diagmat(sqrt(eigval)) * z;
	// return rs;
// }


arma::vec mvrnormCpp(const arma::vec& mu, const arma::mat& sigma){
	unsigned int p = mu.n_elem;
	arma::vec eigval;
	arma::mat eigvec;
	arma::eig_sym(eigval, eigvec, sigma);
	if(arma::any(eigval <= 0)){
		std::string error = std::string("covariance matrix is not positive definite");
		throw std::runtime_error(error);
	}
	arma::vec z = arma::randn(p);
	arma::vec rs = mu + eigvec * arma::diagmat(arma::sqrt(eigval)) * z;
	return rs;
}

