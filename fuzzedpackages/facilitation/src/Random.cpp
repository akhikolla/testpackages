#include"Random.h"
#include<Rcpp.h>
#include<cmath>


double Random(double max){
	return Rcpp::as<double>(Rcpp::runif(1,0,max));
}

bool Bernoulli(double p){
	return Random(1) < p ;
}

double Exponential(double r){
	Rcpp::NumericVector s = Rcpp::rexp(1,r);
	return Rcpp::as<double>(s);
}

double Normal(double m, double v){
	Rcpp::NumericVector s = Rcpp::rnorm(1,m,v);
	return Rcpp::as<double>(s);
}

short RandomSign(){
	return Bernoulli(0.5)?-1:1;
}

Position RandomDirection(){
	double teta = Random(2*PI);
	return Position(cos(teta),sin(teta));
}
