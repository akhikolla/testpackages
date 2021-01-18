#include <RcppArmadillo.h>

double _Soft(double, double);



 double cycle_gaussian_lasso (const arma::vec &, const arma::mat &,  const arma::mat & , const arma::vec &,  
		arma::vec &,  const arma::uvec &, const arma::vec &, double,double,bool,int);


 double cycle_gaussian_MCP (const arma::vec &,const arma::mat &, const arma::mat &, const arma::vec &,
		arma::vec & , const arma::uvec &,const arma::vec &,double, double, double,bool,int);



