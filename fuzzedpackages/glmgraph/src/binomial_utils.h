#include <RcppArmadillo.h>

double _Soft(double, double);


 double cycle_binomial_lasso (const arma::vec &, const arma::mat &, const arma::mat & , const arma::vec &,
		arma::vec &, arma::vec &, const arma::uvec &, const arma::vec &, double,double,int);


 double cycle_binomial_MCP (const arma::vec &,const arma::mat &, const arma::mat &,  const arma::vec &,
		arma::vec & ,arma::vec &,  const arma::uvec &,const arma::vec &,double, double, double,const std::string &,int);



