#ifndef __series_weibull_h__
#define __series_weibull_h__

#include <RcppArmadillo.h>

arma::mat alphagen(double cc , unsigned jrow, unsigned ncol);

arma::vec dWeibullCount_mat(arma::Col<unsigned> x,
			    double shape, double scale,
			    double time, bool logFlag, unsigned jmax);

double dWeibullCount_mat_loglik(arma::Col<unsigned> x,
				arma::vec shape, arma::vec scale,
				double time, unsigned jmax);

arma::vec dWeibullgammaCount_mat(arma::Col<unsigned> x, double shape,
				 double r, double alpha,
				 double time, bool logFlag,
				 unsigned jmax);

double dWeibullgammaCount_mat_loglik(arma::Col<unsigned> x,
				     arma::vec shape,
				     arma::vec r, arma::vec alpha,
				     double time, unsigned jmax);

arma::vec dWeibullCount_acc(arma::Col<unsigned> x, double shape, double scale, 
			    double time, bool logFlag, unsigned jmax,
			    int nmax, double eps, bool printa);

double dWeibullCount_acc_loglik(arma::Col<unsigned> x,
				arma::vec shape, arma::vec scale, 
				double time, bool logFlag, unsigned jmax,
				int nmax, double eps, bool printa);

arma::vec dWeibullgammaCount_acc(arma::Col<unsigned> x, double shape,
				 double r, double alpha, double time,
				 bool logFlag, unsigned jmax,
				 int nmax, double eps,
				 bool printa);

double dWeibullgammaCount_acc_loglik(arma::Col<unsigned> x, arma::vec shape,
				     arma::vec r, arma::vec alpha,
				     double time,
				     bool logFlag, unsigned jmax,
				     int nmax, double eps,
				     bool printa);

#endif 
