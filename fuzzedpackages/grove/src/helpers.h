#ifndef HELPERS_H
#define HELPERS_H

#include "RcppArmadillo.h"
// #include <omp.h>

using namespace Rcpp;
using namespace arma;


double log_exp_x_plus_exp_y(double x, double y);

arma::vec dmvnrm_arma(arma::mat x,  
                      arma::rowvec mean,  
                      arma::mat sigma, 
                      bool logd = true);

arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma);

int sampling(int Kf, arma::vec probs);

arma::uvec integerToBinary(int s, int l = 16);





#endif
