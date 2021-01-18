#ifndef ARMADILLO_MANIPULATIONS_H
#define ARMADILLO_MANIPULATIONS_H
arma::mat direct_sum(arma::field<arma::mat> x);

// faster
arma::mat center_matrix(const arma::mat &x);

#endif
