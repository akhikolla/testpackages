#ifndef rowcolSumsiter_H
#define rowcolSumsiter_H

#include <RcppArmadillo.h>

arma::vec colSumsiter(const arma::sp_mat& x);
arma::vec rowSumsiter(const arma::sp_mat& x);

#endif
