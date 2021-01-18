#ifndef _fdmvnorm_H
#define _fdmvnorm_H
#include <RcppArmadillo.h>
RcppExport SEXP fast_mvnorm_density(SEXP X_mu_diff, SEXP Sigma, SEXP Miss_group_unique, SEXP Miss_group_counts);
#endif
