#include <RcppArmadillo.h>

arma::vec lrpieabs(const arma::vec&, const arma::mat&, const arma::vec&, int, double, const arma::mat&);

arma::vec lraen(arma::vec b, double lam1, double lam2, const arma::vec& w, const arma::mat& x, const arma::vec& y, const arma::vec& dl, int maxiter, double cri);

arma::vec lral(arma::vec b, double lam1, const arma::vec& w, const arma::mat& x, const arma::vec& y, int maxiter, double cri);
