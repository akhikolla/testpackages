#include <RcppArmadillo.h>


double auc(int, const arma::vec&, const arma::vec&);

double logi1(double);

arma::vec logipieabs(const arma::vec&, const arma::mat&, double, const arma::vec&, int, double, const arma::mat&);

arma::vec logiaen(double b0, arma::vec b, double lam1, double lam2, const arma::vec& w, const arma::mat& x, const arma::vec& y, const arma::vec& dl, bool intercept, int maxiter, double cri);

arma::vec logial(double b0, arma::vec b, double lam1, const arma::vec& w, const arma::mat& x, const arma::vec& y, bool intercept, int maxiter, double cri);
