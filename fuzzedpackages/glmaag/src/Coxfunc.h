#include <RcppArmadillo.h>

int tsign(double, int, double, int);

arma::vec invcumsum(int, arma::vec&);

double cidx(int, arma::vec, arma::vec&, arma::vec&);

void tieup(arma::vec&, arma::vec&, arma::uvec, arma::uvec);

void tiedown(arma::vec&, arma::vec&, arma::uvec&, arma::uvec&);

arma::vec getdd(arma::vec, arma::uvec&, arma::uvec&, arma::uvec&);

arma::vec Coxpieabs(bool, int, arma::uvec&, arma::uvec&, arma::vec&, arma::vec&, arma::mat&, arma::vec&, int, double, arma::mat&);

double devianceCox(int, arma::vec, bool, int, arma::uvec&, arma::uvec&, arma::vec&);

arma::vec Coxnet0(int n, int p, int ntie, arma::vec b, arma::mat& x, bool tie, arma::uvec& tie1, arma::uvec& tie2, arma::vec& del, arma::vec& dd, arma::mat l, int maxiter, double cri);

arma::vec Coxaen(int n, int p, int ntie, arma::vec b, double lam1, double lam2, arma::vec w, arma::mat x, bool tie, arma::uvec& tie1, arma::uvec& tie2, arma::vec& del, arma::vec&dd, arma::vec& dl, int maxiter, double cri);

arma::vec Coxal(int n, int p, int ntie, arma::vec b, double lam, arma::vec w, arma::mat x, bool tie, arma::uvec& tie1, arma::uvec& tie2, arma::vec& del, arma::vec& dd, int maxiter, double cri);

arma::vec Coxaagg0(int n, int p, int ntie, arma::vec b, double lam1, double lam2, arma::vec w, arma::mat x, bool tie, arma::uvec& tie1, arma::uvec& tie2, arma::vec& del, arma::vec&dd, arma::mat l, arma::vec dl, int maxiter, double cri);
