#include <RcppArmadillo.h>

RcppExport SEXP mdr(SEXP X, SEXP fold, SEXP status, SEXP t, SEXP cv, SEXP cvp, SEXP top, SEXP na, SEXP fix);

RcppExport SEXP mdrEnsemble(SEXP X, SEXP fold, SEXP status, SEXP t, SEXP top, SEXP oldstatus, SEXP oldX);

arma::vec classOne(const arma::vec a,const arma::vec status, double T);
arma::mat classTwo(const arma::mat a, const arma::vec status, double T);
arma::cube classThree(const arma::mat a, const arma::vec status, double T);
arma::cube classFour(const arma::mat a, const arma::vec status, double T);
arma::vec evalClassOne(const arma::vec a, const arma::vec one, const arma::vec status);
arma::vec evalClassTwo(const arma::mat a, const arma::mat two, const arma::vec status);
arma::vec evalClassThree(const arma::mat a, const arma::cube three, const arma::vec status);
arma::vec evalClassFour(const arma::mat a, const arma::cube four, const arma::vec status);

arma::vec classifyOne(const arma::vec a, const arma::vec one);
arma::vec classifyTwo(const arma::mat a, const arma::mat one);
arma::vec classifyThree(const arma::mat a, const arma::cube one);
arma::vec classifyFour(const arma::mat a, const arma::cube one);
