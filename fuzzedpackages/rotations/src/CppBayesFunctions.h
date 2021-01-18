#ifndef CPPBAYESFUNCTIONS_H
#define CPPBAYESFUNCTIONS_H

#include <RcppArmadillo.h>
// [[Rcpp::interfaces(r, cpp)]]

//////////////////////////////////////////////////////
// Generate matrix Fisher random deviates using C++ //
//////////////////////////////////////////////////////

// [[Rcpp::export]]
Rcpp::NumericVector rcayleyCpp(unsigned int n, double kappa);

//////////////////////////////////////////////////////////
// Generate Maxwell Boltzmann random deviates using C++ //
//////////////////////////////////////////////////////////

double dmbCpp(double r, double kappa);
double arsample_mb_unifCpp(double M, double kappa);
Rcpp::NumericVector rar_mb_Cpp(unsigned int n, double kappa, double M);

// [[Rcpp::export]]
Rcpp::NumericVector rmbCpp(unsigned int n, double kappa);

//////////////////////////////////////////////////////
// Generate matrix Fisher random deviates using C++ //
//////////////////////////////////////////////////////

double dfisherCpp(double r, double kappa);
double arsample_unifCpp(double M, double kappa);
Rcpp::NumericVector rarCpp(unsigned int n, double kappa, double M);

// [[Rcpp::export]]
Rcpp::NumericVector rfisherCpp(unsigned int n, double kappa);

//////////////////////////////////////////////////
// Generate von Mises random deviates using C++ //
//////////////////////////////////////////////////

int sign(double x);

// [[Rcpp::export]]
Rcpp::NumericVector rvmisesCPP(unsigned int n, double kappa);

// [[Rcpp::export]]
arma::mat centerCpp(const arma::mat &Rs, const arma::mat &S);

//////////////////////////////
// log likelihood functions //
//////////////////////////////

// [[Rcpp::export]]
double lpvmises(const arma::mat &Rs, const arma::mat &S, double kappa);

// [[Rcpp::export]]
double lpfisher(const arma::mat &Rs, const arma::mat &S, double kappa);

// [[Rcpp::export]]
double lpcayley(const arma::mat &Rs, const arma::mat &S, double kappa);

arma::mat genrC(const arma::mat &S, double r);

////////////////////////////
// Actual Bayes functions //
////////////////////////////

// [[Rcpp::export]]
arma::mat S_MCMC_CPP(
    const arma::mat &Rs,
    const arma::mat &oldS,
    double rho,
    double kappa,
    int Dist
);

// [[Rcpp::export]]
double kap_MCMC_CPP(
    const arma::mat &Rs,
    double oldKappa,
    double sigma,
    const arma::mat &S,
    int Dist
);

// [[Rcpp::export]]
arma::rowvec afun_CPP(const arma::mat &R1, const arma::mat &R2);

// [[Rcpp::export]]
Rcpp::List both_MCMC_CPP(
    const arma::mat &Rs,
    arma::mat S0,
    double kappa0,
    double rho,
    double sigma,
    int burnin,
    int B,
    int Dist
);

#endif /* CPPBAYESFUNCTIONS_H */
