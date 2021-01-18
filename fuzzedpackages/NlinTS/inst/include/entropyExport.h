#ifndef ENTROPY_H
#define ENTROPY_H

#include <Rcpp.h>
// Discrete variables
double entropy_disc (Rcpp::IntegerVector  & ts_, std::string log = "log");
double joinEntropy_disc (Rcpp::DataFrame & Df, std::string log = "log");
double mutualInformation_disc_u (Rcpp::IntegerVector & X,  Rcpp::IntegerVector & Y, std::string log = "log", bool normalize = false);
double mutualInformation_disc (Rcpp::DataFrame & Df, std::string log = "log", bool normalize = false);
double transferEntropy_disc (Rcpp::IntegerVector & X, Rcpp::IntegerVector & Y, int p, int q, std::string log = "log", bool normalize = false);

// Continuous variables
double entropy_cont (Rcpp::NumericVector & I, int k = 3, std::string log = "loge");
double mutualInformation_cont ( Rcpp::DataFrame & Df, int k = 3, std::string alg = "ksg1", bool normalize = false);
double transferEntropy_cont ( Rcpp::NumericVector & I,  Rcpp::NumericVector & J, int p = 1, int q = 1, int k = 3, bool normalize = false);

#endif
