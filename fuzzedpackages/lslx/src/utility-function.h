// utility functions for package lslx
// written by Po-Hsien Huang psyphh@gmail.com

#include <RcppEigen.h>
#include <limits>
#include <cmath>
#include <algorithm>
#include <cfloat>
#include <string.h>
using namespace Rcpp;
using namespace Eigen;


// [[Rcpp::depends(RcppEigen)]]
// slice columns
Eigen::MatrixXd slice_col(Eigen::MatrixXd x, Rcpp::IntegerVector col_idx);

// slice rows
Eigen::MatrixXd slice_row(Eigen::MatrixXd x, Rcpp::IntegerVector row_idx);

// slice both rows and columns
Eigen::MatrixXd slice_both(Eigen::MatrixXd x, 
                           Rcpp::IntegerVector row_idx, 
                           Rcpp::IntegerVector col_idx);

// expand columns
Eigen::MatrixXd expand_col(Eigen::MatrixXd x, Rcpp::IntegerVector col_idx, int n_col);

// expand both rows and columns
Eigen::MatrixXd expand_both(Eigen::MatrixXd x, 
                            Rcpp::IntegerVector row_idx, 
                            Rcpp::IntegerVector col_idx,
                            int n_row,
                            int n_col);

// vec operator
Eigen::MatrixXd vec(Eigen::MatrixXd x);

// vech operator
Eigen::MatrixXd vech(Eigen::MatrixXd x);

// vech operator for only non-diagonal elements
Eigen::MatrixXd vech_small(Eigen::MatrixXd x);

// create vech idx
Rcpp::IntegerVector create_idx_vech(int n, bool diag);

// create transposed vech idx
Rcpp::IntegerVector create_idx_tvech(int n, bool diag);


// find index of intersected elements
Rcpp::IntegerVector find_idx_match(Rcpp::IntegerVector x, Rcpp::IntegerVector y);

// method for deduplifying both sides
Eigen::MatrixXd deduplify_both(Eigen::MatrixXd x, 
                               Rcpp::IntegerVector idx_vech,
                               Rcpp::IntegerVector idx_tvech,
                               Rcpp::IntegerVector idx_vech_match);

// method for deduplifying left side
Eigen::MatrixXd deduplify_left(Eigen::MatrixXd x, 
                               Rcpp::IntegerVector idx_vech,
                               Rcpp::IntegerVector idx_tvech,
                               Rcpp::IntegerVector idx_vech_match);

// method for deduplifying right side
Eigen::MatrixXd deduplify_right(Eigen::MatrixXd x, 
                                Rcpp::IntegerVector idx_vech,
                                Rcpp::IntegerVector idx_tvech,
                                Rcpp::IntegerVector idx_vech_match);

// method for creating commutation matrix
Eigen::MatrixXd create_commutation(int n);

// create duplication matrix
Eigen::MatrixXd create_duplication(int n);

// method for which function
Rcpp::IntegerVector which(Rcpp::LogicalVector x);

// method for sign function
int sign(double x);
