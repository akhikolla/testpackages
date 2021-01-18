#ifndef ITPPEXT_H_
#define ITPPEXT_H_

// Use BLAS and LAPACK
#define ARMA_USE_BLAS
#define ARMA_USE_LAPACK

/**
 * Various generic functions required by PSGP
 */
#include <string>
#include <RcppArmadillo.h>

#define COLUMN_ORDER 0
#define ROW_ORDER 1

typedef arma::mat mat;
typedef arma::vec vec;
typedef arma::uvec uvec;

namespace psgp_arma
{

vec ltr_vec(mat M);     // Vector of lower triangular elements
vec utr_vec(mat M);     // Vector of upper triangular elements
mat ltr_mat(vec v);     // Lower triangular matrix
mat utr_mat(vec v);     // Upper triangular matrix

double cond(mat M, int p=2); // Condition number for matrix p-norm (1 or 2)

uvec randperm(int n);  // Random permutation of numbers between 0 and N-1
uvec sequence(int from, int to);

vec min(vec u, vec v); // Minimum elements from 2 vectors of equal length

mat concat_cols(mat X, vec y); // Concatenate matrix and vector
mat concat_cols(mat X, mat Y); // Concatenate matrix and matrix

vec mean_rows(mat X);  // vector of column means
vec mean_cols(mat X);  // vector of row means
mat cov(mat X, vec &xmean); // covariance of rows, also returns the mean
mat cov(mat X); // covariance of rows

void normalise(mat &X);
void normalise(mat &X, vec &mean, vec &covdiag);
void denormalise(mat &X, vec mean, vec covdiag);

mat zeros(int m, int n);
mat ones(int m, int n);
vec zeros(int m);
vec ones(int m);

double norm();
double sign(double x);

} // END OF namespace armaEXT


#endif /*ITPPEXT_H_*/
