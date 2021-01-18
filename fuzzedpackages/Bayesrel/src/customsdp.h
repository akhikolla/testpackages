

#ifndef customsdp_h
#define customsdp_h
#include <RcppArmadillo.h>


arma::ivec int_vector_csdp2RArma(int n, int *y);

arma::dvec double_vector_csdp2RArma(int n, double *y);

int * int_vector_R2csdpArma(int n, const arma::ivec& y);

double * double_vector_R2csdpArma(int n, const arma::dvec& y);

struct blockmatrix blkmatrix_R2csdpArma(const Rcpp::List& X);

Rcpp::List blkmatrix_csdp2RArma(const blockmatrix& X);

struct constraintmatrix *constraints_R2csdpArma(const Rcpp::List& A);

void initArma(int n,
              int k,
              struct blockmatrix C,
              double *a,
              struct constraintmatrix *constraints,
              struct blockmatrix *pX0,
              double **py0,
              struct blockmatrix *pZ0);

int custom_sdpCpp(int n,
    int k,
    blockmatrix& C,
    double *a,
    struct constraintmatrix *constraints,
    double constant_offset,
    double *ppobj,
    double *pdobj,
    const arma::cube& car, arma::dvec& out,
    const int printlevel = 0);


#endif /* customsdp_h */
