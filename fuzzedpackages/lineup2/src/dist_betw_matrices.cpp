#include "dist_betw_matrices.h"
#include <Rcpp.h>

using namespace Rcpp;

// calculate the distance between columns of matrix x and columns of matrix y
// d(i,j) = sqrt{ sum_k [x(k,i) - y(k,j)]^2 }
//
// [[Rcpp::export()]]
NumericMatrix rmsd_betw_matrices(const NumericMatrix& x,
                                 const NumericMatrix& y)
{
    int n = x.cols();
    int m = y.cols();
    int p = x.rows();
    if(y.rows() != p)
        throw std::invalid_argument("nrow(x) != nrow(y)");

    NumericMatrix result(n,m);

    for(int i=0; i<n; i++) {
        for(int j=0; j<m; j++) {
            checkUserInterrupt();  // check for ^C from user

            double value = 0.0;
            int count = 0;
            for(int k=0; k<p; k++) {
                if(R_finite(x(k,i)) && R_finite(y(k,j))) {
                    value += (x(k,i) - y(k,j))*(x(k,i) - y(k,j));
                    count++;
                }
            }
            if(count > 0) result(i,j) = sqrt(value/(double)count);
            else result(i,j) = NA_REAL;
        }
    }

    return result;
}

// like rmsd_betw_matrices but using the average absolute difference
//
// [[Rcpp::export()]]
NumericMatrix mad_betw_matrices(const NumericMatrix& x,
                                const NumericMatrix& y)
{
    int n = x.cols();
    int m = y.cols();
    int p = x.rows();
    if(y.rows() != p)
        throw std::invalid_argument("nrow(x) != nrow(y)");

    NumericMatrix result(n,m);

    for(int i=0; i<n; i++) {
        for(int j=0; j<m; j++) {
            checkUserInterrupt();  // check for ^C from user

            double value = 0.0;
            int count=0;
            for(int k=0; k<p; k++) {
                if(R_finite(x(k,i)) && R_finite(y(k,j))) {
                    value += fabs(x(k,i) - y(k,j));
                    count++;
                }
            }
            if(count > 0) result(i,j) = value / (double)count;
            else result(i,j) = NA_REAL;
        }
    }

    return result;
}
