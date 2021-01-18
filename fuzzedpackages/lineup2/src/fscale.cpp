// standardize a vector (subtract mean and divide by SD)

#include "fscale.h"
#include <Rcpp.h>
using namespace Rcpp;

// fscale: standardize the vector x using only the values where
//         both x and y are not missing
// [[Rcpp::export()]]
NumericVector fscale(const NumericVector& x,
                     const NumericVector& y)
{
    const int n = x.size();
    if(y.size() != n)
        throw std::invalid_argument("length(x) != length(y)");

    // result matrix, filled with missing values
    NumericVector result(n);

    double sum=0.0, sumsq=0.0, diff, first=NA_REAL;
    int count=0;

    // get mean and sd
    for(int i=0; i<n; i++) {
        if(R_finite(x[i]) && R_finite(y[i])) {
            count++;
            if(!R_finite(first)) first = x[i]; // first non-missing value
            else {
                // sum(x) and sum(x*x) with x centered at first non-missing value
                sum += (diff=(x[i]-first));
                sumsq += (diff*diff);
            }
            result[i] = x[i];
        }
        else result[i] = NA_REAL; // non-finite values -> NA
    }

    // center and scale the column
    if(count > 1) { /* if count < 2, do nothing */
        sumsq = sqrt((sumsq - (sum*sum)/(double)count)/(double)(count-1));
        sum /= (double)count;
        for(int i=0; i<n; i++) {
            if(R_finite(x[i]) && R_finite(y[i])) {
                result[i] = (x[i] - sum - first)/(sumsq);
            }
        }
    }

    return result;
}
