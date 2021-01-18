#include <Rcpp.h>
#include <cstdint>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
double WCX_Cpp(const NumericVector& X) {
    int n = X.size();
    double TS_sum = 0;
    int i,j;
    for(j = 0; j < n; j++) {
        for(i = 0; i < j; i++) {
            TS_sum += (X[i] + X[j] > 0);
        }
    }

    double TS_value = (double)TS_sum / Rf_choose(n, 2) - 0.5;

    return sqrt((double)n) * TS_value;
}

