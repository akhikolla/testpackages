#include <Rcpp.h>
#include <cstdint>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
double KS_Cpp(const NumericVector& X) {
    int n = X.size();

    NumericVector mX = -X;
    NumericVector w(2*n);

    // concatenate X and -X
    std::copy(X.begin(), X.end(), w.begin());
    std::copy(mX.begin(), mX.end(), w.begin() + n);

    // Get the equivalent of order(w)
    NumericVector sorted = clone(w).sort();
    IntegerVector order = match(sorted, w);

    double TS_sum = 0;
    double TS_value = 0;
    int i;
    for(i = 0; i < n; i++) {
        TS_sum += (order[i] <= n) ? 1.0/n : -1.0/n;
        if (std::abs(TS_sum) > TS_value) {
            TS_value = std::abs(TS_sum);
        }
    }

    return TS_value;
}

