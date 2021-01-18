#include <Rcpp.h>
#include <cstdint>
#include <cmath>
#include <unordered_set>
#include <helpers.h>
using namespace Rcpp;


// [[Rcpp::export]]
double K2_Cpp(const NumericVector& X) {
    int n = X.size();
    int i, j, count = 0;

    int nC2 = Rf_choose(n, 2);
    NumericVector aX = abs(X);

    NumericVector minus_points(nC2);
    NumericVector plus_points(nC2);

    for(i = 0; i < n; i++) {
        for(j = 0; j < i; j++) {
            minus_points[count] = std::abs(X[i] - X[j]);
            plus_points[count] = std::abs(X[i] + X[j]);
            count++;
        }
    }

    // We will sort the potential points to speed up comparison
    aX.sort();
    minus_points.sort();
    plus_points.sort();

    // Iterate over all potential points for maximum and see the
    // sum value and comapre to the max

    double TS_max = 0;

    double TS_sum;

    NumericVector aX2 = 2 * aX;
    for (int x = 0; x < n; x++) {
        TS_sum = std::abs(
          n - count_smaller(aX2, aX2[x]) +
              2 * (count_smaller(minus_points, aX2[x]) -
                   count_smaller(plus_points, aX2[x]))
        );

        if (TS_sum > TS_max) {
            TS_max = TS_sum;
        }
    }


    for (int m = 0; m < nC2; m++) {
        TS_sum = std::abs(
            n - count_smaller(aX2, minus_points[m]) +
                2 * (count_smaller(minus_points, minus_points[m]) -
                count_smaller(plus_points, minus_points[m]))
        );

        if (TS_sum > TS_max) {
            TS_max = TS_sum;
        }
    }


    for (int p = 0; p < nC2; p++) {
        TS_sum = std::abs(
            n - count_smaller(aX2, plus_points[p]) +
                2 * (count_smaller(minus_points, plus_points[p]) -
                count_smaller(plus_points, plus_points[p]))
        );

        if (TS_sum > TS_max) {
            TS_max = TS_sum;
        }

    }

    return TS_max / n; // / n * n;
}

