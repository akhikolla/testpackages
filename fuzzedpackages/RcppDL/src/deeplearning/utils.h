#ifndef UTILS_H_
#define UTILS_H_

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <algorithm>

#include <Rcpp.h>

inline double uniform(double min, double max) {
    Rcpp::RNGScope scope;
    return R::runif(min, max);
}

inline int binomial(int n, double p) {
    Rcpp::RNGScope scope;
    if (p < 0 || p > 1)
        return 0;

    int c = 0;
    double r;

    for (int i = 0; i < n; i++) {
        r = R::runif(0, RAND_MAX) / (RAND_MAX + 1.0);
        if (r < p)
            c++;
    }

    return c;
}

inline double sigmoid(double x) {
    return 1.0 / (1.0 + exp(-x));
}

#endif
