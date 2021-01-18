#ifndef CCHI2
#define CCHI2

#include "gamma.h"
#include <math.h>

/*
Computed the survival function of the chi2 distribution function for x and k
*/
inline double Chi2_sf(double x, double k) {
    return complementedIncompleteGamma(0.5*x,0.5*k);
}

/*
Computed the cumulative distribution function of the chi2 distribution function for x and k
*/
inline double Chi2_cdf(double x, double k) {
    if (k==2.0) {
        return 1.0 - exp(-0.5*x);
    } else {
        return regularizedLowerIncompleteGamma(0.5*x,0.5*k);
    }
}

/*
Probability density function of Chi2
*/
inline double Chi2_pdf(double x, double k) {
    if (x<0.0) return 0.0;
    return pow(x,0.5*k-1.0)*exp(-0.5*x)/(pow(2.0,0.5*k)*tgamma(0.5*k));
}

#endif
