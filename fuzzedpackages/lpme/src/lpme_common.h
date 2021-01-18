#ifndef _lpme_common_h
#define _lpme_common_h

#include <RcppArmadillo.h>

#define ESMALL 1e-10  /* small number */
#define ELARGE 1e+10 /* large number */
#define sqrt2pi 2.5066282746310002416 /* sqrt(2*pi)*/
typedef Rcpp::NumericMatrix::iterator mat_iterator;
using namespace Rcpp;

// Fourier transform for error U
Rcpp::NumericVector FuLap(Rcpp::NumericVector t, double sigU);
Rcpp::NumericVector FuGau(Rcpp::NumericVector t, double sigU);
Rcpp::NumericVector FuLapinv(Rcpp::NumericVector t, double sigU);
Rcpp::NumericVector FuGauinv(Rcpp::NumericVector t, double sigU);
// Fourier transform for Kernel K
Rcpp::NumericVector FK(Rcpp::NumericVector t);

// first derivative for Fourier transform of Kernel K
Rcpp::NumericVector FK1(Rcpp::NumericVector t);

// second derivative for Fourier transform of Kernel K
Rcpp::NumericVector FK2(Rcpp::NumericVector t);

// second-order Kernel K
double K_sec_order(double x);
// Fourier transform for Kernel K
Rcpp::NumericVector FK_sec_order(Rcpp::NumericVector t);
// first derivative for Fourier transform of Kernel K
Rcpp::NumericVector FK1_sec_order(Rcpp::NumericVector t);
// second derivative for Fourier transform of Kernel K
Rcpp::NumericVector FK2_sec_order(Rcpp::NumericVector t);

// function to generate subvectors a[-(ind1:ind2)] and b[-(ind1:ind2)] and save them to w and y respectively
// void allows more than two values can be returned. For example, w and y are returned here.
void subvecij(const Rcpp::NumericVector& a, const Rcpp::NumericVector& b, int ind1, int ind2, Rcpp::NumericVector& w, Rcpp::NumericVector& y);

// function to estimate ghat of x using JASA when U is Laplace
void gjasaLap(Rcpp::NumericVector& res, const Rcpp::NumericVector& x, const Rcpp::NumericVector& t, double dt, const Rcpp::NumericVector& W, 
      const Rcpp::NumericVector& Y, double sigU, double h);

// function to estimate ghat of x using JASA when U is Gaussian
void gjasaGau(Rcpp::NumericVector& res, const Rcpp::NumericVector& x, const Rcpp::NumericVector& t, double dt, const Rcpp::NumericVector& W, 
      const Rcpp::NumericVector& Y, double sigU, double h);

// function to estimate ghat of x using OURS
void gnewLap(Rcpp::NumericVector& res, const Rcpp::NumericVector& x, const Rcpp::NumericVector& input, const Rcpp::NumericVector& output, 
        double beta, double beta2, const Rcpp::NumericVector& mconst, const Rcpp::NumericVector& Kinput, 
        const Rcpp::NumericVector& W, const Rcpp::NumericVector& Y, double sigU, double h);

// function to estimate ghat of x using OURS when error is Gaussian
void gnewGau(Rcpp::NumericVector& ghatofx, const Rcpp::NumericVector& x, const Rcpp::NumericVector& input, const Rcpp::NumericVector& output, 
        double beta, double beta2, const Rcpp::NumericVector& mconst, const Rcpp::NumericVector& Kinput, 
        const Rcpp::NumericVector& W, const Rcpp::NumericVector& Y, double sigU, double h);

#endif
