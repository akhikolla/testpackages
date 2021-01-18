#ifndef DDC_H
#define DDC_H

#ifndef ARMA_DONT_PRINT_ERRORS
#define ARMA_DONT_PRINT_ERRORS
#endif

#ifndef  ARMA_USE_CXX11
#define ARMA_USE_CXX11
#endif

#include "LocScaleEstimators.h"

namespace DDC
{

struct fastRobCorout {
arma::umat ngbrs;
arma::mat robcorrs;
};

struct kbestcorr {
arma::uvec selected;
arma::vec corrs;
};

  
arma::uvec vinter(const arma::uvec &first, const arma::uvec &second);

arma::uvec vdiff(const arma::uvec &first, const arma::uvec &second);


arma::uvec col2cell(const arma::uvec &colNrs, const int n);

arma::uvec row2cell(const arma::uvec &rowNrs, const int n, const int d);

double weightedMedian(arma::vec x, arma::vec weights);


arma::vec predictCol(const arma::vec &colj, const arma::mat &U, const int coln,
                     const arma::umat &ngbrs, const arma::mat &corrweight,
                     const arma::mat &robslopes, const int combinRule);

double slopeMedWLS(const arma::vec &xcol, const arma::vec &colj,
                   double qRegr, double precScale);

arma::vec compSlopes(const arma::vec &colj, arma::uvec ngbrs, const arma::mat &U,
                     double qRegr, double precScale);

double corrGKWLS(arma::vec xcol, double qCorr, arma::vec colj, double precScale);

arma::vec limitFilt(arma::vec v, double qCut);

arma::vec rawEquiGYfilt(const arma::vec &v, double qCut);

arma::vec equiGYfilt(const arma::vec &v, double qCut, const int miter);

kbestcorr kBestCorr(const arma::vec &colj, const arma::mat &U, const int coln,
                    const unsigned int k, double qCorr, double precScale);

double deShrink(const arma::vec &colj, const arma::mat &Z, const int coln,
                   double qRegr, double precScale);


// Code for Fastrobcor

void get_NN_2Set(double *data, double *query, int *D, int *ND, int *NQ, int *K, double *EPS,
                 int *SEARCHTYPE, int *USEBDTREE, double *SQRAD, int *nn_index, double *distances);


arma::vec transClassic( arma::vec y, const double precScale = 1e-12);




fastRobCorout FastRobCorActual(const arma::mat &X, const arma::vec &locX, const arma::vec &scaleX,
                               unsigned int k, int qdim, const unsigned int nCorr,
                               int absCorr = 1,
                               int transFun = 0,  double precScale = 1e-12, int treetype = 0,
                               int searchtype = 1,  double radius = 0,
                               double eps = 0, int includeSelf = 0);

}

#endif
