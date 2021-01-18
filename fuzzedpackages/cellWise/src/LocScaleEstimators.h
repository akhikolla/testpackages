#ifndef LocScaleEstimators_H
#define LocScaleEstimators_H

#ifndef ARMA_DONT_PRINT_ERRORS
#define ARMA_DONT_PRINT_ERRORS
#endif

#ifndef  ARMA_USE_CXX11
#define ARMA_USE_CXX11
#endif

#include "Rmath.h"
#include "RcppArmadillo.h"

namespace LocScaleEstimators // namespace for all univariate location/scale estimators
{
struct locscale {
  double loc;
  double scale;
  double rawloc;
  double rawscale;
  double cfac1;
  double cfac2;
  arma::uvec weights;
};

struct Xlocscale {
  arma::vec loc;
  arma::vec scale;
};



arma::uvec sample(const arma::uvec & sampleFrom,
                  const unsigned int nbSamples,
                  const bool replace); 

//////////////
// Location //
//////////////

void locWeightBiweight(arma::vec &x);
void locWeightHuber15(arma::vec &x);
void locWeightTanh154(arma::vec &x);
double loc1StepM(const arma::vec &x, std::function<void (arma::vec&)> weightFunction, 
                 double initLoc = arma::datum::nan,
                 double initScale = arma::datum::nan, double precScale = 1e-12);

///////////
// SCALE //
///////////

void psiTanh (arma::vec &x, double b = 1.5, double c = 4, double k = 4.1517212, double A = 0.7532528,
                   double B = 0.8430849);
void rhoTanh154(arma::vec &x);
void rhoHuber25(arma::vec &x);
void rhoHuber15(arma::vec &x);
double scale1StepM(const arma::vec &x, std::function<void (arma::vec&)> rhoFunction,
                   double initScale = arma::datum::nan, double precScale = 1e-12);

//////////////
// LOCSCALE //
//////////////

locscale uniMcd(arma::vec y, double alpha = 0.5);

/////////////////////////
// Vectorized LOCSCALE //
/////////////////////////

Xlocscale estLocScale(const arma::mat &X,
                      unsigned int nLocScale, int type, double precScale,
                      const int center, double alpha = 0.5);

///////////
// Ranks //
///////////

arma::vec rank(arma::vec& v);

}

#endif
