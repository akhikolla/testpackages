
#include <Rcpp.h>

#include <algorithm>
#include "choleskyDecomposition.h"
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>

#include "Filter.h"
#include "FilterBessel.h"

using namespace Rcpp;

// [[Rcpp::export(name = ".deconvolveJump")]]
double deconvolveJump(const NumericVector &grid, const NumericVector &observations, const NumericVector &time,
                      const double &leftValue, const double &rightValue,
                      const int &typeFilter, const List &inputFilter, const NumericVector &covariances) {
  Filter* filter;
  
  switch (typeFilter) {
  case 0: 
    filter = new FilterBessel(inputFilter);
    break;
  default:
    stop("unknown filter type");
  }
  
  unsigned int const n = observations.size();
  const int size = observations.size();
  
  double* A = choleskyDecomposition(size, covariances);
  
  const char uplo = 'U';
  const char trans = 'T';
  const char diag = 'N';
  const int incx = 1;
  const int bands = std::min<int>(covariances.size() - 1, size - 1);
  const int ldA = bands + 1;
  
  double ad, cp, optCp = grid[0u] - 1, costs, optCosts = R_PosInf;
  double* centeredObservations = new double[n];
  
  for(unsigned int cpIndex = 0u; cpIndex < static_cast<unsigned int>(grid.size()); ++cpIndex) {
    checkUserInterrupt();
    
    cp = grid[cpIndex];
    
    for (unsigned int i = 0u; i < n; ++i) {
      ad = filter -> antiderivative(time[i] - cp);
      centeredObservations[i] = observations[i] - (leftValue * (1.0 - ad) + rightValue * ad);
    }
    
    F77_CALL(dtbsv)(&uplo, &trans, &diag, &size, &bands, A, &ldA, centeredObservations, &incx);
    
    costs = 0.0;
    for (unsigned int i = 0u; i < n; ++i) {
      costs += centeredObservations[i] * centeredObservations[i];
    }
    
    if (costs < optCosts) {
      optCp = cp;
      optCosts = costs;
    }
  }
  
  delete[] centeredObservations;
  delete[] A;
  delete filter;
  
  return optCp;
}

// [[Rcpp::export(name = ".deconvolvePeak")]]
List deconvolvePeak(const NumericVector &gridLeft, const NumericVector &gridRight, 
                    const NumericVector &observations, const NumericVector &time,
                    const double &leftValue, const double &rightValue,
                    const int &typeFilter, const List &inputFilter, const NumericVector &covariances,
                    const double &tolerance) {
  Filter* filter;
  
  switch (typeFilter) {
  case 0: 
    filter = new FilterBessel(inputFilter);
    break;
  default:
    stop("unknown filter type");
  }
  
  unsigned int const n = observations.size();
  const int size = observations.size();
  
  double* A = choleskyDecomposition(size, covariances);
  
  const char uplo = 'U';
  const char trans = 'T';
  const char diag = 'N';
  const int incx = 1;
  const int bands = std::min<int>(covariances.size() - 1, size - 1);
  const int ldA = bands + 1;
  
  double adLeft, adRight, left, right, value, optLeft = gridLeft[0u] - 1, optRight = gridRight[0u] - 1,
    optValue = 0, costs, optCosts = R_PosInf;
  double* centeredObservations = new double[n];
  double* weights = new double[n];
  
  for(unsigned int leftIndex = 0u; leftIndex < static_cast<unsigned int>(gridLeft.size()); ++leftIndex) {
  for(unsigned int rightIndex = 0u; rightIndex < static_cast<unsigned int>(gridRight.size()); ++rightIndex) {
      checkUserInterrupt();
      
      left = gridLeft[leftIndex];
      right = gridRight[rightIndex];
      
      if (left < right - tolerance) {
        for (unsigned int i = 0u; i < n; ++i) {
          adLeft = filter -> antiderivative(time[i] - left);
          adRight = filter -> antiderivative(time[i] - right);
          
          centeredObservations[i] = observations[i] - (leftValue * (1.0 - adLeft) + rightValue * adRight);
          weights[i] = adLeft - adRight;
        }
        
        F77_CALL(dtbsv)(&uplo, &trans, &diag, &size, &bands, A, &ldA, centeredObservations, &incx);
        F77_CALL(dtbsv)(&uplo, &trans, &diag, &size, &bands, A, &ldA, weights, &incx);
        
        double numerator = 0.0, denominator = 0.0, sumCenteredObservationsSquared = 0.0;
        for (unsigned int i = 0u; i < n; ++i) {
          sumCenteredObservationsSquared += centeredObservations[i] * centeredObservations[i];
          numerator += weights[i] * centeredObservations[i];
          denominator += weights[i] * weights[i];
        }
        value = numerator / denominator;
        
        costs = sumCenteredObservationsSquared - 2.0 * value * numerator + value * value * denominator;
        
        if (costs < optCosts) {
          optLeft = left;
          optRight = right;
          optValue = value;
          optCosts = costs;
        }
      }
    }
  }
  
  delete[] centeredObservations;
  delete[] weights;
  delete[] A;
  delete filter;
  
  return List::create(_["left"] = optLeft, _["right"] = optRight, _["value"] = optValue);
}
