#include "StepGaussCut.h"

/***************
* class StepGauss
* implements virtual class Step for Gaussian data, i.e. (weighted) l2-loss
* Thomas Hotz, 2007-2010
***************/

/*************
* constructor for n data points, with pointers to array of that length
****************/
StepGaussCut::StepGaussCut(unsigned int n, double* xbcs, double* xbcss, double* xbcsv, double* xacs, double* xacss, double* xacsv, int cutBefore, int cutAfter) : Step(n), cbefore(cutBefore), cafter(cutAfter), bcs(xbcs), bcss(xbcss), bcsv(xbcsv), acs(xacs), acss(xacss), acsv(xacsv) {}

/*************
* cost
* calculate cost of a block
*
* in:
* startIndex : the first index in the block, 0 <= startIndex <= endIndex
* endIndex : the last index in the block, startIndex <= endIndex <= N - 1
* 
* out:
* the cost functional of that block
* infinite if jumps get too close, else the inherited cost functional
****************/
double StepGaussCut::cost(unsigned int startIndex, unsigned int endIndex) const {
  if(startIndex == 0) { 
    if(bcsv[endIndex] == NA_REAL)
      return R_PosInf;
    else
      return bcss[endIndex] - bcs[endIndex] * bcs[endIndex] / bcsv[endIndex];
  }
  if(endIndex == N - 1) {
    if(acsv[startIndex - 1] == NA_REAL)
      return R_PosInf; 
    else 
      return bcss[endIndex] - acss[startIndex - 1] - ( bcs[endIndex] - acs[startIndex - 1] ) * ( bcs[endIndex] - acs[startIndex - 1] ) / ( bcsv[endIndex] - acsv[startIndex - 1] ); 
  }
  if(bcsv[endIndex] - acsv[startIndex - 1] < 1) 
    return R_PosInf; 
  else 
    return bcss[endIndex] - acss[startIndex - 1] - ( bcs[endIndex] - acs[startIndex - 1] ) * ( bcs[endIndex] - acs[startIndex - 1] ) / ( bcsv[endIndex] - acsv[startIndex - 1] );
}

// C wrapper
extern "C" {

/*************
* forwardGauss
* function to be called from R
* computes forward selection for block signals, wrapper for StepGauss::forward
*
* in:
* cumSum : vector of cumulative sums
* cumSumSq : vector of cumulative sums of squares
* cumSumVar : vector of cumulative sums of variances, i.e. the expected css
* maxBlocks : a single integer giving the maximal number of blocks
*
* out:
* data.frame comprising rightIndex, number, depth and improve of the candidates
****************/
SEXP forwardGaussCut(SEXP bcumSum, SEXP bcumSumSq, SEXP bcumSumVar, SEXP acumSum, SEXP acumSumSq, SEXP acumSumVar, SEXP maxBlocks, SEXP cbefore, SEXP cafter) {
  // initialise object
  StepGaussCut data = StepGaussCut(length(bcumSum), REAL(bcumSum), REAL(bcumSumSq), REAL(bcumSumVar), REAL(acumSum), REAL(acumSumSq), REAL(acumSumVar), asInteger(cbefore), asInteger(cafter));
  
  // check lengths
  if(data.N < 1) error("bcumSum must have at least one element");
  if(length(bcumSumSq) != (int) data.N) error("bcumSumSq must have same length as bcumSum");
  if(length(bcumSumVar) != (int) data.N) error("bcumSumVar must have same length as bcumSum");
  if(length(acumSum) != (int) data.N) error("acumSum must have same length as bcumSum");
  if(length(acumSumSq) != (int) data.N) error("acumSumSq must have same length as bcumSum");
  if(length(acumSumVar) != (int) data.N) error("acumSumVar must have same length as bcumSum");
  if(length(maxBlocks) != 1) error("maxBlocks must be a single integer");
  
  // run algorithm
  return data.forward(asInteger(maxBlocks));
}

/*************
* pathGauss
* function to be called from R
* computes path of Potts minimisers, wrapper for StepGauss::path
*
* in:
* cumSum : a numeric vector, the cumulative sums of the signal
* cumSumSq : a numeric vector, the cumulative sum of squares of the signal
* cumSumVar : an numeric vector, the cumulative sums of variances, i.e. the expected css
* maxBlocks : a single integer giving the maximal number of blocks
*
* out:
* a list conprimising the right indices of the solution path, and the associated costs
****************/
SEXP pathGaussCut(SEXP bcumSum, SEXP bcumSumSq, SEXP bcumSumVar, SEXP acumSum, SEXP acumSumSq, SEXP acumSumVar, SEXP maxBlocks, SEXP cbefore, SEXP cafter) {
  // initialise object
  StepGaussCut data = StepGaussCut(length(bcumSum), REAL(bcumSum), REAL(bcumSumSq), REAL(bcumSumVar), REAL(acumSum), REAL(acumSumSq), REAL(acumSumVar), asInteger(cbefore), asInteger(cafter));
  
  // check lengths
  if(data.N < 1) error("cumSum must have at least one element");
  if(length(bcumSumSq) != (int) data.N) error("bcumSumSq must have same length as bcumSum");
  if(length(bcumSumVar) != (int) data.N) error("bcumSumVar must have same length as bcumSum");
  if(length(acumSum) != (int) data.N) error("acumSum must have same length as bcumSum");
  if(length(acumSumSq) != (int) data.N) error("acumSumSq must have same length as bcumSum");
  if(length(acumSumVar) != (int) data.N) error("acumSumVar must have same length as bcumSum");
  if(length(maxBlocks) != 1) error("maxBlocks must be a single integer");
  
  // run algorithm
  return data.path(asInteger(maxBlocks)); // the solution path, i.e. p[i, k] is the (i+1)th jump in the solution having k+1 jumps
}

} // end C wrapper
