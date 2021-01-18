#include "StepGaussInhibit.h"

/***************
* class StepGauss
* implements virtual class Step for Gaussian data, i.e. (weighted) l2-loss
* Thomas Hotz, 2007-2010
***************/

/*************
* constructor for n data points, with pointers to array of that length
****************/
StepGaussInhibit::StepGaussInhibit(unsigned int n, double* xcs, double* xcss, double* xcsv, int inhibitStart, int inhibitMiddle, int inhibitEnd) : StepGauss(n, xcs, xcss, xcsv), istart(inhibitStart), imiddle(inhibitMiddle), iend(inhibitEnd) {}

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
double StepGaussInhibit::cost(unsigned int startIndex, unsigned int endIndex) const {
  if(startIndex == 0) { if(csv[endIndex] < istart) return R_PosInf; else return StepGauss::cost(startIndex, endIndex); }
  if(endIndex == N - 1) { if(csv[N - 1] - csv[startIndex - 1] < iend) return R_PosInf; else return StepGauss::cost(startIndex, endIndex); }
  if(csv[endIndex] - csv[startIndex - 1] < imiddle) return R_PosInf; else return StepGauss::cost(startIndex, endIndex);
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
SEXP forwardGaussInhibit(SEXP cumSum, SEXP cumSumSq, SEXP cumSumVar, SEXP maxBlocks, SEXP istart, SEXP imiddle, SEXP iend) {
  // initialise object
  StepGaussInhibit data = StepGaussInhibit(length(cumSum), REAL(cumSum), REAL(cumSumSq), REAL(cumSumVar), asInteger(istart), asInteger(imiddle), asInteger(iend));
  
  // check lengths
  if(data.N < 1) error("cumSum must have at least one element");
  if(length(cumSumSq) != (int) data.N) error("cumSumSq must have same length as cumSum");
  if(length(cumSumVar) != (int) data.N) error("cumSumVar must have same length as cumSum");
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
SEXP pathGaussInhibit(SEXP cumSum, SEXP cumSumSq, SEXP cumSumVar, SEXP maxBlocks, SEXP istart, SEXP imiddle, SEXP iend) {
  // initialise object
  StepGaussInhibit data = StepGaussInhibit(length(cumSum), REAL(cumSum), REAL(cumSumSq), REAL(cumSumVar), asInteger(istart), asInteger(imiddle), asInteger(iend));
  
  // check lengths
  if(data.N <= 1) error("there must be more than one block");
  if(length(cumSumSq) != (int) data.N) error("length of cumSumSq must match cumSum's");
  if(length(cumSumVar) != (int) data.N) error("cumSumVar of rightEnd must match cumSum's");
  if(length(maxBlocks) != 1) error("maxBlocks must be a single integer");
  
  // run algorithm
  return data.path(asInteger(maxBlocks)); // the solution path, i.e. p[i, k] is the (i+1)th jump in the solution having k+1 jumps
}

} // end C wrapper
