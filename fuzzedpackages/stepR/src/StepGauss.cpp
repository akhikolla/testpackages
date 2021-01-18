#include "StepGauss.h"

/***************
* class StepGauss
* implements virtual class Step for Gaussian data, i.e. (weighted) l2-loss
* Thomas Hotz, 2007-2011
***************/

/*************
* constructor for n data points, with pointers to array of that length
****************/
StepGauss::StepGauss(unsigned int n, double* xcs, double* xcss, double* xcsv) : Step(n), cs(xcs), css(xcss), csv(xcsv) {}

/*************
* constructor for n data points and bounds
****************/
StepGauss::StepGauss(unsigned int n, double* xcs, double* xcss, double* xcsv, double* xlb, double* xub) : Step(n, xlb, xub), cs(xcs), css(xcss), csv(xcsv) {}

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
* for Gaussian data: the residual sum of squares after fitting with a mean, i.e. css - cs^2 / csv on this block
****************/
double StepGauss::cost(unsigned int startIndex, unsigned int endIndex) const {
  if(startIndex == endIndex) return 0; else
  if(startIndex == 0) {
    return css[endIndex] - cs[endIndex] * cs[endIndex] / csv[endIndex]; // all data
  } else {
    return css[endIndex] - css[startIndex - 1] - ( cs[endIndex] - cs[startIndex - 1] ) * ( cs[endIndex] - cs[startIndex - 1] ) / ( csv[endIndex] - csv[startIndex - 1] );
  }
}

/*************
* costBound
* calculate cost of a block given bounds
*
* in:
* startIndex : the first index in the block, 0 <= startIndex <= endIndex
* endIndex : the last index in the block, startIndex <= endIndex <= N - 1
* bound : the lower and upper bound for the mean in this interval
* 
* out:
* the cost functional of that block
* for Gaussian data: the residual sum of squares after fitting with a bounded mean on this block
****************/
double StepGauss::costBound(unsigned int startIndex, unsigned int endIndex, const LUBound& bound) const {
  if(bound.lower > bound.upper) return R_PosInf; else
  if(startIndex == 0) {
    double mean = cs[endIndex] / csv[endIndex];
    double bias = fmax2( fmin2( mean, bound.upper ), bound.lower ) - mean;
    return css[endIndex] + ( bias * bias - mean * mean ) * csv[endIndex] ; // all data
  } else {
    double weight = csv[endIndex] - csv[startIndex - 1];
    double mean = ( cs[endIndex] - cs[startIndex - 1] ) / weight;
    double bias = fmax2( fmin2( mean, bound.upper ), bound.lower ) - mean;
    return css[endIndex] - css[startIndex - 1] + ( bias * bias - mean * mean ) * weight;
  }
}

/*************
* estBound
* calculate estimate on a block given bounds
*
* in:
* startIndex : the first index in the block, 0 <= startIndex <= endIndex
* endIndex : the last index in the block, startIndex <= endIndex <= N - 1
* bound : the lower and upper bound for the mean in this interval
* 
* out:
* the estimate for that block
* for Gaussian data: the mean on this block if it isn't outside a boundary, then the latter
****************/
double StepGauss::estBound(unsigned int startIndex, unsigned int endIndex, const LUBound& bound) const {
  if(bound.lower > bound.upper) return R_NaN; else
  if(startIndex == 0) {
    double mean = cs[endIndex] / csv[endIndex];
    #ifdef DEBUGbounded
    Rprintf("  si = %d, ei = %d, mean = %4.2e, lower = %4.2e, upper = %4.2e\n", startIndex, endIndex, mean, bound.lower, bound.upper);
    #endif
    return fmax2( fmin2( mean, bound.upper ), bound.lower );
  } else {
    double weight = csv[endIndex] - csv[startIndex - 1];
    double mean = ( cs[endIndex] - cs[startIndex - 1] ) / weight;
    #ifdef DEBUGbounded
    Rprintf("  si = %d, ei = %d, mean = %4.2e, lower = %4.2e, upper = %4.2e\n", startIndex, endIndex, mean, bound.lower, bound.upper);
    #endif
    return fmax2( fmin2( mean, bound.upper ), bound.lower );
  }
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
SEXP forwardGauss(SEXP cumSum, SEXP cumSumSq, SEXP cumSumVar, SEXP maxBlocks) {
  // initialise object
  StepGauss data = StepGauss(length(cumSum), REAL(cumSum), REAL(cumSumSq), REAL(cumSumVar));
  
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
SEXP pathGauss(SEXP cumSum, SEXP cumSumSq, SEXP cumSumVar, SEXP maxBlocks) {
  // initialise object
  StepGauss data = StepGauss(length(cumSum), REAL(cumSum), REAL(cumSumSq), REAL(cumSumVar));
  
  // check lengths
  if(data.N <= 1) error("there must be more than one block");
  if(length(cumSumSq) != (int) data.N) error("length of cumSumSq must match cumSum's");
  if(length(cumSumVar) != (int) data.N) error("length of cumSumVar must match cumSum's");
  if(length(maxBlocks) != 1) error("maxBlocks must be a single integer");
  
  // run algorithm
  return data.path(asInteger(maxBlocks)); // the solution path, i.e. p[i, k] is the (i+1)th jump in the solution having k+1 jumps
}

/*************
* boundedGauss
* function to be called from R
* computes bounded solution, wrapper for StepGauss::bounded
*
* in:
* cumSum : a numeric vector, the cumulative sums of the signal
* cumSumSq : a numeric vector, the cumulative sum of squares of the signal
* cumSumVar : an numeric vector, the cumulative sums of variances, i.e. the expected css
* start : for every possible left index where intervals with this left index start in the list of intervals (increasing in left indices), NA if none
* rightIndex : right indices of the intervals in the list, increasing for each left index
* lower : the lower bounds for the estimator at the respective interval
* upper : the upper bounds for the estimator at the respective interval
*
* out:
* a list conprimising the right indices and values of the solution, and the associated cost
****************/
SEXP boundedGauss(SEXP cumSum, SEXP cumSumSq, SEXP cumSumVar, SEXP start, SEXP rightIndex, SEXP lower, SEXP upper) {
  // initialise object
  StepGauss data = StepGauss(length(cumSum), REAL(cumSum), REAL(cumSumSq), REAL(cumSumVar), REAL(lower), REAL(upper));
  
  // check lengths
  if(data.N <= 1) error("there must be more than one block");
  if(length(cumSumSq) != (int) data.N) error("length of cumSumSq must match cumSum's");
  if(length(cumSumVar) != (int) data.N) error("length of cumSumVar must match cumSum's");
  if(length(start) != (int) data.N) error("length of start must match cumSum's");
  if(length(lower) != length(upper)) error("lower must have same length as upper");
  if(length(upper) != length(rightIndex)) error("upper must have same length as rightIndex");
  
  Bounds B = Bounds(data.N, INTEGER(start), length(lower), INTEGER(rightIndex), REAL(lower), REAL(upper));

  // run algorithm
  return data.bounded(B); // the optimal feasible solution using minimal number of jumps
}

} // end C wrapper
