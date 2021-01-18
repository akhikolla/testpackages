#include "StepBinom.h"
#include <cmath>

/***************
* class StepBinom
* implements virtual class Step for Binomial data
* Thomas Hotz, 2007-2011
***************/

/*************
* constructor for n data points following Binomial distribution with size trials, with pointers to array of that length
****************/
StepBinom::StepBinom(unsigned int size, unsigned int n, int* xcs, double* xcw) : Step(n), size(size), cs(xcs), cw(xcw) {}

/*************
* constructor for n data points and bounds
****************/
StepBinom::StepBinom(unsigned int size, unsigned int n, int* xcs, double* xcw, double* xlb, double* xub) : Step(n, xlb, xub), size(size), cs(xcs), cw(xcw) {}

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
* for binomial data: the negative log-likelihood of the mean on that block,
* i.e. (cw * size) * log(cw * size) - cs * log(cs) - (cw * size - cs) * log(cw * size - cs), up to a constant -sum_i log choose(size, y[i]) depending only on the data
****************/
double StepBinom::cost(unsigned int startIndex, unsigned int endIndex) const {
  if(startIndex == 0) { // all data
    if(cs[endIndex] == 0 || cw[endIndex] * size - cs[endIndex] < 0.5 ) {
      return 0; // likelihood is 1 since all data are identically 0 or identically size
    } else {
      return cw[endIndex] * size * std::log( cw[endIndex] * size ) - cs[endIndex] * std::log(double( cs[endIndex] )) - ( cw[endIndex] * size - cs[endIndex] ) * std::log( cw[endIndex] * size - cs[endIndex] );
    }
  } else {
    if(cs[endIndex] - cs[startIndex - 1] == 0 || (cw[endIndex] - cw[startIndex - 1]) * size - (cs[endIndex] - cs[startIndex - 1]) < 0.5 ) {
//       Rprintf("    si = %d, ei = %d, cssi1 = %d, csei = %d\n", startIndex, endIndex, cs[startIndex - 1], cs[endIndex]);
      return 0; // likelihood is 1 since all data are identically 0 or identically size
    } else {
      return (cw[endIndex] - cw[startIndex - 1]) * size * std::log( (cw[endIndex] - cw[startIndex - 1]) * size ) - (cs[endIndex] - cs[startIndex - 1]) * std::log(double( (cs[endIndex] - cs[startIndex - 1]) )) - ( (cw[endIndex] - cw[startIndex - 1]) * size - (cs[endIndex] - cs[startIndex - 1]) ) * std::log( (cw[endIndex] - cw[startIndex - 1]) * size - (cs[endIndex] - cs[startIndex - 1]) );
    }
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
* for binomial data: the negative log-likelihood of the bounded mean on that block, i.e. cs * ( log(p) - log(1-p) ) + cw * size * log(1-p), up to a constant -sum_i log choose(size, y[i]) depending only on the data
****************/
double StepBinom::costBound(unsigned int startIndex, unsigned int endIndex, const LUBound& bound) const {
  if(bound.lower > bound.upper || bound.upper < 0 || bound.lower > 1) return R_PosInf; else {
    int S, NS;
    double weight;
    if(startIndex == 0) {
      S = cs[endIndex];
      weight = cw[endIndex];
    } else {
      S = cs[endIndex] - cs[startIndex - 1];
      weight = cw[endIndex] - cw[startIndex - 1];
    }
    NS = weight * size - S;
    if( ( (S == 0) & (bound.lower <= 0) ) || ( (NS < 0.5) & (bound.upper >= 1) ) ) {
      return 0; // likelihood is 1
    } else if( ( (S != 0) & (bound.upper == 0) ) || ( (NS >= 0.5) & (bound.lower == 1) ) ) {
      return R_PosInf; // impossible
    } else {
      double p = fmax2( fmin2( S / weight / size, bound.upper ), bound.lower );
      return -S * std::log(p) - NS * std::log(1-p);
    }
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
* for binomial data: the mean / size on this block if it isn't outside a boundary, then the latter
****************/
double StepBinom::estBound(unsigned int startIndex, unsigned int endIndex, const LUBound& bound) const {
  if(bound.lower > bound.upper || bound.upper < 0 || bound.lower > 1) return R_NaN; else
    if(startIndex == 0) {
      double p = cs[endIndex] / cw[endIndex] / size;
      #ifdef DEBUGbounded
      Rprintf("  si = %d, ei = %d, p = %4.2e, lower = %4.2e, upper = %4.2e\n", startIndex, endIndex, p, bound.lower, bound.upper);
      #endif
      return fmax2( fmin2( p, bound.upper ), bound.lower );
    } else {
      double weight = cw[endIndex] - cw[startIndex - 1];
      double p = ( cs[endIndex] - cs[startIndex - 1] ) / weight / size;
      #ifdef DEBUGbounded
      Rprintf("  si = %d, ei = %d, p = %4.2e, lower = %4.2e, upper = %4.2e\n", startIndex, endIndex, p, bound.lower, bound.upper);
      #endif
      return fmax2( fmin2( p, bound.upper ), bound.lower );
    }
}

// C wrapper
extern "C" {

/*************
* forwardBinom
* function to be called from R
* computes forward selection for block signals, wrapper for StepBinom::forward
*
* in:
* size : a single integer giving the number of trials, specifying the Binomial distribution
* cumSum : vector of cumulative sums
* cumSumWe : vector of cumulative sums of weights
* maxBlocks : a single integer giving the maximal number of blocks
*
* out:
* data.frame comprising rightIndex, number, depth and improve of the candidates
****************/
SEXP forwardBinom(SEXP size, SEXP cumSum, SEXP cumSumWe, SEXP maxBlocks) {
  // initialise object
  StepBinom data = StepBinom(asInteger(size), length(cumSum), INTEGER(cumSum), REAL(cumSumWe));
  
  // check lengths
  if(data.N < 1) error("cumSum must have at least one element");
  if(length(cumSumWe) != (int) data.N) error("cumSumWe must have same length as cumSum");
  if(length(maxBlocks) != 1) error("maxBlocks must be a single integer");
  
  // run algorithm
  return data.forward(asInteger(maxBlocks));
}

/*************
* pathBinom
* function to be called from R
* computes path of Potts minimisers, wrapper for StepBinom::path
*
* in:
* size : a single integer giving the number of trials, specifying the Binomial distribution
* cumSum : a numeric vector, the cumulative sums of the signal
* cumSumWe : vector of cumulative sums of weights
* maxBlocks : a single integer giving the maximal number of blocks
*
* out:
* a list conprimising the right indices of the solution path, and the associated costs
****************/
SEXP pathBinom(SEXP size, SEXP cumSum, SEXP cumSumWe, SEXP maxBlocks) {
  // initialise object
  StepBinom data = StepBinom(asInteger(size), length(cumSum), INTEGER(cumSum), REAL(cumSumWe));
  
  // check lengths
  if(data.N <= 1) error("there must be more than one block");
  if(length(cumSumWe) != (int) data.N) error("cumSumWe must have same length as cumSum");
  if(length(maxBlocks) != 1) error("maxBlocks must be a single integer");
  
  // run algorithm
  return data.path(asInteger(maxBlocks)); // the solution path, i.e. p[i, k] is the (i+1)th jump in the solution having k+1 jumps
}

/*************
* boundedBinom
* function to be called from R
* computes bounded solution, wrapper for StepBinom::bounded
*
* in:
* size : a single integer giving the number of trials, specifying the Binomial distribution
* cumSum : a numeric vector, the cumulative sums of the signal
* cumSumWe : vector of cumulative sums of weights
* start : for every possible left index where intervals with this left index start in the list of intervals (increasing in left indices), NA if none
* rightIndex : right indices of the intervals in the list, increasing for each left index
* lower : the lower bounds for the estimator at the respective interval
* upper : the upper bounds for the estimator at the respective interval
*
* out:
* a list conprimising the right indices and values of the solution, and the associated cost
****************/
SEXP boundedBinom(SEXP size, SEXP cumSum, SEXP cumSumWe, SEXP start, SEXP rightIndex, SEXP lower, SEXP upper) {
  // initialise object
  StepBinom data = StepBinom(asInteger(size), length(cumSum), INTEGER(cumSum), REAL(cumSumWe), REAL(lower), REAL(upper));
  
  // check lengths
  if(data.N <= 1) error("there must be more than one block");
  if(length(cumSumWe) != (int) data.N) error("length of cumSumWe must match cumSum's");
  if(length(start) != (int) data.N) error("length of start must match cumSum's");
  if(length(lower) != length(upper)) error("lower must have same length as upper");
  if(length(upper) != length(rightIndex)) error("upper must have same length as rightIndex");
  
  Bounds B = Bounds(data.N, INTEGER(start), length(lower), INTEGER(rightIndex), REAL(lower), REAL(upper));
  
  // run algorithm
  return data.bounded(B); // the optimal feasible solution using minimal number of jumps
}

} // end C wrapper
