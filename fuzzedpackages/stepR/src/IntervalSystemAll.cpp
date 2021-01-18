#include "IntervalSystemAll.h"
#include "Compute.h"

/***************
* class IntervalSystemAll
* implements the abstract class IntervalSystem for the interval system of all intervals
* Florian Pein, 2015
***************/

IntervalSystemAll::IntervalSystemAll(const unsigned int &n) : IntervalSystem(n * (n + 1u) / 2u) {}
                                    
bool IntervalSystemAll::isInIntervalSystem(const unsigned int &leftIndex,
                                           const unsigned int &rightIndex) const {
  return true;
}

NumericVector IntervalSystemAll::computeMultiscaleStatisticNull(Data * const data) {  
  ComputeStatNull computeStatNull = ComputeStatNull(data -> getN());
  dynamicProgramAll(data, computeStatNull);
  
  return computeStatNull.stat();
} 

NumericVector IntervalSystemAll::computeMultiscaleStatistic(Data * const data, const List &input) {
  ComputeStat computeStat = ComputeStat(data -> getN());
  IntegerVector left  = input["leftIndex"];
  IntegerVector right  = input["rightIndex"];
  NumericVector value  = input["value"];
  
  for (unsigned int k = 0; k < static_cast<unsigned int>(value.size()); ++k) {
    for (unsigned int j = left[k]; j <= static_cast<unsigned int>(right[k]); ++j) {
      checkUserInterrupt();
      data -> reset();
      
      unsigned int i = j + 1u;
      while (i > static_cast<unsigned int>(left[k])) {
        --i;
        data -> addLeft(i);
        computeStat.compute(data, i, j, value[k]);
      }
    }
  }

  return computeStat.stat();
}

List IntervalSystemAll::computeBounds(Data * const data) {
  ComputeBounds computeBounds = ComputeBounds(numberOfIntervals_);
  dynamicProgramAll(data, computeBounds);
  
  return List::create(_["li"] = computeBounds.leftIndex(), _["ri"] = computeBounds.rightIndex(), 
                      _["lower"] = computeBounds.lowerBound(), _["upper"] = computeBounds.upperBound());
}
