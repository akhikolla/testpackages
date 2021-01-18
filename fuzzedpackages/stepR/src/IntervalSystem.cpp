#include "IntervalSystem.h"
#include "Compute.h"

/***************
* class IntervalSystem
* abstract class maintaining the interval system
* Florian Pein, 2015
***************/

IntervalSystem::IntervalSystem(const unsigned int &numberOfIntervals) : 
  numberOfIntervals_(numberOfIntervals) {}

IntervalSystem::~IntervalSystem() {}

NumericVector IntervalSystem::computeMultiscaleStatisticNull(Data * const data) {
  ComputeStatNull computeStatNull = ComputeStatNull(data -> getN());
  dynamicProgram(data, computeStatNull);
  
  return computeStatNull.stat();
}

NumericVector IntervalSystem::computeMultiscaleStatistic(Data * const data, const List &input) {
  ComputeStat computeStat = ComputeStat(data -> getN());
  IntegerVector left  = input["leftIndex"];
  IntegerVector right  = input["rightIndex"];
  NumericVector value  = input["value"];
  
  for (unsigned int k = 0u; k < static_cast<unsigned int>(value.size()); ++k) {
    for (unsigned int j = left[k]; j <= static_cast<unsigned int>(right[k]); ++j) {
      checkUserInterrupt();
      data -> reset();
      
      unsigned int i = j + 1u;
      while (i > static_cast<unsigned int>(left[k])) {
        --i;
        data -> addLeft(i);
        if (isInIntervalSystem(i, j)) {
          computeStat.compute(data, i, j, value[k]);
        }
      }
    }
  }
  
  return computeStat.stat();
}

List IntervalSystem::computeBounds(Data * const data) {
  ComputeBounds computeBounds = ComputeBounds(numberOfIntervals_);
  dynamicProgram(data, computeBounds);
  
  return List::create(_["li"] = computeBounds.leftIndex(), _["ri"] = computeBounds.rightIndex(), 
                      _["lower"] = computeBounds.lowerBound(), _["upper"] = computeBounds.upperBound());
}
