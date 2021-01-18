#include "IntervalSystemDyaLen.h"
#include "Compute.h"

/***************
* class IntervalSystemDyaLen
* implements the abstract class IntervalSystem for the interval system of all intervals for given lengths
* Florian Pein, 2015
***************/
IntervalSystemDyaLen::IntervalSystemDyaLen(const unsigned int &n) : 
IntervalSystem(0u), contained_(std::vector<bool>(n, false)) {
  for (unsigned int i = 1u; i <= n; i *= 2u) {
    contained_[i - 1u] = true;
    numberOfIntervals_ += n - i + 1u;
  }
}
                                    
bool IntervalSystemDyaLen::isInIntervalSystem(const unsigned int &leftIndex,
                                              const unsigned int &rightIndex) const {
  return contained_[rightIndex - leftIndex];
}

NumericVector IntervalSystemDyaLen::computeMultiscaleStatisticNull(Data * const data) {  
  ComputeStatNull computeStatNull = ComputeStatNull(data -> getN());
  dynamicProgramDyaLen(data, computeStatNull);
  
  return computeStatNull.stat();
}

List IntervalSystemDyaLen::computeBounds(Data * const data) {
  ComputeBounds computeBounds = ComputeBounds(numberOfIntervals_);
  dynamicProgramDyaLen(data, computeBounds);
  
  return List::create(_["li"] = computeBounds.leftIndex(), _["ri"] = computeBounds.rightIndex(), 
                      _["lower"] = computeBounds.lowerBound(), _["upper"] = computeBounds.upperBound());
}
