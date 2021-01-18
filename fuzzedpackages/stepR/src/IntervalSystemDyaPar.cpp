#include "IntervalSystemDyaPar.h"
#include "Compute.h"

/***************
* class IntervalSystemDyaPar
* implements the abstract class IntervalSystem for the interval system dyadic partition
* Florian Pein, 2015
***************/
IntervalSystemDyaPar::IntervalSystemDyaPar(const unsigned int &n) : 
IntervalSystem(0u), contained_(std::vector<bool>(n, false)) {
  for (unsigned int i = 1u; i <= n; i *= 2u) {
    contained_[i - 1u] = true;
    numberOfIntervals_ += n / i;
  }
}
                                    
bool IntervalSystemDyaPar::isInIntervalSystem(const unsigned int &leftIndex,
                                              const unsigned int &rightIndex) const {
  return contained_[rightIndex - leftIndex] && (leftIndex % (rightIndex - leftIndex + 1u) == 0u);
}

NumericVector IntervalSystemDyaPar::computeMultiscaleStatisticNull(Data * const data) {  
  ComputeStatNull computeStatNull = ComputeStatNull(data -> getN());
  dynamicProgramDyaPar(data, computeStatNull);
  
  return computeStatNull.stat();
}

List IntervalSystemDyaPar::computeBounds(Data * const data) {
  ComputeBounds computeBounds = ComputeBounds(numberOfIntervals_);
  dynamicProgramDyaPar(data, computeBounds);
  
  return List::create(_["li"] = computeBounds.leftIndex(), _["ri"] = computeBounds.rightIndex(), 
                      _["lower"] = computeBounds.lowerBound(), _["upper"] = computeBounds.upperBound()); 
}
