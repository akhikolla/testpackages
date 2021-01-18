#include "IntervalSystemDyaParLengths.h"
#include "Compute.h"

/***************
* class IntervalSystemDyaParLengths
* implements the abstract class IntervalSystem for the interval system dyadic partition with given lengths
* Florian Pein, 2015
***************/
IntervalSystemDyaParLengths::IntervalSystemDyaParLengths(const unsigned int &n, const List &input) : 
IntervalSystem(0u), contained_(std::vector<bool>(n, false)) {
  IntegerVector lengths = input["lengths"];
  for (unsigned int i = 0u; i < static_cast<unsigned int>(lengths.size()); ++i) {
    contained_[lengths[i] - 1] = true;
    numberOfIntervals_ += n / lengths[i];
  }
}
                                    
bool IntervalSystemDyaParLengths::isInIntervalSystem(const unsigned int &leftIndex,
                                                     const unsigned int &rightIndex) const {
  return contained_[rightIndex - leftIndex] && (leftIndex % (rightIndex - leftIndex + 1u) == 0u);
} 

NumericVector IntervalSystemDyaParLengths::computeMultiscaleStatisticNull(Data * const data) {  
  ComputeStatNull computeStatNull = ComputeStatNull(data -> getN());
  dynamicProgramDyaParLengths(data, computeStatNull);
  
  return computeStatNull.stat();
} 

List IntervalSystemDyaParLengths::computeBounds(Data * const data) {
  ComputeBounds computeBounds = ComputeBounds(numberOfIntervals_);
  dynamicProgramDyaParLengths(data, computeBounds);
  
  return List::create(_["li"] = computeBounds.leftIndex(), _["ri"] = computeBounds.rightIndex(), 
                      _["lower"] = computeBounds.lowerBound(), _["upper"] = computeBounds.upperBound());
}
