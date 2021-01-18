#include "IntervalSystemDyaLenLengths.h"
#include "Compute.h"

/***************
* class IntervalSystemDyaLenLengths
* implements the abstract class IntervalSystem for the interval system dyadic lengths with given lengths
* Florian Pein, 2015
***************/
IntervalSystemDyaLenLengths::IntervalSystemDyaLenLengths(const unsigned int &n, const List &input) : 
IntervalSystem(0u), contained_(std::vector<bool>(n, false)) {
  IntegerVector lengths = input["lengths"];
  for (unsigned int i = 0u; i < static_cast<unsigned int>(lengths.size()); ++i) {
    contained_[lengths[i] - 1] = true;
    numberOfIntervals_ += n - lengths[i] + 1u;
  }
}
                                    
bool IntervalSystemDyaLenLengths::isInIntervalSystem(const unsigned int &leftIndex,
                                                     const unsigned int &rightIndex) const {
  return contained_[rightIndex - leftIndex];
}

NumericVector IntervalSystemDyaLenLengths::computeMultiscaleStatisticNull(Data * const data) {  
  ComputeStatNull computeStatNull = ComputeStatNull(data -> getN());
  dynamicProgramDyaLenLengths(data, computeStatNull);
  
  return computeStatNull.stat();
}

List IntervalSystemDyaLenLengths::computeBounds(Data * const data) {
  ComputeBounds computeBounds = ComputeBounds(numberOfIntervals_);
  dynamicProgramDyaLenLengths(data, computeBounds);
  
  return List::create(_["li"] = computeBounds.leftIndex(), _["ri"] = computeBounds.rightIndex(), 
                      _["lower"] = computeBounds.lowerBound(), _["upper"] = computeBounds.upperBound());
}
