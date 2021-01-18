#include "IntervalSystemAllLengths.h"
#include "Compute.h"

/***************
* class IntervalSystemAllLengths
* implements the abstract class IntervalSystem for the interval system of all intervals with given lengths
* Florian Pein, 2015
***************/
IntervalSystemAllLengths::IntervalSystemAllLengths(const unsigned int &n, const List &input) : 
IntervalSystem(0u), contained_(std::vector<bool>(n, false)) {
  LogicalVector lengths = input["lengths"];
  for (unsigned int i = 0u; i < n; ++i) {
    if (lengths[i] == TRUE) {
      contained_[i] = true;
      numberOfIntervals_ += n - i;
    }
  }
}
                                    
bool IntervalSystemAllLengths::isInIntervalSystem(const unsigned int &leftIndex,
                                                  const unsigned int &rightIndex) const {
  return contained_[rightIndex - leftIndex];
}          

NumericVector IntervalSystemAllLengths::computeMultiscaleStatisticNull(Data * const data) {  
  ComputeStatNull computeStatNull = ComputeStatNull(data -> getN());
  dynamicProgramAllLengths(data, computeStatNull);
  
  return computeStatNull.stat();
}

List IntervalSystemAllLengths::computeBounds(Data * const data) {
  ComputeBounds computeBounds = ComputeBounds(numberOfIntervals_);
  dynamicProgramAllLengths(data, computeBounds);
  
  return List::create(_["li"] = computeBounds.leftIndex(), _["ri"] = computeBounds.rightIndex(), 
                      _["lower"] = computeBounds.lowerBound(), _["upper"] = computeBounds.upperBound());
}
