#include "Compute.h"
#include <algorithm>

// class computeStatNull
ComputeStatNull::ComputeStatNull(const unsigned int &n) : stat_(NumericVector(n, R_NegInf)) {}

const NumericVector& ComputeStatNull::stat() const {
  return stat_;
}

void ComputeStatNull::compute(Data * const data, const unsigned int &li, const unsigned int &ri) {
  stat_[ri - li] = std::max(stat_[ri - li], data -> computeSingleStatNull());
}


// class computeStat
ComputeStat::ComputeStat(const unsigned int &n) : stat_(NumericVector(n, R_NegInf)) {}

const NumericVector& ComputeStat::stat() const {
  return stat_;
}

void ComputeStat::compute(Data * const data, const unsigned int &li, const unsigned int &ri,
                          const double &value) {
  stat_[ri - li] = std::max(stat_[ri - li], data -> computeSingleStat(value));
}

// class computeBounds
ComputeBounds::ComputeBounds(const unsigned int &numberOfIntervals) : 
  leftIndex_(IntegerVector(numberOfIntervals)), rightIndex_(IntegerVector(numberOfIntervals)),
  lowerBound_(NumericVector(numberOfIntervals)), upperBound_(NumericVector(numberOfIntervals)),
  index_(0u) {}

const IntegerVector& ComputeBounds::leftIndex() const {
  return leftIndex_;
}
const IntegerVector& ComputeBounds::rightIndex() const {
  return rightIndex_;
}
const NumericVector& ComputeBounds::lowerBound() const {
  return lowerBound_;
}
const NumericVector& ComputeBounds::upperBound() const {
  return upperBound_;
}

void ComputeBounds::compute(Data * const data, const unsigned int &li, const unsigned int &ri) {
  leftIndex_[index_] = li + 1u;  // R-style
  rightIndex_[index_] = ri + 1u; // R-style
  SingleBounds currentBounds = data -> computeSingleBounds();
  lowerBound_[index_] = currentBounds.lower();
  upperBound_[index_] = currentBounds.upper();
  ++index_;
}
