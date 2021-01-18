#include "DataJsmurf.h"

#include <algorithm>
#include <cmath>

/***************
* class DataJsmurf
* implements abstract class Data for filtered homogeneous gaussian observations using JSMURF idea
* Florian Pein, 2017
***************/
NumericVector DataJsmurf::data_;
double DataJsmurf::standardDeviation_;
unsigned int DataJsmurf::filterLength_;

DataJsmurf::DataJsmurf() : cumulatedSum_(0.0), shortendCumulatedSum_(0.0),
                           intervalLength_(0u), shortendIntervalLength_(0u) {}

void DataJsmurf::cleanUpStaticVariables() {}

void DataJsmurf::setData(const RObject &data, const List &input) {
  data_ = data;
  standardDeviation_ = input["sd"];
  filterLength_ = input["filterLength"];
}

unsigned int DataJsmurf::getN() const {
  return data_.size();
}
                                                  
Data* DataJsmurf::newObject() const {
  return new DataJsmurf();
}

SingleBounds DataJsmurf::computeSingleBounds() const {
  double mean = shortendCumulatedSum_ / shortendIntervalLength_;  
  double help = std::sqrt(2.0 * criticalValues_[intervalLength_ - 1u] / shortendIntervalLength_) *
    standardDeviation_;
  return SingleBounds(mean - help, mean + help);
}

LocalOptimum DataJsmurf::computeLocalOptimum(const unsigned int &leftIndex, const unsigned int &rightIndex, 
                                            const SingleBounds &localBounds,
                                            LocalOptimum * const lastSegment) const {
  double mean = cumulatedSum_ / intervalLength_;
  double est = std::max(localBounds.lower(), std::min(localBounds.upper(), mean));
  double bias = est - mean;
  double costs = (bias * bias - mean * mean) * intervalLength_;
  if (leftIndex == 0u) {
    return LocalOptimum(leftIndex, rightIndex, est, costs, lastSegment);
  }
  
  return LocalOptimum(leftIndex, rightIndex, est, lastSegment -> costs() + costs, lastSegment);
}

double DataJsmurf::computeSingleStatNull() const {
  double help = shortendCumulatedSum_ / standardDeviation_;
  return help * help / 2.0 / shortendIntervalLength_;
}

double DataJsmurf::computeSingleStat(const double &value) const {
  double help = (shortendCumulatedSum_ - shortendIntervalLength_ * value) / standardDeviation_;
  return help * help / 2.0 / shortendIntervalLength_;
}

void DataJsmurf::addLeft(const unsigned int &index) {
  left_ = index;
  if (intervalLength_ == 0u) {
    right_ = index;
  }
  
  intervalLength_++;
  cumulatedSum_ += data_[index];
  
  if (intervalLength_ > filterLength_) {
    shortendIntervalLength_++;
    shortendCumulatedSum_ += data_[index + filterLength_];
  }
}

void DataJsmurf::addRight(const unsigned int &index) {
  right_ = index;
  if (intervalLength_ == 0u) {
    left_ = index;
  }
  
  intervalLength_++;
  cumulatedSum_ += data_[index];
  
  if (intervalLength_ > filterLength_) {
    shortendIntervalLength_++;
    shortendCumulatedSum_ += data_[index];
  }
}

void DataJsmurf::add(Data * const data) {
  DataJsmurf* x = dynamic_cast<DataJsmurf *>(data);
  right_ = x -> right_;
  
  if (intervalLength_ >= filterLength_) {
    shortendCumulatedSum_ += x -> cumulatedSum_;
    shortendIntervalLength_ += x -> intervalLength_;
  } else {
    if (intervalLength_ + x -> intervalLength_ > filterLength_) {
      for (unsigned int index = left_ + filterLength_; index <= right_; ++index) {
        shortendIntervalLength_++;
        shortendCumulatedSum_ += data_[index];
      }
    }
  }
  cumulatedSum_ += x -> cumulatedSum_;
  intervalLength_ += x -> intervalLength_;
}

void DataJsmurf::reset() {
  left_ = 0u;
  right_ = 0u;
  cumulatedSum_ = 0.0;
  shortendCumulatedSum_ = 0.0;
  intervalLength_ = 0u;
  shortendIntervalLength_ = 0u;
}
