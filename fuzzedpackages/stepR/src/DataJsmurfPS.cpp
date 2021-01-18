#include "DataJsmurfPS.h"

#include <algorithm>
#include <cmath>

/***************
* class DataJsmurfPS
* implements abstract class Data for filtered homogeneous gaussian observations using JSMURF idea
* and partial sum test
* Florian Pein, 2017
***************/
NumericVector DataJsmurfPS::data_;
std::vector<double> DataJsmurfPS::varianceSum_;
unsigned int DataJsmurfPS::filterLength_;

DataJsmurfPS::DataJsmurfPS() : cumulatedSum_(0.0), shortendCumulatedSum_(0.0),
                           intervalLength_(0u), shortendIntervalLength_(0u) {}

void DataJsmurfPS::cleanUpStaticVariables() {
  std::vector<double> tmpVarianceSum;
  varianceSum_.swap(tmpVarianceSum);
}

void DataJsmurfPS::setData(const RObject &data, const List &input) {
  data_ = data;
  NumericVector covariances = input["covariances"];
  
  varianceSum_.reserve(data_.size());
  unsigned int m = covariances.size() - 1u;
  
  double sum;
  for (unsigned int i = 1u; i <= static_cast<unsigned int>(data_.size()); ++i) {
    sum = i * covariances[0u];
    for (unsigned int j = 1u; j <= std::min(m, i); ++j) {
      sum += 2.0 * (i - j) * covariances[j];
    }
    varianceSum_.push_back(sum);
  }
  
  filterLength_ = input["filterLength"];
}

unsigned int DataJsmurfPS::getN() const {
  return data_.size();
}
                                                  
Data* DataJsmurfPS::newObject() const {
  return new DataJsmurfPS();
}

SingleBounds DataJsmurfPS::computeSingleBounds() const {
  double mean = shortendCumulatedSum_ / shortendIntervalLength_; 
  double help = std::sqrt(2.0 * criticalValues_[intervalLength_ - 1u] * 
                          varianceSum_[shortendIntervalLength_ - 1u]) / shortendIntervalLength_;
  return SingleBounds(mean - help, mean + help);
}

LocalOptimum DataJsmurfPS::computeLocalOptimum(const unsigned int &leftIndex, const unsigned int &rightIndex, 
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

double DataJsmurfPS::computeSingleStatNull() const {
  return shortendCumulatedSum_ * shortendCumulatedSum_ / varianceSum_[shortendIntervalLength_ - 1u] / 2.0;
}

double DataJsmurfPS::computeSingleStat(const double &value) const {
  double help = shortendCumulatedSum_ - shortendIntervalLength_ * value;
  return help * help / varianceSum_[shortendIntervalLength_ - 1u] / 2.0;
}

void DataJsmurfPS::addLeft(const unsigned int &index) {
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

void DataJsmurfPS::addRight(const unsigned int &index) {
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

void DataJsmurfPS::add(Data * const data) {
  DataJsmurfPS* x = dynamic_cast<DataJsmurfPS *>(data);
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

void DataJsmurfPS::reset() {
  left_ = 0u;
  right_ = 0u;
  cumulatedSum_ = 0.0;
  shortendCumulatedSum_ = 0.0;
  intervalLength_ = 0u;
  shortendIntervalLength_ = 0u;
}
