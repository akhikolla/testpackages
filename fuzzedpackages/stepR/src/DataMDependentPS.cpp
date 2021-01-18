#include "DataMDependentPS.h"

#include <algorithm>
#include <cmath>

/***************
* class DataMDependentPS
* implements abstract class Data for m-dependent homogeneous observations
* using a partial sum test statistic
* Florian Pein, 2015
***************/
NumericVector DataMDependentPS::data_;
std::vector<double> DataMDependentPS::varianceSum_;

DataMDependentPS::DataMDependentPS() : cumulatedSum_(0.0), intervalLength_(0u) {}

void DataMDependentPS::cleanUpStaticVariables() {
  std::vector<double> tmpVarianceSum;
  varianceSum_.swap(tmpVarianceSum);
}

void DataMDependentPS::setData(const RObject &data, const List &input) {
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
}

unsigned int DataMDependentPS::getN() const {
  return data_.size();
}
                                                  
Data* DataMDependentPS::newObject() const {
  return new DataMDependentPS();
}             

SingleBounds DataMDependentPS::computeSingleBounds() const {
  double mean = cumulatedSum_ / intervalLength_;  
  double help = std::sqrt(2.0 * criticalValues_[intervalLength_ - 1u] * varianceSum_[intervalLength_ - 1u]) / 
    intervalLength_;
  return SingleBounds(mean - help, mean + help);
}

LocalOptimum DataMDependentPS::computeLocalOptimum(const unsigned int &leftIndex,
                                                   const unsigned int &rightIndex, 
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

double DataMDependentPS::computeSingleStatNull() const {
  return cumulatedSum_ * cumulatedSum_ / varianceSum_[intervalLength_ - 1u] / 2.0;
}

double DataMDependentPS::computeSingleStat(const double &value) const {
  double help = cumulatedSum_ - intervalLength_ * value;
  return help * help / varianceSum_[intervalLength_ - 1u] / 2.0;
}

void DataMDependentPS::addLeft(const unsigned int &index) {
  cumulatedSum_ += data_[index];
  ++intervalLength_;
}

void DataMDependentPS::addRight(const unsigned int &index) {
  cumulatedSum_ += data_[index];
  ++intervalLength_;
}

void DataMDependentPS::add(Data * const data) {
  DataMDependentPS* x = dynamic_cast<DataMDependentPS *>(data);
  cumulatedSum_ += x -> cumulatedSum_;
  intervalLength_ += x -> intervalLength_;
}

void DataMDependentPS::reset() {
  cumulatedSum_ = 0.0;
  intervalLength_ = 0u;
}
