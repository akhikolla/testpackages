#include "DataHsmuce.h"

#include <algorithm>
#include <cmath>

/***************
* class DataHsmuce
* implements abstract class Data for independent, heterogeneous gaussian observations
* with the mean as parameter of interest and the variance as nuisance parameter
* Florian Pein, 2015
***************/
NumericVector DataHsmuce::data_;

DataHsmuce::DataHsmuce() : cumulatedSum_(0.0), cumulatedSumSquared_(0.0), intervalLength_(0u) {}

void DataHsmuce::cleanUpStaticVariables() {}

void DataHsmuce::setData(const RObject &data) {
  data_ = data;
}

unsigned int DataHsmuce::getN() const {
  return data_.size();
}
                                                  
Data* DataHsmuce::newObject() const {
  return new DataHsmuce();
}          

SingleBounds DataHsmuce::computeSingleBounds() const {
  double mean = cumulatedSum_ / intervalLength_;
  double varEst = (cumulatedSumSquared_ - cumulatedSum_ * cumulatedSum_ / intervalLength_) /
    (intervalLength_ - 1u);
  double help = std::sqrt(2.0 * criticalValues_[intervalLength_ - 1u] / intervalLength_ * varEst);
  return SingleBounds(mean - help, mean + help);
}

LocalOptimum DataHsmuce::computeLocalOptimum(const unsigned int &leftIndex, const unsigned int &rightIndex, 
                                             const SingleBounds &localBounds,
                                             LocalOptimum * const lastSegment) const {
  double mean = cumulatedSum_ / intervalLength_;
  double est = std::max(localBounds.lower(), std::min(localBounds.upper(), mean));
  double bias = est - mean;
  double costs;
  if (intervalLength_ == 1u) {
//     Rcout.precision(17);
//     Rcout << "Pi: " << std::fixed << - std::log(2.0 * std::acos(-1.)) - 1.0 << std::endl;
    costs = -2.83787706640934534;
  } else {
    costs = std::log(cumulatedSumSquared_ / intervalLength_ + bias * bias - mean * mean) * intervalLength_;
  }
  
  if (leftIndex == 0u) {
    return LocalOptimum(leftIndex, rightIndex, est, costs, lastSegment);
  }
  
  return LocalOptimum(leftIndex, rightIndex, est, lastSegment -> costs() + costs, lastSegment);
}

double DataHsmuce::computeSingleStatNull() const {
  double varEst = (cumulatedSumSquared_ - cumulatedSum_ * cumulatedSum_ / intervalLength_) /
    (intervalLength_ - 1u);
  return cumulatedSum_ * cumulatedSum_ / intervalLength_ / varEst / 2.0;
}

double DataHsmuce::computeSingleStat(const double &value) const {
  double varEst = (cumulatedSumSquared_ - cumulatedSum_ * cumulatedSum_ / intervalLength_) /
    (intervalLength_ - 1u);
  double help = cumulatedSum_ - intervalLength_ * value;
  return help * help / intervalLength_ / varEst / 2.0;
}

void DataHsmuce::addLeft(const unsigned int &index) {
  cumulatedSum_ += data_[index];
  cumulatedSumSquared_ += data_[index] * data_[index];
  ++intervalLength_;
}

void DataHsmuce::addRight(const unsigned int &index) {
  cumulatedSum_ += data_[index];
  cumulatedSumSquared_ += data_[index] * data_[index];
  ++intervalLength_;
}

void DataHsmuce::add(Data * const data) {
  DataHsmuce* x = dynamic_cast<DataHsmuce *>(data);
  cumulatedSum_ += x -> cumulatedSum_;
  cumulatedSumSquared_ += x -> cumulatedSumSquared_;
  intervalLength_ += x -> intervalLength_;
}

void DataHsmuce::reset() {
  cumulatedSum_ = 0.0;
  cumulatedSumSquared_ = 0.0;
  intervalLength_ = 0u;
}
