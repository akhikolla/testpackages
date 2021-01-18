#include "DataGauss.h"

#include <algorithm>
#include <cmath>

/***************
* class DataGauss
* implements abstract class Data for independent, homogeneous gaussian observations
* Florian Pein, 2015
***************/
NumericVector DataGauss::data_;
double DataGauss::standardDeviation_;

DataGauss::DataGauss() : cumulatedSum_(0.0), intervalLength_(0u) {}

void DataGauss::cleanUpStaticVariables() {}

void DataGauss::setData(const RObject &data, const List &input) {
  data_ = data;
  standardDeviation_ = input["sd"];
}

unsigned int DataGauss::getN() const {
  return data_.size();
}
                                                  
Data* DataGauss::newObject() const {
  return new DataGauss();
}

SingleBounds DataGauss::computeSingleBounds() const {
  double mean = cumulatedSum_ / intervalLength_;  
  double help = std::sqrt(2.0 * criticalValues_[intervalLength_ - 1u] / intervalLength_) * standardDeviation_;
  return SingleBounds(mean - help, mean + help);
}

LocalOptimum DataGauss::computeLocalOptimum(const unsigned int &leftIndex, const unsigned int &rightIndex, 
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

double DataGauss::computeSingleStatNull() const {
  double help = cumulatedSum_ / standardDeviation_;
  return help * help / 2.0 / intervalLength_;
}

double DataGauss::computeSingleStat(const double &value) const {
  double help = (cumulatedSum_ - intervalLength_ * value) / standardDeviation_;
  return help * help / 2.0 / intervalLength_;
}

void DataGauss::addLeft(const unsigned int &index) {
  cumulatedSum_ += data_[index];
  ++intervalLength_;
}

void DataGauss::addRight(const unsigned int &index) {
  cumulatedSum_ += data_[index];
  ++intervalLength_;
}

void DataGauss::add(Data * const data) {
  DataGauss* x = dynamic_cast<DataGauss *>(data);
  cumulatedSum_ += x -> cumulatedSum_;
  intervalLength_ += x -> intervalLength_;
}

void DataGauss::reset() {
  cumulatedSum_ = 0.0;
  intervalLength_ = 0u;
}
