#include "DataHjsmurfSPS.h"

#include <algorithm>
#include <cmath>

/***************
* class DataHjsmurfSPS
* implements abstract class Data for filtered heterogeneous gaussian observations using JSMURF idea
* and studentized partial sum test
* Florian Pein, 2017
***************/
NumericVector DataHjsmurfSPS::data_;
unsigned int DataHjsmurfSPS::filterLength_;
std::vector<double> DataHjsmurfSPS::correlationSum_;

DataHjsmurfSPS::DataHjsmurfSPS() : cumulatedSum_(0.0), shortendCumulatedSum_(0.0),
                             cumulatedSumSquared_(0.0), shortendCumulatedSumSquared_(0.0),
                             intervalLength_(0u), shortendIntervalLength_(0u) {}

void DataHjsmurfSPS::cleanUpStaticVariables() {
  std::vector<double> tmpCorrelationSum;
  correlationSum_.swap(tmpCorrelationSum);
}

void DataHjsmurfSPS::setData(const RObject &data, const List &input) {
  data_ = data;
  filterLength_ = input["filterLength"];
  NumericVector correlations = input["correlations"];
  
  correlationSum_.reserve(data_.size());
  unsigned int m = correlations.size() - 1u;
  
  double sum;
  for (unsigned int i = 1u; i <= static_cast<unsigned int>(data_.size()); ++i) {
    sum = i * correlations[0u];
    for (unsigned int j = 1u; j <= std::min(m, i); ++j) {
      sum += 2.0 * (i - j) * correlations[j];
    }
    correlationSum_.push_back(sum);
  }
}

unsigned int DataHjsmurfSPS::getN() const {
  return data_.size();
}
                                                  
Data* DataHjsmurfSPS::newObject() const {
  return new DataHjsmurfSPS();
}

SingleBounds DataHjsmurfSPS::computeSingleBounds() const {
  double mean = shortendCumulatedSum_ / shortendIntervalLength_;
  double varEst = (shortendCumulatedSumSquared_ - shortendCumulatedSum_ * shortendCumulatedSum_ / 
                   shortendIntervalLength_) / shortendIntervalLength_ / 
                   (1.0 - correlationSum_[shortendIntervalLength_ - 1u] / shortendIntervalLength_ /
                     shortendIntervalLength_);
  double help = std::sqrt(2.0 * criticalValues_[intervalLength_ - 1u] *
                          correlationSum_[shortendIntervalLength_ - 1u] * varEst) / shortendIntervalLength_;
  return SingleBounds(mean - help, mean + help);
}

LocalOptimum DataHjsmurfSPS::computeLocalOptimum(const unsigned int &leftIndex,
                                                 const unsigned int &rightIndex, 
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

double DataHjsmurfSPS::computeSingleStatNull() const {
  double varEst = (shortendCumulatedSumSquared_ - shortendCumulatedSum_ * shortendCumulatedSum_ / 
                   shortendIntervalLength_) / shortendIntervalLength_ / 
                   (1.0 - correlationSum_[shortendIntervalLength_ - 1u] / shortendIntervalLength_ /
                     shortendIntervalLength_);
  return shortendCumulatedSum_ * shortendCumulatedSum_ / correlationSum_[shortendIntervalLength_ - 1u] /
    varEst / 2.0;
}

double DataHjsmurfSPS::computeSingleStat(const double &value) const {
  double varEst = (shortendCumulatedSumSquared_ - shortendCumulatedSum_ * shortendCumulatedSum_ / 
                   shortendIntervalLength_) / shortendIntervalLength_ / 
                   (1.0 - correlationSum_[shortendIntervalLength_ - 1u] / shortendIntervalLength_ /
                     shortendIntervalLength_);
  double help = shortendCumulatedSum_ - shortendIntervalLength_ * value;
  return help * help / correlationSum_[shortendIntervalLength_ - 1u] / varEst / 2.0;
}

void DataHjsmurfSPS::addLeft(const unsigned int &index) {
  left_ = index;
  if (intervalLength_ == 0u) {
    right_ = index;
  }
  
  intervalLength_++;
  cumulatedSum_ += data_[index];
  cumulatedSumSquared_ += data_[index] * data_[index];
  
  if (intervalLength_ > filterLength_) {
    shortendIntervalLength_++;
    shortendCumulatedSum_ += data_[index + filterLength_];
    shortendCumulatedSumSquared_ += data_[index + filterLength_] * data_[index + filterLength_];
  }
}

void DataHjsmurfSPS::addRight(const unsigned int &index) {
  right_ = index;
  if (intervalLength_ == 0u) {
    left_ = index;
  }
  
  intervalLength_++;
  cumulatedSum_ += data_[index];
  cumulatedSumSquared_ += data_[index] * data_[index];
  
  if (intervalLength_ > filterLength_) {
    shortendIntervalLength_++;
    shortendCumulatedSum_ += data_[index];
    shortendCumulatedSumSquared_ += data_[index] * data_[index];
  }
}

void DataHjsmurfSPS::add(Data * const data) {
  DataHjsmurfSPS* x = dynamic_cast<DataHjsmurfSPS *>(data);
  right_ = x -> right_;
  
  if (intervalLength_ >= filterLength_) {
    shortendCumulatedSum_ += x -> cumulatedSum_;
    shortendCumulatedSumSquared_ += x -> cumulatedSumSquared_;
    shortendIntervalLength_ += x -> intervalLength_;
  } else {
    if (intervalLength_ + x -> intervalLength_ > filterLength_) {
      for (unsigned int index = left_ + filterLength_; index <= right_; ++index) {
        shortendIntervalLength_++;
        shortendCumulatedSum_ += data_[index];
        shortendCumulatedSumSquared_ += data_[index] * data_[index];
      }
    }
  }
  cumulatedSum_ += x -> cumulatedSum_;
  cumulatedSumSquared_ += x -> cumulatedSumSquared_;
  intervalLength_ += x -> intervalLength_;
}

void DataHjsmurfSPS::reset() {
  left_ = 0u;
  right_ = 0u;
  cumulatedSum_ = 0.0;
  shortendCumulatedSum_ = 0.0;
  cumulatedSumSquared_ = 0.0;
  shortendCumulatedSumSquared_ = 0.0;
  intervalLength_ = 0u;
  shortendIntervalLength_ = 0u;
}
