#include "DataJsmurfLR.h"

#include <algorithm>
#include <cmath>

#include "choleskyDecomposition.h"
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>

/***************
* class DataJsmurfLR
* implements abstract class Data for filtered homogeneous gaussian observations using JSMURF idea
* and the likelihood ratio test
* Florian Pein, 2017
***************/
NumericVector DataJsmurfLR::data_;
unsigned int DataJsmurfLR::filterLength_;

int DataJsmurfLR::m_ = 0;
NumericVector DataJsmurfLR::covariances_;
std::vector<bool> DataJsmurfLR::isComputed_;
std::vector<double*> DataJsmurfLR::sigmaInverseOne_;
std::vector<double> DataJsmurfLR::denominator_;

const char DataJsmurfLR::uplo_ = 'U';
const char DataJsmurfLR::trans1_ = 'T';
const char DataJsmurfLR::trans2_ = 'N';
const char DataJsmurfLR::diag_ = 'N';
const int DataJsmurfLR::incx_ = 1;

void DataJsmurfLR::compute(const int &size) {
  const int bands = std::min<int>(DataJsmurfLR::m_ - 1, size - 1);
  
  const int ldA = bands + 1;
  
  double* A = choleskyDecomposition(size, covariances_);
  
  double* x;
  x = new double[size];
  
  std::fill(x, x + size, 1.0);

  F77_CALL(dtbsv)(&uplo_, &trans1_, &diag_, &size, &bands, A, &ldA, x, &incx_);
  F77_CALL(dtbsv)(&uplo_, &trans2_, &diag_, &size, &bands, A, &ldA, x, &incx_);
  
  DataJsmurfLR::sigmaInverseOne_[size - 1u] = x;
  for (int i = 0; i < size; ++i) {
    denominator_[size - 1] += x[i];
  }
  
  isComputed_[size - 1] = true;

  delete[] A;
} 

DataJsmurfLR::DataJsmurfLR() : cumulatedSum_(0.0), shortendCumulatedSum_(0.0),
  intervalLength_(0u), shortendIntervalLength_(0u) {}

void DataJsmurfLR::cleanUpStaticVariables() {
  for (unsigned int i = 0u; i < data_.size(); ++i) {
    if (isComputed_[i]) {
      delete[] sigmaInverseOne_[i]; 
    }
  }
  
  std::vector<bool> tmpIsComputed;
  isComputed_.swap(tmpIsComputed);
  
  std::vector<double*> tmpSigmaInverseOne;
  sigmaInverseOne_.swap(tmpSigmaInverseOne);
  
  std::vector<double> tmpDenominator;
  denominator_.swap(tmpDenominator);
}

void DataJsmurfLR::setData(const RObject &data, const List &input) {
  data_ = data;
  filterLength_ = input["filterLength"];
  covariances_ = input["covariances"];
  m_ = covariances_.size();
  isComputed_.resize(data_.size(), false);
  sigmaInverseOne_.resize(data_.size());
  denominator_.resize(data_.size(), 0.0);
}

unsigned int DataJsmurfLR::getN() const {
  return data_.size();
}
                                                  
Data* DataJsmurfLR::newObject() const {
  return new DataJsmurfLR();
}

SingleBounds DataJsmurfLR::computeSingleBounds() const {
  if (!isComputed_[shortendIntervalLength_ - 1u]) {
    compute(shortendIntervalLength_);
  }
  
  double numerator = 0.0;
  for (unsigned int i = 0u; i < shortendIntervalLength_; i++) {
    numerator += sigmaInverseOne_[shortendIntervalLength_ - 1u][i] * data_[left_ + filterLength_ + i];
  }
  
  double help = std::sqrt(2.0 * criticalValues_[intervalLength_ - 1u] /
                          denominator_[shortendIntervalLength_ - 1u]);
  return SingleBounds(numerator / denominator_[shortendIntervalLength_ - 1u] - help,
                      numerator / denominator_[shortendIntervalLength_ - 1u] + help);
}

LocalOptimum DataJsmurfLR::computeLocalOptimum(const unsigned int &leftIndex, const unsigned int &rightIndex, 
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

double DataJsmurfLR::computeSingleStatNull() const {
  if (!isComputed_[shortendIntervalLength_ - 1u]) {
    compute(shortendIntervalLength_);
  }
  
  double numerator = 0.0;
  for (unsigned int i = 0u; i < shortendIntervalLength_; i++) {
    numerator += sigmaInverseOne_[shortendIntervalLength_ - 1u][i] * data_[left_ + filterLength_ + i];
  }
  
  return numerator * numerator / denominator_[shortendIntervalLength_ - 1u] / 2.0;
}

double DataJsmurfLR::computeSingleStat(const double &value) const {
  if (!isComputed_[shortendIntervalLength_ - 1u]) {
    compute(shortendIntervalLength_);
  }
  
  double numerator = 0.0;
  for (unsigned int i = 0u; i < shortendIntervalLength_; i++) {
    numerator += sigmaInverseOne_[shortendIntervalLength_ - 1u][i] *
      (data_[left_ + filterLength_ + i] - value);
  }
  
  return numerator * numerator / denominator_[shortendIntervalLength_ - 1u] / 2.0;
}

void DataJsmurfLR::addLeft(const unsigned int &index) {
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

void DataJsmurfLR::addRight(const unsigned int &index) {
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

void DataJsmurfLR::add(Data * const data) {
  DataJsmurfLR* x = dynamic_cast<DataJsmurfLR *>(data);
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

void DataJsmurfLR::reset() {
  left_ = 0u;
  right_ = 0u;
  cumulatedSum_ = 0.0;
  shortendCumulatedSum_ = 0.0;
  intervalLength_ = 0u;
  shortendIntervalLength_ = 0u;
}
