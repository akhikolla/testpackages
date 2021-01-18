#include "DataHjsmurfLR.h"

#include <algorithm>
#include <cmath>

#include "choleskyDecomposition.h"
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>

/***************
* class DataHjsmurfLR
* implements abstract class Data for filtered heterogeneous gaussian observations using JSMURF idea
* and likelihood ratio test
* Florian Pein, 2017
***************/
NumericVector DataHjsmurfLR::data_;
unsigned int DataHjsmurfLR::filterLength_;
int DataHjsmurfLR::m_ = 0;
NumericVector DataHjsmurfLR::correlations_;
std::vector<bool> DataHjsmurfLR::isComputed_;
std::vector<double*> DataHjsmurfLR::sigmaInverseOne_;
std::vector<double*> DataHjsmurfLR::cholesky_;
std::vector<double> DataHjsmurfLR::oneSigmaInverseOne_;

const char DataHjsmurfLR::uplo_ = 'U';
const char DataHjsmurfLR::trans1_ = 'T';
const char DataHjsmurfLR::trans2_ = 'N';
const char DataHjsmurfLR::diag_ = 'N';
const int DataHjsmurfLR::incx_ = 1;

void DataHjsmurfLR::compute(const int &size) {
  const int bands = std::min<int>(DataHjsmurfLR::m_ - 1, size - 1);
  
  const int ldA = bands + 1;

  double* A = choleskyDecomposition(size, correlations_);
  cholesky_[size - 1u] = A;
  
  double* x;
  x = new double[size];
  
  std::fill(x, x + size, 1.0);
  
  F77_CALL(dtbsv)(&uplo_, &trans1_, &diag_, &size, &bands, A, &ldA, x, &incx_);
  F77_CALL(dtbsv)(&uplo_, &trans2_, &diag_, &size, &bands, A, &ldA, x, &incx_);
  
  sigmaInverseOne_[size - 1u] = x;
  for (int i = 0; i < size; i++) {
    oneSigmaInverseOne_[size - 1] += x[i];
  }
  
  isComputed_[size - 1] = true;
}

DataHjsmurfLR::DataHjsmurfLR() : cumulatedSum_(0.0), shortendCumulatedSum_(0.0),
                                 cumulatedSumSquared_(0.0), shortendCumulatedSumSquared_(0.0),
                                 intervalLength_(0u), shortendIntervalLength_(0u) {}

void DataHjsmurfLR::cleanUpStaticVariables() {
  for (unsigned int i = 0u; i < data_.size(); ++i) {
    if (isComputed_[i]) {
      delete[] sigmaInverseOne_[i];
      delete[] cholesky_[i];
    }
  }
  
  std::vector<bool> tmpIsComputed;
  isComputed_.swap(tmpIsComputed);
  
  std::vector<double*> tmpSigmaInverseOne;
  sigmaInverseOne_.swap(tmpSigmaInverseOne);
  
  std::vector<double*> tmpCholesky;
  cholesky_.swap(tmpCholesky);
  
  std::vector<double> tmpOneSigmaInverseOne;
  oneSigmaInverseOne_.swap(tmpOneSigmaInverseOne);
}

void DataHjsmurfLR::setData(const RObject &data, const List &input) {
  data_ = data;
  filterLength_ = input["filterLength"];
  correlations_ = input["correlations"];
  m_ = correlations_.size();
  isComputed_.resize(data_.size(), false);
  sigmaInverseOne_.resize(data_.size());
  cholesky_.resize(data_.size());
  oneSigmaInverseOne_.resize(data_.size(), 0.0);
}

unsigned int DataHjsmurfLR::getN() const {
  return data_.size();
}
                                                  
Data* DataHjsmurfLR::newObject() const {
  return new DataHjsmurfLR();
}

SingleBounds DataHjsmurfLR::computeSingleBounds() const {
  if (!isComputed_[shortendIntervalLength_ - 1u]) {
    compute(shortendIntervalLength_);
  }
  
  double* x;
  x = new double[shortendIntervalLength_];
  for (unsigned int i = 0; i < shortendIntervalLength_; i++) {
    x[i] = data_[left_ + filterLength_ + i];
  }
  
  const int size = shortendIntervalLength_;
  const int bands = std::min<int>(m_ - 1, size - 1);
  const int ldA = bands + 1;
  
  F77_CALL(dtbsv)(&uplo_, &trans1_, &diag_, &size, &bands,
           cholesky_[shortendIntervalLength_ - 1u], &ldA, x, &incx_);
  
  double ySigmaInverseY = 0;
  double oneSigmaInverseY = 0;
  
  for (unsigned int i = 0; i < shortendIntervalLength_; i++) {
    ySigmaInverseY += x[i] * x[i];
    oneSigmaInverseY += sigmaInverseOne_[shortendIntervalLength_ - 1u][i] * data_[left_ + filterLength_ + i];
  }
  
  delete[] x;
  
  double mean = shortendCumulatedSum_ / shortendIntervalLength_;
  double help = std::sqrt(oneSigmaInverseY * oneSigmaInverseY -
                          oneSigmaInverseOne_[shortendIntervalLength_ - 1u] * 
                          (ySigmaInverseY - 2.0 * criticalValues_[intervalLength_ - 1u] *
                          (ySigmaInverseY - 2.0 * mean * oneSigmaInverseY +
                          mean * mean * oneSigmaInverseOne_[shortendIntervalLength_ - 1u])));
  return SingleBounds((oneSigmaInverseY - help) / oneSigmaInverseOne_[shortendIntervalLength_ - 1u],
                      (oneSigmaInverseY + help) / oneSigmaInverseOne_[shortendIntervalLength_ - 1u]);
}

LocalOptimum DataHjsmurfLR::computeLocalOptimum(const unsigned int &leftIndex, const unsigned int &rightIndex, 
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

double DataHjsmurfLR::computeSingleStatNull() const {
  if (!isComputed_[shortendIntervalLength_ - 1u]) {
    compute(shortendIntervalLength_);
  }
  
  double* x;
  x = new double[shortendIntervalLength_];
  for (unsigned int i = 0; i < shortendIntervalLength_; i++) {
    x[i] = data_[left_ + filterLength_ + i];
  }
  
  const int size = shortendIntervalLength_;
  const int bands = std::min<int>(m_ - 1, size - 1);
  const int ldA = bands + 1;
  
  F77_CALL(dtbsv)(&uplo_, &trans1_, &diag_, &size, &bands,
           cholesky_[shortendIntervalLength_ - 1u], &ldA, x, &incx_);
  
  double ySigmaInverseY = 0;
  double oneSigmaInverseY = 0;
  
  for (unsigned int i = 0; i < shortendIntervalLength_; i++) {
    ySigmaInverseY += x[i] * x[i];
    oneSigmaInverseY += sigmaInverseOne_[shortendIntervalLength_ - 1u][i] * data_[left_ + filterLength_ + i];
  }
  
  delete[] x;
  
  double mean = shortendCumulatedSum_ / shortendIntervalLength_;
  return ySigmaInverseY / 2.0 / (ySigmaInverseY - 2.0 * mean * oneSigmaInverseY +
                                 mean * mean * oneSigmaInverseOne_[shortendIntervalLength_ - 1u]);
}

double DataHjsmurfLR::computeSingleStat(const double &value) const {
  if (!isComputed_[shortendIntervalLength_ - 1u]) {
    compute(shortendIntervalLength_);
  }
  
  double* x;
  x = new double[shortendIntervalLength_];
  for (unsigned int i = 0; i < shortendIntervalLength_; i++) {
    x[i] = data_[left_ + filterLength_ + i];
  }
  
  const int size = shortendIntervalLength_;
  const int bands = std::min<int>(m_ - 1, size - 1);
  const int ldA = bands + 1;
  
  F77_CALL(dtbsv)(&uplo_, &trans1_, &diag_, &size, &bands,
           cholesky_[shortendIntervalLength_ - 1u], &ldA, x, &incx_);
  
  double ySigmaInverseY = 0;
  double oneSigmaInverseY = 0;
  
  for (unsigned int i = 0; i < shortendIntervalLength_; i++) {
    ySigmaInverseY += x[i] * x[i];
    oneSigmaInverseY += sigmaInverseOne_[shortendIntervalLength_ - 1u][i] * data_[left_ + filterLength_ + i];
  }
  
  delete[] x;

  double mean = shortendCumulatedSum_ / shortendIntervalLength_;
  return (ySigmaInverseY - 2.0 * value * oneSigmaInverseY + 
          value * value * oneSigmaInverseOne_[shortendIntervalLength_ - 1u]) / 2.0 /
            (ySigmaInverseY - 2.0 * mean * oneSigmaInverseY +
              mean * mean * oneSigmaInverseOne_[shortendIntervalLength_ - 1u]);
}

void DataHjsmurfLR::addLeft(const unsigned int &index) {
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

void DataHjsmurfLR::addRight(const unsigned int &index) {
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

void DataHjsmurfLR::add(Data * const data) {
  DataHjsmurfLR* x = dynamic_cast<DataHjsmurfLR *>(data);
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

void DataHjsmurfLR::reset() {
  left_ = 0u;
  right_ = 0u;
  cumulatedSum_ = 0.0;
  shortendCumulatedSum_ = 0.0;
  cumulatedSumSquared_ = 0.0;
  shortendCumulatedSumSquared_ = 0.0;
  intervalLength_ = 0u;
  shortendIntervalLength_ = 0u;
}
