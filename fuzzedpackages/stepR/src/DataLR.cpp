#include "DataLR.h"

#include "choleskyDecomposition.h"
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>

#include <algorithm>

/***************
* class DataLR
* implements abstract class Data for the local likelihood ratio test for given fit
* Florian Pein, 2018
***************/
NumericVector DataLR::obs_;
NumericVector DataLR::obs0_;
NumericVector DataLR::value_;
NumericVector DataLR::covariances_;
unsigned int DataLR::filterLength_;

const char DataLR::uplo_ = 'U';
const char DataLR::trans_ = 'T';
const char DataLR::diag_ = 'N';
const int DataLR::incx_ = 1;

DataLR::DataLR() : len_(0u), Fleft_(NumericVector(0u)), Fright_(NumericVector(0u)),
                   v_(NumericVector(0u)), sol_(NumericVector(0u)), vtAv_(0.0) {}

void DataLR::cleanUpStaticVariables() {}

void DataLR::setData(const RObject &data, const List &input) {
  obs_ = input["obs"];
  obs0_ = input["obs0"];
  value_ = input["value"];
  covariances_ = input["covariances"];
  filterLength_ = input["filterLength"];
}

unsigned int DataLR::getN() const {
  stop("getN() is not implemented for parametric family LR");
}
                                                  
Data* DataLR::newObject() const {
  return new DataLR();
}             

SingleBounds DataLR::computeSingleBounds() const {
  stop("computeSingleBounds() is not implemented for parametric family LR");
}

LocalOptimum DataLR::computeLocalOptimum(const unsigned int &leftIndex,
                                                   const unsigned int &rightIndex, 
                                                   const SingleBounds &localBounds,
                                                   LocalOptimum * const lastSegment) const {
  stop("computeLocalOptimum(...) is not implemented for parametric family LR");
}

double DataLR::computeSingleStatNull() const {
  stop("computeSingleStatNull() is not implemented for parametric family LR");
}

double DataLR::computeSingleStat(const double &value) const {
  stop("computeSingleStat() is not implemented for parametric family LR");
}

void DataLR::addLeft(const unsigned int &index) {
  stop("addLeft() is not implemented for parametric family LR");
}

void DataLR::addRight(const unsigned int &index) {
  stop("addRight() is not implemented for parametric family LR");
}

void DataLR::add(Data * const data) {
  stop("add() is not implemented for parametric family LR");
}

void DataLR::reset() {
  stop("reset() is not implemented for parametric family LR");
}

void DataLR::cleanUpLocalVariables() {
  delete [] chol_;
}

void DataLR::setLocal(const List &input) {
  len_ = input["len"];
  Fleft_ = input["Fleft"];
  Fright_ = input["Fright"];
  v_ = input["v"];
  sol_ = input["sol"];
  vtAv_ = input["vtAv"];
  
  chol_ = choleskyDecomposition(len_ + filterLength_ - 1u, covariances_);
}

double DataLR::computeSingleStat(unsigned int startIndex, unsigned int leftSegment, unsigned int rightSegment) const {
  NumericVector obs = NumericVector(len_ + filterLength_ - 1u);
  double* obsEst = new double[len_ + filterLength_ - 1u];
  double* obs0 = new double[len_ + filterLength_ - 1u];
  
  const int size = len_ + filterLength_ - 1;
  const int bands = std::min<int>(covariances_.size() - 1, len_ + filterLength_ - 2);
  const int ldA = bands + 1;
  
  unsigned int indexObs = startIndex + 1u;

  for (unsigned int i = 0u; i < len_ + filterLength_ - 1u; ++i) {
    obs[i] = obs_[indexObs] - value_[leftSegment] * Fleft_[i] - value_[rightSegment] * Fright_[i];
    obs0[i] = obs0_[indexObs];
    indexObs++;
  }
  
  double estMean = 0.0;
  for (unsigned int i = 0u; i < len_ + filterLength_ - 1u; ++i) {
    estMean += obs[i] * sol_[i];
  }
  estMean /= vtAv_;
  
  for (unsigned int i = 0u; i < len_ + filterLength_ - 1u; ++i) {
    obsEst[i] = obs[i] - v_[i] * estMean;
  }
  
  F77_CALL(dtbsv)(&uplo_, &trans_, &diag_, &size, &bands, chol_, &ldA, obs0, &incx_);
  F77_CALL(dtbsv)(&uplo_, &trans_, &diag_, &size, &bands, chol_, &ldA, obsEst, &incx_);
  
  double stat = 0.0;
  for (unsigned int i = 0u; i < len_ + filterLength_ - 1u; ++i) {
    stat += obs0[i] * obs0[i] - obsEst[i] * obsEst[i];
  }
  
  delete [] obs0;
  delete [] obsEst;
  
  return stat;
}
