#include "Data2Param.h"

/***************
* class Data2Param
* implements abstract class Data for the local two parameter test for given fit
* Florian Pein, 2018
***************/
NumericVector Data2Param::obs_;
NumericVector Data2Param::T0_;
NumericVector Data2Param::Tobs_;
NumericVector Data2Param::value_;
NumericVector Data2Param::var_;
unsigned int Data2Param::filterLength_;

Data2Param::Data2Param() : len_(0u), Fleft_(NumericVector(0u)), Fright_(NumericVector(0u)),
                           v_(NumericVector(0u)), sumv2_(0.0), sumSigmaL_(0.0),
                           sumSigmaR_(0.0), sumW_(0.0), w_(NumericVector(0u)),
                           sigmaL_(NumericVector(0u)), sigmaR_(NumericVector(0u)) {}

void Data2Param::cleanUpStaticVariables() {}

void Data2Param::setData(const RObject &data, const List &input) {
  obs_ = input["obs"];
  T0_ = input["T0"];
  Tobs_ = input["Tobs"];
  value_ = input["value"];
  var_ = input["var"];
  filterLength_ = input["filterLength"];
}

unsigned int Data2Param::getN() const {
  stop("getN() is not implemented for parametric family 2Param");
}
                                                  
Data* Data2Param::newObject() const {
  return new Data2Param();
}             

SingleBounds Data2Param::computeSingleBounds() const {
  stop("computeSingleBounds() is not implemented for parametric family 2Param");
}

LocalOptimum Data2Param::computeLocalOptimum(const unsigned int &leftIndex,
                                                   const unsigned int &rightIndex, 
                                                   const SingleBounds &localBounds,
                                                   LocalOptimum * const lastSegment) const {
  stop("computeLocalOptimum(...) is not implemented for parametric family 2Param");
}

double Data2Param::computeSingleStatNull() const {
  stop("computeSingleStatNull() is not implemented for parametric family 2Param");
}

double Data2Param::computeSingleStat(const double &value) const {
  stop("computeSingleStat() is not implemented for parametric family 2Param");
}

void Data2Param::addLeft(const unsigned int &index) {
  stop("addLeft() is not implemented for parametric family 2Param");
}

void Data2Param::addRight(const unsigned int &index) {
  stop("addRight() is not implemented for parametric family 2Param");
}

void Data2Param::add(Data * const data) {
  stop("add() is not implemented for parametric family 2Param");
}

void Data2Param::reset() {
  stop("reset() is not implemented for parametric family 2Param");
}

void Data2Param::cleanUpLocalVariables() {}

void Data2Param::setLocal(const List &input) {
  len_ = input["len"];
  Fleft_ = input["Fleft"];
  Fright_ = input["Fright"];
  v_ = input["v"];
  sumv2_ = input["sumv2"];
  sumSigmaL_ = input["sumSigmaL"];
  sumSigmaR_ = input["sumSigmaR"];
  sumW_ = input["sumW"];
  w_ = input["w"];
  sigmaL_ = input["sigmaL"];
  sigmaR_ = input["sigmaR"];
}

double Data2Param::computeSingleStat(unsigned int startIndex, unsigned int leftSegment, unsigned int rightSegment) const {
  NumericVector obs = NumericVector(len_ + filterLength_ - 1u);
  unsigned int indexObs = startIndex + 1u;
  for (unsigned int i = 0u; i < len_ + filterLength_ - 1u; ++i) {
    obs[i] = obs_[indexObs] - value_[leftSegment] * Fleft_[i] - value_[rightSegment] * Fright_[i];
    indexObs++;
  }
  
  double estMean = 0.0;
  for (unsigned int i = 0u; i < len_ + filterLength_ - 1u; ++i) {
    estMean += obs[i] * v_[i];
  }
  estMean /= sumv2_;
  
  double estSigma2 = 0.0;
  for (unsigned int i = 0u; i < len_ + filterLength_ - 1u; ++i) {
    estSigma2 += w_[i] * (obs[i] - v_[i] * estMean) * (obs[i] - v_[i] * estMean);
  }
  estSigma2 = (estSigma2 - sumSigmaL_ * var_[leftSegment] - sumSigmaR_ * var_[rightSegment]) / sumW_;
  
  if (estSigma2 < 0.0) {
    // return R_NegInf;
    estSigma2 = 0.0;
  }
  
  double stat = 0.0;
  indexObs = startIndex + 1u;
  for (unsigned int i = 0u; i < len_ + filterLength_ - 1u; ++i) {
    double Tvar = w_[i] * estSigma2 + sigmaL_[i] * var_[leftSegment] + sigmaR_[i] * var_[rightSegment];
    double Test = (obs[i] - v_[i] * estMean) * (obs[i] - v_[i] * estMean);
    
    stat += log(T0_[indexObs] / Tvar) + Tobs_[indexObs] / T0_[indexObs] - Test / Tvar;
    indexObs++; 
  }
  
  return stat;
}
