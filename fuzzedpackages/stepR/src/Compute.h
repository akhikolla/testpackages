#ifndef STEPR_H_COMPUTE
#define STEPR_H_COMPUTE

#include <Rcpp.h>
#include "Data.h"

using namespace Rcpp;

class ComputeStatNull {
  private:
    NumericVector stat_;
  
  public:
    const NumericVector& stat() const;
    
    ComputeStatNull(const unsigned int &n);
    void compute(Data * const data, const unsigned int &li, const unsigned int &ri);  
};

class ComputeStat {
  private:
    NumericVector stat_;
  
  public:
    const NumericVector& stat() const;
    
    ComputeStat(const unsigned int &n);
    void compute(Data * const data, const unsigned int &li, const unsigned int &ri, const double &value);  
};

class ComputeBounds {
  private:
    IntegerVector leftIndex_;
    IntegerVector rightIndex_;
    NumericVector lowerBound_;
    NumericVector upperBound_;
    unsigned int index_;
    
  public:
    const IntegerVector& leftIndex() const;
    const IntegerVector& rightIndex() const;
    const NumericVector& lowerBound() const;
    const NumericVector& upperBound() const;
    ComputeBounds(const unsigned int &numberOfIntervals);
    void compute(Data * const data, const unsigned int &li, const unsigned int &ri);  
};

#endif
