#ifndef STEPR_H_LOCALOPTIMUM
#define STEPR_H_LOCALOPTIMUM

#include <Rcpp.h>

using namespace Rcpp;

/***************
* class LocalOptimum
* maintains local solutions
* Florian Pein, 2015
***************/
class LocalOptimum {
  private:
    unsigned int leftIndex_;
    unsigned int rightIndex_;
    double estimatedValue_;
    double costs_;
    LocalOptimum * lastSegment_;
    
    
  public:
    // constructors for xn data points
    LocalOptimum();
    LocalOptimum(const unsigned int &leftIndex, const unsigned int &rightIndex, const double &estimatedValue,
                 const double &costs, LocalOptimum * const lastSegment);
    
    // get functions
    unsigned int leftIndex() const;
    unsigned int rightIndex() const;
    double estimatedValue() const;
    double costs() const;
    LocalOptimum * lastSegment() const;
    
    void update(const LocalOptimum &newSolution);
};

#endif
