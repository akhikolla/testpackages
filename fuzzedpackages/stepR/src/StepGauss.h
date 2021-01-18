#include "Step.h"

/***************
* class StepGauss
* implements virtual class Step for Gaussian data, i.e. (weighted) l2-loss
* Thomas Hotz, 2007-2011
***************/

class StepGauss : public Step {
  public:
    double* cs; // cumulative sums
    double* css; // cumulative sums of squares
    double* csv; // cumulative sums of variances, i.e. the expected css
                 // equivalently the number of data points leading to the corresponding css if the variance is assumed to be constant
                 // can also be viewed as cumulative weights
    
    // initialization
    StepGauss(unsigned int n, double* xcs, double* xcss, double* xcsv); // constructor for n data points, with pointers to arrays of that length
    StepGauss(unsigned int n, double* xcs, double* xcss, double* xcsv, double* xlb, double* xub); // constructor for n data points and bounds
    
    // needed for all distributions
    virtual double cost(unsigned int startIndex, unsigned int endIndex) const; // calculate cost of a block
    virtual double costBound(unsigned int startIndex, unsigned int endIndex, const LUBound& bound) const; // calculate cost of a block given bounds
    virtual double estBound(unsigned int startIndex, unsigned int endIndex, const LUBound& bound) const; // corresponding estimate
};
