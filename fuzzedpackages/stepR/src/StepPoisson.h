#include "Step.h"

/***************
* class StepPoisson
* implements virtual class Step for Poisson data
* Thomas Hotz, 2007-2011
***************/

class StepPoisson : public Step {
  public:
    int* cs; // cumulative sums
    double* cw; // cumulative weights
                 // equivalently the number of data points leading to the corresponding cs if the wieghts are constant
    
    // initialization
    StepPoisson(unsigned int n, int* xcs, double* xcw); // constructor for n data points, with pointers to arrays of that length
    StepPoisson(unsigned int n, int* xcs, double* xcw, double* xlb, double* xub); // constructor for n data points and bounds
    
    // needed for all distributions
    virtual double cost(unsigned int startIndex, unsigned int endIndex) const; // calculate cost of a block
    virtual double costBound(unsigned int startIndex, unsigned int endIndex, const LUBound& bound) const; // calculate cost of a block given bounds
    virtual double estBound(unsigned int startIndex, unsigned int endIndex, const LUBound& bound) const; // corresponding estimate
};
