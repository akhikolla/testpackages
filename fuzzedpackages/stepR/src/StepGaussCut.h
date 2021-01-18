#include "StepGauss.h"

/***************
* class StepGaussCut
* like class StepGaussInhibit but cut close jumps
* Thomas Hotz, 2007-2010
***************/

class StepGaussCut : public Step {
  public:
    int cbefore; // how much to cut out to the left
    int cafter; // how much to cut out to the right
    double* bcs; // before cumulative sums
    double* bcss; // before cumulative sums of squares
    double* bcsv; // before cumulative sums of variances, i.e. the expected css
    double* acs; // after cumulative sums
    double* acss; // after cumulative sums of squares
    double* acsv; // after cumulative sums of variances, i.e. the expected css
    
    // initialization
    StepGaussCut(unsigned int n, double* xbcs, double* xbcss, double* xbcsv, double* xacs, double* xacss, double* xacsv, int cutBefore, int cutAfter); // constructor
    
    // needed for all distributions
    virtual double cost(unsigned int startIndex, unsigned int endIndex) const; // calculate cost of a block, infinite if jumps get too close
};
