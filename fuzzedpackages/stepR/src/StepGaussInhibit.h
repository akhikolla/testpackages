#include "StepGauss.h"

/***************
* class StepGaussInhibit
* extends class StepGauss to inhibit close jumps
* Thomas Hotz, 2007-2010
***************/

class StepGaussInhibit : public StepGauss {
  public:
    int istart; // the minimal distance to the last jump before the data (at -1)
    int imiddle; // the minimal distance between two jumps
    int iend; // the minimal distance to the last jump (at n)
    
    // initialization
    StepGaussInhibit(unsigned int n, double* xcs, double* xcss, double* xcsv, int inhibitStart, int inhibitMiddle, int inhibitEnd); // constructor
    
    // needed for all distributions
    virtual double cost(unsigned int startIndex, unsigned int endIndex) const; // calculate cost of a block, infinite if jumps get too close
};
