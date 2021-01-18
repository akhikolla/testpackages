#ifndef STEPR_H_SINGLEBOUNDS
#define STEPR_H_SINGLEBOUNDS

#include <Rcpp.h>

using namespace Rcpp;

/***************
* class SingleBounds
* maintains lower and upper bound
* Florian Pein, 2015
***************/
class SingleBounds {
  private:
    double lower_; 
    double upper_;      

  public:  
    // constructors for xn data points
    SingleBounds();
    SingleBounds(const double &lower, const double &upper);

    double lower() const;
    double upper() const;
    bool feasible() const;
    void intersect(const SingleBounds &newBounds);
    void reset();
};

#endif
