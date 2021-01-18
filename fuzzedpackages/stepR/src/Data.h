#ifndef STEPR_H_DATA
#define STEPR_H_DATA

#include <Rcpp.h>
#include "SingleBounds.h"
#include "LocalOptimum.h"

using namespace Rcpp;

/***************
* class Data
* abstract class maintaining the data
* Florian Pein, 2015
***************/
class Data {
  protected:
    static NumericVector criticalValues_;
  
  public:
    virtual ~Data();
    
    // clean up of static variables
    virtual void cleanUpStaticVariables() = 0;
    
    // set criticalValues
    static void setCriticalValues(const List &input);
    
    // return the number of data points
    virtual unsigned int getN() const = 0;
    
    // returns a reference to a new object of the same class, delete has to be called later in the code
    virtual Data* newObject() const = 0;
  
    // computes the test statistic for a single interval under the null hypothesis
    virtual double computeSingleStatNull() const = 0;
    // computes the test statistic for a single interval
    virtual double computeSingleStat(const double &value) const = 0;    
    // computes bounds for a single interval
    virtual SingleBounds computeSingleBounds() const = 0;
    // computes the local optimum for a single interval for given bounds
    virtual LocalOptimum computeLocalOptimum(const unsigned int &leftIndex, const unsigned int &rightIndex, 
                                             const SingleBounds &localBounds,
                                             LocalOptimum * const lastSegment) const = 0;
    
    // update functions for cumulated sums etc.
    virtual void addLeft(const unsigned int &index) = 0;
    virtual void addRight(const unsigned int &index) = 0;
    virtual void add(Data * const data) = 0;
    virtual void reset() = 0;
    
    // required for local tests for given fits
    // clean and set of local variables and computes single stat for given start index and segments
    virtual void cleanUpLocalVariables();
    virtual void setLocal(const List &input);
    virtual double computeSingleStat(unsigned int startIndex, unsigned int leftSegment, unsigned int rightSegment) const;
};

#endif
