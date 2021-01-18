#ifndef STEPR_H_DATA2PARAM
#define STEPR_H_DATA2PARAM

#include "Data.h"
  
/***************
* class Data2Param
* implements abstract class Data for the local two parameter test for given fit
* Florian Pein, 2018
***************/
class Data2Param : public Data {
  private:
    static NumericVector obs_;
    static NumericVector T0_;
    static NumericVector Tobs_;
    static NumericVector value_;
    static NumericVector var_;
    static unsigned int filterLength_;
    
    unsigned int len_;
    NumericVector Fleft_;
    NumericVector Fright_;
    NumericVector v_;
    double sumv2_;
    double sumSigmaL_;
    double sumSigmaR_;
    double sumW_;
    NumericVector w_;
    NumericVector sigmaL_;
    NumericVector sigmaR_;
    
  public:
    Data2Param();
    
    // clean up of static variables
    void cleanUpStaticVariables();
    
    // set data
    static void setData(const RObject &data, const List &input);
    
    // return the number of data points
    unsigned int getN() const;
    
    // returns a reference to a new object of the same class, delete has to be called later in the code
    Data* newObject() const;
    
    // computes the bounds for a single interval
    SingleBounds computeSingleBounds() const;
    // computes the local optimum for a single interval for given bounds
    LocalOptimum computeLocalOptimum(const unsigned int &leftIndex, const unsigned int &rightIndex, 
                                     const SingleBounds &localBounds, LocalOptimum * const lastSegment) const;
    // computes the test statistic for a single interval under the null hypothesis
    double computeSingleStatNull() const;
    // computes the test statistic for a single interval
    double computeSingleStat(const double &value) const;
    
    // update functions for cumulated sums etc.
    void addLeft(const unsigned int &index);
    void addRight(const unsigned int &index);
    void add(Data * const data);    
    void reset();
    
    // required for local tests for given fits
    // clean and set of local variables and computes single stat for given start index and segments
    void cleanUpLocalVariables();
    void setLocal(const List &input);
    double computeSingleStat(unsigned int startIndex, unsigned int leftSegment, unsigned int rightSegment) const;
};

#endif
