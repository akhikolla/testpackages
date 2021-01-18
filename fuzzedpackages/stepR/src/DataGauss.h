#ifndef STEPR_H_DATAGAUSS
#define STEPR_H_DATAGAUSS

#include "Data.h"

/***************
* class DataGauss
* implements abstract class Data for independent, homogeneous gaussian observations
* Florian Pein, 2015
***************/
class DataGauss : public Data {
  private:
    static NumericVector data_;
    static double standardDeviation_; // the global standard deviation
    
    double cumulatedSum_;             // sum of data from current left to current right index
    unsigned int intervalLength_;     // current interval length
                                      // i.e. current right index - current left index + 1
    
  public:
    DataGauss();
    
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
};

#endif
