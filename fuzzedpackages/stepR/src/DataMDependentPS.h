#ifndef STEPR_H_DATAMDEPENDENTPS
#define STEPR_H_DATAMDEPENDENTPS

#include "Data.h"
#include <vector>
  
/***************
* class DataMDependentPS
* implements abstract class Data for m-dependent homogeneous observations
* using a partial sum test statistic
* Florian Pein, 2015
***************/
class DataMDependentPS : public Data {
  private:
    static NumericVector data_;
    static std::vector<double> varianceSum_;
    
    double cumulatedSum_;            // sum of data from current left to current right index
    unsigned int intervalLength_;    // current interval length
                                     // i.e. current right index - current left index + 1
    
  public:
    DataMDependentPS();
    
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
