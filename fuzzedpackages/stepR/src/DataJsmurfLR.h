#ifndef STEPR_H_DATAJSMURFLR
#define STEPR_H_DATAJSMURFLR

#include "Data.h"

/***************
* class DataJsmurf
* implements abstract class Data for filtered homogeneous gaussian observations using JSMURF idea
* and the likelihood ratio test
* Florian Pein, 2017
***************/
class DataJsmurfLR : public Data {
  private:
    static NumericVector data_;
    static unsigned int filterLength_;
    static int m_;
    static NumericVector covariances_;
    
    static std::vector<bool> isComputed_;
    static std::vector<double*> sigmaInverseOne_;
    static std::vector<double> denominator_;
    
    static const char uplo_;
    static const char trans1_;
    static const char trans2_;
    static const char diag_;
    static const int incx_;
    
    static void compute(const int &size);
    
    unsigned int left_;  // current left index
    unsigned int right_; // current right index
    
    double cumulatedSum_;             // sum of data from current left to current right index
    double shortendCumulatedSum_;     // sum of data from current left plus filterLength_ to current right index
    unsigned int intervalLength_;     // current interval length
                                      // i.e. current right index - current left index + 1
    unsigned int shortendIntervalLength_;                                  
    
  public:
    DataJsmurfLR();
    
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
