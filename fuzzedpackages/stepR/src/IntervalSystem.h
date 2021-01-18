#ifndef STEPR_H_INTERVALSYSTEM
#define STEPR_H_INTERVALSYSTEM

#include <Rcpp.h>
#include "Data.h"

using namespace Rcpp;

/***************
* class IntervalSystem
* abstract class maintaining the interval system
* Florian Pein, 2015
***************/
class IntervalSystem {
  protected:
    unsigned int numberOfIntervals_;
  
  public:
    IntervalSystem(const unsigned int &numberOfIntervals);
    virtual ~IntervalSystem();
    
    virtual bool isInIntervalSystem(const unsigned int &leftIndex, const unsigned int &rightIndex) const = 0;
    
    virtual NumericVector computeMultiscaleStatisticNull(Data * const data);
    virtual NumericVector computeMultiscaleStatistic(Data * const data, const List &input);
    virtual List computeBounds(Data * const data);
    
    template <typename T>
    void dynamicProgram(Data * const data, T t) {
      for (unsigned int j = 0u; j < data -> getN(); ++j) {
        checkUserInterrupt();
        data -> reset();
    
        unsigned int i = j + 1u;
        while (i > 0u) {
          --i;
          data -> addLeft(i);
          if (isInIntervalSystem(i, j)) {
            t.compute(data, i, j);
          }
        }
      }
    }
};

#endif
