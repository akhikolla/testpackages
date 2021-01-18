#ifndef STEPR_H_INTERVALSYSTEMALL
#define STEPR_H_INTERVALSYSTEMALL

#include "IntervalSystem.h"

/***************
* class IntervalSystemAll
* implements the abstract class IntervalSystem for the interval system of all intervals
* Florian Pein, 2015
***************/
class IntervalSystemAll : public IntervalSystem {
  public:
    IntervalSystemAll(const unsigned int &n);
    
    bool isInIntervalSystem(const unsigned int &leftIndex,
                            const unsigned int &rightIndex) const;
    NumericVector computeMultiscaleStatisticNull(Data * const data); 
    NumericVector computeMultiscaleStatistic(Data * const data, const List &input);
    List computeBounds(Data * const data);
    
    template <typename T>
    void dynamicProgramAll(Data * const data, T t) {
      for (unsigned int j = 0u; j < data -> getN(); ++j) {
        checkUserInterrupt();
        data -> reset();
    
        unsigned int i = j + 1u;
        while (i > 0u) {
          --i;
          data -> addLeft(i);
          t.compute(data, i, j);
        }
      }
    }
};

#endif
