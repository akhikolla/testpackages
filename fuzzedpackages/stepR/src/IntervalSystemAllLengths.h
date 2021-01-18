#ifndef STEPR_H_INTERVALSYSTEMALLLENGTHS
#define STEPR_H_INTERVALSYSTEMALLLENGTHS

#include "IntervalSystem.h"
#include <vector>

/***************
* class IntervalSystemAllLengths
* implements the abstract class IntervalSystem for the interval system of all intervals with given lengths
* Florian Pein, 2015
***************/
class IntervalSystemAllLengths : public IntervalSystem {
  private:
    std::vector<bool> contained_;
  
  public:
    IntervalSystemAllLengths(const unsigned int &n, const List &input);
    
    bool isInIntervalSystem(const unsigned int &leftIndex, const unsigned int &rightIndex) const;
    
    NumericVector computeMultiscaleStatisticNull(Data * const data); 
    List computeBounds(Data * const data);
    
    template <typename T>
    void dynamicProgramAllLengths(Data * const data, T t) {
      for (unsigned int j = 0u; j < data -> getN(); ++j) {
        checkUserInterrupt();
        data -> reset();
    
        unsigned int i = j + 1u;
        while (i > 0u) {
          --i;
          data -> addLeft(i);
          if (contained_[j - i]) {
            t.compute(data, i, j);
          }
        }
      }
    }
};

#endif
