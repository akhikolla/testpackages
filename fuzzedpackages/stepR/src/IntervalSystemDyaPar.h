#ifndef STEPR_H_INTERVALSYSTEMDYAPAR
#define STEPR_H_INTERVALSYSTEMDYAPAR

#include "IntervalSystem.h"
#include <vector>

/***************
* class IntervalSystemDyaPar
* implements the abstract class IntervalSystem for the interval system dyadic partition
* Florian Pein, 2015
***************/
class IntervalSystemDyaPar : public IntervalSystem {
  private:
    std::vector<bool> contained_;
  
  public:
    IntervalSystemDyaPar(const unsigned int &n);
    
    bool isInIntervalSystem(const unsigned int &leftIndex, const unsigned int &rightIndex) const;
    
    NumericVector computeMultiscaleStatisticNull(Data * const data); 
    List computeBounds(Data * const data);
    
    template <typename T>
    void dynamicProgramDyaPar(Data * const data, T t) {
      unsigned int n = data -> getN();
      std::vector<Data*> dataVector;
      dataVector.reserve(n);
      
      for (unsigned int i = 0u; i < n; ++i) {
        dataVector.push_back(data -> newObject());
        dataVector[i] -> addRight(i);
        t.compute(dataVector[i], i, i);
      }
      
      unsigned int lenDividedBy2 = 1u;
      for (unsigned int len = 2u; len <= n; len *= 2u) {
        checkUserInterrupt();
        for (unsigned int i = 0u, j = len - 1u; j < n; i += len, j += len) {
          dataVector[i] -> add(dataVector[i + lenDividedBy2]);
          t.compute(dataVector[i], i, j);
        }
        lenDividedBy2 = len;
      }
      
      for (unsigned int i = 0u; i < n; ++i) {
        delete dataVector[i];
      }     
    }
};

#endif
