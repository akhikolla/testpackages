#ifndef STEPR_H_INTERVALSYSTEMDYAPARLENGTHS
#define STEPR_H_INTERVALSYSTEMDYAPARLENGTHS

#include "IntervalSystem.h"
#include <vector>

/***************
* class IntervalSystemDyaParLengths
* implements the abstract class IntervalSystem for the interval system dyadic partition with given lengths
* Florian Pein, 2015
***************/
class IntervalSystemDyaParLengths : public IntervalSystem {
  private:
    std::vector<bool> contained_;
  
  public:
    IntervalSystemDyaParLengths(const unsigned int &n, const List &input);
    
    bool isInIntervalSystem(const unsigned int &leftIndex, const unsigned int &rightIndex) const;
    
    NumericVector computeMultiscaleStatisticNull(Data * const data);
    List computeBounds(Data* data);
    
    template <typename T>
    void dynamicProgramDyaParLengths(Data * const data, T t) {
      unsigned int n = data -> getN();
      std::vector<Data*> dataVector;
      dataVector.reserve(n);
      
      if (contained_[0]) {
        for (unsigned int i = 0u; i < n; ++i) {
          dataVector.push_back(data -> newObject());
          dataVector[i] -> addRight(i);
          t.compute(dataVector[i], i, i);
        }
      } else {
        for (unsigned int i = 0u; i < n; ++i) {
          dataVector.push_back(data -> newObject());
          dataVector[i] -> addRight(i);
        }
      }
      
      unsigned int lenDividedBy2 = 1u;
      for (unsigned int len = 2u; len <= n; len *= 2u) {
        checkUserInterrupt();
        if (contained_[len - 1u]) {
          for (unsigned int i = 0u, j = len - 1u; j < n; i += len, j += len) {
            dataVector[i] -> add(dataVector[i + lenDividedBy2]);
            t.compute(dataVector[i], i, j);
          }
        } else {
          for (unsigned int i = 0u, j = len - 1u; j < n; i += len, j += len) {
            dataVector[i] -> add(dataVector[i + lenDividedBy2]);
          }
        }

        lenDividedBy2 = len;
      }
      
      for (unsigned int i = 0u; i < n; ++i) {
        delete dataVector[i];
      }     
    }
};

#endif
