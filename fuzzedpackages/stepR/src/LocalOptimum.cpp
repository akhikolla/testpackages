#include "LocalOptimum.h"

/***************
* class LocalOptimum
* maintains local solutions
* Florian Pein, 2015
***************/
LocalOptimum::LocalOptimum() : leftIndex_(NA_INTEGER), rightIndex_(NA_INTEGER), estimatedValue_(NA_REAL),
                               costs_(R_PosInf), lastSegment_(NULL) {}
                      
LocalOptimum::LocalOptimum(const unsigned int &leftIndex, const unsigned int &rightIndex,
                           const double &estimatedValue, const double &costs,
                           LocalOptimum * const lastSegment) : leftIndex_(leftIndex),
                           rightIndex_(rightIndex), estimatedValue_(estimatedValue), costs_(costs),
                           lastSegment_(lastSegment) {}

unsigned int LocalOptimum::leftIndex() const {
  return leftIndex_;
}

unsigned int LocalOptimum::rightIndex() const {
  return rightIndex_;
}

double LocalOptimum::estimatedValue() const {
  return estimatedValue_;
}

double LocalOptimum::costs() const {
  return costs_;
}

LocalOptimum * LocalOptimum::lastSegment() const {
  return lastSegment_;
}
                          
/*************
* update
* update local solution
****************/
void LocalOptimum::update(const LocalOptimum &newSolution) {
  if (newSolution.costs() < costs_) {
    leftIndex_ = newSolution.leftIndex();
    rightIndex_ = newSolution.rightIndex();
    estimatedValue_ = newSolution.estimatedValue();
    costs_ = newSolution.costs();
    lastSegment_ = newSolution.lastSegment();
  }
}
