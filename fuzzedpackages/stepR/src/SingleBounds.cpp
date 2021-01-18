#include "SingleBounds.h"
#include <algorithm>

/***************
* class SingleBounds
* maintains lower and upper bound
* Florian Pein, 2015
***************/

/*************
* constructor with default bounds set to -+Inf
****************/
SingleBounds::SingleBounds() : lower_(R_NegInf), upper_(R_PosInf) {}

/*************
* constructor given lower and upper bound
****************/
SingleBounds::SingleBounds(const double &lower, const double &upper) : lower_(lower), upper_(upper) {}

// get functions
double SingleBounds::lower() const {
  return lower_;
}
double SingleBounds::upper() const {
  return upper_;
}

/*************
* feasible
* check feasibility of bounds
****************/
bool SingleBounds::feasible() const {
  return lower_ <= upper_;
}

/*************
* intersect
* insects with newBounds 
****************/
void SingleBounds::intersect(const SingleBounds &newBounds) {
  lower_ = std::max<double>(lower_, newBounds.lower());
  upper_ = std::min<double>(upper_, newBounds.upper());
}

/*************
* reset
* resets the bounds 
****************/
void SingleBounds::reset() {
  lower_ = R_NegInf;
  upper_ = R_PosInf;
}
