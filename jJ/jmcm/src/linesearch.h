// Copyright (C) 2015-2016 The University of Manchester
//
// Written by Yi Pan - ypan1988@gmail.com

#ifndef JMCM_LINESEARCH_H_
#define JMCM_LINESEARCH_H_

#include <algorithm>
#include <cmath>
#include <limits>

#include <RcppArmadillo.h>

namespace pan {

template <typename T>
class LineSearch {
 public:
  LineSearch();   // Constructor
  ~LineSearch();  // Destructor

  void GetStep(T &func, arma::vec &x, arma::vec &p, const double kStepMax);

  void set_message(bool message) { message_ = message; }

 protected:
  bool message_;
  bool IsInfOrNaN(double x);
};  // class LineSearch

#include "linesearch_impl.h"

}  // namespace pan

#endif  // JMCM_LINESEARCH_H_
