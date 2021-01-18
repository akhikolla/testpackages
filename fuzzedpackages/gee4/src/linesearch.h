/* Copyright 2016 The University of Manchester.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   Written by Yi Pan - ypan1988@gmail.com
   ==============================================================================*/

#ifndef LINESEARCH_H_  // NOLINT
#define LINESEARCH_H_

#include <RcppArmadillo.h>
#include <cmath>

#include <algorithm>
#include <limits>

namespace dragonwell {

  // Line Searches and Backtracking
  template <typename T>
    class LineSearch {
  public:
    LineSearch() {}   // Constructor
    ~LineSearch() {}  // Destructor

    bool GetStep(T &func, double *f, arma::vec *x, const arma::vec &g,
                 const arma::vec &p, const double stepmax);

    void set_message(bool message) { message_ = message; }

  protected:
    bool message_;
    bool IsInfOrNaN(double x);
  };  // class LineSearch

#include "linesearch_impl.h"  // NOLINT

}  // namespace dragonwell

#endif  // LINESEARCH_H_ NOLINT
