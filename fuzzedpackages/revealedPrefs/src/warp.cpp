////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////  Weak Axiom of Revealed Preferences (WARP) functions  //////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/*  
Copyright 2014 Julien Boelaert.

This file is part of revealedPrefs.

revealedPrefs is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

revealedPrefs is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with revealedPrefs.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <RcppArmadillo.h>
#include <vector>
using namespace Rcpp;

////////////////////////////////////////////////////////////////////////////////
// check WARP
RcppExport SEXP CheckWarp(SEXP x, SEXP p, SEXP afriat) { try {
  // import quantities and prices matrices (x, p)
  NumericMatrix r_x(x), r_p(p);
  unsigned n_obs= r_x.nrow();
  arma::mat mat_x(r_x.begin(), r_x.nrow(), r_x.ncol());
  arma::mat mat_p(r_p.begin(), r_p.nrow(), r_p.ncol());
  double afriat_par= as<double>(afriat);
  
  for (unsigned i_row= 0; i_row < n_obs; i_row++) {
    for (unsigned i_col= i_row+1; i_col < n_obs; i_col++) {
      if (afriat_par * arma::dot(mat_p.row(i_row), mat_x.row(i_row)) >=
            arma::dot(mat_p.row(i_row), mat_x.row(i_col))) {
        if (afriat_par * arma::dot(mat_p.row(i_col), mat_x.row(i_col)) >=
              arma::dot(mat_p.row(i_col), mat_x.row(i_row))) {
          if (arma::accu(arma::abs(mat_x.row(i_col) - mat_x.row(i_row))) != 0)
          { // if p_i x_i >= p_i x_k AND p_k x_k >= p_k x_i AND x_i != x_k
            NumericVector violators(2);
            violators(0)= i_row + 1;
            violators(1)= i_col + 1;
            return List::create(Named("violation", wrap(true)), 
                                Named("violators", violators));
          }
        }
      }
    }
  }
  // if no violation found
  return List::create(Named("violation", wrap(false)));
} catch(std::exception &ex) {  
  forward_exception_to_r(ex);
} catch(...) { 
  ::Rf_error("c++ exception (unknown reason)"); 
}
  // return to avoid CRAN warning:
  return wrap("ok");
}
