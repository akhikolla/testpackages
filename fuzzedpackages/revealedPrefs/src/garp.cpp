////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/////////  Generalized Axiom of Revealed Preferences (GARP) functions  /////////
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
// Auxiliary functions

// Revealed Prefs from matrices x and p, for obs i and k
// afriat_par in (0,1), 1 for standard RP
// return 0 if p_i x_i < p_i x_k,
// return 1 if p_i x_i == p_i x_k,
// return 2 if p_i x_i > p_i x_k,
unsigned CheckRevPref(arma::mat *x, arma::mat *p, unsigned obs_i,
                      unsigned obs_k, double afriat_par) {
  double d_test= afriat_par * arma::dot(p[0].row(obs_i), x[0].row(obs_i)) -
    arma::dot(p[0].row(obs_i), x[0].row(obs_k));
  if (d_test < 0) return 0;
  if (d_test > 0) return 2;
  return 1;
}

// boolean check GARP from an arma::mat argument (used by CpUp)
bool ViolateGarp(arma::mat mat_px, double afriat_par) {
  unsigned n_obs= mat_px.n_rows;
  arma::mat direct_prefs(arma::zeros(n_obs, n_obs)), 
    direct_strict(arma::zeros(n_obs, n_obs));
    
  // Compute direct (and strict) revealed preferences from matrix p*x
  // exit as soon as contradiction between direct preferences (WARP violation)
  // pref(i,j)=1 iff bundle i directly prefered to bundle j
  for (unsigned i_row= 0; i_row < n_obs; i_row++) {
    for (unsigned i_col= i_row; i_col < n_obs; i_col++) {
      // i_row prefered to i_col?
      if (afriat_par * mat_px(i_row, i_row) > mat_px(i_row, i_col)) {
        direct_prefs(i_row, i_col)= 1;
        direct_strict(i_row, i_col)= 1;
      } else if (afriat_par * mat_px(i_row, i_row) == mat_px(i_row, i_col)) {
        direct_prefs(i_row, i_col)= 1;
        direct_strict(i_row, i_col)= 0;
      }

      // i_col prefered to i_row?
      if (afriat_par * mat_px(i_col, i_col) > mat_px(i_col, i_row)) {
        direct_prefs(i_col, i_row)= 1;
        direct_strict(i_col, i_row)= 1;
      } else if (afriat_par * mat_px(i_col, i_col) == mat_px(i_col, i_row)) {
        direct_prefs(i_col, i_row)= 1;
        direct_strict(i_col, i_row)= 0;
      }
      
      // violation in direct preferences?
      if ((direct_prefs(i_row, i_col) + direct_strict(i_col, i_row) == 2) ||
            (direct_prefs(i_col, i_row) + direct_strict(i_row, i_col) == 2)) {
        return true;
      }
    }
  }

  // Variant of Warshall's algorithm to all find indirect preferences 
  // Exit as soon as a contradiction with direct strict preferences is found
  // (ie GARP violation)
  arma::mat indirect(direct_prefs);
  for (unsigned k= 0; k < n_obs; k++) {
    for (unsigned i_row= 0; i_row < n_obs; i_row++) {
      for (unsigned i_col= 0; i_col < n_obs; i_col++) {
        if (indirect(i_row, i_col) == 0) {
          if(indirect(i_row, k) != 0 && indirect(k, i_col) != 0) {
            indirect(i_row, i_col)= 1;
            if (direct_strict(i_col, i_row) != 0) {
              return true;
            }
          }
        }
      }
    }
  }
  
  // if no violation found
  return false;
}

////////////////////////////////////////////////////////////////////////////////
// GARP - Variant of Floyd-Warshall algorithm
RcppExport SEXP CheckGarp(SEXP px, SEXP afriat) { try {
  // import prices*quantities matrix (px)
  NumericMatrix r_px(px);
  unsigned n_obs= r_px.nrow();
  arma::mat mat_px(r_px.begin(), n_obs, n_obs);
  double afriat_par= as<double>(afriat);
  
  // Compute direct (and strict) revealed preferences from matrix p*x
  // exit as soon as contradiction between direct preferences
  // pref(i,j)=1 iff bundle i directly prefered to bundle j
  arma::mat direct_prefs(arma::zeros(n_obs, n_obs)), 
    direct_strict(arma::zeros(n_obs, n_obs));
  for (unsigned i_row= 0; i_row < n_obs; i_row++) {
    for (unsigned i_col= i_row; i_col < n_obs; i_col++) {
      // i_row prefered to i_col?
      if (afriat_par * mat_px(i_row, i_row) > mat_px(i_row, i_col)) {
        direct_prefs(i_row, i_col)= 1;
        direct_strict(i_row, i_col)= 1;
      } else if (afriat_par * mat_px(i_row, i_row) == mat_px(i_row, i_col)) {
        direct_prefs(i_row, i_col)= 1;
        direct_strict(i_row, i_col)= 0;
      }

      // i_col prefered to i_row?
      if (afriat_par * mat_px(i_col, i_col) > mat_px(i_col, i_row)) {
        direct_prefs(i_col, i_row)= 1;
        direct_strict(i_col, i_row)= 1;
      } else if (afriat_par * mat_px(i_col, i_col) == mat_px(i_col, i_row)) {
        direct_prefs(i_col, i_row)= 1;
        direct_strict(i_col, i_row)= 0;
      }
      
      // violation in direct preferences?
      if ((direct_prefs(i_row, i_col) + direct_strict(i_col, i_row) == 2) ||
            (direct_prefs(i_col, i_row) + direct_strict(i_row, i_col) == 2)) {
        NumericVector violators(2), strict(2);
        violators(0)= i_row + 1; violators(1)= i_col + 1;
        strict(0)= direct_strict(i_row, i_col); 
        strict(1)= direct_strict(i_col, i_row); 
        return List::create(Named("violation", wrap(true)),
                            Named("violators", violators), 
                            Named("strict", strict),
                            Named("direct.violation", wrap(true)));
      }
    }
  }
  
  // Variant of Warshall's algorithm to all find indirect preferences 
  // Exit as soon as a contradiction with direct strict preferences is found
  // (ie GARP violation)
  arma::mat indirect(direct_prefs);
  for (unsigned k= 0; k < n_obs; k++) {
    for (unsigned i_row= 0; i_row < n_obs; i_row++) {
      for (unsigned i_col= 0; i_col < n_obs; i_col++) {
        if (indirect(i_row, i_col) == 0) {
          if(indirect(i_row, k) != 0 && indirect(k, i_col) != 0) {
            indirect(i_row, i_col)= 1;
            if (direct_strict(i_col, i_row) != 0) {
              NumericVector violators(2), strict(2);
              violators(0)= i_row + 1; violators(1)= i_col + 1;
              strict(0)= 0; strict(1)= 1;
              return List::create(Named("violation", wrap(true)),
                                  Named("violators", violators),
                                  Named("strict", strict),
                                  Named("direct.violation", wrap(false))) ;
            }
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
