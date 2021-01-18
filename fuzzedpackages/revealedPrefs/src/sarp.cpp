////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/////////  Strong Axiom of Revealed Preferences (SARP) functions  //////////////
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

// Slack Revealed Prefs from matrices x and p, for obs i and k
// afriat_par in (0,1), 1 for standard RP
// return true if afriat_par * p_i x_i >= p_i x_k,
// return false if afriat_par * p_i x_i < p_i x_k,

bool SlackRevPref(arma::mat *x, arma::mat *p, unsigned obs_i, unsigned obs_k, 
                  double afriat_par) {
  if (afriat_par * arma::dot(p[0].row(obs_i), x[0].row(obs_i)) -
        arma::dot(p[0].row(obs_i), x[0].row(obs_k)) >= 0) 
    return true;
  return false;
}

////////////////////////////////////////////////////////////////////////////////
// check SARP with a variant of the Floyd-Warshall algorithm
RcppExport SEXP CheckSarp(SEXP x, SEXP p, SEXP afriat) { try {
  // import quantities and prices matrices (x, p)
  NumericMatrix r_x(x), r_p(p);
  unsigned n_obs= r_x.nrow();
  arma::mat mat_x(r_x.begin(), r_x.nrow(), r_x.ncol());
  arma::mat mat_p(r_p.begin(), r_p.nrow(), r_p.ncol());
  double afriat_par= as<double>(afriat);
  
  // Compute direct revealed preferences from matrix p*x
  // exit as soon as contradiction between direct preferences
  // pref(i,j)=1 iff bundle i directly prefered to bundle j
  arma::mat direct_prefs(arma::zeros(n_obs, n_obs));
  for (unsigned i_row= 0; i_row < n_obs; i_row++) {
    for (unsigned i_col= i_row; i_col < n_obs; i_col++) {
      // i_row prefered to i_col?
      if (afriat_par * arma::dot(mat_p.row(i_row), mat_x.row(i_row)) >=
            arma::dot(mat_p.row(i_row), mat_x.row(i_col)))
        direct_prefs(i_row, i_col)= 1;

      // i_col prefered to i_row?
      if (afriat_par * arma::dot(mat_p.row(i_col), mat_x.row(i_col)) >=
            arma::dot(mat_p.row(i_col), mat_x.row(i_row)))
        direct_prefs(i_col, i_row)= 1;
      
      // SARP violation in direct preferences? (ie. WARP violation)
      if (direct_prefs(i_row, i_col) + direct_prefs(i_col, i_row) == 2) {
        if(arma::accu(arma::abs(mat_x.row(i_row) - mat_x.row(i_col))) != 0) {
          // if p_i x_i >= p_i x_k AND p_k x_k >= p_k x_i AND x_i != x_k
          NumericVector violators(2);
          violators(0)= i_row + 1; violators(1)= i_col + 1;
          return List::create(Named("violation", wrap(true)),
                              Named("violators", violators), 
                              Named("direct.violation", wrap(true)));
        }
      }
    }
  }
  
  // Variant of Warshall's algorithm to all find indirect preferences 
  // Exit as soon as SARP violation found
  arma::mat indirect(direct_prefs);
  for (unsigned k= 0; k < n_obs; k++) {
    for (unsigned i_row= 0; i_row < n_obs; i_row++) {
      for (unsigned i_col= 0; i_col < n_obs; i_col++) {
        if (indirect(i_row, i_col) == 0) {
          if(indirect(i_row, k) != 0 && indirect(k, i_col) != 0) {
            indirect(i_row, i_col)= 1;
            if (indirect(i_col, i_row) == 1) {
              if(arma::accu(arma::abs(mat_x.row(i_row) - mat_x.row(i_col))) != 0) 
              {
                // if p_i x_i >= p_i x_k AND p_k x_k >= p_k x_i AND x_i != x_k
                NumericVector violators(2);
                violators(0)= i_row + 1; violators(1)= i_col + 1;
                return List::create(Named("violation", wrap(true)),
                                    Named("violators", violators),
                                    Named("direct.violation", wrap(false))) ;
              }
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


////////////////////////////////////////////////////////////////////////////////
// SARP - recursive function for depth-first search
// returns empty vector if no violation found
// returns vector of violating path if violation found
// (used by DeepSarp)
std::vector<unsigned> RecSarp(unsigned cur_obs, std::vector<bool> * tabu, 
                              unsigned *n_tabu,
                              std::vector<unsigned> ascendence, 
                              arma::mat *x, arma::mat *p, 
                              double afriat_par) {
  for (unsigned i_search= 0; i_search < x[0].n_rows; i_search++) {
    if (!tabu[0][i_search] && i_search!=cur_obs) {
      if (SlackRevPref(x, p, cur_obs, i_search, afriat_par)) {
        // if cur_obs prefered to i_search, check if i_search is in ascendence
        bool b_asc= false;
        unsigned i_asc;
        for (i_asc= 0; i_asc < ascendence.size(); i_asc++) {
          if (ascendence[i_asc] == i_search) {
            b_asc= true;
            break;
          }
        }
        if (b_asc) {
          // if i_search is in ascendence, check for any x(k_asc) != x(cur_obs)
          bool b_all_equal= true;
          for (unsigned k_asc= i_asc; k_asc < ascendence.size(); k_asc++) {
            // try all ascendents k_asc since i_asc
            if (arma::accu(arma::abs(x[0].row(cur_obs) - 
                                      x[0].row(ascendence[k_asc]))) != 0) {
              b_all_equal= false;
              break;
            }
          }
          if (!b_all_equal) {
            // SARP violation
            ascendence.push_back(cur_obs);
            ascendence.push_back(i_search);
            return(ascendence);
          } 
          // if cycle of equal quantities, no action
        } else {
          // if not in ascendence, pursue depth-first search
          std::vector<unsigned> new_asc= ascendence;
          new_asc.push_back(cur_obs);
          
          std::vector<unsigned> pursue= RecSarp(i_search, tabu, n_tabu,
                                           new_asc, x, p, afriat_par);
          if (pursue.size())
            return pursue; // if violation in descendence, return violation
        }
      }
    }
  }
  
  // once all the available choices have been tried tabu the current observation
  tabu[0][cur_obs]= true;
  (*n_tabu)++;
  return std::vector<unsigned>(0);
}

////////////////////////////////////////////////////////////////////////////////
// Check SARP using depth-first search
RcppExport SEXP DeepSarp(SEXP quanti, SEXP prices, SEXP afriat) { try {
  // import arguments
  NumericMatrix r_quanti(quanti), r_prices(prices);
  unsigned n_obs= r_quanti.nrow();
  arma::mat mat_q(r_quanti.begin(), r_quanti.nrow(), r_quanti.ncol());
  arma::mat mat_p(r_prices.begin(), r_prices.nrow(), r_prices.ncol());
  double afriat_par= as<double>(afriat);
  
  // tabu list to be passed to the recursive search
  std::vector<bool> tabu(n_obs, false);
  unsigned n_tabu= 0;
  
  std::vector<unsigned> df_search;
  bool found_violation= false;
  for (unsigned current_obs=0; !found_violation && n_tabu < n_obs ; current_obs++) {
    if (!tabu[current_obs]) {
      // launch recursive search with a non-tabu observation
      df_search= RecSarp(current_obs, &tabu, &n_tabu,
                         std::vector<unsigned>(0), 
                         &mat_q, &mat_p, afriat_par);
      if (df_search.size())
        found_violation= true;
    }
  }
  
  if (found_violation) {
    return List::create(Named("violation", wrap(true)),
                        Named("path", wrap(df_search)));
  }
  return List::create(Named("violation", wrap(false)));
} catch(std::exception &ex) {  
  forward_exception_to_r(ex);
} catch(...) { 
  ::Rf_error("c++ exception (unknown reason)"); 
}
  // return to avoid CRAN warning:
  return wrap("ok");
}
