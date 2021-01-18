////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////////  Simulation functions //////////////////////////////////
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

// Recursive depth-first SARP 
// (used by SimAxiom, defined in sarp.cpp)
std::vector<unsigned> RecSarp(unsigned cur_obs, std::vector<bool> * tabu, 
                              unsigned *n_tabu,
                              std::vector<unsigned> ascendence, 
                              arma::mat *x, arma::mat *p, 
                              double afriat_par);

// Check revealed prefs (used by RecGarp, defined in garp.cpp)
unsigned CheckRevPref(arma::mat *x, arma::mat *p, unsigned obs_i,
                      unsigned obs_k, double afriat_par);

// Recursive depth-first GARP 
// returns empty vector if no violation found
// returns vector of violating path if violation found
// (used by SimAxiom)
std::vector<unsigned> RecGarp(unsigned cur_obs, std::vector<bool> * tabu, 
                              unsigned *n_tabu, 
                              std::vector<unsigned> *hist_tabu,
                              std::vector<unsigned> ascendence, 
                              std::vector<bool> strict_asc, 
                              arma::mat *x, arma::mat *p, 
                              double afriat_par) {
  for (unsigned i_search= 0; i_search < x[0].n_rows; i_search++) {
    if (!tabu[0][i_search] && i_search!=cur_obs) {
      unsigned rev_pref= CheckRevPref(x, p, cur_obs, i_search, afriat_par);
      if (rev_pref) {
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
          // if i_search is in ascendence, check if there is a strict ascendence
          // (start looking where it left off, at beginning of loop)
          bool b_strict= false;
          for (; i_asc < ascendence.size(); i_asc++) {
            if (strict_asc[i_asc]) {
              b_strict= true;
              break;
            }
          } 
          if (b_strict) { // found strict cycle, exit
            ascendence.push_back(cur_obs);
            ascendence.push_back(i_search);
            return(ascendence);
          }
        } else {
          // if not in ascendence, pursue depth-first search
          std::vector<unsigned> new_asc= ascendence;
          std::vector<bool> new_strict= strict_asc;
          new_asc.push_back(cur_obs);
          if (rev_pref == 2) {
            new_strict.push_back(true);
          } else new_strict.push_back(false);
          
          std::vector<unsigned> pursue= RecGarp(i_search, 
                                                tabu, n_tabu, hist_tabu,
                                                new_asc, new_strict, x, p, 
                                                afriat_par);
          if (pursue.size())
            return pursue; // if violation in descendence, return violation
        }
      }
    }
  }
  // once all the available choices have been tried,
  // if no nonstrict cycle found, tabu the current obs  
  tabu[0][cur_obs]= true;
  (*n_tabu)++;
  hist_tabu[0].push_back(cur_obs);
  return std::vector<unsigned>(0);
}


////////////////////////////////////////////////////////////////////////////////
// Generate simulated axiom-consistent data
RcppExport SEXP SimAxiom(SEXP nobs, SEXP ngoods, SEXP afriat, SEXP maxit, 
                         SEXP pmin, SEXP pmax, SEXP qmin, SEXP qmax, 
                         SEXP axiom) { try {
  // import arguments
  double afriat_par= as<double>(afriat);
  double p_min= as<double>(pmin);
  double p_max= as<double>(pmax);
  double q_min= as<double>(qmin);
  double q_max= as<double>(qmax);
  unsigned n_obs= as<unsigned>(nobs);
  unsigned n_goods= as<unsigned>(ngoods);
  unsigned max_it= as<unsigned>(maxit);  
  CharacterVector r_axiom(axiom);
  char *the_axiom= r_axiom[0];
  
  arma::mat mat_q= q_min + (q_max - q_min) * arma::randu(1, n_goods);
  arma::mat mat_p= p_min + (p_max - p_min) * arma::randu(1, n_goods);
  arma::mat cur_mat_q(mat_q);
  arma::mat cur_mat_p(mat_p);
    
  bool b_violation;
  unsigned cur_length= 1;
  unsigned cur_it= 0;
  // Variables for depth-first search (GARP or SARP)
  unsigned n_tabu= 0;
  std::vector<bool> tabu;
  std::vector<unsigned> hist_tabu; // history of tabu obs, ie. utility ordering
  std::vector<unsigned> df_search;
  while ((cur_length < n_obs) & (cur_it < max_it)) {
    b_violation= false;
    cur_it++;
    
    // Generate a random candidate observation
    cur_mat_q.insert_rows(0, (q_max - q_min) * arma::randu(1, n_goods));
    cur_mat_p.insert_rows(0, (p_max - p_min) * arma::randu(1, n_goods));
    
    if (strcmp(the_axiom, "WARP")==0) {
      // WARP : test new observation against the previous ones
      for (unsigned i_row= 1; i_row < cur_length + 1; i_row++) {
        if (afriat_par * arma::dot(cur_mat_p.row(0), cur_mat_q.row(0)) >=
              arma::dot(cur_mat_p.row(0), cur_mat_q.row(i_row))) {
          if (afriat_par * arma::dot(cur_mat_p.row(i_row), cur_mat_q.row(i_row)) >=
                arma::dot(cur_mat_p.row(i_row), cur_mat_q.row(0))) {
            if (arma::accu(arma::abs(cur_mat_p.row(i_row) - cur_mat_q.row(0))) != 0)
            { // if p_i x_i >= p_i x_k AND p_k x_k >= p_k x_i AND x_i != x_k
              b_violation= true;
              break;
            }
          }
        }
      }
    } else {
      // GARP or SARP : launch recursive search from the generated observation
      tabu= std::vector<bool>(cur_mat_q.n_rows + 1, false);
      n_tabu= 0;
      hist_tabu.clear();
      df_search.clear();
      if (strcmp(the_axiom, "GARP")==0) {
        df_search= RecGarp(0, &tabu, &n_tabu, &hist_tabu,
                           std::vector<unsigned>(0), 
                           std::vector<bool>(0),
                           &cur_mat_q, &cur_mat_p,
                           afriat_par);
      } else {
        df_search= RecSarp(0, &tabu, &n_tabu,
                           std::vector<unsigned>(0), 
                           &cur_mat_q, &cur_mat_p,
                           afriat_par);
      }
      if (df_search.size())  b_violation= true;
    }
    
    if (b_violation) {
      // violation
      cur_mat_q.shed_row(0);
      cur_mat_p.shed_row(0);
    } else {
      mat_q= cur_mat_q;
      mat_p= cur_mat_p;
      cur_length++;
    }
  }
    
  return List::create(Named("x", wrap(mat_q)), 
                      Named("p", wrap(mat_p)),
                      Named("iter", wrap(cur_it)), 
                      Named("nobs", wrap(cur_length)));
} catch(std::exception &ex) {  
  forward_exception_to_r(ex);
} catch(...) { 
  ::Rf_error("c++ exception (unknown reason)"); 
}
  // return to avoid CRAN warning:
  return wrap("ok");
}
