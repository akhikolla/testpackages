////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////////////////////  Crawford-Pendakur functions ///////////////////////////
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

// boolean check GARP from an arma::mat argument 
// (used by CpUp, defined in garp.cpp)
bool ViolateGarp(arma::mat mat_px, double afriat_par);

////////////////////////////////////////////////////////////////////////////////
// Crawford-Pendakur lower bound
RcppExport SEXP CpLow(SEXP x, SEXP p, SEXP samples, SEXP afriat) { try {
  // import quantities and prices matrices (x, p)
  NumericMatrix r_x(x), r_p(p);
  unsigned n_obs= r_x.nrow();
  arma::mat mat_x(r_x.begin(), r_x.nrow(), r_x.ncol());
  arma::mat mat_p(r_p.begin(), r_p.nrow(), r_p.ncol());
  double afriat_par= as<double>(afriat);
  
  // import samples matrix
  NumericMatrix r_samp(samples);
  arma::mat mat_samp(r_samp.begin(), r_samp.nrow(), r_samp.ncol()); 
  unsigned n_tries= r_samp.ncol();
  
  // main loop : "winners" are pairwise GARP incompatible observations
  std::vector<unsigned> best_winners, winners, hist_n_clust;
  unsigned the_obs, the_winner, n_incompatible;

  for (unsigned i_try= 0; i_try < n_tries; i_try++) {
    winners.clear();
    winners.push_back(mat_samp(0, i_try));
    for (unsigned i_obs= 1; i_obs < n_obs; i_obs++) {
      n_incompatible= 0;
      the_obs= mat_samp(i_obs, i_try); // index of current observation
      for (unsigned i_winner= 0; i_winner < winners.size(); i_winner++) {
        the_winner= winners[i_winner]; // index of current winner

        if (afriat_par * arma::dot(mat_p.row(the_obs), mat_x.row(the_obs)) >=  
            arma::dot(mat_p.row(the_obs), mat_x.row(the_winner))) {
          if (afriat_par * arma::dot(mat_p.row(the_winner), mat_x.row(the_winner)) >  
                arma::dot(mat_p.row(the_winner), mat_x.row(the_obs))) {
            n_incompatible++;
          }
        } else if (afriat_par * arma::dot(mat_p.row(the_obs), mat_x.row(the_obs)) >
            arma::dot(mat_p.row(the_obs), mat_x.row(the_winner))) {
          if (afriat_par * arma::dot(mat_p.row(the_winner), mat_x.row(the_winner)) >=
                arma::dot(mat_p.row(the_winner), mat_x.row(the_obs))) {
            n_incompatible++;
          }
        }
      }
      // new winner if pairwise-incompatible with all previous winners
      if (n_incompatible == winners.size()) 
        winners.push_back(the_obs);
    }
    hist_n_clust.push_back(winners.size());
    if (winners.size() > best_winners.size())
      best_winners= winners;
  }
  return List::create(Named("violators", wrap(best_winners)),
                      Named("n.clust", wrap((unsigned) best_winners.size())),
                      Named("hist.n.types", wrap(hist_n_clust)));
} catch(std::exception &ex) {  
  forward_exception_to_r(ex);
} catch(...) { 
  ::Rf_error("c++ exception (unknown reason)"); 
}
  // return to avoid CRAN warning:
  return wrap("ok");
}

////////////////////////////////////////////////////////////////////////////////
// Crawford-Pendakur upper bound classic algorithm 
// (complete Floyd-Warshall GARP check at each stage)
RcppExport SEXP CpUp(SEXP px, SEXP samples, SEXP afriat) { try {
  // import prices*quantities matrix (px)
  NumericMatrix r_px(px);
  unsigned n_obs= r_px.nrow();
  arma::mat mat_px(r_px.begin(), n_obs, n_obs);
  double afriat_par= as<double>(afriat);
  
  // import samples matrix
  NumericMatrix r_samp(samples);
  arma::mat mat_samp(r_samp.begin(), r_samp.nrow(), r_samp.ncol()); 
  unsigned n_tries= r_samp.ncol();
  
  // main loop
  arma::uvec best_clustering(n_obs); // best clustering
  arma::uvec best_clusterpop; // clusterpops (max n_clust=n_obs)
  unsigned best_n_clust= 0; // best found number of clusters

  NumericVector hist_n_clust(r_samp.ncol());
    
  arma::uvec clustering(n_obs); // clustering vector
  arma::uvec clusterpop(n_obs); // cluster population vector (max size n_obs)
  
  arma::uvec current_indices;
  arma::mat current_matrix;
  arma::uvec tmp_clusterpop;
  arma::uvec cluster_order;
  for (unsigned i_try= 0; i_try < n_tries; i_try++) {
    unsigned n_clust= 1;
    clustering.zeros(); // zero = no cluster
    clustering(mat_samp(0, i_try))= 1; // cluster of first sample is 1
    clusterpop.zeros(); // reset clusterpops
    clusterpop(0)= 1;
    
    for (unsigned ind_samp= 1; ind_samp < n_obs; ind_samp++) {
      bool b_found= false;
      unsigned i_obs= mat_samp(ind_samp, i_try); // index of current observation
      
      tmp_clusterpop= clusterpop.rows(0, n_clust - 1);
      cluster_order= arma::stable_sort_index(tmp_clusterpop, "descend"); // descend
      for (unsigned index_clust= 0; index_clust < n_clust; index_clust++) {
        // take the "index_clust"th biggest cluster
        unsigned i_clust= (unsigned)cluster_order(index_clust); 
        
        current_indices= find(clustering == (i_clust + 1));
        current_indices.insert_rows(current_indices.n_rows, 
                                    arma::ones<arma::uvec>(1));
        current_indices(current_indices.n_rows-1)= i_obs;
        
        current_matrix= mat_px(current_indices, current_indices);
        if(!ViolateGarp(current_matrix, afriat_par)) { // no GARP violation
          clustering(i_obs)= i_clust + 1;
          clusterpop(i_clust)= clusterpop(i_clust) + 1;
          b_found= true;
          break;
        }
      }
      if (!b_found) {// if no fitting cluster found 
        clusterpop(n_clust)= 1;
        n_clust++;
        clustering(i_obs)= n_clust;
      }
    }
    
    // check against current best run
    if ( (i_try == 0) || (n_clust < best_n_clust) ) {
      best_clustering= clustering;
      best_clusterpop= clusterpop;
      best_n_clust= n_clust;
    }
    hist_n_clust(i_try)= n_clust;
  }

  return List::create(Named("clustering", wrap(best_clustering)),
                      Named("cluster.pop", wrap(best_clusterpop)),
                      Named("hist.n.types", hist_n_clust));
} catch(std::exception &ex) {  
  forward_exception_to_r(ex);
} catch(...) { 
  ::Rf_error("c++ exception (unknown reason)"); 
}
  // return to avoid CRAN warning:
  return wrap("ok");
}

////////////////////////////////////////////////////////////////////////////////
// Fast Floyd Crawford-Pendakur upper bound: 
// one step of Floyd algorithm for each (obs, cluster) couple
RcppExport SEXP FastUp(SEXP px, SEXP samples, SEXP afriat) { try {
  // import prices*quantities matrix (px)
  NumericMatrix r_px(px);
  unsigned n_obs= r_px.nrow();
  arma::mat mat_px(r_px.begin(), n_obs, n_obs);
  double afriat_par= as<double>(afriat);
  
  // import samples matrix
  NumericMatrix r_samp(samples);
  arma::mat mat_samp(r_samp.begin(), r_samp.nrow(), r_samp.ncol()); 
  unsigned n_tries= r_samp.ncol();
  
  // compute direct preferences (0 : no pref, 1 : equality, 2 : strict preference)
  arma::umat direct= arma::zeros<arma::umat>(n_obs, n_obs);
  for (unsigned i_row= 0; i_row < n_obs; i_row++) {
    for (unsigned i_col= 0; i_col < n_obs; i_col++) {
      if (afriat_par * mat_px(i_row, i_row) == mat_px(i_row, i_col)) {
        direct(i_row, i_col)= 1;
      } else if (afriat_par * mat_px(i_row, i_row) > mat_px(i_row, i_col)) {
        direct(i_row, i_col)= 2;
      }
    }
  }
  
  // main loop
  arma::uvec best_clustering(n_obs); // best clustering
  arma::uvec best_clusterpop; // clusterpops (max n_clust=n_obs)
  unsigned best_n_clust= 0; // best found number of clusters

  NumericVector hist_n_clust(r_samp.ncol());
    
  arma::uvec clustering(n_obs); // clustering vector
  arma::uvec clusterpop(n_obs); // cluster population vector (max size n_obs)
  
  // preference matrices for each cluster (max size n_obs)
  // pref_mat[k](i,j) : 0 - no known preference of i on j (p_i q_i < p_i q_j)
  //                    1 - equality (p_i q_i = p_i q_j) & transitivity
  //                    2 - strict preference (p_i q_i > p_i q_j) & transitivity
  arma::umat *pref_mat= new arma::umat[n_obs]; 
  
  arma::uvec crt_indices; // current indices
  arma::umat crt_pref; // current transitive preferences matrix
  arma::uvec tmp_clusterpop;
  arma::uvec cluster_order;
  bool b_change_margin, b_change_main;
  for (unsigned i_try= 0; i_try < n_tries; i_try++) {
    unsigned n_clust= 1;
    clustering.zeros(); // zero = no cluster
    clustering(mat_samp(0, i_try))= 1; // cluster of first sample is 1
    clusterpop.zeros(); // reset clusterpops
    clusterpop(0)= 1;
    for (unsigned i_prefmat= 0; i_prefmat < n_obs; i_prefmat++)
      pref_mat[i_prefmat].clear(); // clear all preference matrices
    pref_mat[0].ones(1, 1); // first preference matrix gets a single 1
    
    // loop on observations
    for (unsigned ind_samp= 1; ind_samp < n_obs; ind_samp++) {
      bool b_found= false;
      unsigned i_obs= mat_samp(ind_samp, i_try); // index of current observation
      
      tmp_clusterpop= clusterpop.rows(0, n_clust - 1);
      cluster_order= arma::stable_sort_index(tmp_clusterpop, "descend"); // descending
      for (unsigned index_clust= 0; index_clust < n_clust; index_clust++) {
        // take the "index_clust"th biggest cluster
        unsigned i_clust= (unsigned)cluster_order(index_clust); 
        
        // find the indices of observations already included in the cluster
        crt_indices= find(clustering == (i_clust + 1));
        crt_indices.insert_rows(crt_indices.n_rows, 
                                    arma::ones<arma::uvec>(1));
        crt_indices(crt_indices.n_rows-1)= i_obs;
        
        // index of new obs in crt_pref
        unsigned the_k= clusterpop(i_clust); 
        
        // build current preference matrix: previous one and new observation
        crt_pref.clear();
        crt_pref= pref_mat[i_clust];
        crt_pref.insert_rows(crt_pref.n_rows, 
                                 arma::trans(arma::ones<arma::uvec>(crt_pref.n_cols)));
        crt_pref.insert_cols(crt_pref.n_cols, 
                                 arma::ones<arma::uvec>(crt_pref.n_rows));
        for (unsigned i_row= 0; i_row < crt_pref.n_rows ; i_row++) {
          crt_pref(i_row, the_k)= direct(crt_indices(i_row), i_obs);
          crt_pref(the_k, i_row)= direct(i_obs, crt_indices(i_row));
        }
        
        // enforce transitivity on the new preference matrix
        b_change_main= true;
        while (b_change_main) {
          b_change_main= false;
          
          // transitivity on margin
          b_change_margin= true;
          while (b_change_margin) {
            b_change_margin= false;
            
            for (unsigned the_i= 0; the_i <= the_k; the_i++) {
              for (unsigned the_j= 0; the_j < the_k; the_j++) {

                if (crt_pref(the_i, the_k) < 2) {
                  unsigned tmp_ijk= crt_pref(the_i, the_j) * crt_pref(the_j, the_k);
                  if (tmp_ijk > 2) tmp_ijk= 2;
                  if (crt_pref(the_i, the_k) < tmp_ijk) {
                    crt_pref(the_i, the_k)= tmp_ijk;
                    b_change_margin= true;
                  }
                }
                
                if (crt_pref(the_k, the_i) < 2) {
                  unsigned tmp_kji= crt_pref(the_k, the_j) * crt_pref(the_j, the_i);
                  if (tmp_kji > 2) tmp_kji= 2;
                  if (crt_pref(the_k, the_i) < tmp_kji) {
                    crt_pref(the_k, the_i)= tmp_kji;
                    b_change_margin= true;
                  }
                }
              }
            }
          }
          
          // transitivity on previous preferences
          for (unsigned the_i= 0; the_i < the_k; the_i++) {
            for (unsigned the_j= 0; the_j < the_k; the_j++) {
              if (crt_pref(the_i, the_j) < 2) {
                unsigned tmp_ikj= crt_pref(the_i, the_k) * crt_pref(the_k, the_j);
                if (tmp_ikj > 2) tmp_ikj= 2;
                if (crt_pref(the_i, the_j) < tmp_ikj) {
                  crt_pref(the_i, the_j)= tmp_ikj;
                  b_change_main= true;
                }

              }
            }
          }
        }
        
        // check for GARP violation
        b_found= true;
        for (unsigned the_i= 0; the_i <= the_k; the_i++)
          if (crt_pref(the_i, the_i) == 2)
            b_found= false; // GARP violation
        
        if (b_found) { // no GARP violation
          clustering(i_obs)= i_clust + 1;
          clusterpop(i_clust)= clusterpop(i_clust) + 1;
          
          // sort crt_pref with increasing indices and save to pref_mat
          pref_mat[i_clust].clear();
          pref_mat[i_clust]= (arma::umat) crt_pref(arma::sort_index(crt_indices),
                                                   arma::sort_index(crt_indices));
          break;
        }
      }
      if (!b_found) {// if no fitting cluster found 
        clusterpop(n_clust)= 1;
        pref_mat[n_clust]= arma::ones<arma::umat>(1, 1);
        n_clust++;
        clustering(i_obs)= n_clust;
      }
    }
    
    // check against current best run
    if ( (i_try == 0) || (n_clust < best_n_clust) ) {
      best_clustering= clustering;
      best_clusterpop= clusterpop;
      best_n_clust= n_clust;
    }
    hist_n_clust(i_try)= n_clust;
  }
  
  delete[] pref_mat;

  return List::create(Named("clustering", wrap(best_clustering)),
                      Named("cluster.pop", wrap(best_clusterpop)),
                      Named("hist.n.types", hist_n_clust));
} catch(std::exception &ex) {  
  forward_exception_to_r(ex);
} catch(...) { 
  ::Rf_error("c++ exception (unknown reason)"); 
}
  // return to avoid CRAN warning:
  return wrap("ok");
}
