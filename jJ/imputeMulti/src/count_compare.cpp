#include <Rcpp.h>
#include <iostream>
#include <cmath>
#include <string>

using namespace Rcpp;
using namespace std;


// @title count compare (internal only)
// @description Given a dataset and a data.frame of comparison patterns,
// count the number of occurances of each pattern. Internally called by \code{count_levels}.
// @param dat A \code{data.frame}. All variables must be factors
// @param x A \code{data.frame} consisting of all possible patterns for matching
// with \code{dat}.
// @param hasNA A string. Denotes if \code{dat} has complete data or not.
  // \code{"no"} - there are no missing values, count observed patterns
  // \code{"count.obs"} - there are missing values, count the marginally observed patterns
  // \code{"count.miss"} - there are missing values, count the full observed-and-missing patterns

// [[Rcpp::export]]
IntegerVector count_compare (IntegerMatrix& x, IntegerMatrix& dat, const std::string& hasNA) {

  int nr_x = x.nrow(), nr_dat = dat.nrow(), nc_x = x.ncol();
  IntegerVector out(nr_x, 0);

  // Basic strategy: Loop through dat to find match in x. Each j in J has 1 match i in I
  // Once match is found, break and j++
  // match criteria differs by hasNA.
  // many to one matching (many dat.row(j) may match one x.row(i))

  if (hasNA == "no") {
    for (int j = 0; j < nr_dat; j++) {
      for (int i = 0; i < nr_x; i++) {
        if (is_true(all(x.row(i) == dat.row(j)))) {
            ++out[i];
            break; // x.row(i) are unique -- can only have one match
        }
      }
    }
  } else if (hasNA == "count.obs") {
    for (int j = 0; j < nr_dat; j++) {
      for (int i = 0; i < nr_x; i++) {
        // find observed values. Observed values in dat(j,) must match same values in x(i,)
        IntegerVector x_exist, dat_exist;
        for (int k = 0; k < nc_x; k++) {
          if (dat.row(j)[k] != NA_INTEGER) {
            dat_exist.push_back(dat.row(j)[k]);
            x_exist.push_back(x.row(i)[k]);
          }
        }
        if (is_true(all(x_exist == dat_exist))) {
            ++out[i];
            break; // x.row(i) are unique -- can only have one match
        }
      }
    }
  } else if (hasNA == "count.miss") {
    for (int j = 0; j < nr_dat; j++) {
      for (int i = 0; i < nr_x; i++) {
        // find NA values
        // split both rows into two vectors: one of indices of missing values and one of observed values
        // both must match for equivalence
        IntegerVector x_na, x_exist, dat_na, dat_exist;
        for (int k = 0; k < nc_x; k++) {
          if (x.row(i)[k] == NA_INTEGER) {
            x_na.push_back(k);
          } else {
            x_exist.push_back(x.row(i)[k]);
          }
          if (dat.row(j)[k] == NA_INTEGER) {
            dat_na.push_back(k);
          } else {
            dat_exist.push_back(dat.row(j)[k]);
          }
        }
        // if (same number NA, same positions NA, same values existing) ==> match
        if (x_na.size() == dat_na.size() && is_true(all(x_na == dat_na)) &&
            is_true(all(x_exist == dat_exist))) {
          ++out[i];
          break; // x.row(i) are unique -- can only have one match
        }
      }
    }
  } else {
    // cout << "ERROR: hasNA improperly specified" << endl;
    return -1;
  }
  // return counts
  return out;
}
