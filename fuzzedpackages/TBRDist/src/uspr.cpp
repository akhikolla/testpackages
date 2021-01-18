#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

#include <iostream>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <list>
#include <memory>
#include <ctime>
#include <cstdlib>
#include "uspr/utree.h"
#include "uspr/unode.h"
#include "uspr/uforest.h"
#include "uspr/tbr.h"
#include "uspr/uspr.h"

// [[Rcpp::export]]
IntegerVector uspr_dist (const StringVector tree1,
                         const StringVector tree2,
                         const LogicalVector useTbrApproxEstimate,
                         const LogicalVector useTbrEstimate,
                         const LogicalVector useReplugEstimate) {
  /* use[...]Estimate all default to TRUE */

  USE_TBR_APPROX_ESTIMATE = useTbrApproxEstimate[0];
  USE_TBR_ESTIMATE = useTbrEstimate[0];
  USE_REPLUG_ESTIMATE = useReplugEstimate[0];

  // label maps to allow string labels
  map<string, int> label_map = map<string, int>();
  map<int, string> reverse_label_map = map<int, string>();

  // read input trees
  if (tree1.size() != tree2.size()) {
    throw length_error("Number of trees in tree1 and tree2 must match");
  }
  IntegerVector ret(tree1.size());
  for (int i = 0; i != tree1.size(); i++) {
    // load into data structures
    string tr1 = as<string>(tree1(i));
    string tr2 = as<string>(tree2(i));
    uforest F1 = uforest(tr1, &label_map, &reverse_label_map);
    F1.normalize_order();
    uforest F2 = uforest(tr2, &label_map, &reverse_label_map);
    F2.normalize_order();
    ret(i) = uspr_distance(F1, F2);
  }
  return ret;
}

// [[Rcpp::export]]
List tbr_dist (const StringVector tree1,
               const StringVector tree2,
               const LogicalVector printMafs,
               const LogicalVector countMafs,
               const LogicalVector optimize,
               const LogicalVector protectB,
               const LogicalVector exact,
               const LogicalVector approximate) {
  /* optimize, protectB and *Estimate default to TRUE, all others to FALSE */
  bool PRINT_mAFS = printMafs[0];
  bool COUNT_mAFS = countMafs[0];

  bool DEFAULT_OPTIMIZATIONS = optimize[0];

  bool COMPUTE_TBR = exact[0];
  bool COMPUTE_TBR_APPROX = approximate[0];

  if (DEFAULT_OPTIMIZATIONS == false) {
    OPTIMIZE_2B = false;
    OPTIMIZE_PROTECT_A = false;
    OPTIMIZE_PROTECT_B = false;
    OPTIMIZE_BRANCH_AND_BOUND = false;
  }

  OPTIMIZE_PROTECT_B = protectB[0];

  // label maps to allow string labels
  map<string, int> label_map = map<string, int>();
  map<int, string> reverse_label_map = map<int, string>();

  // read input trees
  if (tree1.size() != tree2.size()) {
    throw length_error("Number of trees in tree1 and tree2 must match");
  }
  IntegerVector tbr_exact(tree1.size()),
    tbr_above(tree1.size()),
    tbr_below(tree1.size()),
    n_maf(tree1.size());
  StringVector maf_1(tree1.size()),
    maf_2(tree1.size());
  for (int i = 0; i < tree1.size(); i++) {
    // load into data structures
    string tr1 = as<string>(tree1(i));
    string tr2 = as<string>(tree2(i));
    uforest F1 = uforest(tr1, &label_map, &reverse_label_map);
    F1.normalize_order();
    uforest F2 = uforest(tr2, &label_map, &reverse_label_map);
    F2.normalize_order();
    // compute TBR distance
    if (COMPUTE_TBR_APPROX) {
      tbr_above(i) = tbr_high_lower_bound(F1, F2);
      Rcpp::checkUserInterrupt();
      tbr_below(i) = tbr_low_upper_bound(F1, F2);
      Rcpp::checkUserInterrupt();
    }
    if (COMPUTE_TBR) {
      uforest *MAF1 = NULL;
      uforest *MAF2 = NULL;
      tbr_exact[i] = tbr_distance(F1, F2, /*quiet = */ true, &MAF1, &MAF2);
      Rcpp::checkUserInterrupt();

      if (MAF1 != NULL) {
        maf_1(i) = MAF1->str(false, &reverse_label_map);
        delete MAF1;
      }
      if (MAF2 != NULL) {
        maf_2(i) = MAF2->str(false, &reverse_label_map);
        delete MAF2;
      }
    }
    if (PRINT_mAFS) {
      n_maf(i) = tbr_print_mAFs(F1, F2);
    }
    else if (COUNT_mAFS) {
      n_maf(i) = tbr_count_mAFs(F1, F2);
    }
  }
  List ret = List::create(tbr_exact, tbr_above, tbr_below, n_maf, maf_1, maf_2);
  return (ret);
}

// [[Rcpp::export]]
List replug_dist (StringVector tree1,
                  StringVector tree2) {
  /* opt, replugEstimate defaults to TRUE */

  // label maps to allow string labels
  map<string, int> label_map = map<string, int>();
  map<int, string> reverse_label_map = map<int, string>();

  // read input trees
  if (tree1.size() != tree2.size()) {
    throw length_error("Number of trees in tree1 and tree2 must match");
  }
  IntegerVector replug(tree1.size());
  StringVector maf_1(tree1.size()),
    maf_2(tree1.size());
  for (int i = 0; i < tree1.size(); i++) {
    // load into data structures
    string tr1 = as<string>(tree1(i));
    string tr2 = as<string>(tree2(i));
    uforest F1 = uforest(tr1, &label_map, &reverse_label_map);
    F1.normalize_order();
    uforest F2 = uforest(tr2, &label_map, &reverse_label_map);
    F2.normalize_order();

    uforest *MAF1 = NULL;
    uforest *MAF2 = NULL;
    replug(i) = replug_distance(F1, F2,  /*quiet = */ true, &MAF1, &MAF2);

    if (MAF1 != NULL) {
      maf_1(i) = MAF1->str(false, &reverse_label_map);
      delete MAF1;
    }
    if (MAF2 != NULL) {
      maf_2(i) = MAF2->str(false, &reverse_label_map);
      delete MAF2;
    }

  }

  List ret = List::create(replug, maf_1, maf_2);
  return (ret);
}
