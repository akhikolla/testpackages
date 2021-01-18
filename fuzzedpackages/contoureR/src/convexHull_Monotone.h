#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <Rcpp.h>
#include "structs.h"

using namespace Rcpp;
using namespace std;

#ifndef __CONVEXHULLAM__
#define __CONVEXHULLAM__

vector<Node2> convex_hull(vector<Node2> P,bool includeColinear);
vector<Node2> convexHullAM(NumericVector x, NumericVector y, bool includeColinear);
IntegerVector convexHullAM_Indexes(NumericVector x, NumericVector y, bool includeColinear, bool zeroBased);
vector<int> convexHullAM_IndexesVector(NumericVector x, NumericVector y, bool includeColinear, bool zeroBased);

#endif /* __CONVEXHULLAM__ */
