/*--------------------------------------------------------
        Implementation of Andrew's Monotone Chain 
                2D Convex Hull Algorithm.
           Adapted From: http://goo.gl/krARoZ

Adapted to Rcpp by:
  Author:   Nicholas Hamilton, 
  Date:     23rd August 2015
  Email:    n.hamilton@unsw.edu.au
   
--------------------------------------------------------*/
#include <algorithm>
#include <vector>
#include <Rcpp.h>
#include "structs.h"
#include "convexHull_Monotone.h"

using namespace Rcpp;
using namespace std;

typedef double coord_t;   // coordinate type
typedef double coord2_t;  // must be big enough to hold 2*max(|coordinate|)^2

//Convex Hull via Andrews Monotone -- Return Indexes
// [[Rcpp::export]]
IntegerVector convexHullAM_Indexes(NumericVector x, NumericVector y, bool includeColinear=false, bool zeroBased=true){
  vector<Node2> p = convexHullAM(x,y,includeColinear);
  IntegerVector ret(p.size());
  for(size_t i = 0; i < p.size(); i++)
    ret(i) = p[i].id + (zeroBased ? 0 : 1);
  return ret;
}

vector<int> convexHullAM_IndexesVector(NumericVector x, NumericVector y, bool includeColinear=false,bool zeroBased=true){
  vector<Node2> p = convexHullAM(x,y,includeColinear);
  vector<int> ret;
  for(size_t i = 0; i < p.size(); i++)
    ret.push_back(p[i].id + (zeroBased ? 0 : 1));
  return ret;
}

//Convex Hull via Andrews Monotone -- Return x,y coordinates
// [[Rcpp::export]]
NumericMatrix convexHullAM_Points(NumericVector x, NumericVector y,bool includeColinear=false){
  vector<Node2> p = convexHullAM(x,y,includeColinear); 
  NumericMatrix ret(p.size(),2);
  for(size_t i = 0; i < p.size(); i++){
    ret(i,0) = p[i].x;
    ret(i,1) = p[i].y;
  }
  return ret;
}

vector<Node2> convexHullAM(NumericVector x, NumericVector y,bool includeColinear=false){
  vector<Node2> p;
  if(x.size() != y.size())
    throw std::invalid_argument("x and y must be of the same length");
  for(int i = 0; i < x.size(); i++)
    p.push_back( Node2(i,x(i),y(i)) );
  p = convex_hull(p,includeColinear);
  return p;
}

 
// 2D cross product of OA and OB vectors, i.e. z-component of their 3D cross product.
// Returns a positive value, if OAB makes a counter-clockwise turn,
// negative for clockwise turn, and zero if the points are collinear.
coord2_t cross(const Node2 &O, const Node2 &A, const Node2 &B){
	return (A.x - O.x) * (coord2_t)(B.y - O.y) - (A.y - O.y) * (coord2_t)(B.x - O.x);
}

//Comparator
bool compare_node2(const Node2& a, const Node2& b){
  return a.x < b.x || (a.x == b.x && a.y < b.y);
}
 
// Returns a list of points on the convex hull in counter-clockwise order.
// Note: the last point in the returned list is the same as the first one.
vector<Node2> convex_hull(vector<Node2> P, bool includeColinear){
	int n = P.size(), k = 0;
	vector<Node2> H(2*n);
 
	// Sort points lexicographically
	sort(P.begin(), P.end(),compare_node2);
 
	// Build lower hull
	for (int i = 0; i < n; i++) {
    if(!includeColinear)
		  while (k >= 2 && cross(H[k-2], H[k-1], P[i]) <= 0) k--;
    else
      while (k >= 2 && cross(H[k-2], H[k-1], P[i]) <  0) k--;
		H[k++] = P[i];
	}
 
	// Build upper hull
	for (int i = n-2, t = k+1; i >= 0; i--) {
    if(!includeColinear)
		  while (k >= t && cross(H[k-2], H[k-1], P[i]) <= 0) k--;
    else
      while (k >= t && cross(H[k-2], H[k-1], P[i]) <  0) k--;
		H[k++] = P[i];
	}
 
	H.resize(k);
	return H;
}
