#include <Rcpp.h>
#include <limits>
#include <math.h>
#include <algorithm>
#include <vector>
#include "constants.h"          //Global Constants
#include "functions.h"
#include "structs.h"

//Namespaces
using namespace Rcpp;
using namespace std;

void getVecsByRefZeroBased(const IntegerMatrix& dm, const NumericMatrix& xyz, int ix, Vec3 * out){ 
  for(int j = 0, k = 0; j < 3; j++){
    k = dm(ix,j);
    out[j] = Vec3(xyz(k,0),xyz(k,1),xyz(k,2)); 
  }
}

//Vector Cross Product of plane defined by V0-a  and V0-c
Vec3 crossProd(const Vec3& V0, const Vec3& b, const Vec3& c){
  Vec3 B = b - V0, C = c - V0, result;
  result.x = B.y*C.z - B.z*C.y;
  result.y = B.z*C.x - B.x*C.z;
  result.z = B.x*C.y - B.y*C.x;
  return result;
}

IntegerMatrix operator+=(const IntegerMatrix& a, const int& b){
    IntegerMatrix r = a;
    for(int i = 0; i < a.ncol(); i++){ 
      for(int j = 0; j < a.nrow(); j++){ 
        r(j,i) += b; 
      } 
    }
    return r;
}


bool pointInLine(Vec2 p1, Vec2 p2, Vec2 pTest, bool strictlyBetween = true){
  bool isBetween = (!strictlyBetween) || 
                   (pTest.x >= min(p1.x,p2.x) - D_TOL && 
                    pTest.x <= max(p1.x,p2.x) + D_TOL && 
                    pTest.y >= min(p1.y,p2.y) - D_TOL && 
                    pTest.y <= max(p1.y,p2.y) + D_TOL);
  if(!isBetween)     return false;
  if(abs(p1.x-p2.x) < D_TOL) return abs(p1.x-pTest.x) < D_TOL; //Vertical
  if(abs(p1.y-p2.y) < D_TOL) return abs(p1.y-pTest.y) < D_TOL; //Horizontal
  double slope = (p2.y-p1.y)/(p2.x-p1.x);
  return (abs(slope*(pTest.x - p1.x) - pTest.y) < D_TOL);
}

Vec3 interpolateZVal(Vec3 a, Vec3 b, double z){ 
  Vec3 c = b - a; c *= (z-a.z)/c.z; c += a; c.z = z; return c; 
}








