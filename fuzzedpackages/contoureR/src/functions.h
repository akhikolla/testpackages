#include <Rcpp.h>
#include "templates.h"
#include "structs.h"

//Namespaces
using namespace Rcpp;
using namespace std;

#ifndef __functions__
#define __functions__

//Function declarations
Vec3 crossProd(const Vec3& a, const Vec3& b, const Vec3& c);

//Get Vector Array by Reference
void getVecsByRefZeroBased(const IntegerMatrix& dm, const NumericMatrix& xyz, int ix, Vec3 *out);

//Check Equality w/ Tolerance
inline bool isEqual(double a, double b){  return abs(a-b) <= D_TOL; }

//Operators
IntegerMatrix operator+=(const IntegerMatrix& a, const int& b);

//Check if a point lies in a line
bool pointInLine(Vec2 p1, Vec2 p2, Vec2 pTest,bool strictlyBetween);

//Interpolation
Vec3 interpolateZVal(Vec3 a, Vec3 b, double z);

#endif // __functions__
