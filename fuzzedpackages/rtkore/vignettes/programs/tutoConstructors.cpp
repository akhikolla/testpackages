#include "STKpp.h"
using namespace STK;
int main(int argc, char** argv)
{
  // Array2D constructors
  Array2DDiagonal<Real> a(3); a= 1.; // same as Array2DDiagonal<Real> a(3, 1.);
  Array2DSquare<Real> b(3); b = 0.; // same as Array2DSquare<Real> b(3, 0.);
  Array2DUpperTriangular<Real> c(3, 3); c= 2.; // same as Array2DUpperTriangular<Real> c(3, 3, 2.);
  Array2DLowerTriangular<Real> d(3, 3); d = -2.;// same as Array2DLowerTriangular<Real> d(3, 3, -2.);
  Array2DVector<Real> e(3); e= 5.; // same as Array2DVector<Real> e(3, 5.);
  Array2DPoint<Real> f(3); f = 6.; // same as Array2DPoint<Real> f(3, 6.);
  // CArray constructors
  CArraySquare<Real> g(3); g= 0.; // same as CArraySquare<Real> g(3, 0.);
  CArrayVector<Real> j(3); j= 5.; // same as CArrayVector<Real> j(3, 5.);
  CArrayPoint<Real> k(3); k = 6.; // same as CArrayPoint<Real> k(3, 6.);
  return 0;
}
