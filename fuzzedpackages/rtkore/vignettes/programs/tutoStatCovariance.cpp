#include "STKpp.h"
using namespace STK;
/** @ingroup tutorial */
int main(int argc, char** argv)
{
  // create coovariance matrix and its Cholesky decomposition
  CArraySquare<Real, 3> s; s << 2.0, 0.8, 0.36,
                                0.8, 2.0, 0.8,
                                0.36,0.8, 1.0;
  Array2DLowerTriangular<Real> L; Array2DDiagonal<Real> D;
  cholesky(s, D, L);
  // create correlated data set
  CArray<Real, 100, 3> a;
  a = a.randGauss() * D.sqrt() * L.transpose() + 1;
  stk_cout << "True sigma=\n" << s;
  stk_cout << "Estimated sigma=\n" << Stat::covariance(a);
  return 0;
}
