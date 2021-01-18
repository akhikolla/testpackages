#include "STKpp.h"
using namespace STK;
/** @ingroup tutorial */
int main(int argc, char** argv)
{
  ArraySquareX a(4); a << 1.,2.,3.,4.
                        , 1.,2.,3.,4.
                        , 1.,1.,1.,1.
                        , 1.,1.,1.,1.;
  // create a reference
  ArrayXX b(a.sub(Range(2), Range(3)), true);
  b = -1.;
  std::cout << "Modified a=\n" << a << "\n";

  CArray<Real, 3, 4> c; c << 1.,2.,3.,4.
                           , 1.,2.,3.,4.
                           , 1.,1.,1.,1.;
  // create a reference
  CVector3 d(c.col(1), true);
  d = -1.;
  std::cout << "Modified c=\n" << c;
  return 0;
}
