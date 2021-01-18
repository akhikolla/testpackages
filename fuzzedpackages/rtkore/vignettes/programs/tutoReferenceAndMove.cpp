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

  PointX m;
  m.move(Stat::mean(a));
  std::cout << "m= " << m;
  return 0;
}
