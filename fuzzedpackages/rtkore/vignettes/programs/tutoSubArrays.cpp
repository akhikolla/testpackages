#include "STKpp.h"
using namespace STK;
int main(int argc, char** argv)
{
  VectorX b(3, 0);
  std::cout << "b=\n" << b;
  b.sub(Range(2)) = 1.;
  std::cout << "b=\n" << b << "\n";
  ArrayXX a(3, 4); a << 1.,2.,3.,4.
                      , 1.,2.,3.,4.
                      , 1.,1.,1.,1.;
  std::cout << "a=\n" << a;
  a.col(1) = b;
  a.row(1, Range(2,2)) = 0.;
  a.sub(Range(1,2), Range(2,2)) = 9.;
  std::cout << "a=\n" << a;
  return 0;
}
