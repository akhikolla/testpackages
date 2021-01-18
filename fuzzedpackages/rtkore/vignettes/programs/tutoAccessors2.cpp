#include "STKpp.h"
using namespace STK;
int main(int argc, char** argv)
{
  ArrayXX t(5, 5); // array of size 5x5
  for (int i=t.beginRows(); i<t.endRows(); i++)
  {
    PointX r(t.row(i), true); // create a reference on the i-th row of t
    // fill the i-th row of t with the number i
    r= i; // same as for (int j=r.begin(); j<r.end(); j++) { r[j] = i;}
  }
  std::cout << "t =\n" << t;
  return 0;
}
