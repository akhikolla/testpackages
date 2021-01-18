#include <STKpp.h>
using namespace STK;
/** @ingroup tutorial */
int main(int argc, char** argv)
{
  CArray<Real, 4, 4> a(4,4);
  a << 0, 1, 2, 3,
       2, 3, 4, 5,
       2, 1, 6, 7,
       0, 3,-1, 2;
  stk_cout << "Inverse general matrix:\n";
  stk_cout << "-----------------------\n";
  stk_cout << a*invert(a);
  stk_cout << "\nInverse upper-symmetric matrix:\n";
  stk_cout << "-------------------------------\n";
  stk_cout << a.upperSymmetrize();
  stk_cout << a.upperSymmetrize()*invert(a.upperSymmetrize());
  stk_cout << "\nInverse lower-symmetric matrix:\n";
  stk_cout << "-------------------------------\n";
  stk_cout << a.lowerSymmetrize();
  stk_cout << a.lowerSymmetrize()*invert(a.lowerSymmetrize());
  return 0;
}
