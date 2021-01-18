#include "STKpp.h"
using namespace STK;
int main(int argc, char *argv[])
{
  ArrayXX A(5, 6);
  Law::Normal law(1.,2.);
  A.rand(law);
  stk_cout << _T("A =\n") << A;
  // Remove 2 columns
  A.popBackCols(2);
  // erase two rows in position 1
  A.eraseRows(1,2);
  stk_cout << _T("A =\n") << A;

  PointX B(6); B << 1, 2, 3, 4, 5, 6;
  stk_cout << _T("B =\n") << B;
  // erase three elements at the position 2
  B.erase(2,3);
  stk_cout << _T("B =\n") << B;
}
