#include "STKpp.h"
using namespace STK;
int main(int argc, char *argv[])
{
  ArrayXX A(3, 2), B(3, 3);
  A << 1, 2, 3, 4, 5, 6; B << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  stk_cout << _T("B.isRef() = ") << B.isRef() << _T('\n');
  // merge B with A
  A.merge(B);
  stk_cout << _T("A =\n") << A << _T('\n');
  stk_cout << _T("B.isRef() = ") << B.isRef() << _T('\n');
  stk_cout << _T("B =\n") << B<< _T('\n');
  // merge C with A
  VectorX C(3); C << 4, 5, 6;
  A.merge(C);
  stk_cout << _T("A =\n") << A << _T('\n');
  // B has been invalidate by A.merge(C) and cannot be used anymore
  // stk_cout << _T("B =\n") << B; produce an error at execution
  stk_cout << _T("C.isRef() = ") << C.isRef() << _T('\n');
  stk_cout << _T("C =\n") << C;
}
