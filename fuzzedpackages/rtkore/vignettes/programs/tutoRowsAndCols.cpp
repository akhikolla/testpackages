#include "STKpp.h"
using namespace STK;
int main(int argc, char *argv[])
{
  CArrayXX A(4, 6); CArray2X B(2,4, -1.);
  A << 1, 2, 3, 4, 1, 2,
       4, 3, 1, 3, 2, 4,
       1, 3, 4, 2, 1, 4,
       2, 3, 1, 4, 3, 2;
  stk_cout << _T("A =\n") << A;
  stk_cout << _T("B =\n") << B;
  A.sub(B.rows(), B.cols()) = B;  // copy B in A(0:1,0:3)
  stk_cout << _T("A =\n") << A;
  A.row(3,B.cols()) = meanByCol(B); // copy a row-vector in row 3 of A
  stk_cout << _T("A =\n") << A;
}
