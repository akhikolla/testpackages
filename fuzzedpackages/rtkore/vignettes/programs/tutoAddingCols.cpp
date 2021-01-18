#include "STKpp.h"
using namespace STK;
int main(int argc, char *argv[])
{
  ArrayXX A(4, Range(0,2));
  Law::Normal law(1.,2.);
  A.rand(law);
  stk_cout << _T("A =\n") << A;
  // Adding a column with value 5
  A.pushFrontCols(Const::VectorX(4)*5);
  // insert an uninitialized columns at place 2
  A.insertCols(2, 2);
  // set x^2 and x^3 with x in columns 1
  A.col(2) = A.col(1).square();
  A.col(3) = A.col(1).cube();
  // add an uninitialized columns at the end
  A.pushBackCols(1);
  // set x^2 with x in columns 4
  A.col(5) = A.col(4).square();
  stk_cout << _T("A =\n") << A;
  // do the same with a vector
  VectorX B(4); B << 1, 3, 5, 6;
  stk_cout << _T("B = ") << B.transpose();
  // insert 1 uninitialized element in place 2
  B.insertElt(2); B[2] = 3; // set value to 3
  stk_cout << _T("B = ") << B.transpose();
}
