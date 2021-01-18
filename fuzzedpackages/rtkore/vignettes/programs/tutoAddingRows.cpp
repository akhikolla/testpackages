#include "STKpp.h"
using namespace STK;
int main(int argc, char *argv[])
{
  ArrayXX A(Range(0,2), 4);
  Law::Normal law(1.,2.);
  A.rand(law);
  stk_cout << _T("A =\n") << A;
  // Adding a column with 1
  A.pushFrontRows(Const::PointX(4));
  // Insert 2 uninitialize rows at place 2
  A.insertRows(2, 2);
  // Set x^2 and x^3 with x in row 1
  A.row(2) = A.row(1).square();
  A.row(3) = A.row(1).cube();
  // Add 2 unitialized rows
  A.pushBackRows(2);
  // Set x^2 and x^3 with x in row 4
  A.row(5) = A.row(4).square();
  A.row(6) = A.row(4).cube();
  stk_cout << _T("A =\n") << A;

  PointX B(4); B << 1, 3, 5, 6;
  stk_cout << _T("B =\n") << B;
  // insert an unitialized element at place 2
  B.insertElt(2); B[2] = 3; // set value to 3
  stk_cout << _T("B =\n") << B;
}
