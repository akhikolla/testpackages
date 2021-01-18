#include "STKpp.h"
using namespace STK;
/** @ingroup tutorial */
int main(int argc, char** argv)
{
 CArray22 a;
 a << 0, 1,
      2, 3;
 stk_cout << "a.sum()= " << a.sum()   << _T("\n");
 stk_cout << "a.minElt()= " << a.minElt() << _T("\n");
 stk_cout << "a.maxElt()= " << a.maxElt() << _T("\n");
 stk_cout << "(a > 0).count()= " << (a > 0).count() << _T("\n");
 stk_cout << "(a > 0).all()= " << (a > 0).all() << _T("\n");
 stk_cout << "(a > 0).any()= " << (a > 0).any() << _T("\n");
 int i,j;
 stk_cout << "a.minElt(i,j)= " << a.minElt(i,j) << _T("\n");
 stk_cout << "i= " << i << ", j= " << j << _T("\n");
 return 0;
}
