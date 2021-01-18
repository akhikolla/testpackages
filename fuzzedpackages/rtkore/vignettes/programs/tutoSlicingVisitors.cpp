#include "STKpp.h"
using namespace STK;
/** @ingroup tutorial */
int main(int argc, char** argv)
{
 CArray2X a(2,4);
 a << 0, 1, 2, 3
    , 4, 5, 6, 7;
 stk_cout << "min(a)      = " << min(a);
 stk_cout << "minByRow(a) = " << minByRow(a);
 stk_cout << "mean(a) = "     << mean(a);
 stk_cout << "count(a > 0)= " << count(a > 0);
 stk_cout << "all(a > 0)  = " << all(a > 0);
 stk_cout << "any(a > 0)  = " << any(a > 0);
 a(1,1) = a(1,3) = Arithmetic<Real>::NA();
 stk_cout << "count(a.isNA())= "   << count(a.isNA());
 stk_cout << "countByRow(a.isNA())= "   << countByRow(a.isNA());
 return 0;
}
