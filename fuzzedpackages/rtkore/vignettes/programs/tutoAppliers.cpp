#include "STKpp.h"
using namespace STK;
/** @ingroup tutorial */
int main(int argc, char** argv)
{
 CArray2X a(2, 5);
 stk_cout << "a.randUnif() =\n"   <<  a.randUnif();
 stk_cout << "a.randGauss()=\n" << a.randGauss();
 Law::Gamma law(1, 1.5);
 stk_cout << "a.rand(law)  =\n" << a.rand(law);
 stk_cout << "a.setValue(1)=\n" << a.setValue(1);
 stk_cout << "a+=2 =\n" << (a+=2);
 stk_cout << "a*=2 =\n" << (a*=2);
 stk_cout << "a/=2 =\n" << (a/=2);
 stk_cout << "a-=2 =\n" << (a-=2);
 return 0;
}
