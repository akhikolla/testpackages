#include "STKpp.h"
using namespace STK;
/** @ingroup tutorial */
int main(int argc, char** argv)
{
  Law::Normal l(2, 1);
  stk_cout << "l.pdf(2)= "     << l.pdf(2)     << "\n";
  stk_cout << "l.lpdf(2)= "    << l.lpdf(2)    << "\n";
  stk_cout << "l.cdf(3.96)= "  << l.cdf(3.96)  << "\n";
  stk_cout << "l.cdfc(3.96)= " << l.cdfc(3.96) << "\n";
  stk_cout << "l.lcdf(3.96)= " << l.lcdf(3.96) << "\n";
  stk_cout << "l.lcdfc(3.96)= "<< l.lcdfc(3.96)<< "\n";
  stk_cout << "l.icdf(0.975)= "<< l.icdf(0.975)<< "\n";
  stk_cout << "l.rand()= "     << l.rand()     << "\n";
  CArray<Real, 2, 3> a;  a = 0.5;
  stk_cout << "-----------\n";
  stk_cout << "a=\n"         << a;
  stk_cout << "a.pdf(l)=\n"  << a.pdf(l);
  stk_cout << "a.lpdf(l)=\n" << a.lpdf(l);
  stk_cout << "a.cdf(l)=\n"  << a.cdf(l);
  stk_cout << "a.lcdf(l)=\n" << a.lcdf(l);
  stk_cout << "a.cdfc(l)=\n" << a.cdfc(l);
  stk_cout << "a.lcdfc(l)=\n"<< a.lcdfc(l);
  stk_cout << "a.icdf(l)=\n" << a.icdf(l);
}
