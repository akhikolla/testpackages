#include "STKpp.h"
using namespace STK;
int main(int argc, char** argv)
{
 CArray3X a(3,5);  a.randGauss();
 stk_cout << "Stat::mean(a)=\n"     << Stat::mean(a);
 stk_cout << "Stat::meanByRow(a)=\n"<< Stat::meanByRow(a);
 stk_cout << "Stat::variance(a)=\n" << Stat::variance(a);
 // compute the biased variance (divided by N=3) with mean fixed to 0
 stk_cout << "Stat::varianceWithFixedMean(a,0,false)=\n"
          << Stat::varianceWithFixedMean(a, Const::VectorX(5)*0, false);
 a(1,2) = a(2,2) = Arithmetic<Real>::NA();
 stk_cout << "Stat::meanSafe(a)=\n"      << Stat::meanSafe(a);
 stk_cout << "Stat::meanSafeByRow(a)=\n" << Stat::meanSafeByRow(a);
 return 0;
}
