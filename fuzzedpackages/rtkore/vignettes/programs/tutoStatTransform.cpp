#include "STKpp.h"
using namespace STK;
int main(int argc, char** argv)
{
  CArray<Real, 5, 8> A;
  CArrayPoint<Real, 8> mu, std, mean;
  Law::Normal law(1,2);
  A.rand(law);
  stk_cout << "mean(A)=\n" << (mu=Stat::mean(A));
  stk_cout << "variance(A)=\n"
           << Stat::varianceWithFixedMean(A, mu, false)<<"\n";
  // standardize using empirical mean and standard deviation (std is computed during)
  Stat::standardize(A, mu, std);
  stk_cout << "mean(A)=\n" << (mean=Stat::mean(A));
  stk_cout << "variance(A)=\n"
           << Stat::varianceWithFixedMean(A, mean, false)<<"\n";
  // undo standardization
  Stat::unstandardize(A, mu, std);
  stk_cout << "mean(A)=\n" << (mean=Stat::mean(A));
  stk_cout << "variance(A)=\n"
           << Stat::varianceWithFixedMean(A, mean, false);

  return 0;
}
