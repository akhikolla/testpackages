#include "utils_entropy.h"

double Entropy(double x)// note that this is the negative entropy !!
{
  double res;
  if (x > 0) res = x * log2(x);
  else res = 0;
  return (res);
}
