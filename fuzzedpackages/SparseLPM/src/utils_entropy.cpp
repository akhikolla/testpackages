#include "utils_entropy.h"

double Entropy(double x)
{
  double res;
  if (x > 0) res = x * log(x);
  else res = 0;
  return (-res);
}


