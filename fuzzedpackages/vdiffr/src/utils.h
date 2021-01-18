#ifndef __UTILS__
#define __UTILS__

#include <cmath>

double dbl_format(double x) {
  if (std::abs(x) < 0.01)
    return 0.00;
  else
    return x;
}

#include <freetypeharfbuzz.h>


#endif
