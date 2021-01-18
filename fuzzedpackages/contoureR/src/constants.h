#include <limits>
#include <math.h>

using namespace std;

#ifndef __constants__
#define __constants__

// Define some local scope constants
const double D_MIN  = numeric_limits<double>::min();
const double D_MAX  = numeric_limits<double>::max();
const double D_NAN  = numeric_limits<double>::quiet_NaN();
const double D_INF  = numeric_limits<double>::infinity();
const double D_TOL  = pow(D_MIN,0.5);

//Other Constants
const int N_COL     = 7;

#endif
