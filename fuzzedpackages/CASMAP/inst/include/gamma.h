#ifndef CGAMMAC
#define CGAMMAC

#include "math.h"

double regularizedLowerIncompleteGamma(double x, double alpha);

/*
Function to compute a complemented incomplete gamma function based on a continued fraction
*/
double complementedIncompleteGamma(double x, double alpha);

#endif
