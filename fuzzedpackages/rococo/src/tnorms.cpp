#include "tnorms.h"

// minimum t-norm

double min_tnorm(double in1, double in2)
{
    return MIN(in1, in2);
}

// product t-norm

double prod_tnorm(double in1, double in2)
{
    return in1 * in2;
}
 
// Lukasiewicz t-norm

double lukasiewicz_tnorm(double in1, double in2)
{
    return MAX(0, in1 + in2 - 1);
}
