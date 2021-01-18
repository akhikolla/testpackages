#ifndef _TNORMS_H
#define _TNORMS_H

#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) < (b)) ? (b) : (a))

double min_tnorm(double, double);
double prod_tnorm(double, double);
double lukasiewicz_tnorm(double, double);

#endif
