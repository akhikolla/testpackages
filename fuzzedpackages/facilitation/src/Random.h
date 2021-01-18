#include<cstdlib>
#ifndef _RANDOM_MINE
#define _RANDOM_MINE
#include"Position.h"
double Random(double max);

bool Bernoulli(double p);

double Exponential(double r);

double Normal(double m, double v);

short RandomSign();

Position RandomDirection();
#endif
