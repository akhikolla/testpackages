#include <Rcpp.h>
#include "randgen.h"

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void)
{
	return R::runif(0,1);
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void)
{
	return R::runif(0,1);
}

double genrand_gaussian(void)
{
	return R::rnorm(0,1);
}

