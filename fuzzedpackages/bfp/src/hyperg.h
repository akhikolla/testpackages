#ifndef HYPERG_H_
#define HYPERG_H_

#include <R.h> 
#include <Rmath.h>
#include <Rdefines.h>

#ifdef __cplusplus
extern "C" {
#endif 		

// 	#include "mconf.h"
// 	#include <math.h> 
// 	#include <string.h>
// 	#include <float.h>

double hyp2f1(double a, double b, double c, double x);
	
#ifdef __cplusplus
}
#endif



double logPsi(double b, double c, int n, int p, double R2);
double logBF_hyperg(double R2, int n, int p, double alpha);
double posteriorExpectedg_hyperg(double R2, int n, int p, double alpha, double logBF);
double posteriorExpectedShrinkage_hyperg(double R2, int n, int p, double alpha, double logBF);


#endif /*HYPERG_H_*/
