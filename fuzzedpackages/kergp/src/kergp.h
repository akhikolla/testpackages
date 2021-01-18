#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

SEXP scores_covMan(SEXP fun, SEXP Xt, SEXP par, SEXP weights, SEXP rho);

SEXP covMat_covMan(SEXP fun, SEXP Xt,  SEXP par, SEXP compGrad, SEXP index, 
		   SEXP rho);   

SEXP covMatMat_covMan(SEXP fun, SEXP X1t, SEXP X2t, SEXP par, SEXP compGrad, 
		      SEXP index,  SEXP rho);

SEXP varVec_covMan(SEXP fun, SEXP Xt,  SEXP par, SEXP compGrad, SEXP index, 
		   SEXP rho);   

SEXP covMat_covTS(SEXP fun, SEXP Xt, SEXP par, SEXP parMap, SEXP compGrad, 
		  SEXP index, SEXP rho);
 
SEXP covMatMat_covTS(SEXP fun, SEXP X1t, SEXP X2t, SEXP par, SEXP parMap, 
		     SEXP compGrad, SEXP index, SEXP rho);

SEXP varVec_covTS(SEXP fun, SEXP Xt, SEXP par, SEXP parMap, SEXP compGrad, 
		  SEXP index, SEXP rho);
 
SEXP scores_covTS(SEXP fun, SEXP Xt, SEXP par, SEXP parMap, SEXP weights,  
		  SEXP rho); 

SEXP k1ExpC(SEXP x1, SEXP x2, SEXP par); 

SEXP k1GaussC(SEXP x1, SEXP x2, SEXP par);

SEXP k1PowExpC(SEXP x1, SEXP x2, SEXP par);

SEXP k1Matern3_2C(SEXP x1, SEXP x2, SEXP par);

SEXP k1Matern5_2C(SEXP x1, SEXP x2, SEXP par);

void kern1Gauss(int *n, double *u, double *range,  double *var,  
		double *kern, double *dkern);

void kern1Exp(int *n, double *u, double *range, double *var, 
              double *kern, double *dkern);

void kern1PowExp(int *n, double *u, double *rangeShape, double *var, 
		 double *kern, double *dkern);

double C_covWhiteNoise(const double *x1, const int *n1, 
		       const double *x2, const int *n2, 
                       const int *d, const int *i1, const int *i2,  const double *var); 

void C_covMat1Mat2_WhiteNoise(const double *x1, const int *n1, 
                              const double *x2, const int *n2, 
                              const int *d, const double *var, double *ans);

SEXP corLev_CompSymm(SEXP par, SEXP nlevels, SEXP lowerSQRT, SEXP compGrad);

SEXP corLev_Symm(SEXP par, SEXP nlevels, SEXP lowerSQRT, SEXP compGrad);

SEXP corLev_LowRank(SEXP par, SEXP nlevels, SEXP rank, SEXP lowerSQRT, SEXP compGrad);

SEXP k1FunExpC(SEXP x);

SEXP k1FunMatern3_2C(SEXP x);

SEXP k1FunMatern5_2C(SEXP x);

SEXP k1FunGaussC(SEXP x);

SEXP k1FunPowExpC(SEXP x, SEXP alpha);
