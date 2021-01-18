/* ************************************************************
   
   Various matrix utility functions

   Author: Søren Højsgaard

   ************************************************************ */ 

/* ************************************************************
   ***************************************************************
   These functions *DO* have R counterparts (based on SEXPs)
   ************************************************************
   ************************************************************ */

int C_midx(int *ii, int *jj, int *nrow);

void C_submat(double *X, int *nrX, int *ncX,
	     int *idx1, int *n1, int *idx2, int *n2, double *ans);

void C_transpose(double *X, int *nrX, int *ncX, double *ans);

void C_symmetrize(double *X, int *nrX);

void C_issym(double *X, int *nrX, int *ncX, double *ans);

void C_matadd(double *X, int *nrX, int *ncX,
	     double *Y, int *nrY, int *ncY, double *ans);

void C_matsubt(double *X, int *nrX, int *ncX,
	      double *Y, int *nrY, int *ncY, double *ans);

void C_schursubt(double *X, int *nrX, int *ncX,
		int *idx1, int *n1, int *idx2, int *n2, double *ans);

/* ************************************************************
   ************************************************************
   These functions *DO NOT* have R counterparts (based on SEXPs)
   ************************************************************
   ************************************************************ */

void C_matprod(double *X, int *nrX, int *ncX,
	      double *Y, int *nrY, int *ncY, double *ans);

void C_solve(double *A, int *nrA, double *B, int *ncB);

void C_inverse(double *A, int *nrA);

void C_determinant(double *Ain, int *nrA, double *ans);

void C_traceABsym(double *A, int *nrA, int *ncA,
		 double *B, int *nrB, int *ncB, double *ans);
