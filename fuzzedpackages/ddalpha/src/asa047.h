/*******************************************************************************
Author:
Original FORTRAN77 version by R ONeill; C++ version by John Burkardt.
http://people.sc.fsu.edu/~jburkardt/cpp_src/asa047/asa047.html
*******************************************************************************/

#ifndef ASA047
#define ASA047

void nelmin(double fn(double x[]), int n, double start[], double xmin[],
	double *ynewlo, double reqmin, double step[], int konvge, int kcount,
	int *icount, int *numres, int *ifault);
void timestamp();

#endif

