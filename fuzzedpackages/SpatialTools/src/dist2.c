#include <R.h>
#include <Rmath.h>

void dist2_c(double *x1, double *x2, int *nc, int *nr1, int *nr2, double *d)
{
	int ncol = nc[0];
	int nrow1 = nr1[0];
	int nrow2 = nr2[0];
	int i, j, k;
	
	for(j = 0; j < nrow2; j++)
	{
		for(i = 0; i < nrow1; i++)
		{
			for(k = 0; k < ncol; k++)
			{
				d[j*nrow1 + i] += (x1[nrow1*k + i] - x2[nrow2*k + j])*(x1[nrow1*k + i] - x2[nrow2*k + j]);
			}
			d[j*nrow1 + i] = R_pow(d[j*nrow1 + i], 0.5);
		}
	}
}
