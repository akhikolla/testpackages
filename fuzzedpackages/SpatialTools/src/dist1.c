#include <R.h>
#include <Rmath.h>

void dist1_c(double *x, int *nc, int *nr, double *d)
{
	int ncol = nc[0];
	int nrow = nr[0];
	int i, j, k;
	
	for(i = 0; i < nrow - 1; i++)
	{
		for(j = i+1; j < nrow; j++)
		{
			for(k = 0; k < ncol; k++)
			{
				d[nrow*i + j] += R_pow(x[nrow*k + i] - x[nrow*k + j], 2);
			}
			d[nrow*i + j] = R_pow(d[nrow*i + j], 0.5);
			d[nrow*j + i] = d[nrow*i + j];
		}
	}
}
