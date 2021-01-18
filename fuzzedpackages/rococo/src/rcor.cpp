#include <Rcpp.h>
#include "tnorms.h"
#include "rcor.h"

using namespace Rcpp;

void get_sums (NumericMatrix mat_x, NumericMatrix mat_y, IntegerVector perm,
	       double (*tnorm_fp)(double, double), double *sum, double *sum_t)
{
    int i, j, n = mat_x.nrow();
    *sum = *sum_t = 0.0;
	
    for (i = 0; i < n; i++)
    {
	for (j = i + 1; j < n; j++)
	{
	    *sum += tnorm_fp(mat_x(i, j), mat_y(perm(i), perm(j)));
	    *sum += tnorm_fp(mat_x(j, i), mat_y(perm(j), perm(i)));
			
	    *sum_t += tnorm_fp(mat_x(i, j), mat_y(perm(j), perm(i)));
	    *sum_t += tnorm_fp(mat_x(j, i), mat_y(perm(i), perm(j)));
	}
    }
}

RcppExport SEXP rcor (SEXP matx, SEXP maty, SEXP tnorm)
{
    int i, tnorm_sel = as<int>(tnorm);

    NumericMatrix mat_x(matx); 
    NumericMatrix mat_y(maty); 
	
    IntegerVector perm(mat_x.nrow());	

    double (*tnorm_fp)(double, double), c, d;

    switch(tnorm_sel)
    {
        case  1 : tnorm_fp = min_tnorm; break;
        case  2 : tnorm_fp = prod_tnorm; break;
        case  3 : tnorm_fp = lukasiewicz_tnorm; break;
        default : tnorm_fp = min_tnorm;
    }

    for (i = 0; i < mat_x.nrow(); i++)
	perm[i] = i;

    get_sums(mat_x, mat_y, perm, tnorm_fp, &c, &d);	
	
    return List::create(_["c"] = c, _["d"] = d);
}
