#include <float.h>
#include <Rcpp.h>
#include "tnorms.h"
#include "rcor.h"

using namespace Rcpp;

/* Steinhaus-Johnson-Trotter algorithm for creating all permutations of 1:n
   the following function moves to the next iteration */

bool next_permutation(IntegerVector perm, IntegerVector sign)
{
    int i = 1, lMval = -1, lMindex = -1, n = perm.size();

    if (sign[0] > 0)
    {
	lMval = perm[0];
	lMindex = 0;
    }

    for (; i < (n - 1); i++)
    {
	if (perm[i] > lMval)
	{
	    if (sign[i] > 0)
	    {
		if (perm[i] > perm[i + 1])
		{
		    lMval = perm[i];
		    lMindex = i;
		}
	    }
	    else
	    {
		if (perm[i] > perm[i - 1])
		{
		    lMval = perm[i];
		    lMindex = i;
		}
	    }
	}
    }

    if (sign[i] < 0 && perm[i] > lMval)
    {
	lMval = perm[i];
	lMindex = i;
    }

    if (lMval <= 0)
	return false;

    if (sign[lMindex] < 0)
    {
	perm[lMindex] = perm[lMindex - 1];
	perm[lMindex - 1] = lMval;

	int aux = sign[lMindex];
	sign[lMindex] = sign[lMindex - 1];
	sign[lMindex - 1] = aux;
    }
    else
    {
	perm[lMindex] = perm[lMindex + 1];
	perm[lMindex + 1] = lMval;

	int aux = sign[lMindex];
	sign[lMindex] = sign[lMindex + 1];
	sign[lMindex + 1] = aux;
    }

    for (i = 0; i < n; i++)
	if (perm[i] > lMval)
	    sign[i] = -sign[i];

    return true;
}

/* Rcpp wrapper around next_permutation() for making this function usable from
   R (for exact test in conjunction with user-defined t-norm) */

RcppExport SEXP permNextWrapper(SEXP perm, SEXP sign)
{
    IntegerVector permI(perm);
    IntegerVector signI(sign);

    bool ret = next_permutation(permI, signI);

    if (ret)
	return List::create(_["perm"] = permI, _["sign"] = signI);
    else
	return R_NilValue;
}


RcppExport SEXP rcor_exacttest(SEXP matx, SEXP maty, SEXP tnorm, SEXP tests,
			       SEXP ogamma, SEXP alt, SEXP storeValues)
{
    int i, cnt = 0, oldcnt = 0, tnorm_sel = as<int>(tnorm);
    double c, d;
    NumericVector gamma(1);

    NumericMatrix mat_x(matx); 
    NumericMatrix mat_y(maty); 

    int num_tests = IntegerVector(tests)[0];
    int alternative = IntegerVector(alt)[0];
    // 0 == two.sided, 1 == less, 2 == greater
    double old_gamma = NumericVector(ogamma)[0];

    double (*tnorm_fp)(double, double);
	
    IntegerVector perm(mat_x.nrow());
    IntegerVector sign(mat_x.nrow());

    switch(tnorm_sel)
    {
        case  1 : tnorm_fp = min_tnorm; break;
        case  2 : tnorm_fp = prod_tnorm; break;
        case  3 : tnorm_fp = lukasiewicz_tnorm; break;
        default : tnorm_fp = min_tnorm;
    }
	
    for (i = 0; i < mat_x.nrow(); i++)
    {
	perm[i] = i;
	sign[i] = -1;
    }

    double mean = 0, M2 = 0, delta;

    LogicalVector store(storeValues);

    int numGamma = (store[0] ? num_tests : 0);

    NumericVector permGamma(numGamma);

    i = 0;

    do
    {
	get_sums(mat_x, mat_y, perm, tnorm_fp, &c, &d);
	gamma[0] = (fabs(c + d) <= DBL_EPSILON ? 0 : (c - d) / (c + d));

	if (numGamma) permGamma[i] = gamma[0];

	/* It would, of course, be more efficient to define 'gamma' as a
	   simple double. Unfortunately, this allows the compiler to put
	   'gamma' into a register that could possibly have a different
	   precision, thereby, leading to false results when comparing two
	   numbers that should be considered as exactly equal. */

	switch (alternative)
	{
	    case 0: if (fabs(gamma[0]) >= fabs(old_gamma)) cnt++; break;
	    case 1: if (gamma[0] <= old_gamma)             cnt++; break;
	    case 2: if (gamma[0] >= old_gamma)             cnt++; break;
	}

	oldcnt = cnt;

	i++;

	// computation of mean and variance according to Knuth/Welford

	delta = gamma[0] - mean;
	mean += (delta / i);
	M2 += (delta * (gamma[0] - mean));
    }
    while(next_permutation(perm, sign));

    if (numGamma)
	return List::create(_["cnt"] = cnt, _["ogamma"] = old_gamma,
			    _["H0mu"] = mean,
			    _["H0sd"] = sqrt(M2 / (num_tests - 1)),
	                    _["values"] = permGamma);
    else
	return List::create(_["cnt"] = cnt, _["ogamma"] = old_gamma,
			    _["H0mu"] = mean,
			    _["H0sd"] = sqrt(M2 / (num_tests - 1)));
}
