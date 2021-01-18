#include <float.h>
#include <Rcpp.h>
#include "tnorms.h"
#include "rcor.h"
#include "rcor_permtest.h"

using namespace Rcpp;

/* produce Knuth shuffle of 1:n */

void shuffle_in_place(IntegerVector idx)
{
    int i, j, s;
    int len = idx.size();
	
    for (i = 0; i < len; i++)
    {
	j = (int)floor(R::runif(0, 1) * (i + 1));
	if (j > i) j = 0;

	s = idx[i];
	idx[i] = idx[j];
	idx[j] = s;
    }
}

RcppExport SEXP rcor_permtest(SEXP matx, SEXP maty, SEXP tnorm, SEXP tests,
			      SEXP ogamma, SEXP alt, SEXP storeValues)
{
    int i, cnt = 0, tnorm_sel = as<int>(tnorm);
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
	
    for (i = 0; i < mat_x.nrow(); i++)
	perm[i] = i;

    double mean = 0, M2 = 0, delta;

    LogicalVector store(storeValues);

    int numGamma = (store[0] ? num_tests : 0);

    NumericVector permGamma(numGamma);

    switch(tnorm_sel)
    {
        case  1 : tnorm_fp = min_tnorm; break;
        case  2 : tnorm_fp = prod_tnorm; break;
        case  3 : tnorm_fp = lukasiewicz_tnorm; break;
        default : tnorm_fp = min_tnorm;
    }

    for (i = 0; i < num_tests; i++)
    {
	shuffle_in_place(perm);	
		
	get_sums(mat_x, mat_y, perm, tnorm_fp, &c, &d);
	gamma[0] = (fabs(c + d) <= DBL_EPSILON) ? 0 : (c - d) / (c + d);

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

	// computation of mean and variance according to Knuth/Welford

	delta = gamma[0] - mean;
	mean += (delta / (i + 1));
	M2 += (delta * (gamma[0] - mean));
    }

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
