// Title: Statistique de test de Brys-Hubert-Struyf MC-LR
// Ref. (book or article): Brys, G., Hubert, M. and Struyf, A. (2008), Goodness-of-fit tests based on a robust measure of skewness, 
//						   Computational Statistics, Vol. 23, Issue 3, pp. 429-442. 

#include <R.h>
#include "Rmath.h"

extern "C" {

  void stat16(double *x, int *xlen, double *level, int *nblevel, char **name, int *getname, double *statistic, int *pvalcomp, double *pvalue, double *critvalL, double *critvalR, int *usecrit, int *alter, int *decision, double *paramstat, int *nbparamstat) {

// If the test statistic can only be in category 3 or 4 (see just below), INDICATE the following line accordingly. Else, leave it commented.
// 0: two.sided=bilateral, 1: less=unilateral, 2: greater=unilateral, 3: bilateral test that rejects H0 only for large values of the test statistic, 
// 4: bilateral test that rejects H0 only for small values of the test statistic
    alter[0] = 3;
    
    int i, j=0, n=xlen[0];
    if (getname[0] == 1) {
// Here, INDICATE the name of your statistic
      const char *nom = "$T_{MC-LR}$";
// Here, INDICATE the number of parameters of your statistic
      nbparamstat[0] = 0;
// Here, INDICATE the default values of the parameters:
      if (name[0][0] == '1') { // To prevent writing in a non declared address (because maybe paramstat has not be initialized with a sufficient length since the correct value of nbparamstat[0] may be unkown yet).

     }
// The following 7 lines should NOT be modified
      const char *space = " ";
      while (nom[j] != '\0') {
	name[j][0] = nom[j];
	j++;
      }
      for (i=j;i<50;i++) name[i][0] = space[0];
      return;
    }

    if (n>3) {
// Computation of the value of the test statistic
      double mc_C_d(double *z, int n, double *eps, int *iter);
      void R_rsort (double* x, int n);
      double pchisq(double q, double df, int lower_tail, int log_p);
      double *x1, *x2, *x3;
      double w1, w2, w3;
	  int in1, in2, in3;
	  double *eps;
	  eps = new double[2];
	  eps[0] = 2.220446e-16;
	  eps[1] = 2.225074e-308;
	  int *iter;
	  iter = new int[2];
	  
      double omega1, omega2, omega3, invV11, invV12, invV13, invV21, invV22, invV23, invV31, invV32, invV33, vec1, vec2, vec3, statTMCLR;

      
      in1 = n;
      x1 = new double[n];
      for (i=0;i<=(n-1);i++) x1[i] = x[i];
      R_rsort(x1,n);
      //      double M;
      //      if ((n%2) == 0) M = (x1[n/2]+x1[n/2-1])/2.0; else M = x1[n/2]; // sample median   // PROBLEME ICI !!! M EST UTILISEE NULLE PART !!!

      if ((n%2) == 0) {
		in2 = n/2;
		in3 = n/2;
		x2 = new double[in2];
		x3 = new double[in3];
		for (i=0;i<=(n/2-1);i++) x2[i] = x1[i];
		for (i=(n/2);i<=(n-1);i++) x3[i-(n/2)] = x1[i];
      } else {
		in2 = n/2+1;
		in3 = n/2+1;
		x2 = new double[in2];
		x3 = new double[in3];
		for (i=0;i<=(n/2);i++) x2[i] = x1[i];
		for (i=(n/2);i<=(n-1);i++) x3[i-(n/2)] = x1[i];
      }

      iter[0] = 1000;
      iter[1] = 0;	  
      w1 = mc_C_d(x1,in1,eps,iter);
      iter[0] = 1000;
      iter[1] = 0;	  
      w2 = mc_C_d(x2,in2,eps,iter);
      iter[0] = 1000;
      iter[1] = 0;	  
      w3 = mc_C_d(x3,in3,eps,iter);
	  
	  // w1 = 0.0;
	  // w2 = 0.0;
	  // w3 = 0.0;
	  
      omega1 = 0.0;
      omega2 = 0.198828;
      omega3 = 0.198828;

      vec1 = w1-omega1;
      vec2 = -w2-omega2;
      vec3 = w3-omega3;

      invV11 = 0.8571890822945882;
      invV12 = -0.1051268907484579;
      invV13 = 0.1051268907484580;
      invV21 = -0.1051268907484579;
      invV22 = 0.3944817329840534;
      invV23 = -0.01109532299714422;
      invV31 = 0.1051268907484579;
      invV32 = -0.01109532299714422;
      invV33 = 0.3944817329840535;

      statTMCLR = (double)n*((vec1*invV11+vec2*invV21+vec3*invV31)*vec1 + (vec1*invV12+vec2*invV22+vec3*invV32)*vec2 + (vec1*invV13+vec2*invV23+vec3*invV33)*vec3);
	  
      statistic[0] = statTMCLR; // Here is the test statistic value

if (pvalcomp[0] == 1) {
	// If possible, computation of the p-value.
	#include "pvalues/pvalue16.cpp"
}

// We take the decision to reject or not to reject the null hypothesis H0
    for (i=0;i<=(nblevel[0]-1);i++) {
	if (usecrit[0] == 1) { // We use the provided critical values
	  if (statistic[0] > critvalR[i]) decision[i] = 1; else decision[i] = 0; // two.sided (but in this case only one possible critical value)
	  } else {
	    if (pvalue[0] < level[i]) decision[i] = 1; else decision[i] = 0; // We use the p-value
	  }
    }
    
// If applicable, we free the unused array of pointers
      delete[] x1;
      delete[] x2;
      delete[] x3;
	  delete[] eps;
	  delete[] iter;

    }

// We return
    return;
            
  }



/*
    Algorithm for the skewness estimator medcouple (MC)
    --------------------------------------------------
    ( originally matlabmc.c and also  mc/mcrsoft/spmc.c )
    From R package robustbase
*/

// // #include <stdlib.h>  // Pierre modified this on 22/12/2015
#include <math.h>

#include <inttypes.h>
// -> int64_t

/*
  including

 whimed_i(a,iw,n): the weighted high median of an array a of length n,
		   using the positive integer weights iw[].

 * which is in ./wgt_himed.c_templ
 *		 ~~~~~~~~~~~~~~~~~
*/

/* Includes the auxiliary function

   h_kern(a,b, ai,bi,ab, eps):	the values h(a,b)  needed to compute the mc
*/
static
double h_kern(double a, double b, int ai, int bi, int ab, double eps);

double whimed_i(double *a, int *iw, int n, double *acand, double *a_srt, int *iw_cand);



/* MM:	The tolerance  'eps1' and 'eps2' can now be passed from R;
 *	the original code had only one 'eps' for both and hardcoded
 *	   eps =  0.0000000000001;  (== 1e-13 )
 *
 * MK:  eps1: for (relative) "equality" checks
 *      eps2: used to check for over- and underflow, respectively
 *      therefore I suggest eps1 = DBL_EPS and eps2 = DBL_MIN
 */
double mc_C_d(double *z, int n, double *eps, int *iter)
{
/* NOTE:
    eps	  = c(eps1, eps2)
    iter := c(maxit, trace.lev)  as input
          = c(it, converged)     as output
*/
    int trace_lev = iter[1], it = 0;
    // Rboolean converged = TRUE;
	bool converged = TRUE;
    double medc; // "the" result
    static const double Large = DBL_MAX / 4.;

    if (n < 3) {
	medc = 0.; 
	iter[0] = it; /* to return */
	iter[1] = converged;
	return medc;
    }
    /* copy data before sort()ing in place, also reflecting it -- dealing with +-Inf.
       NOTE: x[0] "empty" so we can use 1-indexing below */
    double *x  = (double *) R_alloc(n+1, sizeof(double));
    x[0] = 0;
    for (int i = 0; i < n; i++) {
	double zi = z[i];
	x[i+1] = - ((zi == R_PosInf) ? Large :
		    (zi == R_NegInf ? -Large : zi));
    }

    R_rsort(&x[1], n); /* full sort */

    double xmed; // := median( x[1:n] ) = - median( z[0:(n-1)] ):
    if (n%2) { /* n even */
	xmed = x[(n/2)+1];
    }
    else { /* n  odd */
	int ind = (n/2);
	xmed = (x[ind] + x[ind+1])/2;
    }

    if (fabs(x[1] - xmed) < eps[0] * (eps[0] + fabs(xmed))) {
	medc = -1.; 
	iter[0] = it; /* to return */
	iter[1] = converged;
	return medc;
    } else if (fabs(x[n] - xmed) < eps[0] * (eps[0] + fabs(xmed))) {
	medc =	1.; 
	iter[0] = it; /* to return */
	iter[1] = converged;
	return medc;
    }
    /* else : median is not at the border ------------------- */

    if(trace_lev)
	Rprintf("mc_C_d(z[1:%d], trace_lev=%d): Median = %g (not at the border)\n",
		n, trace_lev, -xmed);

    int i,j;
    /* center x[] wrt median --> such that then  median( x[1:n] ) == 0 */
    for (i = 1; i <= n; i++)
	x[i] -= xmed;

    /* Now scale to inside [-0.5, 0.5] and flip sign such that afterwards
    *  x[1] >= x[2] >= ... >= x[n] */
    double xden = -2 * fmax2(-x[1], x[n]);
    for (i = 1; i <= n; i++)
	x[i] /= xden;
    xmed /= xden;
    if(trace_lev >= 2)
	Rprintf(" x[] has been rescaled (* 1/s) with s = %g\n", -xden);

    j = 1;
    double x_eps = eps[0] * (eps[0] + fabs(xmed));
    while (j <= n && x[j] > x_eps) { /* test relative to xmed */
	/* x1[j] = x[j]; */
	j++;
    }
    if(trace_lev >= 2)
	Rprintf("   x1[] := {x | x_j > x_eps = %g}    has %d (='j-1') entries\n",
		x_eps, j-1);
    i = 1;
    double *x2 = x+j-1; /* pointer -- corresponding to  x2[i] = x[j]; */
    while (j <= n && x[j] > -x_eps) { /* test relative to xmed */
	/* x1[j] = x[j]; */
        /* x2[i] = x[j]; */
        j++;
        i++;
    }
    /* now  x1[] := {x | x_j > -eps}  also includes the median (0) */
    if(trace_lev >= 2)
        Rprintf("'median-x' {x | -eps < x_i <= eps} has %d (= 'k') entries\n",
		i-1);
    int h1 = j-1, /* == size of x1[] == the sum of those two sizes above */
    /* conceptually,  x2[] := {x | x_j <= eps}   (which includes the median 0) */
	h2 = i + (n-j);// == size of x2[] == maximal size of whimed() arrays

    if(trace_lev)
	Rprintf("  now allocating 2+5 work arrays of size (1+) h2=%d each:\n", h2);
    /* work arrays for whimed_i() :  allocate *once* only !! */
    double *acand  = (double *) R_alloc(h2, sizeof(double)),
	   *a_srt  = (double *) R_alloc(h2, sizeof(double));

    int    *iw_cand= (int *)	R_alloc(h2, sizeof(int)),
    /* work arrays for the fast-median-of-table algorithm:
     *  currently still with  1-indexing */
	*left  = (int *) R_alloc((h2+1), sizeof(int)),
	*right = (int *) R_alloc((h2+1), sizeof(int)),
	*p     = (int *) R_alloc((h2+1), sizeof(int)),
	*q     = (int *) R_alloc((h2+1), sizeof(int));

    for (i = 1; i <= h2; i++) {
	left [i] = 1;
	right[i] = h1;
    }
    int64_t nr = ((int64_t) h1) * ((int64_t) h2), // <-- careful to *NOT* overflow
	knew = nr/2 +1;
    if(trace_lev >= 2)
	Rprintf(" (h1,h2, nr, knew) = (%d,%d, %.0f, %.0f)\n",
		h1,h2, (double)nr, (double)knew);

    double trial = -2./* -Wall */;
    double *work= (double *) R_alloc(n, sizeof(double));
    int	   *iwt = (int *)    R_alloc(n, sizeof(int));
    //    Rboolean IsFound = FALSE;
	bool IsFound = FALSE;
    int nl = 0,
	neq = 0;
    /* MK:  'neq' counts the number of observations in the
     *      inside the tolerance range, i.e., where left > right + 1,
     *      since we would miss those when just using 'nl-nr'.
     *      This is to prevent index overflow in work[] later on.
     *      left might be larger than right + 1 since we are only
     *      testing with accuracy eps_trial and therefore there might
     *      be more than one observation in the `tolerance range`
     *      between < and <=.
     */
    while (!IsFound && (nr-nl+neq > n) && it < iter[0])
    {
	int64_t sum_p, sum_q;
	it++;
	j = 0;
	for (i = 1; i <= h2; i++)
	    if (left[i] <= right[i]) {
		iwt[j] = right[i] - left[i]+1;
		int k = left[i] + (iwt[j]/2);
		work[j] = h_kern(x[k], x2[i], k, i, h1+1, eps[1]);
		j++;
	    }
	if(trace_lev >= 4) {
	    Rprintf(" before whimed(): work and iwt, each [0:(%d-1)]:\n", j);
	    if(j >= 100) {
		for(i=0; i < 90; i++) Rprintf(" %8g", work[i]);
		Rprintf("\n  ... ");
		for(i=j-4; i < j; i++)Rprintf(" %8g", work[i]);
		Rprintf("\n");
		for(i=0; i < 90; i++) Rprintf(" %8d", iwt [i]);
		Rprintf("\n  ... ");
		for(i=j-4; i < j; i++)Rprintf(" %8d", iwt [i]);
		Rprintf("\n");
	    } else { // j <= 99
		for(i=0; i < j; i++) Rprintf(" %8g", work[i]);
		Rprintf("\n");
		for(i=0; i < j; i++) Rprintf(" %8d", iwt [i]);
		Rprintf("\n");
	    }
	}
	trial = whimed_i(work, iwt, j, acand, a_srt, iw_cand);
	double eps_trial = eps[0] * (eps[0] + fabs(trial));
	if(trace_lev >= 3)
	    Rprintf("%2s it=%2d, whimed(*, n=%6d)= %8g ", " ", it, j, trial);

	j = 1;
	for (i = h2; i >= 1; i--) {
	    while (j <= h1 && h_kern(x[j],x2[i],j,i,h1+1,eps[1]) - trial > eps_trial) {
		// while (j <= h1 && h_kern(x[j],x2[i],j,i,h1+1,eps[1]) > trial) {
		if (trace_lev >= 5)
		    Rprintf("\nj=%3d, i=%3d, x[j]=%g, x2[i]=%g, h=%g",
			    j, i, x[j], x2[i],
			    h_kern(x[j],x2[i],j,i,h1+1,eps[1]));
		j++;
	    }
/* 	    for(; j <= h1; j++) { */
/* 		register double h = h_kern(x[j],x2[i],j,i,h1+1,eps[1]); */
/* 		if(h > trial) break; */
/* 	    } */
	    p[i] = j-1;
	}
	j = h1;
	for (i = 1, sum_p=0, sum_q=0; i <= h2; i++) {
	    while (j >= 1 && trial - h_kern(x[j],x2[i],j,i,h1+1,eps[1]) > eps_trial)
		// while (j >= 1 && h_kern(x[j],x2[i],j,i,h1+1,eps[1]) < trial)
		j--;
	    q[i] = j+1;

	    sum_p += p[i];
	    sum_q += j;/* = q[i]-1 */
	}

	if(trace_lev >= 3) {
	    if (trace_lev == 3)
		Rprintf("sum_(p,q)= (%.0f,%.0f)", (double)sum_p, (double)sum_q);
	    else { /* trace_lev >= 4 */
		Rprintf("\n%3s p[1:%d]:", "", h2);
		//		Rboolean lrg = h2 >= 100;
		bool lrg = h2 >= 100;
		int i_m = lrg ? 95 : h2;
		for(i = 1; i <= i_m; i++) Rprintf(" %2d", p[i]);
		if(lrg) Rprintf(" ...");
		Rprintf(" sum=%4.0f\n%3s q[1:%d]:", (double)sum_p, "", h2);
		for(i = 1; i <= i_m; i++) Rprintf(" %2d", q[i]);
		if(lrg) Rprintf(" ...");
		Rprintf(" sum=%4.0f\n", (double)sum_q);
	    }
	}

	if (knew <= sum_p) {
	    if(trace_lev >= 3)
		Rprintf("; sum_p >= kn\n");
	    for (i = 1, neq = 0; i <= h2; i++) {
		right[i] = p[i];
		if (left[i] > right[i]+1) neq += left[i]-right[i]-1;
	    }
	    nr = sum_p;
	}
	else { /* knew > sum_p */
	    IsFound = (knew <= sum_q); /* i.e. sum_p < knew <= sum_q */;

	    if(trace_lev >= 3)
		Rprintf("; s_p < kn ?<=? s_q: %s\n", IsFound ? "TRUE": "no");
	    if(IsFound) {
		medc = trial;
	    } else { /*	 knew > sum_q */
	        for (i = 1; i <= h2; i++) {
		    left[i] = q[i];
		    if (left[i] > right[i]+1) neq += left[i]-right[i]-1;
		}
		nl = sum_q;
	    }
	}
	R_CheckUserInterrupt();

    } /* end while loop */

    converged = IsFound || (nr-nl+neq <= n);
    if(!converged) {
	warning("maximal number of iterations (%d =? %d) reached prematurely\n",
		 iter[0], it);
	/* still: */
	medc = trial;
    }

    if (converged && !IsFound) { /* e.g., for  mc(1:4) : */
	j = 0;
	for (i = 1; i <= h2; i++) {
	    if (left[i] <= right[i]) {
		for (int k = left[i]; k <= right[i]; k++) {
		    work[j] = -h_kern(x[k],x2[i],k,i,h1+1,eps[1]);
		    j++;
		}
	    }
	}
	if(trace_lev)
	    Rprintf("  not found [it=%d,  (nr,nl) = (%d,%d)],"
		    " -> (knew-nl, j) = (%d,%d)\n",
		    it, nr, nl, knew-nl, j);
	/* using rPsort(work, n,k), since we don't need work[] anymore:*/
	rPsort(work, /* n = */ j, /* k = */ knew-nl-1);
	medc = - work[knew-nl-1];
    }

    if(converged && trace_lev >= 2)
	Rprintf("converged in %d iterations\n", it);

    iter[0] = it; /* to return */
    iter[1] = converged;

    return medc;

} /* end{ mc_C_d } */


static
double h_kern(double a, double b, int ai, int bi, int ab, double eps)
{
/*     if (fabs(a-b) <= DBL_MIN) */
    /* check for zero division and positive b */
    if (fabs(a-b) < 2.0*eps || b > 0)
	return sign((double)(ab - (ai+bi)));

    /* else */
    return (a+b)/(a-b);
}




double whimed_i(double *a, int *w, int n, double* a_cand, double *a_srt, int* w_cand)
{

/*
  Algorithm to compute the weighted high median in O(n) time.

  The whimed is defined as the smallest a[j] such that the sum
  of the weights of all a[i] <= a[j] is strictly greater than
  half of the total weight.

  Arguments:

  a: double array containing the observations
  n: number of observations
  w: array of (int/double) weights of the observations.
*/

    int n2, i, kcand;
    /* sum of weights: `int' do overflow when  n ~>= 1e5 */
    int64_t wleft, wmid, wright, w_tot, wrest;

    double trial;

    w_tot = 0;
    for (i = 0; i < n; ++i)
	w_tot += w[i];
    wrest = 0;

/* REPEAT : */
    do {
	for (i = 0; i < n; ++i)
	    a_srt[i] = a[i];
	n2 = n/2;/* =^= n/2 +1 with 0-indexing */
	rPsort(a_srt, n, n2);
	trial = a_srt[n2];

	wleft = 0;    wmid  = 0;    wright= 0;
	for (i = 0; i < n; ++i) {
	    if (a[i] < trial)
		wleft += w[i];
	    else if (a[i] > trial)
		wright += w[i];
	    else
		wmid += w[i];
	}
	/* wleft = sum_{i; a[i]	 < trial}  w[i]
	 * wmid	 = sum_{i; a[i] == trial}  w[i] at least one 'i' since trial is one a[]!
	 * wright= sum_{i; a[i]	 > trial}  w[i]
	 */
	kcand = 0;
	if (2 * (wrest + wleft) > w_tot) {
	    for (i = 0; i < n; ++i) {
		if (a[i] < trial) {
		    a_cand[kcand] = a[i];
		    w_cand[kcand] = w[i];	++kcand;
		}
	    }
	}
	else if (2 * (wrest + wleft + wmid) <= w_tot) {
	    for (i = 0; i < n; ++i) {
		if (a[i] > trial) {
		    a_cand[kcand] = a[i];
		    w_cand[kcand] = w[i];	++kcand;
		}
	    }
	    wrest += wleft + wmid;
	}
	else {
	    return trial;
	    /*==========*/
	}
	n = kcand;
	for (i = 0; i < n; ++i) {
	    a[i] = a_cand[i];
	    w[i] = w_cand[i];
	}
    } while(1);

} /* _WHIMED_ */

 
 
}
