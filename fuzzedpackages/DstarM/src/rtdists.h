#ifndef RTDISTS_H
#define RTDISTS_H

#include <stddef.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>


// Header to tie into R
#include <R.h>

#if __GNUC__ >= 3
#define  jv_noreturn  __attribute__ ((noreturn))
#define  jv_malloc  __attribute__ ((malloc))
#define  jv_fprintf  __attribute__ ((format (printf, 2, 3)))
#else
#define  jv_noreturn  /* no noreturn */
#define  jv_malloc  /* no malloc */
#define  jv_fprintf  /* no printf */
#endif

/* from "xmalloc.c" */

extern  void *xmalloc (size_t size) jv_malloc;
#define xnew(T,N) ((T *)xmalloc((N)*sizeof(T)))
extern  void *xrealloc (void *ptr, size_t newsize);
#define xrenew(T,OLD,N) ((T *)xrealloc(OLD,(N)*sizeof(T)))
extern  void  xfree (void *ptr);
extern  char *xstrdup (const char *s) jv_malloc;
extern  char *xstrndup (const char *s, size_t n) jv_malloc;

/* from "dataset.c" */

enum parameter_index {
	p_a,
	p_v,
	p_t0,
	p_d,
	p_szr,
	p_sv,
	p_st0,
	p_count
};

/* from "precision.c" */

extern double TUNE_PDE_DT_MIN;
extern double TUNE_PDE_DT_MAX;
extern double TUNE_PDE_DT_SCALE;

extern double TUNE_DZ;
extern double TUNE_DV;
extern double TUNE_DT0;

extern double TUNE_INT_T0;
extern double TUNE_INT_Z;

extern int precision_set;

extern void set_precision(double p);

/* from "parameters.c" */

#define LIMIT_a 10.0
#define LIMIT_v 50.0

/* from density.c  */

#define EPSILON 1e-6

#ifndef HAVE_FABS
double
fmax(double a, double b)
{
	return a > b ? a : b;
}
#endif

/* forward declarations from Rfastdm.c */

 // Forward declarations for R-callable functions
void dfastdm   (int *in_numvalues, double *in_params, double *in_RTs, double *in_precision, double *out_densities_a, double *out_densities_0);
void dfastdm_b (int *in_numvalues, double *in_params, double *in_RTs, double *in_precision, int *boundary, double *out_densities);

// Reimplementing these for R compatibility
extern void params_write(double para[p_count], double zr, double precision, int n);
extern void params_check(double para[p_count], double zr);

void params_write(double para[p_count], double zr, double precision, int n)
{
	Rprintf("a = %g\n", para[p_a]);
	Rprintf("zr = %g\n", zr);
	Rprintf("v = %g\n", para[p_v]);
	Rprintf("t0 = %g\n", para[p_t0]);
	Rprintf("d = %g\n", para[p_d]);
	Rprintf("szr = %g\n", para[p_szr]);
	Rprintf("sv = %g\n", para[p_sv]);
	Rprintf("st0 = %g\n", para[p_st0]);
	Rprintf("precision = %g\n", precision);
	if (n > 0) {  Rprintf("n = %d\n", n);}
}

void params_check(double para[p_count], double zr)
{
	if (para[p_a] <= 0)                     { Rf_error ("error: invalid parameter a=%g\n", para[p_a]); }
	if (para[p_szr] < 0 || para[p_szr] > 1) { Rf_error ("error: invalid parameter szr=%g\n", para[p_szr]); }
	if (para[p_st0] < 0)                    { Rf_error ("error: invalid parameter st0=%g\n", para[p_st0]); }
	if (para[p_sv] < 0)                     { Rf_error ("error: invalid parameter sv=%g\n",	para[p_sv]); }
	if (para[p_t0] - fabs(0.5*para[p_d]) - 0.5*para[p_st0] < 0) { Rf_error ("error: invalid parameter combination t0=%g, d=%g, st0=%g\n", para[p_t0], para[p_d], para[p_st0]); }
	if (zr - 0.5*para[p_szr] <= 0)          { Rf_error ("error: invalid parameter combination zr=%g, szr=%g\n", zr, para[p_szr]); }
	if (zr + 0.5*para[p_szr] >= 1)          { Rf_error ("error: invalid parameter combination zr=%g, szr=%g\n",	zr, para[p_szr]); }

#ifdef ADDITIONAL_BOUNDS
    if (para[p_a] > 2)                      { Rf_error ("Additional Bounds error: invalid parameter a=%g\n", para[p_a]); }
    if (para[p_sv] > 2)                     { Rf_error ("Additional Bounds error: invalid parameter sv=%g\n", para[p_sv]); }
	if (para[p_v] < -5 || para[p_v] > 5)    { Rf_error ("Additional Bounds error: invalid parameter v=%g\n", para[p_v]); }
#endif
}


/* from density.c  */
// replaced isinf() with !R_FINITE()

static int
imax(int a, int b)
{
	return a > b ? a : b;
}


struct para {
	double t, a, zr, v, st0, szr, sv;
};

struct function {
	double (*f) (double x, void *data);
	void *data;
};

static double
integrate(struct function *F, double a, double b, double step_width)
{
	double width = b-a;
	int N = imax(4, (int) (width / step_width));
	double step = width / N;
	double x;
	double result = 0;

	for(x = a+0.5*step; x < b; x += step) {
		result += step * F->f(x, F->data);
	}
	return result;
}


static double
g_minus_small_time(double t, double zr, int N)
{
	int i;
	double sum = 0;

	for(i = -N/2; i <= N/2; i++) {
		double d = 2*i + zr;
		sum += exp(-d*d / (2*t)) * d;
	}

	return sum / sqrt(2*M_PI*t*t*t);
}

static double
g_minus_large_time(double t, double zr, int N)
{
	int i;
	double sum = 0;

	for(i = 1; i <= N; i++) {
		double d = i * M_PI;
		sum += exp(-0.5 * d*d * t) * sin(d*zr) * i;
	}

	return sum * M_PI;
}

static double
g_minus_no_var(double t, double a, double zr, double v)
{
	int N_small, N_large;
	double simple, factor, eps;
	double ta = t/(a*a);

	factor = exp(-a*zr*v - 0.5*v*v*t) / (a*a);
	if (!R_FINITE(factor)) {
		return 0;
	}
	eps = EPSILON / factor;

	N_large = (int)ceil(1/(M_PI*sqrt(t)));
	if (M_PI*ta*eps < 1) {
		N_large = imax(N_large,
					   (int)ceil(sqrt(-2*log(M_PI*ta*eps) / (M_PI*M_PI*ta))));
	}

	if (2*sqrt(2*M_PI*ta)*eps < 1) {
		N_small = (int)ceil(fmax(sqrt(ta) + 1,
								 2 + sqrt(-2*ta*log(2*eps*sqrt(2*M_PI*ta)))));
	} else {
		N_small = 2;
	}

	if (N_small < N_large) {
		simple = g_minus_small_time(t/(a*a), zr, N_small);
	} else {
		simple = g_minus_large_time(t/(a*a), zr, N_large);
	}
	return factor * simple;
}

static double
integral_v_g_minus(double zr, void *data)
{
	//struct para *P = data;
	struct para *P = reinterpret_cast<struct para*>(data);
	double t = P->t;
	double a = P->a;
	double v = P->v;
	double sv = P->sv;
	int N_small, N_large;
	double simple, factor, eps;
	double ta = t/(a*a);

	factor = 1 / (a*a * sqrt(t * sv*sv + 1)) * exp(-0.5 * (v*v*t + 2*v*a*zr - a*zr*a*zr*sv*sv) / (t*sv*sv+1));
	if (!R_FINITE(factor)) {
		return 0;
	}
	eps = EPSILON / factor;

	if (P->sv == 0) {
		return g_minus_no_var(P->t, P->a, zr, P->v);
	}

	N_large = (int)ceil(1 / (M_PI*sqrt(t)));
	if (M_PI*ta*eps < 1) {
		N_large = imax(N_large,
					   (int)ceil(sqrt(-2*log(M_PI*ta*eps) / (M_PI*M_PI*ta))));
	}

	if (2*sqrt(2*M_PI*ta)*eps < 1) {
		N_small = (int)ceil(fmax(sqrt(ta)+1,
								 2+sqrt(-2*ta*log(2*eps*sqrt(2*M_PI*ta)))));
	} else {
		N_small = 2;
	}

	if (N_small < N_large) {
		simple = g_minus_small_time(t/(a*a), zr, N_small);
	} else {
		simple = g_minus_large_time(t/(a*a), zr, N_large);
	}
	return factor * simple;
}

static double
integral_z_g_minus(double t, void *data)
{
	//struct para *P = data;
	struct para *P = reinterpret_cast<struct para*>(data);
	double res;

	if (t <= 0) return 0;

	P->t = t;
	if (P->szr == 0) {
		res = integral_v_g_minus(P->zr, P);
	} else {
		struct function F = {
			integral_v_g_minus,
			P
		};
		res = integrate(&F, P->zr - .5*P->szr, P->zr + .5*P->szr,
						TUNE_INT_Z) / P->szr;
	}
	return res;
}

static double
integral_t0_g_minus(double t, void *data)
{
	//struct para *P = data;
	struct para *P = reinterpret_cast<struct para*>(data);
	double res;

	if (P->st0 == 0) {
		res = integral_z_g_minus(t, P);
	} else {
		struct function F = {
			integral_z_g_minus,
			P
		};
		res = integrate(&F, t - .5*P->st0, t + .5*P->st0, TUNE_INT_T0) / P->st0;
	}
	return res;
}

double
g_minus(double t, const double *para)
{
	struct para P;

	P.a = para[p_a];
	P.zr = para[p_count];
	P.v = para[p_v];
	P.szr = para[p_szr];
	P.sv = para[p_sv];
	P.st0 = para[p_st0];

	return integral_t0_g_minus(t - para[p_t0] - 0.5*para[p_d], &P);
}

double
g_plus(double t, const double *para)
{
	struct para P;

	P.a = para[p_a];
	P.zr = 1 - para[p_count];
	P.v = -para[p_v];
	P.szr = para[p_szr];
	P.sv = para[p_sv];
	P.st0 = para[p_st0];

	return integral_t0_g_minus(t - para[p_t0] + 0.5*para[p_d], &P);
}

/* from "precision.c" */

double  TUNE_DZ;
double  TUNE_DV;
double  TUNE_DT0;

double  TUNE_PDE_DT_MIN = 1e-6;
double  TUNE_PDE_DT_MAX = 1e-6;
double  TUNE_PDE_DT_SCALE = 0.0;

double  TUNE_INT_T0;
double  TUNE_INT_Z;

int  precision_set = 0;


void
set_precision (double p)
/* Try to achieve an accuracy of approximately 10^{-p} for the CDF.  */
{
	TUNE_PDE_DT_MIN = pow(10, -0.400825*p-1.422813);
	TUNE_PDE_DT_MAX = pow(10, -0.627224*p+0.492689);
	TUNE_PDE_DT_SCALE = pow(10, -1.012677*p+2.261668);
	TUNE_DZ = pow(10, -0.5*p-0.033403);
	TUNE_DV = pow(10, -1.0*p+1.4);
	TUNE_DT0 = pow(10, -0.5*p-0.323859);

	TUNE_INT_T0 = 0.089045 * exp(-1.037580*p);
	TUNE_INT_Z = 0.508061 * exp(-1.022373*p);

	precision_set = 1;
}

/*from Rfastdm2.c*/
#define MAX_VALUES 1000000

// Order of parameters (note: doesn't include zr..)
// enum parameter_index {
//	p_a,       // Boundary separation
//	p_v,       // Mean of the drift
//	p_t0,      // Non-decision time
//	p_d,       // Difference between boundaries of non-decision time IS THIS NEW AS WELL???
//	p_szr,     // width of zr distribution
//	p_sv,      // standard deviation of v distribution
//	p_st0,     // width of t0 distribution
//	p_count    // Num parameters but often used internally to refer to zr (...!?)
//};
//  zr         // Mean of diffusion starting point relative to boundaries

// Global variables for common stuff (nasty shortcut, used more to ensure correctness between functions)
double g_precision;
double g_params[7];  // An array of 8 parameters is passed in from R, but internally the code uses 7 + zr
double g_zr;
int    g_num_values;

// Internal helper function to set up common elements for each function
void _setup (double num_values, double *params, double precision)
{
#ifdef SHOW_INFO
  #ifdef HAVE_LIBBLAS
    Rprintf ("Note: Using LIBBLAS.\n");
  #else
    Rprintf ("Note: Not using LIBBLAS.\n");
  #endif
#endif

#ifdef FASTDM_DEBUG
    Rprintf ("In _setup.\n");
#endif
    // Set precision
    g_precision = precision;
    set_precision (g_precision);

    // Set up the parameters
    for (int i=0; i<7;i++)
    {
        g_params[i] = params[i];
    }
    g_zr = params[7];

#ifdef SHOW_INFO
    Rprintf ("Writing out parameters.\n");
    params_write(g_params, g_zr, g_precision, 0);
#endif

	params_check(g_params, g_zr);

    // Get number of values we're calculating
    g_num_values = (int)num_values;
    if ((g_num_values <= 0) || (g_num_values > MAX_VALUES))
    {
        Rf_error ("Number of values requested is either <= 0 or exceeds maximum of %d\n", MAX_VALUES);
    }
}

// R-callable probability density function (PDF) for fastdm
void dfastdm   (int *in_numvalues, double *in_params, double *in_RTs, double *in_precision, double *out_densities_a, double *out_densities_0)
{
#ifdef FASTDM_DEBUG
    Rprintf ("In dfastdm.\n");
#endif

    _setup (*in_numvalues, in_params, *in_precision);

    double *RTs = in_RTs;

    // Coercing params into format g_plus/g_minus wants
    double params[8];
    for (int i=0; i<7; i++) { params[i] = g_params[i]; }
    params[7] = g_zr;

    // Loop through each RT and add to the output vector
    for (int i = 0; i < g_num_values; i++)
    {
        out_densities_a[i] = g_plus(RTs[i], params);
        out_densities_0[i] = -g_minus(RTs[i], params);
#ifdef SHOW_INFO
        Rprintf ("Calculating density for RT[%2d] = %3.4g: upper: %.4g, lower: %.4g.\n", i, RTs[i], out_densities_a[i], out_densities_0[i]);
#endif
	}
}

// R-callable probability density function (PDF) for fastdm - pass boundary to retrieve (1 = lower, 2 = upper)
void dfastdm_b (int *in_numvalues, double *in_params, double *in_RTs, double *in_precision, int *boundary, double *out_densities)
{
#ifdef FASTDM_DEBUG
    Rprintf ("In dfastdm.\n");
#endif

    _setup (*in_numvalues, in_params, *in_precision);

    int bound = *boundary;
    if ((bound < 1) || (bound > 2)) Rf_error ("Error: invalid boundary!\n");

    // Coercing params into format g_plus/g_minus wants
    double params[8];
    for (int i=0; i<7; i++) { params[i] = g_params[i]; }
    params[7] = g_zr;

    if (bound == 2)  // Calc upper
    {
        for (int i = 0; i < g_num_values; i++)
        {
            out_densities[i] = g_plus(in_RTs[i], params);
    	}
    }
    else // Calc lower
    {
        for (int i = 0; i < g_num_values; i++)
        {
            out_densities[i] = -g_minus(in_RTs[i], params);
    	}
    }
}



#endif
