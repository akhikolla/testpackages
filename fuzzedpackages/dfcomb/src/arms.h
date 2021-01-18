#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

typedef struct point {    /* a point in the x,y plane */
  double x,y;             /* x and y coordinates */
  double ey;              /* exp(y-ymax+YCEIL) */
  double cum;             /* integral up to x of rejection envelope */
  int f;                  /* is y an evaluated point of log-density */
  struct point *pl,*pr;   /* envelope points to left and right of x */
} POINT;

/* *********************************************************************** */

typedef struct envelope {  /* attributes of the entire rejection envelope */
  int cpoint;              /* number of POINTs in current envelope */
  int npoint;              /* max number of POINTs allowed in envelope */
  int *neval;              /* number of function evaluations performed */
  double ymax;             /* the maximum y-value in the current envelope */
  POINT *p;                /* start of storage of envelope POINTs */
  double *convex;          /* adjustment for convexity */
} ENVELOPE;

/* *********************************************************************** */

typedef struct funbag { /* everything for evaluating log density          */
  void *mydata;      /* user-defined structure holding data for density */
  double (*myfunc)(double x, void *mydata);
                     /* user-defined function evaluating log density at x */
} FUNBAG;

/* *********************************************************************** */

typedef struct metropolis { /* for metropolis step */
  int on;            /* whether metropolis is to be used */
  double xprev;      /* previous Markov chain iterate */
  double yprev;      /* current log density at xprev */
} METROPOLIS;

/* *********************************************************************** */

//#define RAND_MAX 2147483647      /* For Sun4 : remove this for some systems */
#define XEPS  0.00001            /* critical relative x-value difference */
#define YEPS  0.1                /* critical y-value difference */
#define EYEPS 0.001              /* critical relative exp(y) difference */
#define YCEIL 50.

typedef double (*urand)();

/* declarations for functions defined in this file */

int arms_simple (int ninit, double *xl, double *xr,
                 double (*myfunc)(double x, void *mydata), void *mydata,
                 int dometrop, double *xprev, double *xsamp, urand ur);

int arms (double *xinit, int ninit, double *xl, double *xr,
          double (*myfunc)(double x, void *mydata), void *mydata,
          double *convex, int npoint, int dometrop, double *xprev, double *xsamp,
          int nsamp, double *qcent, double *xcent, int ncent,
          int *neval, urand ur);

int initial (double *xinit, int ninit, double xl, double xr, int npoint,
	     FUNBAG *lpdf, ENVELOPE *env, double *convex, int *neval,
             METROPOLIS *metrop);

void sample(ENVELOPE *env, POINT *p, urand ur);

void invert(double prob, ENVELOPE *env, POINT *p);

int test(ENVELOPE *env, POINT *p, FUNBAG *lpdf, METROPOLIS *metrop, urand ur);

int update(ENVELOPE *env, POINT *p, FUNBAG *lpdf, METROPOLIS *metrop);

void cumulate(ENVELOPE *env);

int meet (POINT *q, ENVELOPE *env, METROPOLIS *metrop);

double area(POINT *q);

double expshift(double y, double y0);

double logshift(double y, double y0);

double perfunc(FUNBAG *lpdf, ENVELOPE *env, double x);
