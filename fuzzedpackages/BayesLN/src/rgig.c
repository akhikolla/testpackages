
//Adaption of the C code included in the ghyp R package



/*
 *
 * Ester Pantaleo and Robert B. Gramacy, 2010
 *
 * adapted from the C code in the monomvm package for R.
 */


#include "rgig.h"
#include <math.h>
#include <assert.h>
#include <R.h>
#include <Rmath.h>

#define ZTOL sqrt(DOUBLE_EPS)

/*
 * gig_y_gfn:
 *
 * evaluate the function that we need to find the root
 * of in order to construct the optimal rejection
 * sampler to obtain GIG samples
 */

double gig_y_gfn(double y, double m, double beta, double lambda)
{
  double y2, g;
  y2 = y * y;
  g = 0.5 * beta * y2 * y;
  g -= y2 * (0.5 * beta * m + lambda + 1.0);
  g += y * ((lambda - 1.0) * m - 0.5 * beta) + 0.5 * beta * m;
  return(g);
}

/* THIS FUNCTION HAS BEEN MODIFIED TO DEAL WITH GIG_Y_GFN (extra args) */

double zeroin_gig(ax,bx,f,tol, m, beta, lambda)	/* An estimate to the root  */
double ax;				/* Left border | of the range	*/
double bx;  				/* Right border| the root is seeked*/
/* Function under investigation	*/
double (*f)(double x, double m, double beta, double lambda);
double tol;				/* Acceptable tolerance	*/
double m;                               /* specific to gig_y_gfn */
double beta;                            /* specific to gig_y_gfn */
double lambda;                          /* specific to gig_y_gfn */
{
  double a,b,c;				/* Abscissae, descr. see above	*/
  double fa;				/* f(a)				*/
  double fb;				/* f(b)				*/
  double fc;				/* f(c)				*/

  a = ax;  b = bx;  fa = (*f)(a, m, beta, lambda);  fb = (*f)(b, m, beta, lambda);
  c = a;   fc = fa;

  for(;;)		/* Main iteration loop	*/
  {
    double prev_step = b-a;		/* Distance from the last but one*/
					/* to the last approximation	*/
    double tol_act;			/* Actual tolerance		*/
    double p;      			/* Interpolation step is calcu- */
    double q;      			/* lated in the form p/q; divi- */
  					/* sion operations is delayed   */
 					/* until the last moment	*/
    double new_step;      		/* Step at this iteration       */

    if( fabs(fc) < fabs(fb) )
    {                         		/* Swap data for b to be the 	*/
	a = b;  b = c;  c = a;          /* best approximation		*/
	fa=fb;  fb=fc;  fc=fa;
    }
    tol_act = 2.0*DOUBLE_EPS*fabs(b) + tol/2.0;
    new_step = (c-b)/2.0;

    if( fabs(new_step) <= tol_act || fb == (double)0 )
      return b;				/* Acceptable approx. is found	*/

    			/* Decide if the interpolation can be tried	*/
    if( fabs(prev_step) >= tol_act	/* If prev_step was large enough*/
	&& fabs(fa) > fabs(fb) )	/* and was in true direction,	*/
    {					/* Interpolatiom may be tried	*/
	register double t1,cb,t2;
	cb = c-b;
	if( a==c )			/* If we have only two distinct	*/
	{				/* points linear interpolation 	*/
	  t1 = fb/fa;			/* can only be applied		*/
	  p = cb*t1;
	  q = 1.0 - t1;
 	}
	else				/* Quadric inverse interpolation*/
	{
	  q = fa/fc;  t1 = fb/fc;  t2 = fb/fa;
	  p = t2 * ( cb*q*(q-t1) - (b-a)*(t1-1.0) );
	  q = (q-1.0) * (t1-1.0) * (t2-1.0);
	}
	if( p>(double)0 )		/* p was calculated with the op-*/
	  q = -q;			/* posite sign; make p positive	*/
	else				/* and assign possible minus to	*/
	  p = -p;			/* q				*/

	if( p < (0.75*cb*q-fabs(tol_act*q)/2.0)	/* If b+p/q falls in [b,c]*/
	    && p < fabs(prev_step*q/2.0) )	/* and isn't too large	*/
	  new_step = p/q;			/* it is accepted	*/
					/* If p/q is too large then the	*/
					/* bissection procedure can 	*/
					/* reduce [b,c] range to more	*/
					/* extent			*/
    }

    if( fabs(new_step) < tol_act ) {	/* Adjust the step to be not less*/
      if( new_step > (double)0 )	/* than tolerance		*/
	new_step = tol_act;
      else
	new_step = -tol_act;
    }

    a = b;  fa = fb;			/* Save the previous approx.	*/
    b += new_step;  fb = (*f)(b, m, beta, lambda);  /* Do step to a new approxim. */
    if( (fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0) )
    {                 			/* Adjust c for it to have a sign*/
      c = a;  fc = fa;                  /* opposite to that of b	*/
    }
  }

}


/*
 * rgig:
 *
 * a C implementation of the R code for rgig from
 * the ghyp v_1.5.2 package.
 * Modified to use it in the Gibbs sampler
 */

double rgig( const double lambda, const double chi, const double psi){
double samps;

  /* special case which is basically a gamma distribution */
  if((chi < ZTOL) & (lambda > 0.0)) {
     samps = rgamma(lambda, 2.0/psi);
    return samps;
  }

  /* special cases which is basically an inverse gamma distribution */
  if((psi < ZTOL) & (lambda < 0.0)) {
     samps = 1.0/rgamma(0.0-lambda, 2.0/chi);
    return samps;
  }


  /*
   * begin general purpose rgig code, which was basically
   * translated from the R function rgig in the ghyp package v_1.5.2
   */

  double alpha, beta, beta2, m, m1, lm1, lm12, upper, yM, yP, a, b, c, R1, R2, Y;
  int  need;

  alpha = sqrt(chi/psi);
  beta2 = psi*chi;
  beta = sqrt(psi*chi);
  lm1 = lambda - 1.0;
  lm12 = lm1*lm1;
  m = (lm1 + sqrt(lm12 + beta2))/beta;
  m1 = m + 1.0/m;

  upper = m;
  while (gig_y_gfn(upper, m, beta, lambda) <= 0) { upper *= 2.0; }

  yM = zeroin_gig(0.0, m, gig_y_gfn, ZTOL, m, beta, lambda);
  yP = zeroin_gig(m, upper, gig_y_gfn, ZTOL, m, beta, lambda);

  a = (yP - m) * pow(yP/m, 0.5 * lm1);
  a *=  exp(-0.25 * beta * (yP + 1.0/yP - m1));
  b = (yM - m) * pow(yM/m, 0.5 * lm1);
  b *= exp(-0.25 * beta * (yM + 1/yM - m1));
  c = -0.25 * beta * m1 + 0.5 * lm1 * log(m);


    need = 1;
    while (need) {
      R1 = unif_rand();
      R2 = unif_rand();

      Y = m + a * R2/R1 + b * (1.0 - R2)/R1;
      if (Y > 0.0) {
	  if (-log(R1) >= - 0.5 * lm1 * log(Y) + 0.25 * beta * (Y + 1.0/Y) + c) {
	    need = 0;
	  }
      }
    }
    samps = Y*alpha;


    return samps;
}
