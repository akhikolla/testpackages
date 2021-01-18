/***************************************************************************
error integral is straight from Press et al. and
Chi2-related algorithms are straight from Abramovitz and Stegun's Handbook

inverse cumulative normal is...
//=====================================================
 * ndtri code
"You can use the algorithm, including any of the computer implementations
listed in the section Computer implementations, for whatever purpuse you want,
but please show common courtesy and give credit where credit is due."
 * Lower tail quantile for standard normal distribution function.
 *
 * This function returns an approximation of the inverse cumulative
 * standard normal distribution function.  I.e., given P, it returns
 * an approximation to the X satisfying P = Pr{Z <= X} where Z is a
 * random variable from the standard normal distribution.
 *
 * The algorithm uses a minimax approximation by rational functions
 * and the result has a relative error whose absolute value is less
 * than 1.15e-9.
 *
 * Author:      Peter J. Acklam
 * Time-stamp:  2002-06-09 18:45:44 +0200
 * E-mail:      jacklam@math.uio.no
 * WWW URL:     http://www.math.uio.no/~jacklam
 *
 * C implementation adapted from Peter's Perl version
 *
source says that there is one more accurate algo:
"The URL http://lib.stat.cmu.edu/apstat/241 contains FORTRAN code which was
published in Applied Statistics, vol. 37, pp. 477-484, 1988.

The FORTRAN code contains two algorithms, one which is accurate to 1 part in 10^7
and one which is accurate to 1 part in 10^16."

further

"The Cephes library at the URL http://www.netlib.org/cephes/index.html
contains the C program ndtri.c.

I don't know anything about the accuracy of this algorithm."
====================================================

"@" F. Rousset 1997-2006 for whatever remains (!)
 ***************************************************************************/

#include <cmath> // rajout pour Dev-C++ ...
#include <iostream>
#include <limits>
#ifdef COMPATIBILITYRCPP
#include <Rcpp.h>
#endif
#include "genepop.h"
#include "proba.h"
#include "tools.h"

using namespace std;
double erfcc(double x);
double erfcc(double x) {     //conplementary error integral; erfcc in Press et al
    double t, z, val;
    z = fabs( x );
    t = 1.0 / ( 1.0 + 0.5 * z );
    val = t * exp( -z * z - 1.26551223 + t *
            ( 1.00002368 + t *( 0.37409196 + t *
             ( 0.09678418 + t *( -0.18628806 + t *
              ( 0.27886807 + t *( -1.13520398 + t *
               ( 1.48851587 + t *( -0.82215223 + t *0.17087277)))))))));
if( x >= 0.0 )	return val;
return 2.0 - val;
}

double ndtr(double x) {return 0.5 * erfcc(-x/sqrt(2));}



//double ndtri(double x) {return gsl_cdf_gaussian_Pinv(x,1.);}

static const double a[] =
{
	-3.969683028665376e+01,
	 2.209460984245205e+02,
	-2.759285104469687e+02,
	 1.383577518672690e+02,
	-3.066479806614716e+01,
	 2.506628277459239e+00
};

static const double b[] =
{
	-5.447609879822406e+01,
	 1.615858368580409e+02,
	-1.556989798598866e+02,
	 6.680131188771972e+01,
	-1.328068155288572e+01
};



double ndtri(double p)
{
	double q, r;

  static const double c[] =
    {
    -7.784894002430293e-03,
    -3.223964580411365e-01,
    -2.400758277161838e+00,
    -2.549732539343734e+00,
    4.374664141464968e+00,
    2.938163982698783e+00
    };
  
  static const double d[] =
    {
    7.784695709041462e-03,
    3.224671290700398e-01,
    2.445134137142996e+00,
    3.754408661907416e+00
    };
  
	if (p < 0 || p > 1)
	{
        noR_cout<<"\nCumulative inverse normal (ndtri) argument out of range\n";
		return 0.0;
	}
	else if (p == 0)
	{
		return -numeric_limits<double>::infinity();
	}
	else if (p == 1)
	{
		return numeric_limits<double>::infinity();
	}
	else if (p < 0.02425)
	{
		/* Rational approximation for lower region */
		q = sqrt(-2*log(p));
		return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
			((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
	}
	else if (p > 0.97575)
	{
		/* Rational approximation for upper region */
		q  = sqrt(-2*log(1-p));
		return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
			((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
	}
	else
	{
		/* Rational approximation for central region */
    		q = p - 0.5;
    		r = q*q;
		return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
			(((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
	}
}


//------------------------------------------------------------------------------
//--------------------------- Fonction grosDDL() -------------------------------
//------------------------------------------------------------------------------
void grosDDL(float nu, float chi, float &pchi){

float xx, xy;
double zxx, tt, qxx;
    //AbramovitzS 26.4.14
    xx = float(pow(double(chi)/ nu, double(1.0)/3) - (1. - 2. /(9 * nu))); //correction misplaced parenthesis 06/2010
    xx = xx / sqrt(2 / (9 * nu));
    if (xx > 5) pchi = -1;
    else {
        // approximation by the normal CDF
        xy = 0;
        if (xx < 0) xy = 1;
        xx = fabs(xx);
        zxx = double(exp(-xx * xx / 2)) / sqrt(2 * 3.1415926);
        tt = 1 / (1 + 0.2316419 * xx);
        qxx = zxx * (0.31938153 * tt - 0.356563782 * pow(double(tt), 2) + 1.781477937 * pow(double(tt), 3) - 1.821255978 * pow(double(tt), 4) + 1.330274429 * pow(double(tt), 5));
        if (xy == 1) pchi = 1 - float(qxx);
        else {

            pchi = float(qxx);
        }
    }

}



//------------------------------------------------------------------------------
//---------- chi2 code using algos form Numerical recipes 3rd ed. --------------
//------------------------------------------------------------------------------

/*-----------ln(gamma(xx))----------------------------*/
double gammln(double xx) {
  double x,tmp,y,ser;
  static const double cof[14]={57.1562356658629235,-59.5979603554754912,
                               14.1360979747417471,-0.491913816097620199,.339946499848118887e-4,
                                 .465236289270485756e-4,-.983744753048795646e-4,.158088703224912494e-3,
                               -.210264441724104883e-3,.217439618115212643e-3,-.164318106536763890e-3,
                                .844182239838527433e-4,-.261908384015814087e-4,.368991826595316234e-5};
  
  y=x=xx;
  tmp = x+5.24218750000000000; //Rational 671/128.
  tmp = (x+0.5)*log(tmp)-tmp;
  ser = 0.999999999999997092;
  for (int j=0;j<14;j++) ser += cof[j]/++y;
  return tmp+log(2.5066282746310005*ser/x);
}
/*-----------------------------------------------------------------------------*/

double gser(const double a, const double x) {
  const double EPS = numeric_limits<double>::epsilon();
  double sum,del,ap;
  double gln=gammln(a);
  ap=a;
  del=sum=1.0/a;
  for (;;) {
    ++ap;
    del *= x/ap;
    sum += del;
    if (fabs(del) < fabs(sum)*EPS) {
      return sum*exp(-x+a*log(x)-gln);
    }
  }
}

double gcf(const double a, const double x) {
  const double EPS = numeric_limits<double>::epsilon();
  const double FPMIN = numeric_limits<double>::min()/EPS;
  double an,b,c,d,del,h;
  double gln=gammln(a);
  b=x+1.0-a; 
  c=1.0/FPMIN;
  d=1.0/b;
  h=d;
  for (int i=1;;i++) { 
    an = -i*(i-a);
    b += 2.0;
    d=an*d+b;
      if (fabs(d) < FPMIN) d=FPMIN;
      c=b+an/c;
      if (fabs(c) < FPMIN) c=FPMIN;
      d=1.0/d;
      del=d*c;
      h *= del;
      if (fabs(del-1.0) <= EPS) break;
  }
  return exp(-x+a*log(x)-gln)*h; 
}

void chi2(float &pchi, float nu, float chi) {
  if (chi == 0.0) pchi= 1.0;
  else if (nu >= 100) grosDDL(nu, chi, pchi);
  else if (chi < nu+1.0) pchi= 1.0-gser(0.5*nu,0.5*chi); // Use the series representation.
  else pchi= gcf(0.5*nu,0.5*chi); // Use the continued fraction representation.
}

// legacy version of chi2(), of uncertain origin (tracing back to first BASIC implementation?)
void old_chi2 (float &pchi, float nu, float chi){
  float mm, sum;
  int j, i, cc;
  if (nu > 100) grosDDL(nu, chi, pchi); //Lance la fonction grosDDL()
  else {
    if (chi > 200) pchi = -1;
    else {
      cc = int(nu) / 2;
      mm = chi / 2;
      sum = 1;
      for (j = 1; j <= cc - 1; j++) {
        i = cc - j;
        sum = sum * mm / float(i) + 1;
      }
      pchi = exp(-mm) * sum;
    }
  }
}


