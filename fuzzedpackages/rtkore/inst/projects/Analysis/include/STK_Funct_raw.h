/*--------------------------------------------------------------------*/
/*  Copyright (C) 2004-2016  Serge Iovleff, Université Lille 1, Inria
    
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.
    
    You should have received a copy of the GNU Lesser General Public
    License along with this program; if not, write to the
    Free Software Foundation, Inc.,
    59 Temple Place,
    Suite 330,
    Boston, MA 02111-1307
    USA
    
    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
*/

/*
 * Project:  Analysis
 * Purpose:  raw mathematical functions
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Funct_raw.h
 *  @brief In this file we declare raw the functions.
 * 
 *  raw functions are generic functions that can be used in various part
 *  of the STatistiK project. No test is done about the arguments.
 **/

#ifndef STK_FUNCT_RAW_H
#define STK_FUNCT_RAW_H

#include "STK_Funct_gamma.h"
#include "STK_Algo.h"
#include "STK_Funct_Util.h"

namespace STK
{

namespace Funct
{
// forward declaration
Real beta_pdf_raw(Real const& x, Real const& a, Real const& b);

Real binomial_pdf_raw(Real const& x, Real const& n, Real const& p);
Real binomial_pdf_raw(int x, int n, Real const& p);
Real binomial_lpdf_raw(Real const& x, Real const& n, Real const& p);
Real binomial_lpdf_raw(int x, int n, Real const& p);

Real poisson_pdf_raw(Real const& x, Real const& lambda);
Real poisson_pdf_raw(int x, Real const& lambda);
Real poisson_lpdf_raw(Real const& x, Real const& lambda);
Real poisson_lpdf_raw(int x, Real const& lambda);

Real erf_raw(Real const& a);
Real erfc_raw(Real const& a);
Real normal_cdf_raw(Real const& x);
Real normal_pdf_raw(Real const& x);

Real psi_raw(Real x);

/** @ingroup Analysis
 *  @brief compute the partial deviance
 *  \f$ d_0(a,b) = a\log(a/b)+ b - a \f$.
 **/
inline Real dev0(Real const& a, Real const& b)
{
  // special values
  if (a == 0.) return b;
  if (b == a)  return 0.;
  // ratio
  Real v = (a-b)/(a+b);
  // small v -> use dl
  if (std::abs(v)<0.125)
  {
    Real sum  = (a-b)*v, ej = 2*a*v, term;
    v = v*v;
    int n =0;
    do
    { sum += (term = ((ej *= v) / ( ((++n)<<1) +1 ) ));}
    while (std::abs(term) > std::abs(sum) * Arithmetic<Real>::epsilon());
    return sum;
  }
  // else compute directly the result
  return (a*log(double(a/b))+b-a);
}

/** @ingroup Analysis
 *  Compute the function
 *  \f[ B_1(a,b,x) = \frac{ x^{a} (1-x)^{b}}{B(a,b)} \f]
 *  using  the partial deviance \f$ (a+b) * (p*log(x/p)+q*log((1-x)/q)) \f$.
 *  @param a,b,x parameters of the generalized beta
 *  @param xm1 true if @e x is to be taken as @e 1-x
 **/
inline Real b1(Real const& a, Real const& b, Real const& x, bool xm1)
{
  if (x == 0) return 0;
  if (x == 1) return 0;
  Real s = a+b, sx = s*x, sy = s*(1.-x);
  return ( std::exp(- Const::_LNSQRT2PI_
                    + 0.5 * ( a<b ? std::log(a) + log1p(-a/s) : std::log(b) + log1p(-b/s))
                    + (lgammaStirlingError(s)-lgammaStirlingError(a)-lgammaStirlingError(b))
                    - (xm1 ? dev0(a, sy)+dev0(b, sx) : dev0(a, sx)+dev0(b, sy))
                   )
        );
}

/** @ingroup Analysis
 *  @brief compute the partial deviance \f$g_0(x) = x\log(x)+ 1 - x\f$.
 */
inline Real g0(Real const& x)
{
  // special values
  if (x == 0.) return 1.;
  if (x == 1.) return 0.;
  // general case
  Real v = (x-1.)/(x+1.);
  // if 7/9 < x < 9/7  use Taylor serie of log((1+v)/(1-v))
  if (v < 0.125)
  {
    Real sum  = (x-1)*v, ej = 2*x*v, term;
    v = v*v;
    int n =0;
    do
    { sum += (term = ( (ej *= v) /( ((++n)<<1) +1)));}
    while (std::abs(term) > std::abs(sum) * Arithmetic<Real>::epsilon());
    return sum;
  }
  // else compute directly
  return x*Real(log(double(x)))-x + 1.;
}

/** @ingroup Analysis
 *  @brief compute the function \f$ \log(1+x) \f$.
 *  @param x value to evaluate the function
 **/
inline Real log1p(Real const& x)
{
  // trivial values
  if (x == 0.)  return 0.;
  if (x == -1.) return(-Arithmetic<Real>::infinity());
  // small values
  if (std::abs(x) < .375)
  {
    // create functor and compute the alternate serie
    Serielog1p f(x);
    return x * (1. - x * sumAlternateSerie(f));
  }
  // large values : compute directly the result
  return log(1. + x);
}

/** @ingroup Analysis
 *  @brief compute the function \f$ \exp(x)-1 \f$.
 **/
inline Real expm1(Real const& x)
{
  // small values of x
  if (std::abs(x) < Const::_LN2_)
  {
    // a + 1 != 1 -> use compute Taylor serie else use first order
    // Taylor expansion  S = x + (1/2!) x^2 + (1/3!) x^3 + ...
   if (std::abs(x) > Arithmetic<Real>::epsilon())
    {
      Real term = x, sum =x, n =1.;
      do
      { sum += (term *= x/(++n));}
      while (std::abs(term) > std::abs(sum) * Arithmetic<Real>::epsilon()) ;
      return sum ;
    }
    else return (x / 2. + 1.) * x;
  }
  return exp(double(x))-1.;
}

/** @ingroup Analysis
 *  @brief Compute the beta density function.
 *
 *  Compute the function (beta pdf)
 *  \f[ B(a,b,x) = \frac{ x^{a-1} (1-x)^{b-1}}{B(a,b)}. \f]
 *  @param a,b,x parameters of the beta density
 */
inline Real beta_pdf_raw(Real const& x, Real const& a, Real const& b)
{
  // trivial cases
  if (x < 0 || x > 1) return(0);
  if (x == 0)
  {
    if(a > 1) return(0);
    if(a < 1) return(Arithmetic<Real>::infinity());
    return(b); // a == 1
  }
  if (x == 1)
  {
    if(b > 1) return(0);
    if(b < 1) return(Arithmetic<Real>::infinity());
    return(a); // b == 1
  }
  // limit cases with x \in (0,1)
  if (a == 0 || b == 0) { return(0);}
  // general case
  Real s = a+b, sx = s*x, sy = s*(1.-x);
  return ( std::exp(- Const::_LNSQRT2PI_
                    + 0.5 *( a<b ? std::log(a) + log1p(-a/s) : std::log(b) + log1p(-b/s))
                    + lgammaStirlingError(s)-lgammaStirlingError(a)-lgammaStirlingError(b)
                    - dev0(a, sx) - dev0(b, sy)
                    - std::log(x) - log1p(-x)
                   )
         );
}

/** @ingroup Analysis
 *  @brief Compute the generalized binomial probability mass function.
 *  Compute the function
 *  \f[ B(n,p,x) = p^{x} (1-p)^{n-x} \frac{\Gamma(n+1)}{\Gamma(x+1)\Gamma(n-x+1)} \f].
 *  The function assume that there is no NA or infinite values and
 *  that <em>0<=p<=1, n>=0, 0<=x<=n</em>
 *  @param x,n,p parameters of the binomial density
 */
inline Real binomial_pdf_raw(Real const& x, Real const& n, Real const& p)
{
  // trivial and limit cases
  if (p == 0) return ((x == 0) ? 1 : 0);
  if (p == 1) return ((x == n) ? 1 : 0);
  if (x == 0) return ((n == 0) ? 1 : std::exp(n*log1p(-p)) );
  if (x == n) return std::pow(p,n);
  // other cases
  Real y = n-x;
  return ( std::exp(- Const::_LNSQRT2PI_ - 0.5 *( std::log(x) + log1p(-x/n))
                    + lgammaStirlingError(n)-lgammaStirlingError(x)-lgammaStirlingError(y)
                    - dev0(x, n*p) - dev0(y, n- n*p)
                   )
        );
}

/** @ingroup Analysis
 *  @brief Compute the binomial probability mass function.
 *
 *  Compute the function
 *  \f[ B(n,p,x) = p^{x} (1-p)^{n-x} \binom{n}{x} \f].
 *  The function assume that there is no NA or infinite values and
 *  that <em>0<=p<=1, n>=0, 0<=x<=n</em>
 *  @param x,n,p parameters of the binomial density
 */
inline Real binomial_pdf_raw(int x, int n, Real const& p)
{ return binomial_pdf_raw(Real(x), Real(n), p);}

/** @ingroup Analysis
 *  @brief Compute the generalized binomial log-probability mass function.
 *  Compute the function
 *  \f[ B(n,p,x) = x\log(p) (n-x)\log(1-p)
 *      \log\left(\frac{\Gamma(n+1)}{\Gamma(x+1)\Gamma(n-x+1)}\right)
 *  \f].
 *  The function assume that there is no NA or infinite values and
 *  that <em>0<=p<=1, n>=0, 0<=x<=n</em>
 *  @param x,n,p parameters of the binomial density
 */
inline Real binomial_lpdf_raw(Real const& x, Real const& n, Real const& p)
{
  // trivial and limit cases
  if (p == 0) return((x == 0) ? 0 : -Arithmetic<Real>::infinity());
  if (p == 1) return((x == n) ? 0 : -Arithmetic<Real>::infinity());
  if (x == 0) return((n == 0) ? 0 : n*log1p(-p));
  if (x == n) return(n*std::log(p));
  // other cases
  Real y = n-x;
  return ( - Const::_LNSQRT2PI_ - 0.5 *( std::log(x) + log1p(-x/n))
           + lgammaStirlingError(n)-lgammaStirlingError(x)-lgammaStirlingError(y)
           - dev0(x, n*p) - dev0(y, n- n*p)
        );
}

/** @ingroup Analysis
 *  @brief Compute the binomial log-probability mass function.
 *
 *  Compute the function
 *  \f[ B(n,p,x) = x\log(p) (n-x)\log(1-p)
 *      \log\left(\frac{\Gamma(n+1)}{\Gamma(x+1)\Gamma(n-x+1)}\right)
 *  \f].
 *  The function assume that there is no NA or infinite values and
 *  that <em>0<=p<=1, n>=0, 0<=x<=n</em>
 *  @param x,n,p parameters of the binomial density
 */
inline Real binomial_lpdf_raw(int x, int n, Real const& p)
{
  // trivial and limit cases
  if (p == 0) return((x == 0) ? 0 : -Arithmetic<Real>::infinity());
  if (p == 1) return((x == n) ? 0 : -Arithmetic<Real>::infinity());
  if (x == 0) return((n == 0) ? 0 : n*log1p(-p));
  if (x == n) return(n*std::log(p));
  // other cases
  int y = n-x;
  return( - Const::_LNSQRT2PI_ - 0.5 *( std::log(x) + log1p(-x/(Real)n))
          + lgammaStirlingError(n)-lgammaStirlingError(x)-lgammaStirlingError(y)
          - dev0(x, n*p) - dev0(y, n - n*p)
        );
}

/** @ingroup Analysis
 *  @brief Compute the Poisson density.
 *  Compute the function:
 *  \f[
 *    p(x, \lambda) = e^{-\lambda} \frac{\lambda^x}{\Gamma(x+1)}
 *  \f]
 *  with good accuracy using the partial deviance.
 *  This is the version for x Real.
 *
 *  @see http://www.herine.net/stat/software/dbinom.html
 *
 *  @param x value to evaluate the function
 *  @param lambda value of the parameter
 */
inline Real poisson_pdf_raw(Real const& x, Real const& lambda)
{
  // check trivial values
  if (x<0.) return( 0. );
  // if lambda is 0, we have P(X=0) = 1
  if (lambda==0.) return( (x==0.) ? 1. : 0. );
  // special value
  if (x==0.) return( Real(std::exp(-lambda)) );
  // stirling approximation and deviance
  return( std::exp(-lgammaStirlingError(x)-dev0(x, lambda))/(Const::_SQRT2PI_*std::sqrt(x)));
}

/** @ingroup Analysis
 *  @brief Compute the poisson density with integer value.
 *  Compute the function:
 *  \f[
 *    p(x, \lambda) = e^{-\lambda} \frac{\lambda^x}{x!}
 *  \f]
 *  with good accuracy using the partial deviance.
 *  This is the version for x integer.
 *
 *  @see http://www.herine.net/stat/software/dbinom.html
 *
 *  @param x value to evaluate the function
 *  @param lambda value of the parameter
 */
inline Real poisson_pdf_raw(int x, Real const& lambda)
{ return poisson_pdf_raw(Real(x), lambda);}

/** @ingroup Analysis
 *  @brief Compute the log-poisson density.
 *  Compute the function:
 *  \f[
 *    p(x, \lambda) = -\lambda + x \log(\lambda) -\log(\Gamma(x+1))
 *  \f]
 *  with good accuracy using the partial deviance.
 *  This is the version for x Real.
 *
 *  @see http://www.herine.net/stat/software/dbinom.html
 *
 *  @param x value to evaluate the function
 *  @param lambda value of the parameter
 */
inline Real poisson_lpdf_raw(Real const& x, Real const& lambda)
{
  // check trivial values
  if (x<0.) return( -Arithmetic<Real>::infinity() );
  // if lambda is 0, we have P(X=0) = 1
  if (lambda==0.) return( (x==0.) ? 0 : -Arithmetic<Real>::infinity() );
  // special value
  if (x==0.) return( -lambda );
  // stirling approximation and deviance
  return( -lgammaStirlingError(x)-dev0(x, lambda)-Const::_LNSQRT2PI_-std::log(x)/2.);
}

/** @ingroup Analysis
 *  @brief Compute the log-poisson density with integer value.
 *  Compute the function:
 *  \f[
 *    p(x, \lambda) =  -\lambda + x \log(\lambda) -\log(x!)
 *  \f]
 *  with good accuracy using the partial deviance.
 *  This is the version for x integer.
 *
 *  @see http://www.herine.net/stat/software/dbinom.html
 *
 *  @param x value to evaluate the function
 *  @param lambda value of the parameter
 */
inline Real poisson_lpdf_raw(int x, Real const& lambda)
{ return poisson_lpdf_raw(Real(x), lambda);}

/** @ingroup Analysis
 *  @brief Compute the error function erf(a)
 *  Compute the function
 *   \f[
 *   \mathrm{erf}(a) = \frac{2}{\sqrt{\pi}} \int_0^{a} e^{-t^2} dt
 *  \f]
 *  @param[in] a value to evaluate the function
 */
inline Real erf_raw(Real const& a)
{
  Real T[5] =
  { 9.60497373987051638749E0,
    9.00260197203842689217E1,
    2.23200534594684319226E3,
    7.00332514112805075473E3,
    5.55923013010394962768E4
  };
  Real U[5] =
  {
    /* 1.00000000000000000000E0,*/
    3.35617141647503099647E1,
    5.21357949780152679795E2,
    4.59432382970980127987E3,
    2.26290000613890934246E4,
    4.92673942608635921086E4
  };
  if ( std::abs(a)> 1.0) return ( 1.0 - erfc(a) );
  Real z = a * a;
  return a * evalPolynomial<4>( z, T)
           / evalPolynomial1<5>( z, U);
}

/** @ingroup Analysis
 *  @brief Compute the complementary error function erfc(a)
 *  Compute the function
 *   \f[
 *   \mathrm{erfc}(a) = \frac{2}{\sqrt{\pi}} \int_a^{+\infty} e^{-t^2} dt
 *  \f]
 *  @param[in] a value to evaluate the function
 */
inline Real erfc_raw(Real const& a)
{
  Real P[9] =
  { 2.46196981473530512524E-10,
    5.64189564831068821977E-1,
    7.46321056442269912687E0,
    4.86371970985681366614E1,
    1.96520832956077098242E2,
    5.26445194995477358631E2,
    9.34528527171957607540E2,
    1.02755188689515710272E3,
    5.57535335369399327526E2
  };
  Real Q[8] =
  {
    /* 1.00000000000000000000E0,*/
    1.32281951154744992508E1,
    8.67072140885989742329E1,
    3.54937778887819891062E2,
    9.75708501743205489753E2,
    1.82390916687909736289E3,
    2.24633760818710981792E3,
    1.65666309194161350182E3,
    5.57535340817727675546E2
  };
  Real R[6] =
  { 5.64189583547755073984E-1,
    1.27536670759978104416E0,
    5.01905042251180477414E0,
    6.16021097993053585195E0,
    7.40974269950448939160E0,
    2.97886665372100240670E0
  };
  Real S[6] =
  {
    /* 1.00000000000000000000E0,*/
    2.26052863220117276590E0,
    9.39603524938001434673E0,
    1.20489539808096656605E1,
    1.70814450747565897222E1,
    9.60896809063285878198E0,
    3.36907645100081516050E0
  };
  Real x = std::abs(a);
  if ( x < 1.0)  return ( 1.0 - erf_raw(a) );
  Real z = exp(-a * a);
  Real y = ( x < 8.0) ? (z * evalPolynomial<8>(x, P))/evalPolynomial1<8>(x, Q)
                      : (z * evalPolynomial<5>(x, R))/evalPolynomial1<6>(x, S);
  return (a < 0) ? 2.0 - y : y;
}

/** @ingroup Analysis
 *  @brief Compute the cumulative distribution function of the normal density.
 *  Compute the cumulative distribution function of the normal density
 *   \f[
 *   \mathrm{\Phi}(x)
 *         = \frac{1}{\sqrt{2\pi}} \int_{-\infty}^{x} e^{-t^2/2} dt
 *         =  ( 1 + \mathrm{erf}(z) ) / 2
 *         =  \mathrm{erfc}(z) / 2
 *  \f]
 * where \f$ z = x/sqrt(2) \f$. Computation use functions
 * @code erf @endcode and @code erfc @endcode.
 *
 *  @param x value to evaluate the function
 */
inline Real normal_cdf_raw(Real const& x)
{
  Real t = x * Const::_1_SQRT2_, z = std::abs(t);
  if ( z < Const::_1_SQRT2_) return 0.5 + 0.5 * erf(t);
  Real y = 0.5 * erfc(z);
  return ( t > 0) ? 1. - y : y;
}

/** @ingroup Analysis
 *  @brief compute the probability distribution function of the normal density.
 *  compute the probability density function of the normal density
 *   \f[
 *   \mathrm{\phi}(x) = \frac{1}{\sqrt{2\pi}} e^{-x^2/2}
 *  \f]
 *  @param x value to evaluate the function
 */
inline Real normal_pdf_raw(Real const& x)
{ return Const::_1_SQRT2PI_ * exp(-0.5 * x * x);}

/** @ingroup Analysis
 *  @brief Compute the psi function.
 *  Compute the psi function
 * \f[
 *   \psi(x)  =\frac{d}{dx}(-\ln(\Gamma(x)))
 * \f]
 * the logarithmic derivative of the gamma function.
 *
 * For integer x, we use the formula
 * \f[
 * \psi(n) = -\gamma + \sum_{k=1}^{n-1}  \frac{1}{k}.
 * \f]
 * This formula is used for 0 < n <= 20.
 *
 * If x is negative, it is transformed to a positive argument by the reflection
 * formula \f$ \psi(1-x) = \psi(x) + \pi \cot(\pi x) \f$.
 *
 * For general positive x, the argument is made greater than 10
 * using the recurrence  \f$ \psi(x+1) = \psi(x) + 1/x \f$.
 * Then the following asymptotic expansion is applied:
 * \f[
 * \psi(x) = \log(x) - \frac{1}{2} x - \sum_{k=1}^{\infty} \frac{B_{2k}}{2k x^{2k}}
 * \f]
 * where the \f$ B_{2k} \f$ are Bernoulli numbers.
 **/
inline Real psi_raw( Real x)
{
  // special value
  if (x==1.) return -Const::_EULER_;
  Real y, nz = 0.0, p = std::floor(x);
  bool negative = false;
  if( x < 0.0 ) // use reflection formula
  {
    negative = true;
    // Remove the zeros of tan(pi x) by subtracting the nearest integer from x
    if( (nz = x - p) != 0.5 )
    {
      if( nz > 0.5 ) { nz = x - p - 1.;}
      nz = Const::_PI_/std::tan(Const::_PI_*nz);
    }
    x = 1.0 - x;
  }
  // check for positive integer up to 20
  if( (x <= 20.0) && (x == p) )
  {
    y = 0.0;
    for( int i=int(x)-1; i>1; i-- ) { y += 1.0/(Real)i; }
    y += -Const::_EULER_ + 1.;
  }
  else // not integer or greater than 20
  {
    Real w = 0.0;
    while( x < 10.0 ) // shift to
    {
      w += 1.0/x;
      x += 1.0;
    }
    Real z = 1.0/(x * x);
    y = std::log(x) - (0.5/x)
      - z * evalPolynomial<6>( z, Const::bernouilliNumbersArrayDivBy2K) -  w;
  }
  return (negative) ? y - nz : y;
}

} // namespace Funct

} // namespace STK

#endif // STK_FUNCT_RAW_H
