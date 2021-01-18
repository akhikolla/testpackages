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
 * Purpose:  Declaration of functions around the gamma function
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Funct_gamma.h
 *  @brief In this file we declare functions around the gamma function.
 **/

#ifndef STK_FUNCT_GAMMA_H
#define STK_FUNCT_GAMMA_H

#include <vector>

#include <STKernel/include/STK_Misc.h>
#include <Sdk/include/STK_Macros.h>

#include "STK_Const_Math.h"
#include "STK_Const_Sequences.h"

namespace STK
{

namespace Funct
{
// forward declaration
Real factorial( int);
Real lfactorial( int);

Real factorial( Real const&);
Real lfactorial( Real const&);

Real gamma( Real const&);
Real lgamma( Real const&);

// raw versions
Real factorial_raw( int);
Real lfactorial_raw( int);

Real factorial_raw( Real const&);
Real lfactorial_raw( Real const&);

Real gamma_raw( Real const&);
Real lgamma_raw( Real const&);

/** @ingroup Analysis
 *  @brief Compute the Lanzcos correction series for the gamma function
 *  with n = 21 terms.
 *  @param z given value for the lanzcos Series
 **/
inline Real lanczosSerie(Real const& z)
{
  Real sum = 0.0;
  for (int k=20; k>=0; k--)
    sum += Const::lanczosCoefArray[k]/(z+k);
  return 2.0240434640140357514731512432760e-10 + sum;
}

/** @ingroup Analysis
 *  @brief Compute the gamma function using the Lanzcos expansion
 *  using n = 21 terms and r= 22.618910.
 *  @param z given value for the gamma function
 **/
inline Real gammaLanczos(Real const& z)
{
  // 2 * sqrt(e/pi) = 1.86038273...
  return 1.8603827342052657173362492472666631120594218414085774528900013
       * exp((z-0.5)*(log(z+22.118910)-1.))
       * lanczosSerie(z);
}

/** @ingroup Analysis
 *  @brief Compute the Stirling's series for the lgamma function.
 *  @param z given value for the stirling Series
 **/
inline double stirlingSerie(Real const& z)
{
  const Real z2 = z * z;
  return (z <= 50) ? ( Const::stirlingCoefArray[0]
                       + ( Const::stirlingCoefArray[1]
                         + ( Const::stirlingCoefArray[2]
                           + Const::stirlingCoefArray[3]/z2
                           )/z2
                         )/z2
                     )/z
                   : ( Const::stirlingCoefArray[0]
                       + ( Const::stirlingCoefArray[1]
                         + Const::stirlingCoefArray[2]/z2
                         )/z2
                     )/z
                   ;
}

/** @ingroup Analysis
 *  @brief This function computes the gamma function using the
 *  Stirling approximation.
 *
 *  This approximation is valid for large values of z.
 *  @param z given value for the gamma function
 */
inline Real gammaStirling( Real const& z)
{ return Const::_SQRT2PI_ * exp((z-0.5)*(log(z)-1.)+stirlingSerie(z)-0.5);}

/** @ingroup Analysis
 *  @brief This function computes the log gamma function using the
 *  Stirling approximation.
 *
 *  This approximation is valid for large values of z.
 *  @param z given value for log gamma function
 */
inline Real lgammaStirling( Real const& z)
{ return ( Const::_LNSQRT2PI_ + (z-0.5)*log(z) - z + stirlingSerie(z) );}

/** @ingroup Analysis
 *  @brief Compute the error when we compute  \f$ \ln(\Gamma(z)) \f$
 *  using the Stirling's formula.
 *  Computes the ln of the error term in Stirling's formula.
 *  For z <100, integers or half-integers, use stored values.
 *  For z >= 100, uses the stirling serie
 *  \f$ 1/12n - 1/360n^3 + ... \f$
 *  For other z < 100, uses lgamma directly (don't use this to write lgamma!)
 * \f[
 *  \ln(SE(z)) = \ln(\Gamma(z)) - (z - 1/2) \ln(z) + z - \ln(2\pi)/2
 * \f]
 *  @param z given value for the gamma function
 **/
inline Real lgammaStirlingError(Real const& z)
{
  int n = (int)std::floor(z);
  // z is a discrete value
  if (z == n)
  { return (n<100) ? Const::lgammaStirlingErrorArray[n] : stirlingSerie(z);}
  // z is a discrete value halves
  if (z== n + 0.5)
  { return (n<100) ? Const::lgammaStirlingErrorHalvesArray[n] : stirlingSerie(z);}
  // other cases
  return (z > 16) ? stirlingSerie(z) : lgamma(z) - (z-0.5)*log(z) + z - Const::_LNSQRT2PI_;
}

/** @ingroup Analysis
 *  @brief Compute the error when we compute  \f$ \ln(\Gamma(z)) \f$
 *  using the Stirling's formula and z is an integer.
 *  Computes the ln of the error term in Stirling's formula.
 *  For z <100, integers or half-integers, use stored values.
 *  For z >= 100, uses the stirling series
 *  \f$ 1/12n - 1/360n^3 + ... \f$
 *  For other z < 100, uses lgamma directly (don't use this to write lgamma!)
 * \f[
 *  \ln(SE(z)) = \ln(\Gamma(z)) - (z - 1/2) \ln(z) + z - \ln(2\pi)/2
 * \f]
 *  @param n given value for the gamma function
 **/
inline Real lgammaStirlingError(int n)
{ return (n<100.0) ? Const::lgammaStirlingErrorArray[n] : stirlingSerie(n);}


/** @ingroup Analysis
 *  @brief This function computes \f$ n! \f$ for integer argument.
 *  Compute \f$ n!=1\times 2\times \ldots \times n \f$
 *  using the values stored in @c factorial³Array for n<51
 *  and using the @c gamma function for n>50.
 *  @param n given value for the factorial function
 **/
inline Real factorial( int n)
{
 // Check if z is Available and finite
 if (!Arithmetic<int>::isFinite(n)) return Arithmetic<Real>::NA();
 if (n < 0) { STKDOMAIN_ERROR_1ARG(factorial,n,"Negative argument");}
 return factorial_raw(n);
}
inline Real factorial_raw( int n)
{ return (n < 51) ? Const::factorialArray[n] : Funct::gamma_raw(n + 1.);}

/** @ingroup Analysis
 *  @brief This function computes \f$ z! \f$ when z is an integer in
 *  a Real format.
 *  Compute \f$ z!=1\times 2\times \ldots \times z \f$
 *  using the values stored in @c factorialArray for n<51
 *  and using the @c gamma function for n>50.
 *  @param z given value for the factorial function
**/
inline Real factorial( Real const& z)
{
  // Check if z is Available and finite
  if (!Arithmetic<Real>::isFinite(z)) return z;
  // Negative integers or reals arguments not allowed
  if (z < 0 ||(z != std::floor(z))) { STKDOMAIN_ERROR_1ARG(factorial,z,"Negative or not discrete argument");}
  return factorial_raw(z);
}
inline Real factorial_raw( Real const& z)
{
  const int n = std::floor(z);
  return (n < 51) ? Const::factorialArray[n] : Funct::gamma_raw(z + 1.);
}

/** @ingroup Analysis
 *  @brief This function computes \f$ \ln(n!) \f$ for integer argument.
 *  Compute \f$ \ln(n!)=\ln(2)+ \ldots \ln(n) \f$
 *  using the values stored in @c factorialLnArray for n<51
 *  and using the @c gamma function for n>50.
 *  @param n given value for the factorial function
**/
inline Real lfactorial(int n)
{
  // Check if z is Available and finite
  if (!Arithmetic<int>::isFinite(n)) return n;
  // Negative integers or reals arguments not allowed
  if (n < 0) { STKDOMAIN_ERROR_1ARG(lfactorial,n,"Negative integer argument");}
  // if n is a small number return a constant
  return lfactorial_raw(n);
}
inline Real lfactorial_raw(int n)
{ return (n < 51) ? Const::factorialLnArray[n] : lgamma( n + 1.);}

/** @ingroup Analysis
 *  @brief This function computes \f$ \ln(z!) \f$ when z is an integer
 *  in Real format.
 *  Compute \f$ \ln(z!)=\ln(2)+ \ldots \ln(z) \f$
 *  using the values stored in factorialLnArray for n<51
 *  and using the lgamma function for n>50.
 *  @param z given value for the factorial function
**/
inline Real lfactorial(Real const& z)
{
  // Check if z is Available and finite
  if (!Arithmetic<Real>::isFinite(z)) return z;
  // Negative integers or real arguments not allowed
  if ((z < 0)||(z != std::floor(z)))
  { STKDOMAIN_ERROR_1ARG(lfactorial,z,"Negative integer or decimal argument");}
  // if n is a small number return a constant
  return lfactorial_raw(z);
}
inline Real lfactorial_raw(Real const& z)
{
  const int n = (int)std::floor(z);
  return (n < 51) ? Const::factorialLnArray[n] : lgamma(z + 1.);
}

/** @ingroup Analysis
 *  @brief This function computes the function \f$ \Gamma(z) \f$.
 *  The gamma function is valid when z is non zero nor a negative
 *  integer. The negative part is computed using the reflection formula
 *  \f[
 *   \Gamma(z) \Gamma(1-z) = \frac{\pi}{\sin(\pi z)}.
 *  \f]
 * if |z| <8 we use the gamma Lanczos method, else we use the Stirling
 * approximation method.
 * @param z given value for the gamma function
**/
inline Real gamma(Real const& z)
{
  // Check if z is Available and finite
  if (!Arithmetic<Real>::isFinite(z)) return z;
  // Negative integer argument not allowed
  if ( z<=0 && z == std::floor(z))
  { STKDOMAIN_ERROR_1ARG(Funct::gamma,z,"Negative integer or zero argument");}
  return gamma_raw(z);
}

inline Real gamma_raw(Real const& z)
{
  int n = (int)std::floor(z);
  // If z is an integer (z integer and z<0 is not possible)
  if (z == n)
  { return (z < 51) ? Const::factorialArray[(int)z-1] : gammaStirling(z);}
  // compute the absolute value of x -> y and compute the sign of the gamma function
  Real y = std::abs(z);
  int ny = std::floor(y);
  Real signgam = (z<0 && isEven(ny)) ? -1 : 1, value;
  // if y is an integer and a half, use reflection formula
  if (y == ny+0.5)
  { value = (ny<50) ?  Const::factorialHalvesArray[ny] : gammaStirling(y);}
  else // normal case
  {
    if (y <= 8)
    {
      Real r = y-ny;            // r in [0,1]
      value = gammaLanczos(r);  // use Lanzcos approximation
      // scale result
      for (int i=0; i<ny; i++) value *= (r+i);
    }
    else // shift value
    {
      if (y <=16)
      {
        value = gammaStirling(y+6.);
        for (int i=5; i>=0; --i) value /= (y+i);
      }
      else
        value = gammaStirling(y);
    }
  }
  // z >0 terminate
  if (z>0) return value;
  // z <0  -> use reflection formula and check overflow
  Real den = y*sin(Const::_PI_ *y)*value;
  return (den == 0.0) ?  signgam * Arithmetic<Real>::infinity() : -Const::_PI_/den;
}

/** @ingroup Analysis
 *  @brief This function computes the function \f$ \ln(\Gamma(z)) \f$.
 *  The log gamma function is valid when z is non zero nor a negative
 *  integer.
 *  if |z| <16 we use the gamma Lanczos method, else we use the Stirling
 *  approxiamtion method.
 *  @param z given value for the gamma function
**/
inline Real lgamma(Real const& z)
{
  // Check if z is Available and finite
  if (!Arithmetic<Real>::isFinite(z)) return z;
  // Negative integer argument not allowed
  if ( z<=0 && z == std::floor(z))
  { STKDOMAIN_ERROR_1ARG(Funct::lgamma,z,"Negative integer or zero argument");}
  return lgamma_raw(z);
}
inline Real lgamma_raw(Real const& z)
{
  // compute the absolute value of x -> y
  Real y = std::abs(z), value;
  int ny = std::floor(y);
  // If x is an integer
  if (y == ny)
  { value = (ny<=50) ? Const::factorialLnArray[ny-1] : lgammaStirling(y);}
  else
  {
    // If x is an integer and a half
    if (y == ny+0.5)
    { value = (ny<50) ? Const::factorialLnHalvesArray[ny] : lgammaStirling(y);}
    else
    {
      // small values -> use gamma function
      if (y <= 16)  return log(std::abs(Funct::gamma_raw(z)));
      else          value = lgammaStirling(y);
    }
  }
  // z >0 terminate
  if (z>0) return value;
  // z <0  -> use reflectiong formula
  Real sinpiy = std::abs(sin(Const::_PI_ *y));
  // overflow
  return (sinpiy == 0.0) ? -Arithmetic<Real>::infinity()
                         : Const::_LNSQRTPI_2_+(z-0.5)*log(y)-z-log(sinpiy)+stirlingSerie(z);
}

} // namespace Funct

} // namespace STK

#endif // STK_FUNCT_GAMMA_H
