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
 * Purpose:  Usual mathematical functions
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Funct_Util.h
 *  @brief In this file we declare usual Real functions.
 **/

#ifndef STK_FUNCT_UTIL_H
#define STK_FUNCT_UTIL_H

#include "STK_ISerie.h"

namespace STK
{

namespace Funct
{
/** @ingroup Analysis
 *  This Series computes
 *  \f[ \frac{x^n}{n+2} \f]
 *  for @e n >= 0.
 */
class Serielog1p: public ISerie<Serielog1p>
{
  public:
    inline Serielog1p( Real const& x): x_(x), n_(2), xpown_(1.) {}
    inline Real firstImpl() const { return 1./n_;}
    inline Real nextImpl() const { return (xpown_ *= x_)/(++n_);}

  private:
    const   Real& x_;
    mutable Real n_;
    mutable Real xpown_;
};

/** @ingroup Analysis
 *  @brief Polynomial evaluator.
 *  Evaluate the quantity
 *  \f[
 *      P(x) = P[0] x^n  +  P[1] x^(n-1)  + \ldots  +  P[n]
 *  \f]
 *  @note Coefficients are stored in reverse order.
 * There is no checks for out of bounds.
 */
template<int N>
inline Real evalPolynomial( Real x, const Real* P)
{ return P[N] +  x * evalPolynomial<N-1>(x,P);}
template<>
inline Real evalPolynomial<0>( Real x, const Real* P)
{ return P[0];}

/** @ingroup Analysis
 *  @brief Polynomial evaluator.
 *  Evaluate the quantity
 *  \f[
 *      P(x) = 1. * x^n  +  P[0] x^(n-1)  + \ldots  +  P[n-1]
 *  \f]
 *  @note Coefficients are stored in reverse order.
 * There is no checks for out of bounds.
 * @tparam N degree of the polynomial
 */
template<int N>
inline Real evalPolynomial1( Real x, const Real* P)
{ return P[N-1] +  x * evalPolynomial1<N-1>(x,P);}
template<>
inline Real evalPolynomial1<0>( Real x, const Real* P)
{ return 1.;}

} // namespace Funct

} // namespace STK

#endif // STK_FUNCT_UTIL_H
