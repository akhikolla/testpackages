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
 * Purpose:  Declaration of functions around the gamma function rqtio
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Funct_gammaRatio.h
 *  @brief In this file we declare functions around the gamma
 *  function ratio.
 **/

#ifndef STK_FUNCT_GAMMARATIO_H
#define STK_FUNCT_GAMMARATIO_H

#include "STKernel/include/STK_Real.h"

namespace STK
{

namespace Funct
{

#ifdef IS_RTKPP_LIBRARY

inline Real gammaRatio(Real const& a, Real const& x, bool lower_tail)
{ return Rf_pgamma(x, a, 1., lower_tail, false);}

#else

Real gammaRatio_raw(Real const& a, Real const& x, bool lower_tail);

inline Real gammaRatioQ_raw(Real const& a, Real const& x)
{ return gammaRatio_raw(a, x, false);}

inline Real gammaRatioP_raw(Real const& a, Real const& x)
{ return gammaRatio_raw(a, x, true);}

#endif

/** @ingroup Analysis
 *  @brief Compute the incomplete gamma functions ratio.
 * Compute the incomplete gamma function ratio P(a,x)
 *  \f[ P(a, x) = \frac{1}{\Gamma(a)}
 *      \int_0^x e^{-t} t^{a-1} dt
 *  \f]
 *  @param a parameter of the gamma ratio function
 *  @param x value to evaluate the gamma ratio function
 *  @param lower_tail @c true if we want the lower tail, @c false otherwise
 **/
Real gammaRatio(Real const& a, Real const& x, bool lower_tail);

/** @ingroup Analysis
 *  @brief Compute the incomplete gamma function ratio Q(a,x).
 *  Compute the incomplete gamma function ratio Q(a,x)
 *  \f[ Q(a, x) = \frac{1}{\Gamma(a)}
 *      \int_x^\infty e^{-t} t^{a-1} dt
 *  \f]
 *  @param a parameter of the gamma ratio function
 *  @param x value to evaluate the gamma ratio function
 **/
inline Real gammaRatioQ(Real const& a, Real const& x)
{ return gammaRatio(a, x, false);}

/** @ingroup Analysis
 *  @brief Compute the incomplete gamma function ratio P(a,x).
 *  Compute the incomplete gamma function ratio P(a,x)
 *  \f[ P(a, x) = \frac{1}{\Gamma(a)}
 *      \int_0^x e^{-t} t^{a-1} dt
 *  \f]
 *  @param a parameter of the gamma ratio function
 *  @param x value to evaluate the gamma ratio function
 **/
inline Real gammaRatioP(Real const& a, Real const& x)
{ return gammaRatio(a, x, true);}

} // namespace Funct

} // namespace STK
	
#endif /*STK_FUNCT_GAMMARATIO_H*/
