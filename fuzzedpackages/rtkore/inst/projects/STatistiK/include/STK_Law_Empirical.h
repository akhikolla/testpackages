/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2016  Serge Iovleff, Université Lille 1, Inria

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
 * Project:  stkpp::STatistiK::Law
 * Purpose:  Empirical probability distribution.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Law_Cauchy.h
 *  @brief In this file we define the Empirical probability distribution.
 **/

#ifndef STK_LAW_CAUCHY_H
#define STK_LAW_CAUCHY_H

#include "STK_Law_IUnivLaw.h"
#include <Sdk/include/STK_Macros.h>
#include <STKernel/include/STK_Real.h>

namespace STK
{

namespace Law
{
/** @ingroup Laws
 *  @brief Empirical distribution law.
 *
 *  In statistics, an empirical distribution function is the distribution
 *  function associated with the empirical measure of a sample. This cumulative
 *  distribution function is a step function that jumps up by \f$ 1/n \f$ at
 *  each of the @e n data points. Its value at any specified value of the
 *  measured variable is the fraction of observations of the measured variable
 *  that are less than or equal to the specified value.
 *
 *  The empirical distribution function is an estimate of the cumulative
 *  distribution function that generated the points in the sample. It converges
 *  with probability 1 to that underlying distribution, according to the
 *  Glivenko–Cantelli theorem.
 **/
class Empirical: public IUnivLaw<Real>
{
  public:
    typedef IUnivLaw<Real> Base;

    /** Destructor. **/
    inline virtual ~Empirical() {}

  protected:
};

#ifdef IS_RTKPP_LIB

/
#endif

} // namespace Law

} // namespace STK

#endif /*STK_LAWCAUCHY_H*/
