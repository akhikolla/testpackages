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
 * created on: 23 janv. 2013
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Law_Exponential.h
 *  @brief In this file we implement the exponential law.
 **/

#ifndef STK_LAW_EXPONENTIAL_H
#define STK_LAW_EXPONENTIAL_H

#include "STK_Law_IUnivLaw.h"
#include <Sdk/include/STK_Macros.h>
#include <STKernel/include/STK_Real.h>

namespace STK
{

namespace Law
{
/** @ingroup Laws
 *  @brief Exponential distribution law.
 *
 *  In probability theory and statistics, the <em>exponential distribution</em>
 *  (a.k.a. negative exponential distribution) is the probability distribution
 *  that describes the time between events in a Poisson process, i.e. a process
 *  in which events occur continuously and independently at a constant average
 *  rate. It is a particular case of gamma distribution. It is the continuous
 *  analogue of the geometric distribution, and it has the key property of being
 *  memoryless. In addition to being used for the analysis of Poisson processes,
 *  it is found in various other contexts.
 *
 *  The probability density function (pdf) of an exponential distribution is
 *  \f[
 *    f(x; \lambda) = \frac{1}{\lambda} e^{- x/\lambda} 1_{x\geq 0}
 *  \f]
 *  where \f$\lambda>0\f$ is the scale (inverse rate) parameter.
**/
class Exponential: public IUnivLaw<Real>
{
  public:
    typedef IUnivLaw<Real> Base;
    /** constructor. */
    inline Exponential( Real const& scale=1)
                     : Base(String(_T("Exponential"))), scale_(scale)
    {
      // check parameters
      if ( !Arithmetic<Real>::isFinite(scale) || scale <= 0 )
        STKDOMAIN_ERROR_1ARG(Exponential::Exponential,scale,invalid argument);
    }
    /** destructor. */
	  inline virtual ~Exponential() {}

    /** Generate a pseudo Exponential random variate. */
    virtual Real rand() const;
    /** Give the value of the pdf at x.
     *  @param x a real value
     **/
    virtual Real pdf( Real const& x) const;
    /** Give the value of the log-pdf at x.
     *  @param x a real value
     **/
    virtual Real lpdf( Real const& x) const;
    /** The cumulative distribution function is
     * \f[
     *  F(t; \lambda)= 1- e^{- t/\lambda}
     * \f]
     *  @param t a real value
     **/
    virtual Real cdf( Real const& t) const;
    /** The inverse cumulative distribution function is
     * \f[
     * F^{-1}(p; \lambda) = - \lambda\log(1-p).
     * \f]
     *  @param p a probability
     **/
    virtual Real icdf( Real const& p) const;

    /** Generate a pseudo Exponential random variate with the specified
     *  parameter.
     *  @param scale the scale of the distribution
     **/
    static Real rand( Real const& scale);
    /** Give the value of the pdf at x.
     *  @param x a real value
     *  @param scale the scale of the distribution
     **/
    static Real pdf( Real const& x, Real const& scale);
    /** Give the value of the log-pdf at x.
     *  @param x a real value
     *  @param scale the scale of the distribution
     **/
    static Real lpdf( Real const& x, Real const& scale);
    /** Compute he cumulative distribution function
     *  @param t a real value
     *  @param scale the scale of the distribution
     **/
    static Real cdf( Real const& t, Real const& scale);
    /** Compute rhe inverse cumulative distribution function
     *  @param p a probability
     *  @param scale the scale of the distribution
     **/
    static Real icdf( Real const& p, Real const& scale);

#ifdef IS_RTKPP_LIB
    /** @return log-cumulative distribution function */
    virtual Real lcdf( Real const& t) const;
    /** @return complement of cumulative distribution function */
    virtual Real cdfc( Real const& t) const;
    /** @return log-complement of cumulative distribution function */
    virtual Real lcdfc( Real const& t) const;
#endif

  protected:
    /** The scale parameter. */
    Real scale_;
};

#ifdef IS_RTKPP_LIB

inline Real Exponential::rand() const
{
#ifdef _OPENMP
//#pragma omp critical
#endif
GetRNGstate(); Real s = Rf_rexp(scale_); PutRNGstate(); return s;
}

inline Real Exponential::pdf( Real const& x) const {   return Rf_dexp(x, scale_, false);}
inline Real Exponential::lpdf( Real const& x) const {   return Rf_dexp(x, scale_, true);}
inline Real Exponential::cdf( Real const& t) const { return Rf_pexp(t, scale_, true, false);}
inline Real Exponential::lcdf( Real const& t) const { return Rf_pexp(t, scale_, true, true);}
inline Real Exponential::cdfc( Real const& t) const { return Rf_pexp(t, scale_, false, false);}
inline Real Exponential::lcdfc( Real const& t) const { return Rf_pexp(t, scale_, false, true);}
inline Real Exponential::icdf( Real const& p) const { return Rf_qexp(p , scale_, true, false);}

inline Real Exponential::rand( Real const& scale)
{
#ifdef _OPENMP
//#pragma omp critical
#endif
GetRNGstate(); Real s = Rf_rexp(scale); PutRNGstate(); return s;
}

inline Real Exponential::pdf(Real const& x, Real const& scale)
{ return Rf_dexp(x, scale, false);}
inline Real Exponential::lpdf(Real const& x, Real const& scale)
{ return Rf_dexp(x, scale, true);}
inline Real Exponential::cdf(Real const& t, Real const& scale)
{ return Rf_pexp(t, scale, true, false);}
inline Real Exponential::icdf(Real const& p, Real const& scale)
{ return Rf_qexp(p, scale, true, false);}

#endif
} // namespace Law

} // namespace STK

#endif /*STK_LAW_EXPONENTIAL_H*/
