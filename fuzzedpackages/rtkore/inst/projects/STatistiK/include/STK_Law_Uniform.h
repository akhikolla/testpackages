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

/** @file STK_Law_Uniform.h
 *  @brief In this file we implement the (continuous) uniform distribution law.
 **/

#ifndef STK_LAW_UNIFORM_H
#define STK_LAW_UNIFORM_H

#include "STK_Law_IUnivLaw.h"
#include "../include/STK_Law_Util.h"
#include <STKernel/include/STK_Real.h>
#include <Sdk/include/STK_Macros.h>

namespace STK
{

namespace Law
{
/** @ingroup Laws
 *  @brief class for the Uniform law distribution.
 *
 *  In probability theory and statistics, the <em>continuous uniform distribution</em>
 *  or rectangular distribution is a family of symmetric probability distributions
 *  such that for each member of the family, all intervals of the same length on
 *  the distribution's support are equally probable. The support is defined by
 *  the two parameters, @e a and @e b, which are its minimum and maximum values.
 *
 *  The probability density function of the continuous uniform distribution is:
 *  \f[
 *    f(x; a, b) = \frac{1}{b-a} 1_{ a \leq x \leq b}.
 *  \f]
**/
class Uniform: public IUnivLaw<Real>
{
  public:
    typedef IUnivLaw<Real> Base;
    /** constructor.
     *  @param a,b the lower and upper bounds
     **/
    inline Uniform( Real const& a =0., Real const& b =1.)
                 : Base(_T("Uniform")), a_(a), b_(b), range_(b_ - a_)
    {
      if (range_ <= 0.)
        STKINVALIDARGUMENT_ERROR_2ARG(Uniform::Uniform, a_, b_,invalid parameters);
    }
    /** copy constructor.
     *  @param law the law to copy
     **/
    inline Uniform( Uniform const& law): Base(law), a_(law.a_), b_(law.b_), range_(law.range_)
    {};
    /** destructor. */
	  inline virtual ~Uniform() {}
    /** @return the lower bound */
    inline Real const& a() const { return a_;}
    /** @return the upper bound */
    inline Real const& b() const { return b_;}
    /** @return the value b-a */
    inline Real const& range() const { return range_;}
    /** @param a set the lower bound */
    inline void setA(Real const& a) { a_ =a;}
    /** @param b set the upper bound */
    inline void setB(Real const& b){ b_ =b;}

    /** Generate a pseudo Uniform random variate. */
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
     *  F(t; a,b)= \frac{t - a}{b-a}
     * \f]
     *  @param t a real value
     **/
    virtual Real cdf( Real const& t) const;
    /** The inverse cumulative distribution function is
     * \f[
     * F^{-1}(p; \lambda) = p (b-a) + a.
     * \f]
     *  @param p a probability
     **/
    virtual Real icdf( Real const& p) const;

    /** Generate a pseudo Uniform random variate.
     *  @param a,b the lower and upper bounds
     **/
    static Real rand( Real const& a, Real const& b);
    /** Give the value of the pdf at x.
     *  @param x a real value
     *  @param a,b the lower and upper bounds
     **/
    static Real pdf( Real const& x, Real const& a, Real const& b);
    /** Give the value of the log-pdf at x.
     *  @param p a probablility
     *  @param a,b the lower and upper bounds
     **/
    static Real lpdf( Real const& p, Real const& a, Real const& b);
    /** Give the value of the cdf at t.
     *  @param t a real value
     *  @param a,b the lower and upper bounds
     **/
    static Real cdf( Real const& t, Real const& a, Real const& b);
    /** Give the value of the quantile at @e p.
     *  @param p a probability
     *  @param a,b the lower and upper bounds
     **/
    static Real icdf( Real const& p, Real const& a, Real const& b);

#ifdef IS_RTKPP_LIB
    /** @return log-cumulative distribution function */
    virtual Real lcdf( Real const& t) const;
    /** @return complement of cumulative distribution function */
    virtual Real cdfc( Real const& t) const;
    /** @return log-complement of cumulative distribution function */
    virtual Real lcdfc( Real const& t) const;
#endif

  protected:
    /** The lower bound. */
    Real a_;
    /** The upper bound. */
    Real b_;

  private:
    Real range_;
};

#ifdef IS_RTKPP_LIB

inline Real Uniform::rand() const
{
#ifdef _OPENMP
//#pragma omp critical
#endif
  GetRNGstate(); Real s = Rf_runif(a_, b_); PutRNGstate(); return s;
}

inline Real Uniform::pdf( Real const& x)  const { return Rf_dunif(x, a_, b_, false);}
inline Real Uniform::lpdf( Real const& x) const { return Rf_dunif(x, a_, b_, true);}
inline Real Uniform::cdf( Real const& t)  const { return Rf_punif(t, a_, b_, true, false);}
inline Real Uniform::lcdf( Real const& t)  const { return Rf_punif(t, a_, b_, true, true);}
inline Real Uniform::cdfc( Real const& t)  const { return Rf_punif(t, a_, b_, false, false);}
inline Real Uniform::lcdfc( Real const& t)  const { return Rf_punif(t, a_, b_, false, true);}
inline Real Uniform::icdf( Real const& p) const { return Rf_qunif(p , a_, b_, true, false);}

inline Real Uniform::rand( Real const& a, Real const& b)
{
#ifdef _OPENMP
//#pragma omp critical
#endif
  GetRNGstate(); Real s = Rf_runif(a, b); PutRNGstate(); return s;
}

inline Real Uniform::pdf(const Real& x, const Real& a, const Real& b)
{ return Rf_dunif(x,a, b, false);}
inline Real Uniform::lpdf(const Real& x, const Real& a, const Real& b)
{ return Rf_dunif(x, a, b, true);}
inline Real Uniform::cdf(const Real& t, const Real& a, const Real& b)
{ return Rf_punif(t, a, b, true, false);}
inline Real Uniform::icdf(const Real& p, const Real& a, const Real& b)
{ return Rf_qunif(p , a, b, true, false);}

#endif

} // namespace Law

} // namespace STK

#endif /*STK_LAW_UNIFORM_H*/
