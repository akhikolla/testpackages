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

/** @file STK_Law_UniformDiscrete.h
 *  @brief In this file we implement the uniform (discrete) law.
 **/

#ifndef STK_LAW_UNIFORMDISCRETE_H
#define STK_LAW_UNIFORMDISCRETE_H

#include "STK_Law_IUnivLaw.h"
#include "../include/STK_Law_Util.h"
#include <STKernel/include/STK_Integer.h>
#include <Sdk/include/STK_Macros.h>

namespace STK
{

namespace Law
{
/** @ingroup Laws
 *  @brief class for the Uniform law distribution.
 *  In probability theory and statistics, the <em> discrete uniform distribution </em>
 *  is a probability distribution whereby a finite number of values are equally
 *  likely to be observed; every one of @e n values has equal probability <em>1/n</em>.
 *  Another way of saying "discrete uniform distribution" would be "a known,
 *  finite number of outcomes equally likely to happen".
 *
 *  The probability density function of the discrete uniform distribution is:
 *  \f[
 *    f(x; a, b) = \frac{1}{b-a+1} 1_{ a \leq x \leq b}, \quad a,b,x\in\mathbb{N}.
 *  \f]
**/
class UniformDiscrete: public IUnivLaw<int>
{
  public:
    typedef IUnivLaw<int> Base;
    /** constructor.
     *  @param a,b the lower and upper bounds
     **/
    inline UniformDiscrete( int a, int b): Base(_T("UniformDiscrete")), a_(a), b_(b), n_(b_ - a_ + 1)
    {
      if (n_ < 0.)
      { STKINVALIDARGUMENT_ERROR_2ARG(UniformDiscrete::UniformDiscrete, a_, b_,invalid parameters);}
    }
    /** copy constructor.
     *  @param law the law to copy
     **/
    inline UniformDiscrete( UniformDiscrete const& law)
                         : Base(law), a_(law.a_), b_(law.b_), n_(law.n_)
    {}
    /** destructor. */
	  inline virtual ~UniformDiscrete() {}
    /** @return the lower bound */
    inline int const& a() const { return a_;}
    /** @return the upper bound */
    inline int const& b() const { return b_;}
    /** @return the value b-a+1 */
    inline Real const& n() const { return n_;}
    /** @param a set the lower bound */
    inline void setA(int a) { a_ =a; n_ = b_-a_+1;}
    /** @param b set the upper bound */
    inline void setB(int b){ b_ =b; n_ = b_ - a_ + 1;}

    /** Generate a pseudo Uniform random variate. */
    virtual int rand() const;
    /** Give the value of the pdf at x.
     *  @param x a real value
     **/
    virtual Real pdf( int const& x) const;
    /** Give the value of the log-pdf at x.
     *  @param x a real value
     **/
    virtual Real lpdf( int const& x) const;
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
    virtual int icdf( Real const& p) const;

    /** Generate a pseudo Uniform random variate.
     *  @param a,b the lower and upper bounds
     **/
    static int rand( int a, int b);
    /** Give the value of the pdf at x.
     *  @param x a real value
     *  @param a,b the lower and upper bounds
     **/
    static Real pdf( Real const& x, int a, int b);
    /** Give the value of the log-pdf at x.
     *  @param p a probablility
     *  @param a,b the lower and upper bounds
     **/
    static Real lpdf( Real const& p, int a, int b);
    /** Give the value of the cdf at t.
     *  @param t a real value
     *  @param a,b the lower and upper bounds
     **/
    static Real cdf( Real const& t, int a, int b);
    /** Give the value of the quantile at @e p.
     *  @param p a probability
     *  @param a,b the lower and upper bounds
     **/
    static int icdf( Real const& p, int a, int b);

  protected:
    /** The lower bound. */
    int a_;
    /** The upper bound. */
    int b_;

  private:
    Real n_;
};

#ifdef IS_RTKPP_LIB

inline int UniformDiscrete::rand() const
{
#ifdef _OPENMP
//#pragma omp critical
#endif
GetRNGstate(); Real s = a_ + Rf_runif(0, double(n_)); PutRNGstate(); return s;
}

inline Real UniformDiscrete::pdf( int const& x) const
{
  if (!Arithmetic<Real>::isFinite(x) ) return x;
  if ((x < a_)||(x > b_)) return 0.;
  return 1./n_;
}
inline Real UniformDiscrete::lpdf( int const& x) const
{
  if (!Arithmetic<Real>::isFinite(x) ) return x;
  if ((x < a_)||(x > b_)) return -Arithmetic<Real>::infinity();
  return -std::log(n_);
}
inline Real UniformDiscrete::cdf( Real const& t) const
{
  if (!Arithmetic<Real>::isFinite(t) ) return t;
  if (t <= a_) return 0.;
  if (t >= b_) return 1.;
  return (b_ - (int)t)/n_;
}
inline int UniformDiscrete::icdf( Real const& p) const
{
  // check parameter
  if ((p > 1.) || (p < 0.))
    STKDOMAIN_ERROR_1ARG(Exponential::icdf,p,invalid argument);

  if (!Arithmetic<Real>::isFinite(p) ) return p;
  if (p == 1.) return b_;
  if (p == 0.) return a_;
  return(int)((1.-p) * a_ + p * b_);
}

inline int UniformDiscrete::rand( int a, int b)
{
#ifdef _OPENMP
//#pragma omp critical
#endif
  GetRNGstate(); Real s = a + Rf_runif(0, double(b - a + 1)); PutRNGstate(); return s;
}

inline Real UniformDiscrete::pdf( Real const& x, int a, int b)
{
  if (!Arithmetic<Real>::isFinite(x) ) return x;
  if ((x < a)||(x > b)) return 0.;
  return 1./Real(b-a+1);
}
inline Real UniformDiscrete::lpdf( Real const& x, int a, int b)
{
  if (!Arithmetic<Real>::isFinite(x) ) return x;
  if ((x < a)||(x > b)) return -Arithmetic<Real>::infinity();
  return -std::log(b-a+1);
}
inline Real UniformDiscrete::cdf(const Real& t, int a, int b)
{ return (b - t)/(b-a+1);}
inline int UniformDiscrete::icdf(const Real& p, int a, int b)
{ return (int)((1.-p) * a + p * b);}

#endif

} // namespace Law

} // namespace STK

#endif /*STK_LAW_UNIFORMDISCRETE_H*/
