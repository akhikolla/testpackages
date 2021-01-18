/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2016 Serge Iovleff

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
 * Project: stkpp::STatistiK::Law
 * Purpose: Normal probability distribution.
 * Author:  Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Law_Normal.h
 *  @brief In this file we define the Normal probability law class.
 **/

#ifndef STK_LAW_NORMAL_H
#define STK_LAW_NORMAL_H

#include "STK_Law_IUnivLaw.h"
#include <Sdk/include/STK_Macros.h>
#include <STKernel/include/STK_Real.h>

namespace STK
{

namespace Law
{
class Normal;
/** A synonymous for the normal law */
typedef Normal Gaussian;

/**
 *  @ingroup Laws
 *  @brief Normal distribution law.
 *
 *  In probability theory, the <em>normal (or Gaussian) distribution</em> is a
 *  very commonly occurring continuous probability distribution. Normal
 *  distributions are extremely important in statistics and are often used in
 *  the natural and social sciences for real-valued random variables whose
 *  distributions are not known.
 *
 *  The normal distribution is immensely useful because of the central limit
 *  theorem, which states that, under mild conditions, the mean of many random
 *  variables independently drawn from the same distribution is distributed
 *  approximately normally, irrespective of the form of the original distribution.
 *
 *  The probability density of normal distribution is:
 *  \f[
 *  f(x;\mu,\sigma) = \frac{1}{\sigma\sqrt{2\pi}}
 *            \exp\left(-\frac{\left(x-\mu\right)^2}{2\sigma^2} \right)
 *  \f]
 * where \f$ \mu \mbox{ and } \sigma\f$ are the mean and the standard deviation.
**/
class Normal: public IUnivLaw<Real>
{
  public:
    typedef IUnivLaw<Real> Base;
    /** Constructor.
     *  @param mu mean of the Normal distribution
     *  @param sigma standard deviation of the Normal distribution
     **/
    Normal( Real const& mu=0., Real const& sigma=1.)
         : Base(String(_T("Normal")))
          , mu_(mu), sigma_(sigma)
    {
      if (!Arithmetic<Real>::isFinite(mu) || !Arithmetic<Real>::isFinite(sigma) || sigma < 0)
      { STKDOMAIN_ERROR_2ARG(Normal::Normal,mu,sigma,invalid argument);}
    }
    /** Destructor. **/
    inline virtual ~Normal() {}
    /** @return the mean */
    inline Real const& mu() const { return mu_;}
    /** @return the standard deviation */
    inline Real const& sigma() const { return sigma_;}
    /** @param mu the mean to set */
    inline void setMu( Real const& mu) { mu_ = mu;}
    /** @param sigma the standard deviation to set */
    inline void setSigma( Real const& sigma)
    {
      if (sigma<0) STKDOMAIN_ERROR_1ARG(Normal::setSigma,sigma,sigma must be >= 0);
       sigma_ = sigma;
    }
    /** @brief Generate a pseudo Normal random variate.
     *  Generate a pseudo Normal random variate
     *  with location parameter @c mu_ and standard deviation @c sigma_.
     *  @return a pseudo normal random variate
     **/
    Real rand() const;
    /** @param x a real value
     *  @return the value of the normal pdf at @c x
     **/
    virtual Real pdf( Real const& x) const;
    /** @return Give the value of the log-pdf at x.
     *  @param x a real value
     **/
    virtual Real lpdf( Real const& x) const;
    /** @brief Compute the cumulative distribution function at t of
     *  the standard normal distribution.
     *
     *  @author  W. J. Cody
     *  @see     http://www.netlib.org/specfun/erf
     *
     *  This is the erfc() routine only, adapted by the
     *  transform cdf(u)=erfc(-u/sqrt(2))/2
     *
     *  @param t a real value
     *  @return the cumulative distribution function value at t
     **/
    virtual Real cdf( Real const& t) const;
    /** @brief Compute the inverse cumulative distribution function at p
     *  of the standard normal distribution.
     *
     *  @author Peter J. Acklam <pjacklam@online.no>
     *  @see    http://home.online.no/~pjacklam/notes/invnorm/index.html
     *
     *  This function is based on the MATLAB code from the address above.
     *
     *  @param p a probability number.
     *  @return the inverse cumulative distribution function value at p.
     **/
    virtual Real icdf( Real const& p) const;

    /** @brief Generate a pseudo Normal random variate.
     *
     *  Generate a pseudo Normal random variate with location @c mu and
     *  standard deviation @c sigma parameters.
     *  @param mu, sigma mean and standard deviation of the Normal distribution
     *  @return a pseudo normal random variate, centered in @c mu and with
     *  standard deviation @c sigma
     **/
    static Real rand( Real const& mu, Real const& sigma);
    /** @param x a real value
     *  @param mu, sigma mean and standard deviation of the Normal distribution
     *  @return the value of the normal pdf at @c x
     **/
    static Real pdf( Real const& x, Real const& mu, Real const& sigma);
    /** @return Give the value of the log-pdf at x.
     *  @param x a real value
     *  @param mu, sigma mean and standard deviation of the Normal distribution
     **/
    static Real lpdf( Real const& x, Real const& mu, Real const& sigma);
    /** @brief Compute the cumulative distribution function at t of
     *  the standard normal distribution.
     *  @param t a real value
     *  @param mu, sigma mean and standard deviation of the Normal distribution
     *  @return the cumulative distribution function value at t
     **/
    static Real cdf( Real const& t, Real const& mu, Real const& sigma);
    /** @brief Compute the inverse cumulative distribution function at p
     *  of the standard normal distribution.
     *  @param p a probability number.
     *  @param mu, sigma mean and standard deviation of the Normal distribution
     *  @return the inverse cumulative distribution function value at p.
     **/
    static Real icdf( Real const& p, Real const& mu, Real const& sigma);

#ifdef IS_RTKPP_LIB
    /** @return log-cumulative distribution function */
    virtual Real lcdf( Real const& t) const;
    /** @return complement of cumulative distribution function */
    virtual Real cdfc( Real const& t) const;
    /** @return log-complement of cumulative distribution function */
    virtual Real lcdfc( Real const& t) const;
#endif

  protected:
    /** The mu parameter. **/
    Real mu_;
    /** The sigma parameter. **/
    Real sigma_;
};

#ifdef IS_RTKPP_LIB

/*  Generate a pseudo Normal random variate. */
inline Real Normal::rand() const
{
#ifdef _OPENMP
//#pragma omp critical
#endif
GetRNGstate(); Real s = Rf_rnorm(mu_, sigma_); PutRNGstate(); return s;
}

inline Real Normal::pdf( Real const& x) const {   return Rf_dnorm4(x, mu_, sigma_, false);}
inline Real Normal::lpdf( Real const& x) const {   return Rf_dnorm4(x, mu_, sigma_, true);}
inline Real Normal::cdf( Real const& t) const { return Rf_pnorm5(t, mu_, sigma_, true, false);}
inline Real Normal::lcdf( Real const& t) const { return Rf_pnorm5(t, mu_, sigma_, true, true);}
inline Real Normal::cdfc( Real const& t) const { return Rf_pnorm5(t, mu_, sigma_, false, false);}
inline Real Normal::lcdfc( Real const& t) const { return Rf_pnorm5(t, mu_, sigma_, false, true);}
inline Real Normal::icdf( Real const& p) const { return Rf_qnorm5(p , mu_, sigma_, true, false);}


inline Real Normal::rand( Real const& mu, Real const& scale)
{
#ifdef _OPENMP
//#pragma omp critical
#endif
GetRNGstate(); Real s = Rf_rnorm(mu, scale); PutRNGstate(); return s;
}

inline Real Normal::pdf(Real const& x, Real const& mu, Real const& scale)
{ return Rf_dnorm4(x,mu, scale, false);}
inline Real Normal::lpdf(Real const& x, Real const& mu, Real const& scale)
{ return Rf_dnorm4(x,mu, scale, true);}
inline Real Normal::cdf(Real const& t, Real const& mu, Real const& scale)
{ return Rf_pnorm5(t, mu, scale, true, false);}
inline Real Normal::icdf(Real const& p, Real const& mu, Real const& scale)
{ return Rf_qnorm5(p , mu, scale, true, false);}

#endif

} // namespace Law

} // namespace STK

#endif /*STK_LAW_NORMAL_H*/
