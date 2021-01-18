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
 * Purpose: LogNormal probability distribution.
 * Author:  Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Law_LogNormal.h
 *  @brief In this file we define the LogNormal probability law class.
 **/

#ifndef STK_LAW_LOGNORMAL_H
#define STK_LAW_LOGNORMAL_H

#include "STK_Law_IUnivLaw.h"
#include <Sdk/include/STK_Macros.h>

namespace STK
{
namespace Law
{
/**
 *  @ingroup Laws
 *  @brief LogNormal distribution law.
 *
 * In probability theory, a <em>log-normal (or loglog-normal) distribution</em> is a
 * continuous probability distribution of a random variable whose logarithm is
 * log-normally distributed. Thus, if the random variable @e X is log-normally
 * distributed, then <em>Y = log(X)</em> has a log-normal distribution. Likewise,
 * if @e Y has a log-normal distribution, then <em>X = exp(Y)</em> has a log-normal
 * distribution. A random variable which is log-normally distributed takes only
 * positive real values.
 *
 *  A variable might be modeled as log-normal if it can be thought of as the
 *  multiplicative product of many independent random variables, each of which
 *  is positive.
 *
 * The probability density function of a log-normal distribution is:
 * \f[
 *   f_X(x;\mu,\sigma) = \frac{1}{ x\sigma \sqrt{2 \pi}}\,
 *   e^{-\frac{(\ln x - \mu)^2}{2\sigma^2}},\ \ x>0
 * \f]
 * where \f$ \mu  \mbox{ and } \sigma\f$ are the location and the
 * scale parameters.
**/
class LogNormal: public IUnivLaw<Real>
{
  public:
    typedef IUnivLaw<Real> Base;
    /** Constructor.
     *  @param mu location of the LogNormal distribution
     *  @param sigma scale of the LogNormal distribution
     **/
    inline LogNormal( Real const& mu=0., Real const& sigma=1.)
                   : Base(_T("LogNormal")), mu_(mu), sigma_(sigma) {}

    /** Destructor. **/
    inline virtual ~LogNormal() {}
    /** @return mu */
    inline Real const& mu() const { return mu_;}
    /** @return sigma */
    inline Real const& sigma() const { return sigma_;}
    /** @param mu the value to set to mu */
    inline void setMu( Real const& mu) { mu_ = mu;}
    /** @param sigma the value to set to sigma */
    inline void setSigma( Real const& sigma)
    {
      if (sigma<0) STKDOMAIN_ERROR_1ARG(LogNormal::setSigma,sigma,sigma must be >= 0);
       sigma_ = sigma;
    }
    /** @brief Generate a pseudo log-normalized LogNormal random variate.
     *
     *  Generate a pseudo log-normalized LogNormal random variate
     *  with location parameter @c mu_ and scale @c sigma_.
     *  @return a pseudo log-normal random variate
     **/
    Real rand() const;
    /** @param x a real value
     *  @return the value of the log-normal pdf at @c x
     **/
    virtual Real pdf( Real const& x) const;
    /** @return Give the value of the log-pdf at x.
     *  @param x a real value
     **/
    virtual Real lpdf( Real const& x) const;
    /** @brief Compute the cumulative distribution function at t of
     *  the standard log-normal distribution.
     *  @param t a real value
     *  @return the cumulative distribution function value at t
     **/
    virtual Real cdf( Real const& t) const;
    /** @brief Compute the inverse cumulative distribution function at p
     *  of the standard log-normal distribution.
     *
     *  @param p a probability number.
     *  @return the inverse cumulative distribution function value at p.
     **/
    virtual Real icdf( Real const& p) const;

    /** @brief Generate a pseudo LogNormal random variate.
     *
     *  Generate a pseudo LogNormal random variate with location @c mu and
     *  scale @c sigma parameters.
     *  @param mu, sigma location and scale of the log-normal law
     *  @return a pseudo log-normal random variate, centered in @c mu and with
     *  scale @c sigma
     **/
    static Real rand( Real const& mu, Real const& sigma);
    /** @param x a real value
     *  @param mu, sigma location and scale of the log-normal law
     *  @return the value of the log-normal pdf at @c x
     **/
    static Real pdf( Real const& x, Real const& mu, Real const& sigma);
    /** @return Give the value of the log-pdf at x.
     *  @param x a real value
     *  @param mu, sigma location and scale of the log-normal law
     **/
    static Real lpdf( Real const& x, Real const& mu, Real const& sigma);
    /** @brief Compute the cumulative distribution function at t of
     *  the standard log-normal distribution.
     *  @param t a real value
     *  @param mu, sigma location and scale of the log-normal law
     *  @return the cumulative distribution function value at t
     **/
    static Real cdf( Real const& t, Real const& mu, Real const& sigma);
    /** @brief Compute the inverse cumulative distribution function at p
     *  of the standard log-normal distribution.
     *
     *  @param p a probability number.
     *  @param mu, sigma location and scale of the log-normal law
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
    /** The location parameter. **/
    Real mu_;
    /** The scale parameter. **/
    Real sigma_;
};

#ifdef IS_RTKPP_LIB

inline Real LogNormal::rand() const
{
#ifdef _OPENMP
//#pragma omp critical
#endif
GetRNGstate(); Real s = Rf_rlnorm(mu_, sigma_); PutRNGstate(); return s;
}

inline Real LogNormal::pdf( Real const& x) const { return Rf_dlnorm(x, mu_, sigma_, false);}
inline Real LogNormal::lpdf( Real const& x) const { return Rf_dlnorm(x, mu_, sigma_, true);}
inline Real LogNormal::cdf( Real const& t) const { return Rf_plnorm(t, mu_, sigma_, true, false);}
inline Real LogNormal::lcdf( Real const& t) const { return Rf_plnorm(t, mu_, sigma_, true, true);}
inline Real LogNormal::cdfc( Real const& t) const { return Rf_plnorm(t, mu_, sigma_, false, false);}
inline Real LogNormal::lcdfc( Real const& t) const { return Rf_plnorm(t, mu_, sigma_, false, true);}
inline Real LogNormal::icdf( Real const& p) const { return Rf_qlnorm(p, mu_, sigma_, true, false);}

inline Real LogNormal::rand( Real const& mu, Real const& sigma)
{
#ifdef _OPENMP
//#pragma omp critical
#endif
GetRNGstate(); Real s = Rf_rlnorm(mu, sigma); PutRNGstate(); return s;
}

inline Real LogNormal::pdf( Real const& x, Real const& mu, Real const& sigma)
{ return Rf_dlnorm(x, mu, sigma, false);}
inline Real LogNormal::lpdf( Real const& x, Real const& mu, Real const& sigma)
{ return Rf_dlnorm(x, mu, sigma, true);}
inline Real LogNormal::cdf( Real const& t, Real const& mu, Real const& sigma)
{ return Rf_plnorm(t, mu, sigma, true, false);}
inline Real LogNormal::icdf( Real const& p, Real const& mu, Real const& sigma)
{ return Rf_qlnorm(p, mu, sigma, true, false);}

#endif

} // namespace Law

} // namespace STK

#endif /*STK_LAW_LOGNORMAL_H*/
