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
 * Purpose: Logistic probability distribution.
 * Author:  Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Law_Logistic.h
 *  @brief In this file we define the Logistic probability law class.
 **/

#ifndef STK_LAW_LOGISTIC_H
#define STK_LAW_LOGISTIC_H

#include "STK_Law_IUnivLaw.h"
#include <Sdk/include/STK_Macros.h>

namespace STK
{

namespace Law
{
/**
 *  @ingroup Laws
 *  @brief Logistic distribution law.
 *
 * In probability theory and statistics, the <em>logistic distribution</em> is a
 * continuous probability distribution. Its cumulative distribution function is
 * the logistic function, which appears in logistic regression and feedforward
 * neural networks. It resembles the logistic distribution in shape but has heavier
 * tails (higher kurtosis).
 *
 * The Logistic distribution with <em>location = m</em> and <em>scale = s>0</em>
 * has distribution function
 * \f[
 *   F(x) = 1 / (1 + exp(-(x-m)/s))
 * \f]
 * and density
 * \f[
 * f(x) = \frac{1}{s} \frac{exp\left(\frac{x-m}{s}\right)}
 *                         {\left(1 + exp\left(\frac{x-m}{s}\right)\right)^2}.
 * \f]
 * It is a long-tailed distribution with mean @e m and variance \f$ \frac{π^2s^2}{3} \f$.
**/
class Logistic: public IUnivLaw<Real>
{
  public:
    typedef IUnivLaw<Real> Base;
    /** Constructor.
     *  @param mu mean of the Logistic distribution
     *  @param scale scale of the Logistic distribution
     **/
    inline Logistic( Real const& mu=0., Real const& scale=1.)
                  : Base(_T("Logistic")), mu_(mu), scale_(scale) {}
    /** Destructor. **/
    inline virtual ~Logistic() {}
    /** @return mu */
    inline Real const& mu() const { return mu_;}
    /** @return mu */
    inline Real const& scale() const { return scale_;}
    /** @param mu the value to set to mu */
    inline void setMu( Real const& mu) { mu_ = mu;}
    /** @param scale the value to set to scale */
    inline void setScale( Real const& scale)
    {
      if (scale<=0) STKDOMAIN_ERROR_1ARG(Logistic::setScale,scale,scale must be > 0);
       scale_ = scale;
    }
    /** @brief Generate a pseudo Logistic random variate.
     *
     *  Generate a pseudo Logistic random variate
     *  with location parameter @c mu_ and scale @c scale_.
     *  @return a pseudo logistic random variate
     **/
    Real rand() const;
    /** @param x a real value
     *  @return the value of the logistic pdf at @c x
     **/
    virtual Real pdf( Real const& x) const;
    /** @return Give the value of the log-pdf at x.
     *  @param x a real value
     **/
    virtual Real lpdf( Real const& x) const;
    /** @brief Compute the cumulative distribution function at t of
     *  the standard logistic distribution.
     *
     *  The cumulative distribution function of the logistic distribution is
     *  also a scaled version of the Hyperbolic function.
     *  \f[
     *   F(t; \mu, s) = \frac{1}{1+e^{-\frac{t-\mu}{s}}}
     *   = \frac{1}{2} + \frac{1}{2} \;\operatorname{tanh}\!\left(\frac{t-\mu}{2s}\right).
     *   \f]
     *
     *  @param t a real value
     *  @return the cumulative distribution function value at t
     **/
    virtual Real cdf( Real const& t) const;
    /** @brief Compute the inverse cumulative distribution function at p
     *  of the standard logistic distribution.
     *
     *  The inverse cumulative distribution function (quantile function) of the
     *  logistic distribution is a generalization of the logit function.
     *  It is defined as follows:
     *  \f[
     *      Q(p;\mu,s) = \mu + s\,\ln\left(\frac{p}{1-p}\right).
     *  \f]
     *  @param p a probability number.
     *  @return the inverse cumulative distribution function value at p.
     **/
    virtual Real icdf( Real const& p) const;

    /** @brief Generate a pseudo Logistic random variate  with location @c mu and
     *  scale @c scale parameters.  scale @c scale parameters.
     *  @param mu, scale location and scale of the Logistic distribution
     *  @return a pseudo logistic random variate, centered in @c mu and with
     *  scale @c scale
     **/
    static Real rand( Real const& mu, Real const& scale);
    /** @param x a real value
     *  @param mu, scale location and scale of the Logistic distribution
     *  @return the value of the logistic pdf at @c x
     **/
    static Real pdf( Real const& x, Real const& mu, Real const& scale);
    /** @return Give the value of the log-pdf at x.
     *  @param x a real value
     *  @param mu, scale location and scale of the Logistic distribution
     **/
    static Real lpdf( Real const& x, Real const& mu, Real const& scale);
    /** @return the cumulative distribution function value at t.
     *  @param t a real value
     *  @param mu, scale location and scale of the Logistic distribution
     **/
    static Real cdf( Real const& t, Real const& mu, Real const& scale);
    /** @return the inverse cumulative distribution function value at p.
     *  @param p a probability number.
     *  @param mu, scale location and scale of the Logistic distribution
     **/
    static Real icdf( Real const& p, Real const& mu, Real const& scale);

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
    /** The scale parameter. **/
    Real scale_;
};

#ifdef IS_RTKPP_LIB

inline Real Logistic::rand() const
{
#ifdef _OPENMP
//#pragma omp critical
#endif
GetRNGstate(); Real s = Rf_rlogis(mu_, scale_); PutRNGstate(); return s;
}

inline Real Logistic::pdf( Real const& x) const { return Rf_dlogis(x, mu_, scale_, false);}
inline Real Logistic::lpdf( Real const& x) const { return Rf_dlogis(x, mu_, scale_, true);}
inline Real Logistic::cdf( Real const& t) const { return Rf_plogis(t, mu_, scale_, true, false);}
inline Real Logistic::lcdf( Real const& t) const { return Rf_plogis(t, mu_, scale_, true, true);}
inline Real Logistic::cdfc( Real const& t) const { return Rf_plogis(t, mu_, scale_, false, false);}
inline Real Logistic::lcdfc( Real const& t) const { return Rf_plogis(t, mu_, scale_, false, true);}
inline Real Logistic::icdf( Real const& p) const { return Rf_qlogis(p, mu_, scale_, true, false);}

inline Real Logistic::rand( Real const& mu, Real const& scale)
{
#ifdef _OPENMP
//#pragma omp critical
#endif
GetRNGstate(); Real s = Rf_rlogis(mu, scale); PutRNGstate(); return s;
}

inline Real Logistic::pdf( Real const& x, Real const& mu, Real const& scale)
{ return Rf_dlogis(x, mu, scale, false);}
inline Real Logistic::lpdf( Real const& x, Real const& mu, Real const& scale)
{ return Rf_dlogis(x, mu, scale, true);}
inline Real Logistic::cdf( Real const& t, Real const& mu, Real const& scale)
{ return Rf_plogis(t, mu, scale, true, false);}
inline Real Logistic::icdf( Real const& p, Real const& mu, Real const& scale)
{ return Rf_qlogis(p, mu, scale, true, false);}

#endif

} // namespace Law

} // namespace STK

#endif /*STK_LAW_LOGISTIC_H*/
