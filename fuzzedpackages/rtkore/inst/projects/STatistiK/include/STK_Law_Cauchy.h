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
 * Purpose:  Cauchy probability distribution.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Law_Cauchy.h
 *  @brief In this file we define the Cauchy probability distribution.
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
 *  @brief Cauchy distribution law.
 *
 *  The <em>Cauchy distribution</em>, named after Augustin Cauchy, is a
 *  continuous probability distribution. The Cauchy distribution has the
 *  probability density function
 *  \f[
 *  f(x; \mu,\gamma) =
 *  \frac{1}{\pi\gamma \left[1 + \left(\frac{x - \mu}{\gamma}\right)^2\right]}
 *   = \frac{ 1}{\pi \gamma }\left[ \frac{\gamma^2}{(x - \mu)^2 + \gamma^2} \right].
 *  \f]
 *  where \f$ \mu  \f$ is the location parameter, specifying the location
 *  of the peak of the distribution, and \f$ \gamma>0 \f$ is the scale
 *  parameter.
 *
 *  The simplest Cauchy distribution is called the standard Cauchy distribution.
 *  It is the distribution of a random variable that is the ratio of two
 *  independent standard normal variables and has the probability density
 *  function \f$  f(x;0,1) \f$.
 **/
class Cauchy: public IUnivLaw<Real>
{
  public:
    typedef IUnivLaw<Real> Base;
    /** Default constructor.
     *  @param mu, scale location and scale of the Cauchy distribution
     **/
    inline Cauchy( Real const& mu=0, Real const& scale=1)
                : Base(String(_T("Cauchy")))
                 , mu_(mu)
                 , scale_(scale)
    {
      // check parameters
      if (  !Arithmetic<Real>::isFinite(mu) || !Arithmetic<Real>::isFinite(scale) || scale <= 0)
      STKDOMAIN_ERROR_2ARG(Cauchy::Cauchy,mu, scale,argument error);
    }

    /** Destructor. **/
    inline virtual ~Cauchy() {}
    /** @return the mu parameter */
    inline Real const& mu() const { return mu_;}
    /** @return the scale parameter */
    inline Real const& scale() const { return scale_;}
    /** @param mu the mu parameter */
    inline void setMu( Real const& mu) { mu_ = mu;}
    /** @param scale the scale parameter */
    inline void setScale( Real const& scale)
    {
      if (scale<0) STKRUNTIME_ERROR_1ARG(Cauchy::setScale,scale,scale must be >= 0);
       scale_ = scale;
    }
    /** Generate a pseudo Cauchy random variate.
     * @return a cauchy random variable **/
    virtual Real rand() const;
    /** @param x a real value
     *  @return the pdf of the cauchy distribution at x
     **/
    virtual Real pdf( Real const& x) const;
    /** @param x a real value
     *  @return the log-pdf of the cauchy distribution at x
     **/
    virtual Real lpdf( Real const& x) const;
    /** The cumulative distribution function of the Cauchy distribution at t is
     * \f[
     *  F(t; \mu,\gamma)=
     * \frac{1}{\pi} \arctan\left(\frac{t-\mu}{\gamma}\right)+\frac{1}{2}
     * \f]
     * @see http://www.faqs.org/faqs/fr/maths/maths-faq-3/
     * for the computation of arctan.
     *
     *  @param t a real value
     *  @return the cdf of the Cauchy distribution at t
     **/
    virtual Real cdf( Real const& t) const;
    /** The inverse cumulative distribution function at p is
     * \f[
     * F^{-1}(p; \mu,\gamma) = \mu + \gamma \tan(\pi (p-1/2)).
     * \f]
     *
     *  @param p a probability number
     *  @return the inverse cdf of the Cauchy distribution at p
     **/
    virtual Real icdf( Real const& p) const;

    /** Generate a pseudo Cauchy random variate with the specified
     *  parameters.
     *  @param mu, scale location and scale of the Cauchy distribution
     *  @return a cauchy random variable
     **/
    static Real rand( Real const& mu, Real const& scale);
    /** @param x a real value
     *  @param mu, scale location and scale of the Cauchy distribution
     *  @return the pdf of the cauchy distribution at x
     **/
    static Real pdf( Real const& x, Real const& mu, Real const& scale);
    /** @param x a real value
     *  @param mu mu of the Cauchy distribution
     *  @param scale scale of the Cauchy distribution
     *  @return the log-pdf of the cauchy distribution at x
     **/
    static Real lpdf( Real const& x, Real const& mu, Real const& scale);
    /** The cumulative distribution function of the Cauchy distribution at t is
     * \f[
     *  F(t; \mu,\gamma)=
     * \frac{1}{\pi} \arctan\left(\frac{t-\mu}{\gamma}\right)+\frac{1}{2}
     * \f]
     * @see http://www.faqs.org/faqs/fr/maths/maths-faq-3/
     * for the computation of arctan.
     *
     *  @param t a real value
     *  @param mu, scale location and scale of the Cauchy distribution
     *  @return the cdf of the Cauchy distribution at t
     **/
    static Real cdf( Real const& t, Real const& mu, Real const& scale);
    /** The inverse cumulative distribution function at p is
     * \f[
     * F^{-1}(p; \mu,\gamma) = \mu + \gamma \tan(\pi (p-1/2)).
     * \f]
     *
     *  @param p a probability number
     *  @param mu, scale location and scale of the Cauchy distribution
     *  @return the inverse cdf of the Cauchy distribution at p
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
    /** The mu parameter. */
    Real mu_;
    /** The scale parameter. */
    Real scale_;
};

#ifdef IS_RTKPP_LIB

/*  Generate a pseudo Cauchy random variate. */
inline Real Cauchy::rand() const
{
#ifdef _OPENMP
//#pragma omp critical
#endif
GetRNGstate(); Real s = Rf_rcauchy(mu_, scale_); PutRNGstate(); return s;
}

inline Real Cauchy::pdf( Real const& x) const { return Rf_dcauchy(x,mu_, scale_, false);}
inline Real Cauchy::lpdf( Real const& x) const { return Rf_dcauchy(x,mu_, scale_, true);}
inline Real Cauchy::cdf( Real const& t) const { return Rf_pcauchy(t, mu_, scale_, true, false);}
inline Real Cauchy::lcdf( Real const& t) const { return Rf_pcauchy(t, mu_, scale_, true, true);}
inline Real Cauchy::cdfc( Real const& t) const { return Rf_pcauchy(t, mu_, scale_, false, false);}
inline Real Cauchy::lcdfc( Real const& t) const { return Rf_pcauchy(t, mu_, scale_, false, true);}
inline Real Cauchy::icdf( Real const& p) const { return Rf_qcauchy(p , mu_, scale_, true, false);}

// static
inline Real Cauchy::rand( Real const& mu, Real const& scale)
{
#ifdef _OPENMP
//#pragma omp critical
#endif
GetRNGstate(); Real s = Rf_rcauchy(mu, scale); PutRNGstate(); return s;
}

inline Real Cauchy::pdf(Real const& x, Real const& mu, Real const& scale)
{ return Rf_dcauchy(x,mu, scale, false);}
inline Real Cauchy::lpdf(Real const& x, Real const& mu, Real const& scale)
{ return Rf_dcauchy(x,mu, scale, true);}
inline Real Cauchy::cdf(Real const& t, Real const& mu, Real const& scale)
{ return Rf_pcauchy(t, mu, scale, true, false);}
inline Real Cauchy::icdf(Real const& p, Real const& mu, Real const& scale)
{ return Rf_qcauchy(p , mu, scale, true, false);}

#endif

} // namespace Law

} // namespace STK

#endif /*STK_LAWCAUCHY_H*/
