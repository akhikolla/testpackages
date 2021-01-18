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
 * Purpose:  Weibull probability distribution.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Law_Weibull.h
 *  @brief In this file we define the Weibull probability distribution.
 **/

#ifndef STK_LAW_WEIBULL_H
#define STK_LAW_WEIBULL_H

#include "STK_Law_IUnivLaw.h"
#include <Sdk/include/STK_Macros.h>

namespace STK
{

namespace Law
{
/** @ingroup Laws
 *  @brief Weibull distribution law.
 *
 *  In probability theory and statistics, the <em>Weibull distribution</em> is a
 *  continuous probability distribution. It is named after Waloddi Weibull, who
 *  described it in detail in 1951.
 *
 * The probability density function of a Weibull random variable is
 * \f[
 *  f(x;\lambda,k) =
 *  \frac{k}{\lambda}\left(\frac{x}{\lambda}\right)^{k-1}e^{-(x/\lambda)^{k}} \quad x\geq0 ,
 *  \f]
 *  where <em>k > 0</em> is the shape parameter and \f$ λ > 0\f$ is the scale parameter
 *  of the distribution.
 **/
class Weibull: public IUnivLaw<Real>
{
  public:
    typedef IUnivLaw<Real> Base;
    /** Default constructor.
     *  @param k, lambda shape and scale (dispersion) parameters
     **/
    inline Weibull( Real const& k = 1., Real const& lambda = 1.)
                 : Base(_T("Weibull")), k_(k), lambda_(lambda)
    {}
    /** destructor */
    inline virtual ~Weibull() {}
    /** @return shape parameter */
    inline Real const& k() const { return k_;}
    /** @return scale parameter */
    inline Real const& lambda() const { return lambda_;}
    /** @param k set the shape parameter */
    inline void setK( Real const& k)
    {
      if (k<=0) STKDOMAIN_ERROR_1ARG(Weibull::setShape,k,shape must be > 0);
      k_ = k;
    }
    /** @param lambda set the scale parameter */
    inline void setLambda( Real const& lambda)
    {
      if (lambda<=0) STKDOMAIN_ERROR_1ARG(Weibull::setScale,lambda,scale must be > 0);
      lambda_ = lambda;
    }
    /** @return a pseudo Weibull random variate. */
    virtual Real rand() const;
    /** @return the value of the pdf
     *  @param x a positive real value
     **/
    virtual Real pdf(Real const& x) const;
    /** @return the value of the log-pdf
     *  @param x a positive real value
     **/
    virtual Real lpdf(Real const& x) const;
    /**The cumulative distribution function for the Weibull distribution is
     *  \f$ F(x;k,\lambda) = 1- e^{-(x/\lambda)^k}.\f$
     *  @return the cumulative distribution function
     *  @param t a positive real value
     **/
    virtual Real cdf(Real const& t) const;
    /**The quantile (inverse cumulative distribution) function for the Weibull
     * distribution is \f$ Q(p;k,\lambda) = \lambda {(-\ln(1-p))}^{1/k} \f$
     *  @return the inverse cumulative distribution function
     *  @param p a probability number
     **/
    virtual Real icdf(Real const& p) const;
    /** @return a pseudo Weibull random variate with the specified parameters.
     *  @param k, lambda shape and scale (dispersion) parameters
     **/
    static Real rand( Real const& k, Real const& lambda);
    /** @return the value of the pdf
     *  @param x a positive real value
     *  @param k, lambda shape and scale (dispersion) parameters
     **/
    static Real pdf(Real const& x, Real const& k, Real const& lambda);
    /** @return the value of the log-pdf
     *  @param x a positive real value
     *  @param k, lambda shape and scale (dispersion) parameters
     **/
    static Real lpdf(Real const& x, Real const& k, Real const& lambda);
    /** @return the cumulative distribution function
     *  @param t a positive real value
     *  @param k, lambda shape and scale (dispersion) parameters
     **/
    static Real cdf(Real const& t, Real const& k, Real const& lambda);
    /** @brief Compute the inverse cumulative distribution function at p
     *  of the standard log-normal distribution.
     *
     *  @param p a probability number.
     *  @param k, lambda shape and scale (dispersion) parameters
     *  @return the inverse cumulative distribution function value at p.
     **/
    static Real icdf(Real const& p, Real const& k, Real const& lambda);

#ifdef IS_RTKPP_LIB
    /** @return log-cumulative distribution function */
    virtual Real lcdf( Real const& t) const;
    /** @return complement of cumulative distribution function */
    virtual Real cdfc( Real const& t) const;
    /** @return log-complement of cumulative distribution function */
    virtual Real lcdfc( Real const& t) const;
#endif

  protected:
    /** The shape parameter */
    Real k_;
    /** The scale parameter */
    Real lambda_;
};

#ifdef IS_RTKPP_LIB

inline Real Weibull::rand() const
{
#ifdef _OPENMP
//#pragma omp critical
#endif
GetRNGstate(); Real s = Rf_rweibull(k_, lambda_); PutRNGstate(); return s;
}

inline Real Weibull::pdf(Real const& x) const { return Rf_dweibull(x, k_, lambda_, false);}
inline Real Weibull::lpdf(Real const& x) const { return Rf_dweibull(x, k_, lambda_, true);}
inline Real Weibull::cdf(Real const& t) const { return Rf_pweibull(t, k_, lambda_, true, false);}
inline Real Weibull::lcdf(Real const& t) const { return Rf_pweibull(t, k_, lambda_, true, true);}
inline Real Weibull::cdfc(Real const& t) const { return Rf_pweibull(t, k_, lambda_, false, false);}
inline Real Weibull::lcdfc(Real const& t) const { return Rf_pweibull(t, k_, lambda_, false, true);}
inline Real Weibull::icdf(Real const& p) const { return Rf_qweibull(p, k_, lambda_, true, false);}

inline Real Weibull::rand( Real const& k, Real const& lambda)
{
#ifdef _OPENMP
//#pragma omp critical
#endif
GetRNGstate(); Real s = Rf_rweibull(k, lambda); PutRNGstate(); return s;
}

inline Real Weibull::pdf(Real const& x, Real const& k, Real const& lambda)
{ return Rf_dweibull(x, k, lambda, false);}
inline Real Weibull::lpdf(Real const& x, Real const& k, Real const& lambda)
{ return Rf_dweibull(x, k, lambda, true);}
inline Real Weibull::cdf(Real const& t, Real const& k, Real const& lambda)
{ return Rf_pweibull(t, k, lambda, true, false);}
inline Real Weibull::icdf( Real const& p, Real const& k, Real const& lambda)
{ return Rf_qweibull(p, k, lambda, true, false);}

#endif

} // namespace Law

} // namespace STK

#endif /*STK_LAWWEIBULL_H*/
