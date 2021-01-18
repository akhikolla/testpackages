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
 * Purpose:  Beta probability distribution.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Law_Beta.h
 *  @brief In this file we define the Beta probability distribution.
 **/

#ifndef STK_LAW_BETA_H
#define STK_LAW_BETA_H

#include "STK_Law_IUnivLaw.h"
#include <Sdk/include/STK_Macros.h>

namespace STK
{

namespace Law
{
/** @ingroup Laws
 *  @brief Beta distribution law.
 *
 *  In probability theory and statistics, the <em> beta distribution </em> is a
 *  family of continuous probability distributions defined on the interval [0, 1]
 *  parameterized by two positive shape parameters that appear as exponents of
 *  the random variable and control the shape of the distribution.
 *
 *  The Beta distribution, is a continuous probability distribution with pdf
 *  \f[
 *  f(x) = \frac{x^{\alpha-1} (1-x)^{\beta-1} }
 *              {{\boldsymbol \beta}(\alpha,\beta)},
 *  \quad x \in [0,1],\ \mbox{ et }\ \alpha,\beta>0
 *  \f]
 *  where \f${\boldsymbol \beta}(a,b) = \frac{\Gamma(\alpha)\Gamma(\beta)}{\Gamma(\alpha+\beta)} \f$
 *  represents the beta function.
 *
 *  @sa STK::Funct::betaRatio
**/
class Beta: public IUnivLaw<Real>
{
  public:
    typedef IUnivLaw<Real> Base;
    /** default constructor. */
    inline Beta( Real const& alpha = .5, Real const& beta = .5)
              : Base(String(_T("Beta"))), alpha_(alpha), beta_(beta)
    {
      if ( !isFinite(alpha) || !isFinite(beta) || alpha <= 0.0 || beta <= 0.0)
        STKDOMAIN_ERROR_2ARG("Beta::Beta",alpha,beta,"argument error");
    }
    /** Dtor. */
    inline virtual ~Beta() {}
    /** @return the alpha value */
    inline Real alpha() const {return alpha_;}
    /** @return the beta value */
    inline Real beta() const { return beta_;}
    /** set the alpha value */
    inline void setAlpha(Real alpha) { alpha_ = alpha;}
    /** set the beta value */
    inline void setBeta(Real beta) {beta_ =beta;}

    /** @return a pseudo Beta random variate.
     *  This function use the Gamma::rand() random generator.
     *  TODO : implement the order statistics when a and b
     *  are small integers (Devroye p. 431).
     *  TODO : implement the rejection method from the normal
     *  pdf when a=b (Devroye p. 434).
     *  TODO : Implement Cheng's algorithm BA for beta pdf
     *  (Devroye p. 438)
     **/
    virtual Real rand() const;
    /** @return the value of the pdf at x. */
    virtual Real pdf( Real const& x) const;
    /** @return the value of the log-pdf at x. */
    virtual Real lpdf( Real const& x) const;
    /** @return The cumulative distribution function */
    virtual Real cdf( Real const& t) const;
    /** @return The inverse cumulative distribution */
    virtual Real icdf( Real const& p) const;

    /** @return a pseudo Beta random variate with the specified parameters. */
    static Real rand( Real const& alpha, Real const& beta);
    /** @return the value of the pdf at x. */
    static Real pdf( Real const& x, Real const& alpha, Real const& beta);
    /** @return the value of the log-pdf at x. */
    static Real lpdf( Real const& x, Real const& alpha, Real const& beta);
    /** @return The cumulative distribution function */
    static Real cdf( Real const& t, Real const& alpha, Real const& beta);
    /** @return The inverse cumulative distribution */
    static Real icdf( Real const& p, Real const& alpha, Real const& beta);

#ifdef IS_RTKPP_LIB
    /** @return log-cumulative distribution function */
    virtual Real lcdf( Real const& t) const;
    /** @return complement of cumulative distribution function */
    virtual Real cdfc( Real const& t) const;
    /** @return log-complement of cumulative distribution function */
    virtual Real lcdfc( Real const& t) const;
#endif

  protected:
    /** First parameter. */
    Real alpha_;
    /** Second parameter. */
    Real beta_;
};

#ifdef IS_RTKPP_LIB
inline Real Beta::rand() const
{
#ifdef _OPENMP
//#pragma omp critical
#endif
  GetRNGstate(); Real s = Rf_rbeta(alpha_, beta_); PutRNGstate(); return s;
}

inline Real Beta::pdf( Real const& x) const {  return Rf_dbeta(x,alpha_, beta_, false);}
inline Real Beta::lpdf( Real const& x) const { return Rf_dbeta(x,alpha_, beta_, true);}
inline Real Beta::cdf( Real const& t) const { return Rf_pbeta(t, alpha_, beta_, true, false);}
inline Real Beta::lcdf( Real const& t) const { return Rf_pbeta(t, alpha_, beta_, true, true);}
inline Real Beta::cdfc( Real const& t) const { return Rf_pbeta(t, alpha_, beta_, false, false);}
inline Real Beta::lcdfc( Real const& t) const { return Rf_pbeta(t, alpha_, beta_, false, true);}
inline Real Beta::icdf( Real const& p) const { return Rf_qbeta(p , alpha_, beta_, true, false);}

inline Real Beta::rand( Real const& alpha, Real const& beta)
{
#ifdef _OPENMP
//#pragma omp critical
#endif
  GetRNGstate(); Real s = Rf_rbeta(alpha, beta); PutRNGstate(); return s;
}

inline Real Beta::pdf(const Real& x, const Real& alpha, const Real& beta)
{ return Rf_dbeta(x,alpha, beta, (int)false);}
inline Real Beta::lpdf(const Real& x, const Real& alpha, const Real& beta)
{ return Rf_dbeta(x,alpha, beta, (int)true);}
inline Real Beta::cdf(const Real& t, const Real& alpha, const Real& beta)
{ return Rf_pbeta(t, alpha, beta, (int)true, (int)false);}
inline Real Beta::icdf(const Real& p, const Real& alpha, const Real& beta)
{ return Rf_qbeta(p , alpha, beta, (int)true, (int)false);}

#endif
} // namespace Law

} // namespace STK

#endif /*STK_LAWBETA_H*/
