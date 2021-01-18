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
 * Purpose:  Gamma probability distribution.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Law_Gamma.h
 *  @brief In this file we define the Gamma probability distribution.
 **/

#ifndef STK_LAW_GAMMA_H
#define STK_LAW_GAMMA_H

#include "STK_Law_IUnivLaw.h"
#include <Sdk/include/STK_Macros.h>

namespace STK
{

namespace Law
{
/** @ingroup Laws
 *  @brief Gamma distribution law.
 *
 *  In probability theory and statistics, the <em>gamma distribution</em>
 *  is a two-parameter family of continuous probability distributions.
 *  The common exponential distribution and chi-squared distribution are special
 *  cases of the gamma distribution.
 *
 *  The probability pdf function of the gamma distribution
 *  can be expressed in terms of the @ref Funct::gamma function:
 *  \f[
 *   f(x;a,b) = \left(\frac{x}{b}\right)^{a-1}
 *                   \frac{e^{-x/b}}{b \, \Gamma(a)}
 *   \ \mathrm{for}\ x>0,\ a>0,\ b>0
 *  \f]
 *  where @a a is the shape parameter and @e b is the scale parameter.
 **/
class Gamma: public IUnivLaw<Real>
{
  public:
    typedef IUnivLaw<Real> Base;
    /** Default constructor.
     *  @param shape shape (position) parameter
     *  @param scale scale (dispersion) parameter
     **/
    Gamma( Real const& shape = 1., Real const& scale = 1.)
        : Base(_T("Gamma")), a_(shape), b_(scale)
#ifndef IS_RTKPP_LIB
         , c_(), d_()
#endif
    {
      // check parameters
      if ( !Arithmetic<Real>::isFinite(a_) || !Arithmetic<Real>::isFinite(b_)
         || a_ <= 0 || b_ <= 0
         )
      { STKDOMAIN_ERROR_2ARG(Gamma::Gamma,a_, b_,arguments not valid);}
#ifndef IS_RTKPP_LIB
      computeCtes();
#endif
    }
    /** destructor */
    inline virtual ~Gamma() {}
    /** @return shape */
    inline Real const& shape() const { return a_;}
    /** @return scale */
    inline Real const& scale() const { return b_;}
    /** @param shape the shape parameter */
    inline void setShape( Real const& shape)
    {
      if (shape<=0) STKDOMAIN_ERROR_1ARG(Gamma::setShape,shape,shape must be > 0);
      a_ = shape;
#ifndef IS_RTKPP_LIB
      computeCtes();
#endif
    }
    /** @param scale the scale parameter */
    inline void setScale( Real const& scale)
    {
      if (scale<=0) STKDOMAIN_ERROR_1ARG(Gamma::setScale,scale,scale must be > 0);
      b_ = scale;
    }
    /** @brief generate a gamma random variate using the G.S algorithm
     *  of Ahrens and Dieter (1974) for 0<a_<1
     *  and Marsaglia fast Gamma generator with enhanced squeeze step for a>1.
     *  @return a pseudo Gamma random variate.
     **/
    virtual Real rand() const;
    /** compute
     *  \f[
     *   f(x;\alpha,\beta) = \left(\frac{x}{\beta}\right)^{\alpha-1}
     *                   \frac{e^{-x/\beta}}{\beta \, \Gamma(\alpha)}
     *   \ \mathrm{for}\ x > 0
     *  \f]
     *  where \f$ \alpha >, \mbox{ et } \beta > 0\f$ are the shape and scale
     *  parameters.
     *  @return the value of the pdf
     *  @param x a positive real value
     **/
    virtual Real pdf(Real const& x) const;
    /** Compute
     *  \f[
     *   \ln(f(x;\alpha,\beta)) = - x/\beta + (\alpha-1) \ln(x)
     *                            - \alpha \ln(\beta) + \ln(\Gamma(\alpha))
     *  \f]
     *  @return the value of the log-pdf
     *  @param x a positive real value
     **/
    virtual Real lpdf(Real const& x) const;
    /** @return the cumulative distribution function
     *  @param t a positive real value
     **/
    virtual Real cdf(Real const& t) const;
    /** @return the inverse cumulative distribution function
     *  @param p a probability number
     **/
    virtual Real icdf(Real const& p) const;

    /** @return a pseudo Gamma random variate with the specified parameters.
     *  @param shape, scale shape and scale parameters
     **/
    static Real rand( Real const& shape, Real const& scale);
    /** @return the value of the pdf
     *  @param x a positive real value
     *  @param shape, scale shape and scale parameters
     **/
    static Real pdf(Real const& x, Real const& shape, Real const& scale);
    /** @return the value of the log-pdf
     *  @param x a positive real value
     *  @param shape, scale shape and scale parameters
     **/
    static Real lpdf(Real const& x, Real const& shape, Real const& scale);
    /** @return the cumulative distribution function
     *  @param t a positive real value
     *  @param shape, scale shape and scale parameters
     **/
    static Real cdf(Real const& t, Real const& shape, Real const& scale);
    /** @return the inverse cumulative distribution function
     *  @param p a probability number
     *  @param shape, scale shape and scale parameters
     **/
    static Real icdf(Real const& p, Real const& shape, Real const& scale);

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
    Real a_;
    /** The scale parameter */
    Real b_;

#ifndef IS_RTKPP_LIB
  private:
    /** First and second constants for rand */
    mutable Real c_, d_;
    /** compute c_ and d_ */
    void computeCtes() const;
#endif
};

#ifdef IS_RTKPP_LIB

inline Real Gamma::rand() const
{
#ifdef _OPENMP
//#pragma omp critical
#endif
GetRNGstate(); Real s = Rf_rgamma(a_, b_); PutRNGstate(); return s;
}

inline Real Gamma::pdf( Real const& x) const { return Rf_dgamma(x, a_, b_, false);}
inline Real Gamma::lpdf( Real const& x) const { return Rf_dgamma(x, a_, b_, true);}
inline Real Gamma::cdf( Real const& t) const { return Rf_pgamma(t, a_, b_, true, false);}
inline Real Gamma::lcdf( Real const& t) const { return Rf_pgamma(t, a_, b_, true, true);}
inline Real Gamma::cdfc( Real const& t) const { return Rf_pgamma(t, a_, b_, false, false);}
inline Real Gamma::lcdfc( Real const& t) const { return Rf_pgamma(t, a_, b_, false, true);}
inline Real Gamma::icdf( Real const& p) const { return Rf_qgamma(p, a_, b_, true, false);}

inline Real Gamma::rand( Real const& a, Real const& b)
{
#ifdef _OPENMP
//#pragma omp critical
#endif
GetRNGstate(); Real s = Rf_rgamma(a, b); PutRNGstate(); return s;
}

inline Real Gamma::pdf( Real const& x, Real const& a, Real const& b)
{ return Rf_dgamma(x, a, b, false);}
inline Real Gamma::lpdf( Real const& x, Real const& a, Real const& b)
{ return Rf_dgamma(x, a, b, true);}
inline Real Gamma::cdf(Real const& t, Real const& a, Real const& b)
{ return Rf_pgamma(t, a, b, true, false);}
inline Real Gamma::icdf(Real const& p, Real const& a, Real const& b)
{ return Rf_qgamma(p, a, b, true, false);}

#endif

} // namespace Law

} // namespace STK

#endif /*STK_LAWGAMMA_H*/
