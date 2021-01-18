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
 * Purpose:  ChiSquared probability distribution.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Law_ChiSquared.h
 *  @brief In this file we define the ChiSquared probability distribution.
 **/

#ifndef STK_LAW_CHISQUARED_H
#define STK_LAW_CHISQUARED_H

#include "STK_Law_IUnivLaw.h"
#include <Sdk/include/STK_Macros.h>

namespace STK
{

namespace Law
{
/** @ingroup Laws
 *  @brief ChiSquared distribution law.
 *
 *  In probability theory and statistics, the <em>chi-squared distribution</em>
 *  with k degrees of freedom is the distribution of a sum of the squares of k
 *  independent standard normal random variables. It is a special case of the
 *  gamma distribution and is one of the most widely used probability
 *  distributions in inferential statistics, e.g., in hypothesis testing or in
 *  construction of confidence intervals.
 *
 *  The probability density function (pdf) of the chi-squared distribution with
 *  @e n degree of freedom is
 *  \f[
 *   f(x;n) = \frac{ x^(n/2-1) e^(-x/2)}{2^(n/2) \Gamma(n/2)}
 *   \ \mathrm{ for }\ x > 0\ \mathrm{ and }\ n \geq 0.
 *  \f]
 **/
class ChiSquared: public IUnivLaw<Real>
{
  public:
    typedef IUnivLaw<Real> Base;
    /** Default constructor.
     *  @param df degree of freedom parameter
     **/
    inline ChiSquared( int df = 1): Base(_T("Chi-squared")), df_(df)
    { if (df<=0) STKDOMAIN_ERROR_1ARG(ChiSquared::ChiSquared,df,df must be > 0);}
    /** destructor */
    inline virtual ~ChiSquared(){}
    /** @return the number of degree of freedom */
    inline int df() const { return df_;}
    /** @param df degree of freedom parameter */
    inline void setDf( int df)
    {
      if (df<=0) STKDOMAIN_ERROR_1ARG(ChiSquared::setShape,df,shape must be > 0);
      df_ = df;
    }
    /** @return a pseudo ChiSquared random variate. */
    virtual Real rand() const;
    /** @return the value of the pdf
     *  @param x a positive real value
     **/
    virtual Real pdf(Real const& x) const;
    /** @return the value of the log-pdf
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

    /** @return a pseudo ChiSquared random variate with the specified parameters.
     *  @param df degree of freedom parameter
     **/
    static Real rand( int df);
    /** @return the value of the pdf
     *  @param x a positive real value
     *  @param df degree of freedom parameter
     **/
    static Real pdf(Real const& x, int df);
    /** @return the value of the log-pdf
     *  @param x a positive real value
     *  @param df degree of freedom parameter
     **/
    static Real lpdf(Real const& x, int df);
    /** @return the cumulative distribution function
     *  @param t a positive real value
     *  @param df degree of freedom parameter
     **/
    static Real cdf(Real const& t, int df);
    /** @return the inverse cumulative distribution function
     *  @param p a probability number
     *  @param df degree of freedom parameter
     **/
    static Real icdf(Real const& p, int df);

#ifdef IS_RTKPP_LIB
    /** @return log-cumulative distribution function */
    virtual Real lcdf( Real const& t) const;
    /** @return complement of cumulative distribution function */
    virtual Real cdfc( Real const& t) const;
    /** @return log-complement of cumulative distribution function */
    virtual Real lcdfc( Real const& t) const;
#endif

  protected:
    /** degree of freedom */
    int df_;
};

#ifdef IS_RTKPP_LIB

/* @return a pseudo ChiSquared random variate. */
inline Real ChiSquared::rand() const
{
#ifdef _OPENMP
//#pragma omp critical
#endif
 GetRNGstate(); Real s = Rf_rchisq(df_); PutRNGstate(); return s;
}

inline Real ChiSquared::pdf(const Real& x) const { return Rf_dchisq(x, df_, false);}
inline Real ChiSquared::lpdf(const Real& x) const { return Rf_dchisq(x, df_, true);}
inline Real ChiSquared::cdf(const Real& t) const { return Rf_pchisq(t, df_, true, false);}
inline Real ChiSquared::lcdf(const Real& t) const { return Rf_pchisq(t, df_, true, true);}
inline Real ChiSquared::cdfc(const Real& t) const { return Rf_pchisq(t, df_, false, false);}
inline Real ChiSquared::lcdfc(const Real& t) const { return Rf_pchisq(t, df_, false, true);}
inline Real ChiSquared::icdf(const Real& p) const { return Rf_qchisq(p, df_, true, false);}

inline Real ChiSquared::rand(int df)
{
#ifdef _OPENMP
//#pragma omp critical
#endif
GetRNGstate(); Real s = Rf_rchisq(df); PutRNGstate(); return s;
}

inline Real ChiSquared::pdf(const Real& x, int df) { return Rf_dchisq(x, df, false);}
inline Real ChiSquared::lpdf(const Real& x, int df) { return Rf_dchisq(x, df, true);}
inline Real ChiSquared::cdf(const Real& t, int df) { return Rf_pchisq(t, df, true, false);}
inline Real ChiSquared::icdf(const Real& p, int df) { return Rf_qchisq(p, df, true, false);}

#endif

} // namespace Law

} // namespace STK

#endif /*STK_LAW_CHISQUARED_H*/
