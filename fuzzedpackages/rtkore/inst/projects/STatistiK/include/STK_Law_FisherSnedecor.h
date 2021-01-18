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
 * Purpose:  FisherSnedecor probability distribution.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Law_FisherSnedecor.h
 *  @brief In this file we define the FisherSnedecor probability distribution.
 **/

#ifndef STK_LAW_FISHERSNEDECOR
#define STK_LAW_FISHERSNEDECOR

#include "STK_Law_IUnivLaw.h"
#include <Sdk/include/STK_Macros.h>

namespace STK
{

namespace Law
{
/** @ingroup Laws
 *  @brief FisherSnedecor distribution law.
 *
 *  In probability theory and statistics, the <em>F-distribution</em> is a
 *  continuous probability distribution. It is also known as Snedecor's F
 *  distribution or the Fisher–Snedecor distribution (after R. A. Fisher and
 *  George W. Snedecor). The F-distribution arises frequently as the null
 *  distribution of a test statistic, most notably in the analysis of variance.
 *
 *  If a random variable @e X has an F-distribution with parameters @e d1 and
 *  @e d2, we write \f$ X \sim F(d1, d2)\f$. Then the probability density
 *  function (pdf) for @e X is given by
 *  \f[
 *  f(x; d_1,d_2) = \frac{\sqrt{\frac{(d_1\,x)^{d_1}\,\,d_2^{d_2}} {(d_1\,x+d_2)^{d_1+d_2}}}}
 *                  {x\,\mathrm{B}\!\left(\frac{d_1}{2},\frac{d_2}{2}\right)}
 *  \f]
 **/
class FisherSnedecor: public IUnivLaw<Real>
{
  public:
    typedef IUnivLaw<Real> Base;
    /** Default constructor.
     *  @param df1, df2 degree of freedom parameters
     **/
    inline FisherSnedecor( int df1 = 1., int df2 = 1)
                        : Base(_T("Fisher-Snedecor")), df1_(df1), df2_(df2)
    {};
    /** destructor */
    inline virtual ~FisherSnedecor() {}
    /** @return first degree of freedom parameter */
    inline int df1() const { return df1_;}
    /** @return second degree of freedom parameter */
    inline int df2() const { return df2_;}
    /** @param df1 set first degree of freedom parameter */
    inline void setDf1( int df1)
    {
      if (df1<=0) STKDOMAIN_ERROR_1ARG(FisherSnedecor::setDf1,df1,degree of freedom must be > 0);
      df1_ = df1;
    }
    /** @param df2 set second degree of freedom parameter */
    inline void setDf2( int df2)
    {
      if (df2<=0) STKDOMAIN_ERROR_1ARG(FisherSnedecor::setDf2,df2,degree of freedom must be > 0);
      df2_ = df2;
    }
    /** @return a pseudo FisherSnedecor random variate. */
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

    /** @return a pseudo FisherSnedecor random variate with the specified parameters.
     *  @param df1, df2 degree of freedom parameters
     **/
    static Real rand( int df1, int df2);
    /** @return the value of the pdf
     *  @param x a positive real value
     *  @param df1, df2 degree of freedom parameters
     **/
    static Real pdf(Real const& x, int df1, int df2);
    /** @return the value of the log-pdf
     *  @param x a positive real value
     *  @param df1, df2 degree of freedom parameters
     **/
    static Real lpdf(Real const& x, int df1, int df2);
    /** @return the cumulative distribution function
     *  @param t a positive real value
     *  @param df1, df2 degree of freedom parameters
     **/
    static Real cdf(Real const& t, int df1, int df2);
    /** @return the inverse cumulative distribution function
     *  @param p a probability number
     *  @param df1, df2 degree of freedom parameters
     **/
    static Real icdf(Real const& p, int df1, int df2);

#ifdef IS_RTKPP_LIB
    /** @return log-cumulative distribution function */
    virtual Real lcdf( Real const& t) const;
    /** @return complement of cumulative distribution function */
    virtual Real cdfc( Real const& t) const;
    /** @return log-complement of cumulative distribution function */
    virtual Real lcdfc( Real const& t) const;
#endif

  protected:
    /** First degree of freedom */
    int df1_;
    /** Second degree of freedom */
    int df2_;
};

#ifdef IS_RTKPP_LIB

inline Real FisherSnedecor::rand() const
{
#ifdef _OPENMP
//#pragma omp critical
#endif
GetRNGstate(); Real s = Rf_rf(df1_, df2_); PutRNGstate(); return s;
}

inline Real FisherSnedecor::pdf(Real const& x) const { return Rf_df(x, df1_, df2_, false);}
inline Real FisherSnedecor::lpdf(Real const& x) const { return Rf_df(x, df1_, df2_, true);}
inline Real FisherSnedecor::cdf(Real const& t) const { return Rf_pf(t, df1_, df2_, true, false);}
inline Real FisherSnedecor::lcdf(Real const& t) const { return Rf_pf(t, df1_, df2_, true, true);}
inline Real FisherSnedecor::cdfc(Real const& t) const { return Rf_pf(t, df1_, df2_, false, false);}
inline Real FisherSnedecor::lcdfc(Real const& t) const { return Rf_pf(t, df1_, df2_, false, true);}
inline Real FisherSnedecor::icdf(Real const& p) const { return Rf_qf(p, df1_, df2_, true, false);}

inline Real FisherSnedecor::rand( int df1, int df2)
{
#ifdef _OPENMP
//#pragma omp critical
#endif
GetRNGstate(); Real s = Rf_rf(df1, df2); PutRNGstate(); return s;
}

inline Real FisherSnedecor::pdf(Real const& x, int df1, int df2)
{ return Rf_df(x, df1, df2, false);}
inline Real FisherSnedecor::lpdf(Real const& x, int df1, int df2)
{ return Rf_df(x, df1, df2, true);}
inline Real FisherSnedecor::cdf(Real const& t, int df1, int df2)
{  return Rf_pf(t, df1, df2, true, false);}
inline Real FisherSnedecor::icdf(Real const& p, int df1, int df2)
{  return Rf_qf(p, df1, df2, true, false);}

#endif

} // namespace Law

} // namespace STK

#endif /*STK_LAW_FISHERSNEDECOR*/
