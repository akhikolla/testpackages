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
 * Purpose:  Student probability distribution.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Law_Student.h
 *  @brief In this file we define the Student probability distribution.
 **/

#ifndef STK_LAW_STUDENT_H
#define STK_LAW_STUDENT_H

#include "STK_Law_IUnivLaw.h"
#include <Sdk/include/STK_Macros.h>

namespace STK
{

namespace Law
{
/** @ingroup Laws
 *  @brief Student distribution law.
 *
 *  In probability and statistics, <em>Student's t-distribution</em> (or simply
 *  the <em>t-distribution</em>) is any member of a family of continuous
 *  probability distributions that arises when estimating the mean of a normally
 *  distributed population in situations where the sample size is small and
 *  population standard deviation is unknown. Whereas a normal distribution
 *  describes a full population, t-distributions describe samples drawn from a
 *  full population; accordingly, the t-distribution for each sample size is
 *  different, and the larger the sample, the more the distribution resembles a
 *  normal distribution.
 *
 *  The t-distribution plays a role in a number of widely used statistical
 *  analyses, including the Student's t-test for assessing the statistical
 *  significance of the difference between two sample means, the construction of
 *  confidence intervals for the difference between two population means, and in
 *  linear regression analysis. The Student's t-distribution also arises in the
 *  Bayesian analysis of data from a normal family
 *
 *  Student's t-distribution has the probability density function given by
 *  \f[
 *  f(t) = \frac{\Gamma(\frac{\nu+1}{2})}{\sqrt{\nu\pi}\,\Gamma(\frac{\nu}{2})}
 *  \left(1+\frac{t^2}{\nu} \right)^{-\frac{\nu+1}{2}},\!
 *  \f]
 *  where \f$ \nu \f$ is the number of degrees of freedom.
 **/
class Student: public IUnivLaw<Real>
{
  public:
    typedef IUnivLaw<Real> Base;
    /** Default constructor.
     *  @param df degree of freedom parameter
     **/
    inline Student( int df = 1.): Base(_T("Student")), df_(df) {}
    /** destructor */
    inline virtual ~Student() {}
    /** @return the number of degree of freedom */
    inline int df() const { return df_;}
    /** @param df degree of freedom parameter */
    inline void setDf( int df)
    {
      if (df<=0) STKDOMAIN_ERROR_1ARG(Student::setShape,df,shape must be > 0);
      df_ = df;
    }
    /** @return a pseudo Student random variate. */
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

    /** @return a pseudo Student random variate with the specified parameters.
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

inline Real Student::rand() const
{
#ifdef _OPENMP
//#pragma omp critical
#endif
GetRNGstate(); Real s = Rf_rt(df_); PutRNGstate(); return s;
}

inline Real Student::pdf(Real const& x) const  { return Rf_dt(x, df_, false);}
inline Real Student::lpdf(Real const& x) const { return Rf_dt(x, df_, true);}
inline Real Student::cdf(Real const& t) const  { return Rf_pt(t, df_, true, false);}
inline Real Student::lcdf(Real const& t) const  { return Rf_pt(t, df_, true, true);}
inline Real Student::cdfc(Real const& t) const  { return Rf_pt(t, df_, false, false);}
inline Real Student::lcdfc(Real const& t) const  { return Rf_pt(t, df_, false, true);}
inline Real Student::icdf(Real const& p) const { return Rf_qt(p, df_, true, false);}

inline Real Student::rand( int df)
{
#ifdef _OPENMP
//#pragma omp critical
#endif
GetRNGstate(); Real s = Rf_rt(df); PutRNGstate(); return s;
}

inline Real Student::pdf(Real const& x, int df)  { return Rf_dt(x, df, false);}
inline Real Student::lpdf(Real const& x, int df) { return Rf_dt(x, df, true);}
inline Real Student::cdf(Real const& t, int df)  { return Rf_pt(t, df, true, false);}
inline Real Student::icdf(Real const& p, int df) { return Rf_qt(p, df, true, false);}

#endif

} // namespace Law

} // namespace STK

#endif /*STK_LAW_STUDENT_H*/
