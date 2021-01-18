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
 * created on: 8 dec. 2014
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Law_Poisson.h
 *  @brief In this file we define the Poisson distribution.
 **/

#ifndef STK_LAW_POISSON_H
#define STK_LAW_POISSON_H

#include "STK_Law_IUnivLaw.h"
#include <Sdk/include/STK_Macros.h>

namespace STK
{

namespace Law
{
/** @ingroup Laws
 *  @brief Poisson distribution law.
 *
 *  In probability theory and statistics, the <em>Poisson distribution</em>,
 *  named after French mathematician Siméon Denis Poisson, is a discrete
 *  probability distribution that expresses the probability of a given number
 *  of events occurring in a fixed interval of time and/or space if these events
 *  occur with a known average rate and independently of the time since the last
 *  event.
 *
 * The Poisson distribution can be applied to systems with a large number of
 * possible events, each of which is rare. How many such events will occur
 * during a fixed time interval? Under the right circumstances, this is a
 * random number with a Poisson distribution.
 *
 * A discrete random variable @e X is said to have a Poisson distribution with
 * parameter \f$ \lambda >0 \f$ (the mean) if the probability mass function
 * of @e X is given by
 * \f[
 *   f(k; \lambda) = P(X=k) = \frac{\lambda^k e^{-\lambda}}{k!}, \quad k=0,1,2,\ldots,
 * \f]
 * The positive real number \f$ \lambda \f$ is equal to the expected value of
 * @e X and also to its variance.
 */
class Poisson: public IUnivLaw<int>
{
  public:
    typedef IUnivLaw<int> Base;
    /** constructor
     * @param lambda mean of a Poisson distribution
     **/
    inline Poisson( Real const& lambda = 1.)
                 : Base(_T("Poisson")), lambda_(lambda) {}
    /** destructor */
    inline virtual ~Poisson() {}
    /** @return the mean */
    inline Real const& lambda() const { return lambda_;}
    /** @param lambda mean to set */
    inline void setLambda(Real const& lambda )
    {
      if (lambda<0)
      { STKDOMAIN_ERROR_1ARG(Poisson::setLambda,lambda,lambda must be >= 0);}
      lambda_ = lambda;
    }

    /** @return a Poisson random variate . */
    virtual int rand() const;
    /** @brief compute the probability distribution function.
     *  Give the value of the pdf at the point x.
     *  @param x a binary value
     *  @return the value of the pdf
     **/
    virtual Real pdf(int const& x) const;
    /** @brief compute the log probability distribution function.
     *  Give the value of the log-pdf at the point x.
     *  @param x an integer value
     *  @return the value of the log-pdf
     **/
    virtual Real lpdf(int const& x) const;
    /** @brief compute the cumulative distribution function
     *  Give the probability that a Poisson random variate is less or equal
     *  to t.
     *  @param t a real value
     *  @return the value of the cdf
     **/
    virtual Real cdf(Real const& t) const;
    /** @brief inverse cumulative distribution function
     *  The quantile is defined as the smallest value @e q such that
     *  <em> F(q) >= p </em>, where @e F is the cumulative distribution function.
     *  @param p a probability number
     *  @return the quantile for p
     **/
    virtual int icdf(Real const& p) const;

    /** @param lambda the mean
     *  @return a int random variate.
     **/
    static int rand(Real const& lambda);
    /** @brief compute the probability distribution function
     *  Give the value of the pdf at the point x.
     *  @param x a binary value
     *  @param lambda the mean
     *  @return the value of the pdf
     **/
    static Real pdf(int const& x, Real const& lambda);
    /** @brief compute the log probability distribution function
     *  Give the value of the log-pdf at the point x.
     *  @param x a binary value
     *  @param lambda the mean
     *  @return the value of the log-pdf
     **/
    static Real lpdf(int const& x, Real const& lambda);
    /** @brief compute the cumulative distribution function
     *  Give the probability that a Poisson random variate is less or equal
     *  to t.
     *  @param t a real value
     *  @param lambda the mean
     *  @return the value of the cdf at t
     **/
    static Real cdf(Real const& t, Real const& lambda);
    /** @brief inverse cumulative distribution function
     *  The quantile is defined as the smallest value @e x such that
     *  <em> F(x) >= p </em>, where @e F is the cumulative distribution function.
     *  @param p a probability number
     *  @param lambda the mean
     *  @return the quantile for p
     **/
    static int icdf(Real const& p, Real const& lambda);

#ifdef IS_RTKPP_LIB
    /** @return log-cumulative distribution function */
    virtual Real lcdf( Real const& t) const;
    /** @return complement of cumulative distribution function */
    virtual Real cdfc( Real const& t) const;
    /** @return log-complement of cumulative distribution function */
    virtual Real lcdfc( Real const& t) const;
#endif

  protected:
    /** mean of the Poisson distribution */
    Real lambda_;
};

#ifdef IS_RTKPP_LIB

inline int Poisson::rand() const
{
#ifdef _OPENMP
//#pragma omp critical
#endif
  GetRNGstate(); int s = (int)Rf_rpois(lambda_); PutRNGstate(); return s;
}

inline Real Poisson::pdf(int const& x) const  { return Rf_dpois((double)x, lambda_, false);}
inline Real Poisson::lpdf(int const& x) const { return Rf_dpois((double)x, lambda_, true);}
inline Real Poisson::cdf(Real const& t) const { return Rf_ppois(t, lambda_, true, false);}
inline Real Poisson::lcdf(Real const& t) const { return Rf_ppois(t, lambda_, true, true);}
inline Real Poisson::cdfc(Real const& t) const { return Rf_ppois(t, lambda_, false, false);}
inline Real Poisson::lcdfc(Real const& t) const { return Rf_ppois(t, lambda_, false, true);}
inline int Poisson::icdf(Real const& p) const { return (int)::Rf_qpois(p, lambda_, true, false);}


inline int Poisson::rand(Real const& lambda)
{
#ifdef _OPENMP
//#pragma omp critical
#endif
  GetRNGstate(); int s = (int)Rf_rpois(lambda); PutRNGstate(); return s;
}

inline Real Poisson::pdf(int const& x, Real const& lambda)
{ return Rf_dpois((double)x, lambda, false);}
inline Real Poisson::lpdf(int const& x, Real const& lambda)
{ return Rf_dpois((double)x, lambda, true);}
inline Real Poisson::cdf(Real const& t, Real const& lambda)
{ return Rf_ppois(t, lambda, true, false);}
inline int Poisson::icdf(Real const& p, Real const& lambda)
{ return (int)Rf_qpois(p, lambda, true, false);}

#endif

} // namespace Law

} // namespace STK

#endif /* STK_LAW_POISSON_H */
