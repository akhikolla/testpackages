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
 * created on: 23 janv. 2013
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Law_Binomial.h
 *  @brief In this file we define the Binomial distribution.
 **/

#ifndef STK_LAW_BINOMIAL_H
#define STK_LAW_BINOMIAL_H

#include "STK_Law_IUnivLaw.h"
#include <Sdk/include/STK_Macros.h>
#include <STKernel/include/STK_Integer.h>


namespace STK
{
namespace Law
{
/** @ingroup Laws
 *  @brief Binomial probability law.
 *
 * In probability theory and statistics, the <em>Binomial distribution</em>
 * with parameters @e n and @e p is the discrete probability distribution
 * of the number of successes in a sequence of @e n independent
 * yes/no experiments, each of which yields success with probability @e p.
 * A success/failure experiment is also called a Bernoulli experiment or
 * Bernoulli trial. When @e n = 1, the binomial distribution is a
 * <em>Bernoulli distribution</em>
 *
 * In general, if the random variable @e X follows the binomial distribution
 * with parameters @e n and @e p, we write \f$ X \sim B(n;p)\f$.
 * The probability of getting exactly @e k successes in @e n trials is given by
 * the probability mass function:
 * \f[ f(k;n,p) = \mathbb{P}(X = k) = {n\choose k}p^k(1-p)^{n-k} \f]
 *
 *  The binomial distribution is the basis for the popular binomial test of
 *  statistical significance.
 */
class Binomial: public IUnivLaw<Integer>
{
  public:
    typedef IUnivLaw<Integer> Base;
    /** constructor
     *  @param prob probability of success in a Binomial trial
     *  @param n the number of trials
     **/
    inline Binomial( int n =1, Real const& prob =0.5)
                  : Base(_T("Binomial")), n_(n), prob_(prob)
    {
      if (prob<0) STKDOMAIN_ERROR_2ARG(Binomial::Binomial,prob,n,prob must be >= 0);
      if (prob>1) STKDOMAIN_ERROR_2ARG(Binomial::Binomial,prob,n,prob must be <= 1);
      if (n<0) STKDOMAIN_ERROR_2ARG(Binomial::Binomial,prob,n,n must be >= 0);
    }
    /** destructor */
    inline virtual ~Binomial() {}
    /** @return the number of trials */
    inline int n() const { return n_;}
    /** @return the probability of success */
    inline Real const& prob() const { return prob_;}
    /** @param n the number of trials to set */
    inline void setN(Integer n)
    {
      if (n<0) STKDOMAIN_ERROR_1ARG(Binomial::setN,n,n must be >= 0);
      n_ = n;
    }
    /** @param prob the probability of success to set */
    inline void setProb(Real const& prob )
    {
      if (prob<0) STKDOMAIN_ERROR_1ARG(Binomial::setProb,prob,prob must be >= 0);
      if (prob>1) STKDOMAIN_ERROR_1ARG(Binomial::setProb,prob,prob must be <= 1);
      prob_ = prob;
    }

    /** @return a Integer random variate . */
    virtual Integer rand() const;
    /** @brief compute the probability distribution function (density)
     *  Give the value of the pdf at the point x.
     *  @param x a binary value
     *  @return the value of the pdf
     **/
    virtual Real pdf(Integer const& x) const;
    /** @brief compute the log probability distribution function
     *  Give the value of the log-pdf at the point x.
     *  @param x a binary value
     *  @return the value of the log-pdf
     **/
    virtual Real lpdf(Integer const& x) const;
    /** @brief compute the cumulative distribution function
     *  Give the probability that a Binomial random variate is less or equal to t.
     *  @param t a real value
     *  @return the value of the cdf
     **/
    virtual Real cdf(Real const& t) const;
    /** @brief inverse cumulative distribution function
     *  The quantile is defined as the smallest value @e x such that
     *  <em> F(x) >= p </em>, where @e F is the cumulative distribution function.
     *  @param p a probability number
     **/
    virtual Integer icdf(Real const& p) const;

    /**
     *  @param n,prob number of trial and probability of success
     *  @return a Integer random variate.
     **/
    static Integer rand(int n, Real const& prob);
    /** @brief compute the probability distribution function (density)
     *  Give the value of the pdf at the point x.
     *  @param x a binary value
     *  @param n,prob number of trial and probability of success
     *  @return the value of the pdf
     **/
    static Real pdf(Integer x, int n, Real const& prob);
    /** @brief compute the log probability distribution function
     *  Give the value of the log-pdf at the point x.
     *  @param x a binary value
     *  @param n,prob number of trial and probability of success
     *  @return the value of the log-pdf
     **/
    static Real lpdf(Integer x, int n, Real const& prob);
    /** @brief compute the cumulative distribution function
     *  Give the probability that a Binomial random variate is less or equal
     *  to t.
     *  @param t a real value
     *  @param n,prob number of trial and probability of success
     *  @return the value of the cdf
     **/
    static Real cdf(Real const& t, int n, Real const& prob);
    /** @brief inverse cumulative distribution function
     *  The quantile is defined as the smallest value @e x such that
     *  <em> F(x) >= p </em>, where @e F is the cumulative distribution function.
     *  @param p a probability number
     *  @param n,prob number of trial and probability of success
     **/
    static Integer icdf(Real const& p, int n, Real const& prob);

#ifdef IS_RTKPP_LIB
    /** @return log-cumulative distribution function */
    virtual Real lcdf( Real const& t) const;
    /** @return complement of cumulative distribution function */
    virtual Real cdfc( Real const& t) const;
    /** @return log-complement of cumulative distribution function */
    virtual Real lcdfc( Real const& t) const;
#endif

  protected:
    /** number of trials */
    int n_;
    /** probability of success in a Bernoulli trial */
    Real prob_;
};

#ifdef IS_RTKPP_LIB

inline int Binomial::rand() const
{
#ifdef _OPENMP
//#pragma omp critical
#endif
GetRNGstate(); int s = Rf_rbinom(n_, prob_); PutRNGstate(); return s;
}

inline Real Binomial::pdf(int const& x) const { return (Real)Rf_dbinom((double)x, (double)n_, prob_, false);}
inline Real Binomial::lpdf(int const& x) const { return (Real)Rf_dbinom((double)x, (double)n_, prob_, true);}
inline Real Binomial::cdf(Real const& t) const { return (Real)Rf_pbinom((double)t, (double)n_, prob_, true, false);}
inline Real Binomial::lcdf(Real const& t) const { return (Real)Rf_pbinom((double)t, (double)n_, prob_, true, true);}
inline Real Binomial::cdfc(Real const& t) const { return (Real)Rf_pbinom(t, (double)n_, prob_, false, false);}
inline Real Binomial::lcdfc(Real const& t) const { return (Real)Rf_pbinom(t, (double)n_, prob_, false, true);}
inline int Binomial::icdf(Real const& p) const { return Rf_qbinom(p, (double)n_, prob_, true, false);}

inline int Binomial::rand(int n, Real const& prob)
{
#ifdef _OPENMP
//#pragma omp critical
#endif
GetRNGstate(); int s = Rf_rbinom(n, prob); PutRNGstate(); return s;
}

inline Real Binomial::pdf(int x, int n, Real const& prob)
{ return (Real)Rf_dbinom(x, (double)n, prob, false);}
inline Real Binomial::lpdf(int x, int n, Real const& prob)
{ return (Real)Rf_dbinom((double)x, (double)n, prob, true);}
inline Real Binomial::cdf(Real const& t, int n, Real const& prob)
{ return (Real)Rf_pbinom(t, (double)n, prob, true, false);}
inline int Binomial::icdf(Real const& p, int n, Real const& prob)
{ return Rf_qbinom(p, (double)n, prob, true, false);}

#endif

} // namespace Law

} // namespace STK

#endif /* STK_LAW_BINOMIAL_H */
