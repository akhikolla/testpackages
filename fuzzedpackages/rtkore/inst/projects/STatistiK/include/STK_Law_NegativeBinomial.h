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

/** @file STK_Law_NegativeBinomial.h
 *  @brief In this file we define the NegativeBinomial distribution.
 **/

#ifndef STK_LAW_NEGATIVEBINOMIAL_H
#define STK_LAW_NEGATIVEBINOMIAL_H

#include "STK_Law_IUnivLaw.h"
#include <Sdk/include/STK_Macros.h>
#include <STKernel/include/STK_Integer.h>

namespace STK
{

namespace Law
{

/** @ingroup Laws
 *  @brief NegativeBinomial probability law.
 *
 * In probability theory and statistics, the <em>negative binomial distribution</em>
 * is a discrete probability distribution of the number of successes in a
 * sequence of independent and identically distributed Bernoulli trials before
 * a specified (non-random) number of failures (denoted @e r) occurs.
 *
 * Suppose there is a sequence of independent Bernoulli trials, each trial having
 * two potential outcomes called “success” and “failure”. In each trial the
 * probability of success is @e p and of failure is <em>(1 − p)</em>. We are
 * observing this sequence until a predefined number @e r of failures has
 * occurred. Then the random number of successes we have seen, @e X, will have
 * the negative binomial (or Pascal) distribution:
 * \f$ X\ \sim\ \text{NB}(r; p) \f$
 *
 * The probability mass function of the negative binomial distribution is
 * \f[
 * f(k; r, p) = \Pr(X = k) = \binom{k+r-1}{k} p^k(1-p)^r \quad\text{for }
 * k = 0, 1, 2, \dots
 * \f]
 **/
class NegativeBinomial: public IUnivLaw<Integer>
{
  public:
    typedef IUnivLaw<Integer> Base;
    /** constructor
     *  @param prob, size probability of success and number of successes in a NegativeBinomial trial
     **/
    inline NegativeBinomial( int size =1, Real const& prob =0.5)
                          : Base(_T("Negative Binomial")), size_(size), prob_(prob) {}

    /** destructor */
    inline virtual ~NegativeBinomial() {}
    /** @return the probability of success */
    inline Real const& prob() const { return prob_;}
    /** @return the number of trials */
    inline int size() const { return size_;}
    /** @param prob the probability of success to set */
    inline void setProb(Real const& prob )
    {
      if (prob<0) STKDOMAIN_ERROR_1ARG(NegativeBinomial::setProb,prob,prob must be >= 0);
      if (prob>1) STKDOMAIN_ERROR_1ARG(NegativeBinomial::setProb,prob,prob must be <= 1);
      prob_ = prob;
    }
    /** @param size the number of successes to set */
    inline void setSize(Integer size)
    {
      if (size<0) STKDOMAIN_ERROR_1ARG(NegativeBinomial::setSize,size,size must be >= 0);
      size_ = size;
    }

    /** @return a Integer random variate . */
    virtual Integer rand() const;
    /** @brief compute the cumulative distribution function
     *  Give the probability that a NegativeBinomial random variate is less or equal
     *  to t.
     *  @param t a real value
     *  @return the value of the cdf
     **/
    virtual Real cdf(Real const& t) const;
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
    /** @brief inverse cumulative distribution function
     *  The quantile is defined as the smallest value @e x such that
     *  <em> F(x) >= p </em>, where @e F is the cumulative distribution function.
     *  @param p a probability number
     **/
    virtual Integer icdf(Real const& p) const;

    /** @param prob, size probability of success and number of successes in a NegativeBinomial trial
     *  @return a Integer random variate.
     **/
    static Integer rand( int size, Real const& prob);
    /** @brief compute the probability distribution function (density)
     *  Give the value of the pdf at the point x.
     *  @param x an integer value
     *  @param prob, size probability of success and number of successes in a NegativeBinomial trial
     *  @return the value of the pdf
     **/
    static Real pdf(Integer x,  int size, Real const& prob);
    /** @brief compute the log probability distribution function
     *  Give the value of the log-pdf at the point x.
     *  @param x an integer value
     *  @param prob, size probability of success and number of successes in a NegativeBinomial trial
     *  @return the value of the log-pdf
     **/
    static Real lpdf(Integer x,  int size, Real const& prob);
    /** @brief compute the cumulative distribution function
     *  Give the value of the cdf at the point t.
     *  @param t a real value
     *  @param size, prob probability of success and number of successes in a NegativeBinomial trial
     *  @return the value of the cdf
     **/
    static Real cdf(Real const& t, int size, Real const& prob);
    /** @brief inverse cumulative distribution function
     *  The quantile is defined as the smallest value @e x such that
     *  <em> F(x) >= p </em>, where @e F is the cumulative distribution function.
     *  @param p a probability number
     *  @param prob, size probability of success and number of successes in a NegativeBinomial trial
     **/
    static Integer icdf(Real const& p,  int size, Real const& prob);

#ifdef IS_RTKPP_LIB
    /** @return log-cumulative distribution function */
    virtual Real lcdf( Real const& t) const;
    /** @return complement of cumulative distribution function */
    virtual Real cdfc( Real const& t) const;
    /** @return log-complement of cumulative distribution function */
    virtual Real lcdfc( Real const& t) const;
#endif

  protected:
    /** number of successes */
    int size_;
    /** probability of success in a Bernoulli trial */
    Real prob_;
};

#ifdef IS_RTKPP_LIB

inline int NegativeBinomial::rand() const
{
#ifdef _OPENMP
//#pragma omp critical
#endif
GetRNGstate(); int s = Rf_rnbinom(size_, prob_); PutRNGstate(); return s;
}

inline Real NegativeBinomial::pdf(Integer const& x) const { return Rf_dnbinom((double)x, size_, prob_, false);}
inline Real NegativeBinomial::lpdf(Integer const& x) const { return Rf_dnbinom((double)x, size_, prob_, true);}
inline Real NegativeBinomial::cdf(Real const& t) const { return Rf_pnbinom(t, size_, prob_, true, false);}
inline Real NegativeBinomial::lcdf(Real const& t) const { return Rf_pnbinom(t, size_, prob_, true, true);}
inline Real NegativeBinomial::cdfc(Real const& t) const { return Rf_pnbinom(t, size_, prob_, false, false);}
inline Real NegativeBinomial::lcdfc(Real const& t) const { return Rf_pnbinom(t, size_, prob_, false, true);}
inline int NegativeBinomial::icdf(Real const& p) const { return (int)Rf_qnbinom(p, size_, prob_, true, false);}

inline int NegativeBinomial::rand(int size, Real const& prob)
{
#ifdef _OPENMP
//#pragma omp critical
#endif
GetRNGstate(); int s = Rf_rnbinom(size, prob); PutRNGstate(); return s;
}

inline Real NegativeBinomial::pdf(Integer x, int size, Real const& prob)
{ return Rf_dnbinom((double)x, size, prob, false);}
inline Real NegativeBinomial::lpdf(Integer x, int size, Real const& prob)
{ return Rf_dnbinom((double)x, size, prob, true);}
inline Real NegativeBinomial::cdf(Real const& t, int size, Real const& prob)
{ return Rf_pnbinom(t, size, prob , true, false);}
inline int NegativeBinomial::icdf(Real const& p, int size, Real const& prob)
{ return (int)::Rf_qnbinom(p, size, prob , true, false);}

#endif

} // namespace Law

} // namespace STK

#endif /* STK_LAW_NEGATIVEBINOMIAL_H */
