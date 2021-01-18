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

/** @file STK_Law_Geometric.h
 *  @brief In this file we define the Geometric distribution.
 **/

#ifndef STK_LAW_GEOMETRIC_H
#define STK_LAW_GEOMETRIC_H

#include "STK_Law_IUnivLaw.h"
#include <Sdk/include/STK_Macros.h>
#include <STKernel/include/STK_Integer.h>

namespace STK
{

namespace Law
{

/** @ingroup Laws
 * @brief Geometric probability law.
 *
 * In probability theory and statistics, the <em>Geometric distribution</em>
 * with parameter @e p is the discrete probability distribution
 * of the number of failures of Bernoulli trials before the first success.
 *
 * In general, if the random variable @e X follows the geometric distribution
 * with parameters @e p, we write \f$ X \sim \mathcal{G}(p)\f$.
 * The probability of getting exactly @e k failures befor a success in Bernoulli
 * trials is given by the probability mass function:
 * \f[ f(k;p) = \mathbb{P}(X = k) = p(1-p)^{k} \mbox{ for } k=0,1,\ldots \f]
 */
class Geometric: public IUnivLaw<Integer>
{
  public:
    typedef IUnivLaw<Integer> Base;
    /** constructor
     * @param prob probability of success in a Bernoulli trial
     **/
    inline Geometric(Real const& prob =0.5): Base(_T("Geometric")), prob_(prob)
    {}
    /** destructor */
    inline virtual ~Geometric(){}
    /** @return the probability of success */
    inline Real const& prob() const { return prob_;}
    /** @param prob the probability of success to set */
    inline void setProb(Real const& prob )
    {
      if (prob<0) STKDOMAIN_ERROR_1ARG(Geometric::setProb,prob,prob must be >= 0);
      if (prob>1) STKDOMAIN_ERROR_1ARG(Geometric::setProb,prob,prob must be <= 1);
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
     *  Give the probability that a Geometric random variate is less or equal
     *  to t.
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
     *  @param prob the probability of success in a Bernoulli trial
     *  @return a Integer random variate.
     **/
    static Integer rand(Real const& prob);
    /** @brief compute the probability distribution function (density)
     *  Give the value of the pdf at the point x.
     *  @param x a binary value
     *  @param prob the probability of success in a Bernoulli trial
     *  @return the value of the pdf
     **/
    static Real pdf(Integer x, Real const& prob);
    /** @brief compute the log probability distribution function
     *  Give the value of the log-pdf at the point x.
     *  @param x a binary value
     *  @param prob the probability of success in a Bernoulli trial
     *  @return the value of the log-pdf
     **/
    static Real lpdf(Integer x, Real const& prob);
    /** @brief compute the cumulative distribution function
     *  Give the probability that a Geometric random variate is less or equal
     *  to t.
     *  @param t a real value
     *  @param prob the probability of success in a Bernoulli trial
     *  @return the value of the cdf
     **/
    static Real cdf(Real const& t, Real const& prob);
    /** @brief inverse cumulative distribution function
     *  The quantile is defined as the smallest value @e x such that
     *  <em> F(x) >= p </em>, where @e F is the cumulative distribution function.
     *  @param p a probability number
     *  @param prob the probability of success in a Bernoulli trial
     **/
    static Integer icdf(Real const& p, Real const& prob);

#ifdef IS_RTKPP_LIB
    /** @return log-cumulative distribution function */
    virtual Real lcdf( Real const& t) const;
    /** @return complement of cumulative distribution function */
    virtual Real cdfc( Real const& t) const;
    /** @return log-complement of cumulative distribution function */
    virtual Real lcdfc( Real const& t) const;
#endif

  protected:
    /** probability of success in a Bernoulli trial */
    Real prob_;
};

#ifdef IS_RTKPP_LIB


inline int Geometric::rand() const
{
#ifdef _OPENMP
//#pragma omp critical
#endif
GetRNGstate(); int s = Rf_rgeom(prob_); PutRNGstate(); return s;
}

inline Real Geometric::pdf(Integer const& x) const { return Rf_dgeom((double)x, prob_, false);}
inline Real Geometric::lpdf(Integer const& x) const { return Rf_dgeom((double)x, prob_, true);}
inline Real Geometric::cdf(Real const& t) const { return Rf_pgeom(t, prob_, true, false);}
inline Real Geometric::lcdf(Real const& t) const { return Rf_pgeom(t, prob_, true, true);}
inline Real Geometric::cdfc(Real const& t) const { return Rf_pgeom(t, prob_, false, false);}
inline Real Geometric::lcdfc(Real const& t) const { return Rf_pgeom(t, prob_, false, true);}
inline int Geometric::icdf(Real const& p) const { return Rf_qgeom(p, prob_, true, false);}

inline int Geometric::rand(Real const& prob)
{
#ifdef _OPENMP
//#pragma omp critical
#endif
GetRNGstate(); int s = Rf_rgeom(prob); PutRNGstate(); return s;
}

inline Real Geometric::pdf(Integer x, Real const& prob)
{ return Rf_dgeom((double)x, prob, false);}
inline Real Geometric::lpdf(Integer x, Real const& prob)
{ return Rf_dgeom((double)x, prob, true);}
inline Real Geometric::cdf(Real const& t, Real const& prob)
{ return Rf_pgeom(t, prob, true, false);}
inline Integer Geometric::icdf(Real const& p, Real const& prob)
{ return (Integer)::Rf_qgeom(p, prob, true, false);}

#endif
} // namespace Law

} // namespace STK

#endif /* STK_LAW_GEOMETRIC_H */
