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

/** @file STK_Law_HyperGeometric.h
 *  @brief In this file we define the HyperGeometric distribution.
 **/

#ifndef STK_LAW_HYPERGEOMETRIC_H
#define STK_LAW_HYPERGEOMETRIC_H

#include "STK_Law_IUnivLaw.h"
#include <Sdk/include/STK_Macros.h>
#include <STKernel/include/STK_Integer.h>

namespace STK
{

namespace Law
{

/** @ingroup Laws
 * @brief HyperGeometric probability law.
 *
 * In probability theory and statistics, the <em>hypergeometric distribution</em>
 * is a discrete probability distribution that describes the probability of @e x
 * successes in @e k draws, without replacement, from a finite population of
 * size <em> m+n </em> that contains exactly @e m successes, wherein each draw
 * is either a success or a failure. In contrast, the binomial distribution
 * describes the probability of @e k successes in@e n draws with replacement.
 *
 * The hypergeometric distribution is used for sampling without replacement.
 * The density of this distribution with parameters <em>m, n</em> and @e k
 * is given by
 * \f[
 *    p(x;m,n,k) = \frac{\binom{m}{x} \binom{n}{k-x}}{\binom{m+n}{k}}
 *    \quad x=0,\ldots,k, \quad m\geq 0,\quad n\geq 0,\quad 0\leq k \leq m+n.
 * \f]
 */
class HyperGeometric: public IUnivLaw<Integer>
{
  public:
    typedef IUnivLaw<Integer> Base;
    /** constructor
     *  @param nbSuccesses, nbFailures, nbDraws number of successes, failures, draws
     **/
    inline HyperGeometric( int nbSuccesses, int nbFailures, int nbDraws)
                        : Base(_T("HyperGeometric"))
                         , nbSuccesses_(nbSuccesses)
                         , nbFailures_(nbFailures)
                         , nbDraws_(nbDraws)
    {}
    /** destructor */
    inline virtual ~HyperGeometric() {}
    /** @return the number of successes */
    inline int nbSuccesses() const { return nbSuccesses_;}
    /** @return the number of failures */
    inline int nbFailures() const { return nbFailures_;}
    /** @return the number of draws */
    inline int nbDraws() const { return nbDraws_;}
    /** @param nbSuccesses the number of successes to set */
    inline void setNbSuccesses(int nbSuccesses )
    {
      if (nbSuccesses<0)
      { STKDOMAIN_ERROR_1ARG(HyperGeometric::setNbSuccesses,nbSuccesses,nbSuccesses must be >= 0);}
      nbSuccesses_ = nbSuccesses;
    }
    /** @param nbFailures the number of failures to set */
    inline void setNbFailures(int nbFailures )
    {
      if (nbFailures<0)
      { STKDOMAIN_ERROR_1ARG(HyperGeometric::setNbFailures,nbFailures,nbFailures must be >= 0);}
      nbFailures_ = nbFailures;
    }
    /** @param nbDraws the number of draws to set */
    inline void setNbDraws(int nbDraws )
    {
      if (nbDraws<0)
      {  STKDOMAIN_ERROR_1ARG(HyperGeometric::setNbDraws,nbDraws,nbDraws must be >= 0);}
      if (nbDraws>nbFailures_+nbSuccesses_)
      { STKDOMAIN_ERROR_1ARG(HyperGeometric::setNbDraws,nbDraws,nbDraws must be <= nbFailures+nbSuccesses);}
      nbDraws_ = nbDraws;
    }
    /** @return a random hypergeometric variate. */
    virtual Integer rand() const;
    /** @brief compute the probability distribution function (density)
     *  Give the value of the pdf at the point x.
     *  @param x an Integer value
     *  @return the value of the pdf
     **/
    virtual Real pdf(Integer const& x) const;
    /** @brief compute the log probability distribution function
     *  Give the value of the log-pdf at the point x.
     *  @param x an Integer value
     *  @return the value of the log-pdf
     **/
    virtual Real lpdf(Integer const& x) const;
    /** @brief compute the cumulative distribution function
     *  Give the probability that a HyperGeometric random variate is less or equal
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

    /** @brief random hypergeometric variate generation.
     *  @param nbSuccesses, nbFailures, nbDraws number of successes, failures, draws
     *  @return a Integer random variate.
     **/
    static Integer rand( int nbSuccesses, int nbFailures, int nbDraws);
    /** @brief compute the probability distribution function.
     *  Give the value of the pdf at the point x.
     *  @param x an Integer value
     *  @param nbSuccesses, nbFailures, nbDraws number of successes, failures, draws
     *  @return the value of the pdf
     **/
    static Real pdf(Integer x, int nbSuccesses, int nbFailures, int nbDraws);
    /** @brief compute the log probability distribution function.
     *  @param x an Integer value
     *  @param nbSuccesses, nbFailures, nbDraws number of successes, failures, draws
     *  @return the value of the log-pdf
     **/
    static Real lpdf(Integer x, int nbSuccesses, int nbFailures, int nbDraws);
    /** @brief compute the cumulative distribution function
     *  Give the probability that a HyperGeometric random variate is less or equal
     *  to t.
     *  @param t a real value
     *  @param nbSuccesses, nbFailures, nbDraws number of successes, failures, draws
     *  @return the value of the cdf
     **/
    static Real cdf(Real const& t, int nbSuccesses, int nbFailures, int nbDraws);
    /** @brief inverse cumulative distribution function
     *  The quantile is defined as the smallest value @e x such that
     *  <em> F(x) >= p </em>, where @e F is the cumulative distribution function.
     *  @param p a probability number
     *  @param nbSuccesses, nbFailures, nbDraws number of successes, failures, draws
     **/
    static Integer icdf(Real const& p, int nbSuccesses, int nbFailures, int nbDraws);

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
    int nbSuccesses_;
    /** number of failures */
    int nbFailures_;
    /** number of draws */
    int nbDraws_;
};

#ifdef IS_RTKPP_LIB

inline int HyperGeometric::rand() const
{
#ifdef _OPENMP
//#pragma omp critical
#endif
GetRNGstate(); int s = Rf_rhyper(nbSuccesses_, nbFailures_, nbDraws_); PutRNGstate(); return s;
}

inline Real HyperGeometric::pdf(Integer const& x) const
{ return Rf_dhyper((double)x, nbSuccesses_, nbFailures_, nbDraws_, false);}
inline Real HyperGeometric::lpdf(Integer const& x) const
{ return Rf_dhyper((double)x, nbSuccesses_, nbFailures_, nbDraws_, true);}
inline Real HyperGeometric::cdf(Real const& t) const
{ return Rf_phyper(t, nbSuccesses_, nbFailures_, nbDraws_, true, false);}
inline Real HyperGeometric::lcdf(Real const& t) const
{ return Rf_phyper(t, nbSuccesses_, nbFailures_, nbDraws_, true, true);}
inline Real HyperGeometric::cdfc(Real const& t) const
{ return Rf_phyper(t, nbSuccesses_, nbFailures_, nbDraws_, false, false);}
inline Real HyperGeometric::lcdfc(Real const& t) const
{ return Rf_phyper(t, nbSuccesses_, nbFailures_, nbDraws_, false, true);}
inline int HyperGeometric::icdf(Real const& p) const
{ return (int)::Rf_qhyper(p, nbSuccesses_, nbFailures_, nbDraws_, true, false);}

inline Integer HyperGeometric::rand( int nbSuccesses, int nbFailures, int nbDraws)
{
#ifdef _OPENMP
//#pragma omp critical
#endif
  GetRNGstate(); int s = Rf_rhyper(nbSuccesses, nbFailures, nbDraws); PutRNGstate(); return s;
}

inline Real HyperGeometric::pdf(Integer x, int nbSuccesses, int nbFailures, int nbDraws)
{ return Rf_dhyper((double)x, nbSuccesses, nbFailures, nbDraws, false);}
inline Real HyperGeometric::lpdf(Integer x, int nbSuccesses, int nbFailures, int nbDraws)
{ return Rf_dhyper((double)x, nbSuccesses, nbFailures, nbDraws, true);}
inline Real HyperGeometric::cdf(Real const& t, int nbSuccesses, int nbFailures, int nbDraws)
{ return Rf_phyper(t, nbSuccesses, nbFailures, nbDraws, true, false);}
inline int HyperGeometric::icdf(Real const& p, int nbSuccesses, int nbFailures, int nbDraws)
{ return (int)::Rf_qhyper(p, nbSuccesses, nbFailures, nbDraws, true, false);}

#endif

} // namespace Law

} // namespace STK

#endif /* STK_LAW_HYPERGEOMETRIC_H */
