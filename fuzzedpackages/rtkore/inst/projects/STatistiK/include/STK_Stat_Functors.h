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
 * Project:  stkpp::STatistiK::StatDesc
 * Purpose:  Compute elementary 1D statistics for all variables.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Stat_Functors.h
 *  @brief This file contain the functors computings statistics.
 **/

#ifndef STK_STAT_FUNCTORS_H
#define STK_STAT_FUNCTORS_H

#include <Arrays/include/STK_ExprBaseFunctors.h>

#define FUNCTION_ON_ARRAY(funct, Op) \
template< class Derived> \
typename hidden::FunctorTraits<Derived, Op >::Row \
funct(Derived const& A) \
{ return typename hidden::FunctorTraits<Derived, Op>::ColOp(A)();} \
template< class Derived, class Weights> \
typename hidden::FunctorTraits<Derived, Op >::Row \
funct(Derived const& A, Weights const& w) \
{ return typename hidden::FunctorTraits<Derived, Op>::ColWeightedOp(A)(w);} \
template< class Derived> \
typename hidden::FunctorTraits<Derived, Op >::Row \
funct##ByCol(Derived const& A) \
{ return typename hidden::FunctorTraits<Derived, Op>::ColOp(A)();} \
template< class Derived, class Weights> \
typename hidden::FunctorTraits<Derived, Op >::Row \
funct##ByCol(Derived const& A, Weights const& w) \
{ return typename hidden::FunctorTraits<Derived, Op>::ColWeightedOp(A)(w);} \
template< class Derived> \
typename hidden::FunctorTraits<Derived, Op >::Col \
funct##ByRow(Derived const& A) \
{ return typename hidden::FunctorTraits<Derived, Op>::RowOp(A)();} \
template< class Derived, class Weights> \
typename hidden::FunctorTraits<Derived, Op >::Col \
funct##ByRow(Derived const& A, Weights const& w) \
{ return typename hidden::FunctorTraits<Derived, Op>::RowWeightedOp(A)(w);}


namespace STK
{


namespace Stat
{
/** @ingroup Stat
 *  Compute the minimal value of the variable V
 **/
template<class Derived>
struct MinOp
{
    typedef typename Derived::Type Type;
    /** constructor */
    inline MinOp( ExprBase<Derived> const&  V):  V_(V.asDerived())
    { STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Derived);}
    /**  @return the minimal value of the variable V
     *  \f[ \min_{i=1}^n v_i \f]
     **/
    Type const operator()() const
    {
      Type min  = Arithmetic<Type>::max();
      for (int i=V_.begin(); i<V_.end(); i++)
      { min = std::min(min, V_[i]);}
      return min;
    }
    /** @return the minimal value of the variable V
     *  \f[ \min_{i=1}^n w_i v_i \f]
     *  @param w the weights
     **/
    template< class Weights>
    Type const operator()( ExprBase<Weights> const&  w) const
    {
      STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Weights);
      Type min  = Arithmetic<Type>::max();
      for (int i=V_.begin(); i<V_.end(); i++)
      { min = std::min(min, w[i]*V_[i]);}
      return min;
    }
  protected:
    Derived const& V_;
};

/** @ingroup Stat
 *  Compute safely the minimal value of the variable V
 **/
template<class Derived>
struct MinSafeOp
{
    typedef typename Derived::Type Type;
    /** constructor */
    inline MinSafeOp( ExprBase<Derived> const&  V) :  V_(V.asDerived())
    { STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Derived);}
    /**  @return the safely computed minimal value of the variable V */
    Type const operator()() const
    {
      Type min  = Arithmetic<Type>::max();
      for (int i=V_.begin(); i<V_.end(); i++)
      { if (Arithmetic<Type>::isFinite(V_[i])) min = std::min(min, V_[i]);}
      return min;
    }
    /** @return the safely computed weighted minimal value of the variable V
     *  @param w the weights
     **/
    template< class Weights>
    Type const operator()( ExprBase<Weights> const&  w) const
    {
      STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Weights);
      Type min  = Arithmetic<Type>::max();
      for (int i=V_.begin(); i<V_.end(); i++)
      { if (!Arithmetic<Type>::isNA(V_[i])) min = std::min(min, w[i]*V_[i]);}
      return min;
    }
  protected:
    Derived const& V_;
};

/** @ingroup Stat
 *  Compute the maximal value of the variable V
 **/
template<class Derived>
struct MaxOp
{
    typedef typename Derived::Type Type;
    /** constructor */
    inline MaxOp( ExprBase<Derived> const&  V) :  V_(V.asDerived())
    { STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Derived);}
    /** @return the maximal value of the variable V
     *  \f[ \max_{i=1}^n v_i \f]
     **/
    Type const operator()() const
    {
      Type max  = -Arithmetic<Type>::max();
      for (int i=V_.begin(); i<V_.end(); i++)
      { max = std::max(max, V_[i]);}
      return max;
    }
    /** @return the weighted maximal value of the variable V
     *  \f[ \max_{i=1}^n w_i v_i \f]
     *  @param w the weights
     **/
    template< class Weights>
    Type const operator()( ExprBase<Weights> const&  w) const
    {
      STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Weights);
      Type max  = -Arithmetic<Type>::max();
      for (int i=V_.begin(); i<V_.end(); i++)
      { max = std::max(max, w[i]*V_[i]);}
      return max;
    }
  protected:
    Derived const& V_;
};

/** @ingroup Stat
 *  Compute safely the maximal value of the variable V
 **/
template<class Derived>
struct MaxSafeOp
{
    typedef typename Derived::Type Type;
    /** constructor */
    inline MaxSafeOp( ExprBase<Derived> const&  V) :  V_(V.asDerived())
    { STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Derived);}
    /** @return the safely computed maximal value of the variable V */
    Type const operator()() const
    {
      Type max  = -Arithmetic<Type>::max();
      for (int i=V_.begin(); i<V_.end(); i++)
      { if (!Arithmetic<Type>::isNA(V_[i])) max = std::max(max, V_[i]);}
      return max;
    }
    /** @return the safely computed weighted maximal value of the variable V
     *  @param w the weights
     **/
    template< class Weights>
    Type const operator()( ExprBase<Weights> const&  w) const
    {
      STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Weights);
      Type max  = -Arithmetic<Type>::max();
      for (int i=V_.begin(); i<V_.end(); i++)
      { if (!Arithmetic<Type>::isNA(V_[i])) max = std::max(max, w[i]*V_[i]);}
      return max;
    }
  protected:
    Derived const& V_;
};

/** @ingroup Stat
 *  Compute the sum of the variable V.
 **/
template<class Derived>
struct SumOp
{
    typedef typename Derived::Type Type;
    /** constructor */
    inline SumOp( ExprBase<Derived> const&  V) :  V_(V.asDerived())
    { STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Derived);}
    /** @return the sum of the variable V
     *  \f[ \sum_{i=1}^n V(i) \f]
     **/
    Type const operator()() const
    {
      Type sum  = 0.;
      for (int i=V_.begin(); i<V_.end(); i++) { sum += V_[i];}
      return sum;
    }
    /** @return the weighted sum of the variable V
     *  \f[ \hat{\mu} = \sum_{i=1}^n w(i) V(i). \f]
     *  @param w the weights
     **/
    template< class Weights>
    Type const operator()( ExprBase<Weights> const&  w) const
    {
      STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Weights);
      // sum the weighted samples
      Type sum  = 0.0;
      for (int i=V_.begin(); i<V_.end(); i++) { sum += w[i] * V_[i];}
      return sum;
    }
  protected:
    Derived const& V_;
};

/** @ingroup Stat
 *  Compute safely the sum of the variable V.
 **/
template<class Derived>
struct SumSafeOp
{
    typedef typename Derived::Type Type;
    /** constructor */
    inline SumSafeOp( ExprBase<Derived> const&  V) :  V_(V.asDerived())
    { STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Derived);}
    /** @return the safely computed sum of the variable V discarding all missing values. */
    Type const operator()() const
    {
      // sum the samples
      Type sum  = 0.0;
      for (int i=V_.begin(); i<V_.end(); i++)
      { if (!Arithmetic<Type>::isNA(V_[i])) sum += V_[i];}
      return sum;
    }
    /** @return the safely computed weighted sum of the variable V
     *  @param w the weights
     **/
    template< class Weights>
    inline Type const operator()( ExprBase<Weights> const&  w) const
    {
      STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Weights);
#ifdef STK_BOUNDS_CHECK
      if (V_.range() != w.range())
        STKRUNTIME_ERROR_NO_ARG(wmeanSafe,V.range()!=w.range());
#endif
      // sum the weighted samples
      Type sum  = 0.0;
      for (int i=V_.begin(); i<V_.end(); i++)
      { if ( (!Arithmetic<Type>::isNA(V_[i])) && (!Arithmetic<Type>::isNA(w[i])))
        { sum += std::abs(w[i]) * V_[i];}
      }
      return sum;
    }
  protected:
    Derived const& V_;
};

/** @ingroup Stat
 *  Compute the mean of the variable V.
 **/
template<class Derived>
struct MeanOp
{
    typedef typename Derived::Type Type;
    /** constructor */
    inline MeanOp( ExprBase<Derived> const&  V) :  V_(V.asDerived())
    { STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Derived);}
    /** @return the mean of the variable V
     *  \f[ \hat{\mu} = \frac{1}{n} \sum_{i=1}^n V(i) \f]
     **/
    Type const operator()() const
    {
      // no samples
      if (V_.empty()) { return Type(0);}
      // sum the samples
      Type sum  = 0.;
      for (int i=V_.begin(); i<V_.end(); i++) { sum += V_[i];}
      return sum /= (Type)(V_.size());
    }
    /** @return the weighted mean value of the variable V
     *  \f[ \hat{\mu} = \frac{1}{\sum_{i=1}^n w(i)} \sum_{i=1}^n w(i) V(i). \f]
     *  @param w the weights
     **/
    template< class Weights>
    Type const operator()( ExprBase<Weights> const&  w) const
    {
      STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Weights);
      // no samples
      if (V_.empty()) { return Type(0);}
      // sum the weighted samples
      Type sum  = 0.0, sumweights= 0.0;
      for (int i=V_.begin(); i<V_.end(); i++)
      {
        Type const weight = std::abs(w[i]);
        sumweights += weight;
        sum        += weight * V_[i];
      }
      // compute the weighted mean. If all weights are 0, we get 0
      return (sumweights) ? sum /= sumweights: Type(0);
    }
  protected:
    Derived const& V_;
};

/** @ingroup Stat
 *  Compute safely mean of the variable V.
 **/
template<class Derived>
struct MeanSafeOp
{
    typedef typename Derived::Type Type;
    /** constructor */
    inline MeanSafeOp( ExprBase<Derived> const&  V) :  V_(V.asDerived())
    { STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Derived);}
    /** @return the safely computed mean of the variable V discarding all missing values. */
    Type const operator()() const
    {
      // no samples
      if (V_.empty()) { return Type(0);}
      // get dimensions
      int nobs = V_.size();
      // sum the samples
      Type sum  = 0.0;
      for (int i=V_.begin(); i<V_.end(); i++)
      { (!Arithmetic<Type>::isNA(V_[i])) ? sum += V_[i] : nobs--;}
      // compute the mean
      return nobs ? (sum /= Type(nobs)) : Type(0);
    }
    /** @return the safely computed weighted mean of the variable V
     *  @param w the weights
     **/
    template< class Weights>
    inline Type const operator()( ExprBase<Weights> const&  w) const
    {
      STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Weights);
      // no samples
      if (V_.empty()) { return Type(0);}
#ifdef STK_BOUNDS_CHECK
      if (V_.range() != w.range())
        STKRUNTIME_ERROR_NO_ARG(wmeanSafe,V.range()!=w.range());
#endif
      // sum the weighted samples
      Type sum  = 0.0, sumweights= 0.0;
      for (int i=V_.begin(); i<V_.end(); i++)
      { if ( (!Arithmetic<Type>::isNA(V_[i])) && (!Arithmetic<Type>::isNA(w[i])))
        {
          Type weight  = std::abs(w[i]);
          sumweights += weight;
          sum        += weight * V_[i];
        }
      }
      // compute the weighted mean. If all weights are 0, we get 0
      return (sumweights) ? (sum / sumweights) : Type(0);
    }
  protected:
    Derived const& V_;
};

/** @ingroup Stat
 *  Compute the variance of the variable V.
 **/
template<class Derived>
struct VarianceOp
{
    typedef typename Derived::Type Type;
    /** constructor */
    inline VarianceOp( ExprBase<Derived> const&  V) :  V_(V.asDerived())
    { STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Derived);}
    /** @return the variance of the variable V.
     *  \f[ \hat{\sigma}^2 = \frac{1}{n} \sum_{i=1}^n (V(i)-\hat{\mu})^2. \f]
     *  @param unbiased @c true if we want an unbiased estimate of the variance,
     *  @c false otherwise
     **/
    inline Type const operator()(bool unbiased) const
    {
      // no samples
      if (V_.empty()) { return Type(0);}

      int nobs = V_.size();
      // Compute the mean and sum
      Type mu = MeanOp<Derived>(V_)();
      Type sum  = 0.0, var  = 0.0, dev;
      for (int i=V_.begin(); i<V_.end(); i++)
      {
        sum += (dev = V_[i] - mu); // deviation from the mean
        var += (dev*dev);         // squared value
      }
      // compute the variance
      return (unbiased) ?
      (nobs > 1) ? (var - (sum*sum)/Type(nobs))/Type(nobs -1)
                 : Arithmetic<Type>::infinity()
      :
      (nobs > 0) ? (var - (sum*sum)/(Type)nobs)/(Type)(nobs)
                 : Type(0);
    }
    /** @return the weighted variance of the variable V.
     *  \f[ \hat{\sigma}^2
     *  = \frac{\sum_{i=1}^n w(i)}{\left( \sum_{i=1}^n w(i))\right)^2-\sum_{i=1}^n w(i)^2}
     *    \sum_{i=1}^n w(i) (V(i)-\hat{\mu})^2.
     *  \f]
     * If there is no weights, this definition reduces to the usual
     * definition of the variance with factor 1/(n-1).
     *  @param w weights
     *  @param unbiased @c true if we want an unbiased estimate of the variance,
     *  @c false otherwise
     **/
    template< class Weights>
    inline Type const operator()( ExprBase<Weights> const&  w, bool unbiased) const
    {
      STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Weights);
      // no samples
      if (V_.empty()) { return Type(0);}
      // Compute the mean
      Type mu = MeanOp<Derived>(V_)(w);
      // sum the weighted samples
      Type sum = 0.0, sum1weights= 0.0, sum2weights = 0.0;
      for (int i=V_.begin(); i<V_.end(); i++)
      {
        Type const weight = w[i];
        sum1weights += weight;
        sum2weights += weight * weight;
        sum         += weight*(V_[i]-mu)*(V_[i]-mu); // deviation from the mean
      }
      // compute the variance
      return(unbiased) ?
      (sum1weights) ? (sum)/(sum1weights- sum2weights/sum1weights) : Type(0)
      :
      (sum1weights) ? (sum/sum1weights) : Type(0);
    }
  protected:
    Derived const& V_;
};

/** @ingroup Stat
 *  Compute safely the variance of the variable V.
 **/
template<class Derived>
struct VarianceSafeOp
{
    typedef typename Derived::Type Type;
    /** constructor */
    inline VarianceSafeOp( ExprBase<Derived> const&  V) :  V_(V.asDerived())
    { STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Derived);}
    /** @return the variance of the variable V discarding all missing values.
     *  \f[ \hat{\sigma}^2 = \frac{1}{n} \sum_{i=1}^n (V(i)-\hat{\mu})^2. \f]
     *  @param unbiased @c true if we want an unbiased estimate of the variance,
     *  @c false otherwise
     **/
    Type const operator()(bool unbiased) const
    {
      // no samples
      if (V_.empty()) { return Type(0);}
      int nobs = V_.size();
      // Compute the mean
      Type mu = MeanSafeOp<Derived>(V_)();
      // sum
      Type sum  = 0.0, var  = 0.0, dev;
      for (int i=V_.begin(); i<V_.end(); i++)
      {
        if (!Arithmetic<Type>::isNA(V_[i]))
        {
          sum += (dev = V_[i] - mu); // deviation from the mean
          var += (dev*dev);         // squared value
        }
        else nobs--;
      }
      // compute the variance
      return (unbiased) ?
      (nobs > 1) ? (var - (sum*sum)/Type(nobs))/Type(nobs -1)
                 : Arithmetic<Type>::infinity()
      :
      (nobs > 0) ? (var - (sum*sum)/(Type)nobs)/(Type)(nobs)
                 : Type(0);
    }
    /** @return the safely computed weighted variance of the variable V.
     *  \f[ \hat{\sigma}^2
     *  = \frac{\sum_{i=1}^n w(i)}{\left( \sum_{i=1}^n w(i))\right)^2-\sum_{i=1}^n w(i)^2}
     *    \sum_{i=1}^n w(i) (V(i)-\hat{\mu})^2.
     *  \f]
     * If there is no weights, this definition reduces to the usual
     * definition of the variance with factor 1/(n-1).
     *  @param w weights
     *  @param unbiased @c true if we want an unbiased estimate of the variance,
     *  @c false otherwise
    **/
    template< class Weights>
    Type const operator()( ExprBase<Weights> const&  w, bool unbiased) const
    {
      STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Weights);
      // no samples
      if (V_.empty()) { return Type(0);}
      // Compute the mean
      Type mu = MeanSafeOp<Derived>(V_)(w);
#ifdef STK_BOUNDS_CHECK
      if (V_.range() != w.range())
        STKRUNTIME_ERROR_NO_ARG(VarianceSafeOp,V.range()!=w.range());
#endif
      // sum the weighted samples
      Type sum = 0.0, sum1weights= 0.0, sum2weights = 0.0;
      for (int i=V_.begin(); i<V_.end(); i++)
      { if ( !Arithmetic<Type>::isNA(V_[i]) && !Arithmetic<Type>::isNA(w[i]) )
        {
          Type weight = std::abs(w[i]);
          sum1weights += weight;
          sum2weights += weight * weight;
          sum         += weight*(V_[i]-mu)*(V_[i]-mu); // deviation from the mean
        }
      }
      // compute the variance
      return (unbiased) ?
      (sum1weights) ? (sum)/(sum1weights- sum2weights/sum1weights) : 0.
      :
      (sum1weights) ? (sum/sum1weights) : 0.;
    }
  protected:
    Derived const& V_;
};

/** @ingroup Stat
 *  Compute the variance of the variable V when the mean is given.
 **/
template<class Derived>
struct VarianceWithFixedMeanOp
{
  typedef typename Derived::Type Type;
  /** constructor */
  inline VarianceWithFixedMeanOp( ExprBase<Derived> const&  V) :  V_(V.asDerived())
  { STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Derived);}
  /** @return the variance of the variable V with fixed mean.
   *  \f[ \hat{\sigma^2} = \frac{1}{n} \sum_{i=1}^n (V(i) - \mu)^2. \f]
   *  using a compensated algorithm.
   *  @note
   *  Chan, Tony F.; Golub, Gene H.; LeVeque, Randall J. (1983).
   *  Algorithms for Computing the Sample Variance: Analysis and Recommendations.
   *  The American Statistician 37, 242-247.
   *  @param mu the fixed mean
   *  @param unbiased @c true if we want an unbiased estimate of the variance,
   *  @c false otherwise
   **/
    Type const operator()( Type const& mu, bool unbiased) const
    {
      // no samples
      if (V_.empty()) { return Type(0);}
      int nobs = V_.size();
      // sum
      Type sum = 0., var = 0., dev;
      for (int i=V_.begin(); i<V_.end(); i++)
      {
        sum += (dev = V_[i] - mu); // deviation from the mean
        var += (dev*dev);         // squared value
      }
      // unbiased variance
      return (unbiased) ?
      (nobs > 1) ? (var - (sum*sum)/(Type)(nobs))/((Type)nobs -1)
                 : Arithmetic<Type>::infinity()
      :
      (nobs > 0) ? (var - (sum*sum)/Type(nobs))/Type(nobs) : Type(0);
    }
    /** @return the weighted variance of the variable V with fixed mean.
     *  \f[ \hat{\sigma^2} = \frac{1}{\sum_{i=1}^n w(i)}
     *                  \sum_{i=1}^n w(i) (V(i) - \mu)^2
     *  \f]
     *  In the unbiased case, the method use the compensated method described in
     *  http://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Weighted_sample_variance
     *  @param w weights
     *  @param mu the mean
     *  @param unbiased @c true if we want an unbiased estimate of the variance,
     *  @c false otherwise
     **/
    template< class Weights>
    Type const operator()( ExprBase<Weights> const& w, Type const& mu, bool unbiased) const
    {
      STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Weights);
      // no samples
      if (V_.empty()) { return Type(0);}
      // sum the weighted samples
      Type sum = 0.0, sum1weights = 0.0, sum2weights = 0.0;
      for (int i=V_.begin(); i<V_.end(); i++)
      {
        Real weight = std::abs(w[i]); // only positive weights
        sum1weights += weight;
        sum2weights += weight * weight;
        sum         += weight * (V_[i]-mu)*(V_[i]-mu);  // deviation from the mean
      }
      // compute the variance
      return (unbiased) ?
      (sum1weights) ?  (sum)/(sum1weights - sum2weights/sum1weights)
                    : Arithmetic<Type>::infinity()
      :
      (sum1weights) ? (sum/sum1weights) : Type(0);
    }
  protected:
    Derived const& V_;
};

/** @ingroup Stat
 *  Compute safely the variance of the variable V when the mean is given.
 **/
template<class Derived>
struct VarianceWithFixedMeanSafeOp
{
    typedef typename Derived::Type Type;
    /** constructor */
    inline VarianceWithFixedMeanSafeOp( ExprBase<Derived> const&  V) :  V_(V.asDerived())
    { STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Derived);}
    /** @return the safely computed variance of the variable V with fixed mean.
     *  \f[ \hat{\sigma^2} = \frac{1}{n} \sum_{i=1}^n (V(i) - \mu)^2. \f]
     *  using a compensated algorithm and discarding the missing values.
     *  @note
     *  Chan, Tony F.; Golub, Gene H.; LeVeque, Randall J. (1983).
     *  Algorithms for Computing the Sample Variance: Analysis and Recommendations.
     *  The American Statistician 37, 242-247.
     *  @param mu the fixed mean
     *  @param unbiased @c true if we want an unbiased estimate of the variance,
     *  @c false otherwise
     **/
    Type const operator()( Type mu, bool unbiased) const
    {
      // no samples
      if (V_.empty()) { return Type(0);}
      int nobs = V_.size();
      // sum
      Type sum  = 0.0, var  = 0.0, dev;
      for (int i=V_.begin(); i<V_.end(); i++)
      {
        if (!Arithmetic<Type>::isNA(V_[i]))
        {
          sum += (dev = V_[i] - mu); // deviation from the mean
          var += (dev*dev);         // squared value
        }
        else nobs--;
      }
      // compute the variance
      return (unbiased) ?
      (nobs > 1) ? (var - (sum*sum)/Type(nobs))/Type(nobs -1)
                 : Arithmetic<Type>::infinity()
      :
      (nobs > 0) ? (var - (sum*sum)/(Type)nobs)/(Type)(nobs)
                 : Type(0.);
    }

    /** @return the safely computed weighted variance of the variable V with fixed mean.
     *  \f[ \hat{\sigma^2} = \frac{1}{\sum_{i=1}^n w(i)}
     *                  \sum_{i=1}^n w(i) (V(i) - \mu)^2
     *  \f]
     *  @param w weights
     *  @param mu the mean
     *  @param unbiased @c true if we want an unbiased estimate of the variance,
     *  @c false otherwise
     **/
    template< class Weights>
    Type const operator()( ExprBase<Weights> const& w, Type const& mu, bool unbiased) const
    {
      STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Weights);
      // no samples
      if (V_.empty()) { return Type(0);}
      if (V_.range() != w.range())
        STKRUNTIME_ERROR_NO_ARG(wmeanSafe,V.range()!=w.range());
      // sum the weighted samples
      Type sum = 0.0, sum1weights = 0.0, sum2weights = 0.0;
      for (int i=V_.begin(); i<V_.end(); i++)
      {
        if ( !Arithmetic<Type>::isNA(V_[i]) && !Arithmetic<Type>::isNA(w[i]) )
        { Type weight = std::abs(w[i]);
          sum1weights += weight;
          sum2weights += weight * weight;
          sum         += weight*(V_[i]-mu)*(V_[i]-mu); // deviation from the mean
        }
      }
      // compute the variance
      return (unbiased) ?
        (sum1weights) ? (sum)/(sum1weights - sum2weights/sum1weights)
                      : Arithmetic<Type>::infinity()
        :
        (sum1weights) ? (sum/sum1weights) : Type(0);
    }
  protected:
    Derived const& V_;
};

/** @ingroup StatDesc
 *  Compute the minimal(s) value(s) of A. If A is a row-vector or a
 *  column-vector then the function will return the usual min of the vector.
 *  If A is a two-dimensional array, the function will return (by value) an
 *  Array2DPoint with the minimal values of each columns.
 *  @sa STK::Stat::mean, STK::Stat::max, STK::Stat::variance,
 *  STK::Stat::varianceWithFixedMean, STK::Stat::min, STK::Stat::sum.
 *  @param A the array
 *  @return the minimal value(s) of A or NA if there is no available value.
 **/
FUNCTION_ON_ARRAY(min, MinOp)

/** @ingroup StatDesc
 *  Compute safely the minimal(s) [weighted] value(s) of A. If A is a row-vector
 *  or a column-vector then the function will return the usual minimal value of
 *  the vector. If A is a two-dimensional array, the function will return
 *  (by value) an STK::Array2DPoint with the minimal values of each columns.
 *
 *  @sa STK::Stat::mean, STK::Stat::max, STK::Stat::variance, STK::Stat::varianceWithFixedMean,
 *  STK::Stat::min, STK::Stat::sum.
 *  @param A the array
 *  @return the minimal value(s) of A or NA if there is no available
 *  value.
 **/
FUNCTION_ON_ARRAY(minSafe, MinSafeOp)

/** @ingroup StatDesc
 *  Compute the maximal(s) value(s) of A. If A is a row-vector or a
 *  column-vector then the function will return the usual max of the vector.
 *  If A is a two-dimensional array, the function will return (by value) an
 *  Array2DPoint with the maximal values of each columns.
 *  @sa STK::Stat::mean, STK::Stat::max, STK::Stat::variance,
 *  STK::Stat::varianceWithFixedMean, STK::Stat::min, STK::Stat::sum.
 *  @param A the array
 *  @return the maximal value(s) of A or NA if there is no available
 *  value.
 **/
FUNCTION_ON_ARRAY(max, MaxOp)

/** @ingroup StatDesc
 *  Compute safely the maximal(s) value(s) of A. If A is a row-vector or a
 *  column-vector then the function will return the usual max of the vector.
 *  If A is a two-dimensional array, the function will return (by value) an
 *  Array2DPoint with the maximal values of each columns.
 *  @sa STK::Stat::mean, STK::Stat::max, STK::Stat::variance, STK::Stat::varianceWithFixedMean, STK::Stat::min, STK::Stat::sum.
 *  @param A the array
 *  @return the maximal value(s) of A or NA if there is no available
 *  value.
 **/
FUNCTION_ON_ARRAY(maxSafe, MaxSafeOp)

/** @ingroup StatDesc
 *  Compute the sum of A. If A is a row-vector or a
 *  column-vector then the function will return the usual sum of the vector.
 *  If A is a two-dimensional array, the function will return (by value) an
 *  Array2DPoint with the sum of each columns.
 *  @sa STK::Stat::mean, STK::Stat::max, STK::Stat::variance, STK::Stat::varianceWithFixedMean, STK::Stat::min
 *  @param A the data
 *  @return the mean(s) or NA if there is no available value
 **/
FUNCTION_ON_ARRAY(sum, SumOp)

/** @ingroup StatDesc
 *  Compute safely the mean(s) value(s) of A. If A is a row-vector or a
 *  column-vector then the function will return the usual mean of the vector.
 *  If A is a two-dimensional array, the function will return (by value) an
 *  Array2DPoint with the mean values of each columns.
 *  @sa STK::Stat::mean, STK::Stat::max, STK::Stat::variance, STK::Stat::varianceWithFixedMean, STK::Stat::min, STK::Stat::sum.
 *  @param A the data
 *  @return the mean(s) or NA if there is no available value
 **/
FUNCTION_ON_ARRAY(sumSafe, SumSafeOp)

/** @ingroup StatDesc
 *  Compute the mean(s) value(s) of A. If A is a row-vector or a
 *  column-vector then the function will return the usual mean of the vector.
 *  If A is a two-dimensional array, the function will return (by value) an
 *  Array2DPoint with the mean values of each columns.
 *  @sa STK::Stat::mean, STK::Stat::max, STK::Stat::variance,
 *  STK::Stat::varianceWithFixedMean, STK::Stat::min, STK::Stat::sum.
 *  @param A the data
 *  @return the mean(s) or NA if there is no available value
 **/
FUNCTION_ON_ARRAY(mean, MeanOp)

/** @ingroup StatDesc
 *  Compute safely the mean(s) value(s) of A. If A is a row-vector or a
 *  column-vector then the function will return the usual mean of the vector.
 *  If A is a two-dimensional array, the function will return (by value) an
 *  Array2DPoint with the mean values of each columns.
 *  @sa STK::Stat::mean, STK::Stat::max, STK::Stat::variance,
 *  STK::Stat::varianceWithFixedMean, STK::Stat::min, STK::Stat::sum.
 *  @param A the data
 *  @return the mean(s) or NA if there is no available value
 **/
FUNCTION_ON_ARRAY(meanSafe, MeanSafeOp)

/** @ingroup StatDesc
 *  Compute the variance(s) value(s) of A. If A is a row-vector or a
 *  column-vector then the function will return the usual variance of the vector.
 *  If A is a two-dimensional array, the function will return (by value) an
 *  Array2DPoint with the variance values of each columns.
 *  @sa STK::Stat::mean, STK::Stat::max, STK::Stat::variance, STK::Stat::varianceWithFixedMean, STK::Stat::min, STK::Stat::sum.
 *  @param A the data
 *  @param unbiased the unbiased variance(s) if @c true or the Maximum-likelihood
 *  variance otherwise (the default)
 *  @return the variance(s) or NA if there is no available value
 **/
template< class Derived>
typename hidden::FunctorTraits<Derived, VarianceOp >::Row
variance(Derived const& A, bool unbiased = false)
{ return typename hidden::FunctorTraits<Derived, VarianceOp>::ColOp(A)(unbiased);}

template< class Derived, class Weights>
typename hidden::FunctorTraits<Derived, VarianceOp >::Row
variance(Derived const& A, ExprBase<Weights> const& w, bool unbiased = false)
{ return typename hidden::FunctorTraits<Derived, VarianceOp>::ColWeightedOp(A)(w, unbiased);}

template< class Derived>
typename hidden::FunctorTraits<Derived, VarianceOp >::Row
varianceByCol(Derived const& A, bool unbiased = false)
{ return typename hidden::FunctorTraits<Derived, VarianceOp>::ColOp(A)(unbiased);}

template< class Derived, class Weights>
typename hidden::FunctorTraits<Derived, VarianceOp >::Row
varianceByCol(Derived const& A, ExprBase<Weights> const& w, bool unbiased = false)
{ return typename hidden::FunctorTraits<Derived, VarianceOp>::ColWeightedOp(A)(w, unbiased);}

template< class Derived>
typename hidden::FunctorTraits<Derived, VarianceOp >::Col
varianceByRow(Derived const& A, bool unbiased = false)
{ return typename hidden::FunctorTraits<Derived, VarianceOp>::RowOp(A)(unbiased);}


template< class Derived, class Weights>
typename hidden::FunctorTraits<Derived, VarianceOp >::Col
varianceByRow(Derived const& A, ExprBase<Weights> const& w, bool unbiased = false)
{ return typename hidden::FunctorTraits<Derived, VarianceOp>::RowWeightedOp(A)(w, unbiased);}

/** @ingroup StatDesc
 *  Compute safely the variance(s) value(s) of A. If A is a row-vector or a
 *  column-vector then the function will return the usual variance of the vector.
 *  If A is a two-dimensional array, the function will return (by value) an
 *  Array2DPoint with the variance values of each columns.
 *  @sa STK::Stat::mean, STK::Stat::max, STK::Stat::variance,
 *  STK::Stat::varianceWithFixedMean, STK::Stat::min, STK::Stat::sum.
 *  @param A the data
 *  @param unbiased the unbiased variance(s) if @c true or the Maximum-likelihood
 *  variance otherwise (the default)
 *  @return the variance(s) or NA if there is no available value
 **/
template< class Derived>
typename hidden::FunctorTraits<Derived, VarianceSafeOp >::Row
varianceSafe(Derived const& A, bool unbiased = false)
{ return typename hidden::FunctorTraits<Derived, VarianceSafeOp>::ColOp(A)(unbiased);}

template< class Derived, class Weights>
typename hidden::FunctorTraits<Derived, VarianceSafeOp >::Row
varianceSafe(Derived const& A, ExprBase<Weights> const& w, bool unbiased = false)
{ return typename hidden::FunctorTraits<Derived, VarianceSafeOp>::ColWeightedOp(A)(w, unbiased);}

template< class Derived>
typename hidden::FunctorTraits<Derived, VarianceSafeOp >::Row
varianceSafeByCol(Derived const& A, bool unbiased = false)
{ return typename hidden::FunctorTraits<Derived, VarianceSafeOp>::ColOp(A)(unbiased);}

template< class Derived, class Weights>
typename hidden::FunctorTraits<Derived, VarianceSafeOp >::Row
varianceSafeByCol(Derived const& A, ExprBase<Weights> const& w, bool unbiased = false)
{ return typename hidden::FunctorTraits<Derived, VarianceSafeOp>::ColWeightedOp(A)(w, unbiased);}

template< class Derived>
typename hidden::FunctorTraits<Derived, VarianceSafeOp >::Col
varianceSafeByRow(Derived const& A, bool unbiased = false)
{ return typename hidden::FunctorTraits<Derived, VarianceSafeOp>::RowOp(A)(unbiased);}

template< class Derived, class Weights>
typename hidden::FunctorTraits<Derived, VarianceSafeOp >::Col
varianceSafeByRow(Derived const& A, ExprBase<Weights> const& w, bool unbiased = false)
{ return typename hidden::FunctorTraits<Derived, VarianceSafeOp>::RowWeightedOp(A)(w, unbiased);}

/** @ingroup StatDesc
 *  Compute the VarianceWithFixedMean(s) value(s) of A. If A is a row-vector or a
 *  column-vector then the function will return the usual variance of the vector.
 *  If A is a two-dimensional array, the function will return (by value) an
 *  Array2DPoint with the variance values of each columns.
 *  @sa STK::Stat::mean, STK::Stat::max, STK::Stat::variance,
 *  STK::Stat::varianceWithFixedMean, STK::Stat::min, STK::Stat::sum.
 *  @param A the data
 *  @param mean The mean (s) to use
 *  @param unbiased the unbiased variance(s) if @c true or the Maximum-likelihood
 *  variance otherwise (the default)
 *  @return the variance(s) or NA if there is no available value
 **/
template< class Derived, class MeanType>
typename hidden::FunctorTraits<Derived, VarianceWithFixedMeanOp >::Row
varianceWithFixedMean(Derived const& A, MeanType const& mean, bool unbiased)
{ return typename hidden::FunctorTraits<Derived, VarianceWithFixedMeanOp>::ColOp(A)(mean, unbiased);}

template< class Derived, class MeanType, class Weights>
typename hidden::FunctorTraits<Derived, VarianceWithFixedMeanOp >::Row
varianceWithFixedMean(Derived const& A, ExprBase<Weights> const& w, MeanType const& mean, bool unbiased)
{ return typename hidden::FunctorTraits<Derived, VarianceWithFixedMeanOp>::ColWeightedOp(A)(w, mean, unbiased);}

template< class Derived, class MeanType>
typename hidden::FunctorTraits<Derived, VarianceWithFixedMeanOp >::Row
varianceWithFixedMeanByCol(Derived const& A, MeanType const& mean, bool unbiased)
{ return typename hidden::FunctorTraits<Derived, VarianceWithFixedMeanOp>::ColOp(A)(mean, unbiased);}

template< class Derived, class MeanType, class Weights>
typename hidden::FunctorTraits<Derived, VarianceWithFixedMeanOp >::Row
varianceWithFixedMeanByCol(Derived const& A, ExprBase<Weights> const& w, MeanType const& mean, bool unbiased)
{ return typename hidden::FunctorTraits<Derived, VarianceWithFixedMeanOp>::ColWeightedOp(A)(w, mean, unbiased);}

template< class Derived, class MeanType>
typename hidden::FunctorTraits<Derived, VarianceWithFixedMeanOp >::Col
varianceWithFixedMeanByRow(Derived const& A, MeanType const& mean, bool unbiased = false)
{ return typename hidden::FunctorTraits<Derived, VarianceWithFixedMeanOp>::RowOp(A)(mean, unbiased);}

template< class Derived, class MeanType, class Weights>
typename hidden::FunctorTraits<Derived, VarianceWithFixedMeanOp >::Col
varianceWithFixedMeanByRow(Derived const& A, ExprBase<Weights> const& w, MeanType const& mean, bool unbiased)
{ return typename hidden::FunctorTraits<Derived, VarianceWithFixedMeanOp>::RowWeightedOp(A)(w, mean, unbiased);}

/** @ingroup StatDesc
 *  Compute safely the VarianceWithFixedMean(s) value(s) of A.
 *  @sa STK::Stat::mean, STK::Stat::max, STK::Stat::variance, STK::Stat::varianceWithFixedMean, STK::Stat::min, STK::Stat::sum.
 *  @param A the data
 *  @param mean The mean (s) to use
 *  @param unbiased the unbiased variance(s) if @c true or the Maximum-likelihood
 *  variance otherwise (the default)
 *  @return the variance(s) or NA if there is no available value
 **/
template< class Derived, class MeanType>
typename hidden::FunctorTraits<Derived, VarianceWithFixedMeanSafeOp >::Row
varianceWithFixedMeanSafe( Derived const& A, MeanType const& mean, bool unbiased)
{ return typename hidden::FunctorTraits<Derived, VarianceWithFixedMeanSafeOp>::ColOp(A)(mean, unbiased);}

template< class Derived, class MeanType, class Weights>
typename hidden::FunctorTraits<Derived, VarianceWithFixedMeanSafeOp >::Row
varianceWithFixedMeanSafe( Derived const& A, ExprBase<Weights> const& w, MeanType const& mean, bool unbiased)
{ return typename hidden::FunctorTraits<Derived, VarianceWithFixedMeanSafeOp>::ColWeightedOp(A)(w, mean, unbiased);}

template< class Derived, class MeanType>
typename hidden::FunctorTraits<Derived, VarianceWithFixedMeanSafeOp >::Row
varianceWithFixedMeanSafeByCol( Derived const& A, MeanType const& mean, bool unbiased)
{ return typename hidden::FunctorTraits<Derived, VarianceWithFixedMeanSafeOp>::ColOp(A)(mean, unbiased);}

template< class Derived, class MeanType, class Weights>
typename hidden::FunctorTraits<Derived, VarianceWithFixedMeanSafeOp >::Row
varianceWithFixedMeanSafeByCol( Derived const& A, ExprBase<Weights> const& w, MeanType const& mean, bool unbiased)
{ return typename hidden::FunctorTraits<Derived, VarianceWithFixedMeanSafeOp>::ColWeightedOp(A)(w, mean, unbiased);}

template< class Derived, class MeanType>
typename hidden::FunctorTraits<Derived, VarianceWithFixedMeanSafeOp >::Col
varianceWithFixedMeanSafeByRow(Derived const& A, MeanType const& mean, bool unbiased)
{ return typename hidden::FunctorTraits<Derived, VarianceWithFixedMeanSafeOp>::RowOp(A)(mean, unbiased);}

template< class Derived, class MeanType, class Weights>
typename hidden::FunctorTraits<Derived, VarianceWithFixedMeanSafeOp >::Col
varianceWithFixedMeanSafeByRow(Derived const& A, ExprBase<Weights> const& w, MeanType const& mean, bool unbiased)
{ return typename hidden::FunctorTraits<Derived, VarianceWithFixedMeanOp>::RowWeightedOp(A)(w, mean, unbiased);}


}  // namespace Stat

}  // namespace STK

#undef FUNCTION_ON_ARRAY

#endif /*STK_STAT_FUNCTORS_H*/
