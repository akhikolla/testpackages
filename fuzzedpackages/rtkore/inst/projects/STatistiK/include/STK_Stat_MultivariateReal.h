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
 * Project:  stkpp::Stat
 * Purpose:  Compute multivariate elementary statistics for a Real 2D container.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Stat_MultivariateReal.h
 *  @brief In this file we specialize the class Multivariate to Real type.
 **/

#ifndef STK_STAT_MULTIVARIATEREAL_H
#define STK_STAT_MULTIVARIATEREAL_H

#include "STKernel/include/STK_Misc.h"

#include "Arrays/include/STK_Array2DPoint.h"
#include "Arrays/include/STK_Array2DSquare.h"

#include "STK_Stat_Covariance.h"
#include "STK_Stat_Multivariate.h"
#include "STK_Stat_Functors.h"

namespace STK
{
namespace Stat
{
/** @ingroup StatDesc
 *  @brief Computation of the Multivariate Statistics of a 2D Container
 *  of Real.
 *
 *  The class @c Multivariate is just a factory class for computing the mean,
 *  the variance and the covariance square matrix of a p_data set stored in a
 *  @c Array with n rows (the samples) and p columns (the variables).
 **/
template <class Array>
class Multivariate<Array, Real>: public IRunnerUnsupervised< Array, typename hidden::Traits<Array>::Col>
{
  typedef typename hidden::Traits<Array>::Row RowVector;
  typedef typename hidden::Traits<Array>::Col ColVector;
  /** type of runner */
  typedef IRunnerUnsupervised< Array, ColVector> Runner;
  typedef typename Array::Type Type_;
  enum
  {
    value_ = hidden::isSame<Type_, Real>::value_
  };

  public:
    /** Default Constructor. */
    Multivariate(): Runner()
                  , nbSamples_(0), nbVar_(0)
                  , min_(), max_(), mean_(), var_(), cov_()
  { STK_STATIC_ASSERT(value_, YOU_CANNOT_USED_THIS_TYPE_OF_DATA_WITH_THIS_OBJECT);}
    /** Constructor.
     *  @param data a reference on the data set
     **/
    Multivariate( Array const& data)
                : Runner(&data)
                , nbSamples_(data.sizeRows()), nbVar_(data.sizeCols())
                , min_(), max_(), mean_(), var_(), cov_()
    { STK_STATIC_ASSERT(value_, YOU_CANNOT_USED_THIS_TYPE_OF_DATA_WITH_THIS_OBJECT);}
    /** Constructor.
     *  @param p_data a pointer on the data set
     **/
    Multivariate( Array const* p_data)
                : Runner(p_data)
                , nbSamples_((p_data) ? p_data->sizeRows() : 0)
                , nbVar_((p_data) ? p_data->sizeCols() : 0)
                , min_(), max_(), mean_(), var_(), cov_()
    { STK_STATIC_ASSERT(value_, YOU_CANNOT_USED_THIS_TYPE_OF_DATA_WITH_THIS_OBJECT);}
    /** copy constructor.
     *  @param stat the statistics to copy
     **/
    Multivariate( Multivariate const& stat)
                : Runner(stat)
                , nbSamples_(stat.nbSamples_), nbVar_(stat.nbVar_)
                , min_(stat.min_), max_(stat.max_), mean_(stat.mean_)
                , var_(stat.var_), cov_(stat.cov_)
    {}
    /** virtual destructor.*/
    virtual ~Multivariate() { }
    /** clone pattern */
    inline virtual Multivariate* clone() const { return new Multivariate(*this);}
    /** @return the number of variables in the p_data_ set (the number of columns) */
    inline int nbVariable() const {return nbVar_;}
    /** @return the number of samples in the p_data_ set (the number of rows) */
    inline int nbSamples() const {return nbSamples_;}
    /** @return the minimal values of the variables in a RowVector */
    inline PointX const& min() const { return min_;}
    /** @return the maximal values of the variables in a RowVector  */
    inline PointX const& max() const { return max_;}
    /** @return the mean of the variables in a RowVector */
    inline PointX const& mean() const { return mean_;}
    /**@return the variance of the variables in a RowVector */
    inline PointX const& variance() const { return var_;}
    /** @return the covariance of the variables in a square matrix */
    inline ArraySquareX const& covariance() const { return cov_;}

    /** run the estimation of the Multivariate statistics. **/
    virtual bool run()
    {
      if (!this->p_data_)
      { this->msg_error_ = STKERROR_NO_ARG(MultivariateArray::run(),data is not set);
        return false;
      }
      try
      {
        mean_.move(Stat::mean(*this->p_data_));
        min_.move(Stat::min(*this->p_data_));
        max_.move(Stat::max(*this->p_data_));
        var_.move(varianceWithFixedMean(*this->p_data_, mean_, false));

        cov_.resize(this->p_data_->cols());
        for (int j=this->p_data_->beginCols(); j<this->p_data_->endCols(); j++)
        {
          cov_(j, j) = var_[j];
          for (int i=this->p_data_->beginCols(); i<j; i++)
          {
            cov_(i, j) = covarianceWithFixedMean(this->p_data_->col(i), this->p_data_->col(j), mean_[i], mean_[j], false);
            cov_(j, i) = cov_(i, j);
          }
        }
      }
      catch (Exception const& error)
      {
        this->msg_error_ += _T("Error in Multivariate::run():\nWhat: ");
        this->msg_error_ += error.error();
        return false;
      }
      // no error
      return true;
    }
    /** run the estimation of the weighted multivariate statistics.
     * @param weights the weights of the samples
     **/
    virtual bool run( ColVector const& weights)
    {
      if (!this->p_data_)
      { this->msg_error_ = STKERROR_NO_ARG(MultivariateArray::run(weights),data is not set);
        return false;
      }
      if (this->p_data_->rows() != weights.range())
      { this->msg_error_ = STKERROR_NO_ARG(MultivariateArray::run(weights),p_data_->rows() != weights.range());
        return false;
      }
      try
      {
        mean_.move(Stat::mean(*this->p_data_, weights));
        min_.move(Stat::min(*this->p_data_,weights));
        max_.move(Stat::max(*this->p_data_, weights));
        var_.move(varianceWithFixedMean(*this->p_data_, weights, mean_, false));

        cov_.resize(this->p_data_->cols());
        // for each variables
        for (int j= this->p_data_->beginCols(); j< this->p_data_->endCols(); j++)
        {
          cov_(j, j) = var_[j];
          // compute the covariances
          for (int i= this->p_data_->beginCols(); i<j; i++)
          {
            cov_(i, j) = covarianceWithFixedMean(this->p_data_->col(i), this->p_data_->col(j), weights, mean_[i], mean_[j]);
            cov_(j, i) = cov_(i, j);
          }
        }
      }
      catch (Exception const& error)
      {
        this->msg_error_  = _T("Error in Multivariate::run(weights): ");
        this->msg_error_ += error.error();
        return false;
      }
      // no error
      return true;
    }

  protected:
    /** number of samples */
    int nbSamples_;
    /** Number of variables */
    int nbVar_;

    /** Vector of the mean of the Variables */
    PointX min_;
    /** Vector of the mean of the Variables */
    PointX max_;
    /** Vector of the mean of the Variables */
    PointX mean_;
    /** Vector of the variance of the variables */
    PointX var_;
    /** Array of the covariance of the variables */
    ArraySquareX cov_;

    /** udpating method in case we set a new data set */
    virtual void update()
    {
      // if there is no data there is nothing to update
      if (this->p_data_)
      {
        nbSamples_ = this->p_data_->sizeRows();
        nbVar_     = this->p_data_->sizeCols();
      }
    }
};


}  // namespace Stat

}  // namespace STK

#endif /*STK_STAT_MULTIVARIATEREAL_H*/
