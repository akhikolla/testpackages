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
 * Project:  stkpp::Model
 * created on: 22 juil. 2011
 * Purpose: define the class IUnivStatModel.
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_ModelGamma_aj_bj.h
 *  @brief In this file we define the class ModelGamma_aj_bj.
 **/

#ifndef STK_JOINTGAMMAMODEL_H
#define STK_JOINTGAMMAMODEL_H

#include <cmath>

#include "STK_Model_Util.h"
#include "STK_IMultiStatModel.h"
#include <Arrays/include/STK_Array2DPoint.h>

#include <STatistiK/include/STK_Law_Gamma.h>
#include <Analysis/include/STK_Algo_FindZero.h>
#include <Analysis/include/STK_Funct_raw.h>

namespace STK
{
// forward declaration
template <class Array, class WColVector = CVectorX> class ModelGamma_aj_bj;
struct ModelGamma_aj_bjParameters;

namespace hidden
{
/** @ingroup hidden
 *  Partial specialization of the StatModelTraits for the ModelGamma_aj_bj
 *  model
 **/
template<class Data_, class WColVector_>
struct StatModelTraits< ModelGamma_aj_bj<Data_, WColVector_> >
{
  /** Type of the container storing the data */
  typedef DataBridge<Data_> Data;
  /** Type of the array storing the weights of the data */
  typedef WColVector_ WColVector;
  /** Type of the data in the container */
  typedef typename Data::Type Type;
  /** Type of the parameters of the ModelGamma_aj_bj */
  typedef ModelGamma_aj_bjParameters Parameters;
};

} // hidden

/** @ingroup StatModels
 *  Structure encapsulating the parameters of a Joint Gamma model.
 *  @sa GammaComponent, ModelGamma_aj_bj
 */
struct ModelGamma_aj_bjParameters
{
  /** default constructor */
  ModelGamma_aj_bjParameters() : shape_(), scale_(), mean_(), meanLog_(), variance_() {}
  /** default constructor */
  ModelGamma_aj_bjParameters( Range const& range)
                            : shape_(range, 1.), scale_(range, 1.)
                            , mean_(range, 1.), meanLog_(range, 0.), variance_(range, 1.)
  {}
  /** copy constructor. @param param the parameters to copy. */
  ModelGamma_aj_bjParameters( ModelGamma_aj_bjParameters const& param)
                            : shape_(param.shape_), scale_(param.scale_)
                            , mean_(param.mean_), meanLog_(param.meanLog_), variance_(param.variance_)
  {}
  /** destructor */
  ~ModelGamma_aj_bjParameters() {}
  /** @return the vector of the mean of the observations */
  inline CPointX const& shape() const {return shape_;}
  /** vector of the mean log of the observations */
  inline CPointX const& scale() const {return scale_;}

  /** resize the parameters only if the range is modified, otherwise, stay
   *  with the current values.
   *  @param range the range of the parameters (= range of the variables of the model)
   **/
  inline void resize(Range const& range)
  {
    if (range != range_)
    {
      shape_.resize(range); shape_ = 1.;
      scale_.resize(range); scale_ = 1.;
      mean_.resize(range); mean_ = 1.;
      meanLog_.resize(range); meanLog_ = 0.;
      variance_.resize(range); variance_ = 1.;
      range_ = range;
    }
  }
  /** vector of the shape */
  CPointX shape_;
  /** vector of the scale */
  CPointX scale_;
  /** vector of the mean of the observations */
  CPointX mean_;
  /** vector of the mean log of the observations */
  CPointX meanLog_;
  /** vector of the variance of the observations */
  CPointX variance_;
  Range range_;
};


/** @ingroup StatModels
 * A joint Gamma model is a statistical model of the following form
 * \f[
 *     f(\mathbf{x}_i|\theta) =
 *     \prod_{j=1}^p\left(\frac{x_i^j}{b_j}\right)^{a_j-1}
 *                   \frac{e^{-x_i^j/b_j}}{b_j \, \Gamma(a_j)},
 *      \quad x_i^j\in\{0,1\}, \quad j=1,\ldots,p, \quad i=1,\ldots,n.
 * \f]
 **/
template <class Data_, class WColVector_>
class ModelGamma_aj_bj: public IMultiStatModel< ModelGamma_aj_bj<Data_, WColVector_> >
{
  protected:
    class dloglikelihood: public IFunction<dloglikelihood >
    {
      public:
        dloglikelihood( Real const& mean, Real const& meanLog)
                      : delta_(meanLog - Real(std::log(mean))) {}
        /** @return the value of the function at a
         * @param a a positive real value
         **/
        inline Real fImpl(Real a) const
        { return (delta_ + std::log(a) - Funct::psi_raw(a));}
        /** @return the minimal value of the function at x */
        inline Real xminImpl() const { return 0;}
      private:
        Real delta_;
    };

  public:
    /** Type of the container storing the data */
    typedef DataBridge<Data_> Data;
    typedef typename hidden::Traits<Data_>::Row RowVector;
    /** Type of the array storing the weights of the data */
    typedef WColVector_ WColVector;
    /** Type of the data in the container */
    typedef typename Data::Type Type;
    /** Base class */
    typedef IMultiStatModel< ModelGamma_aj_bj<Data_, WColVector_> > Base;
    using Base::p_data;
    using Base::param;

    /** default constructor. */
    ModelGamma_aj_bj(): Base() {}
    /** Constructor with data set. */
    ModelGamma_aj_bj(Data const& data): Base(data){}
    /** Constructor with a ptr on the data set. */
    ModelGamma_aj_bj(Data const* p_data): Base(p_data){}
    /** Copy constructor. */
    ModelGamma_aj_bj(ModelGamma_aj_bj const& model): Base(model) {}
    /** destructor */
    ~ModelGamma_aj_bj(){}

    /** @return the vector of the mean of the observations */
    inline CPointX const& shape() const {return param().shape_;}
    /** @return the vector of the mean of the observations */
    inline CPointX const& scale() const {return param().scale_;}
    /** @return the vector of the mean of the observations */
    inline CPointX const& mean() const {return param().mean_;}
    /** vector of the mean log of the observations */
    inline CPointX const& meanLog() const {return param().meanLog_;}
    /** vector of the variance of the observations */
    inline CPointX const& variance() const {return param().variance_;}

    /** compute the number of free parameters */
    int computeNbFreeParameters() const { return 2*p_data()->dataij().sizeCols();}
    /** compute the log Likelihood of an observation. */
    Real computeLnLikelihood( RowVector const& rowData) const;
    /** compute the parameters */
    void computeParameters();
    /** compute the weighted parameters */
    void computeParameters( WColVector const& weights);
    /** Write the parameters on the output stream os */
    void writeParametersImpl(ostream& os) const;

  private:
    /** @return the vector of the shape of the observations */
    inline CPointX& shape() { return param().shape_;}
    /** @return the vector of the scale of the observations */
    inline CPointX& scale() { return param().scale_;}
    /** @return the vector of the mean of the observations */
    inline CPointX& mean() { return param().mean_;}
    /** vector of the mean log of the observations */
    inline CPointX& meanLog() { return param().meanLog_;}
    /** vector of the variance of the observations */
    inline CPointX& variance() { return param().variance_;}
};

/* compute the log Likelihood of an observation. */
template<class Data_, class WColVector_>
Real ModelGamma_aj_bj<Data_, WColVector_>::computeLnLikelihood( RowVector const& rowData) const
{
  Real sum =0.;
  for (Integer j= rowData.begin(); j < rowData.end(); ++j)
  { sum += Law::Gamma::lpdf(rowData[j], shape()[j], scale()[j]);}
  return sum;
}

/* compute the parameters */
template<class Data_, class WColVector_>
void ModelGamma_aj_bj<Data_, WColVector_>::computeParameters()
{
  for (int j=p_data()->dataij().beginCols(); j < p_data()->dataij().endCols(); ++j)
  {
    mean()[j] =  p_data()->dataij().col(j).meanSafe();
    meanLog()[j] = p_data()->dataij().col(j).safe(1.).log().mean();
    variance()[j] = p_data()->dataij().col(j).safe().variance();
    Real x0 = (mean()[j]*mean()[j]) / variance()[j];
    Real x1 = 0.9*x0 +  0.05/(mean()[j] - meanLog()[j]);
    dloglikelihood funct(mean()[j], meanLog()[j]);
    Real a =  Algo::findZero(funct, x0, x1);
    // replace with moment estimate if needed
    if (!Arithmetic<Real>::isFinite(a)) { a =  mean()[j]*mean()[j]/variance()[j];}
    shape()[j] = a;
    scale()[j] = mean()[j]/a;
  }
}
/* compute the weighted parameters */
template<class Data_, class WColVector_>
void ModelGamma_aj_bj<Data_, WColVector_>::computeParameters( WColVector const& weights)
{
  for (int j=p_data()->dataij().beginCols(); j < p_data()->dataij().endCols(); ++j)
  {
    mean()[j] =  p_data()->dataij().col(j).safe().wmean(weights);
    meanLog()[j] = p_data()->dataij().col(j).safe(1).log().wmean(weights);
    variance()[j] = p_data()->dataij().col(j).safe().wvariance(weights);
    Real x0 = (mean()[j]*mean()[j]) / variance()[j];
    Real x1 = 0.9*x0 +  0.05/(mean()[j] - meanLog()[j]);
    dloglikelihood funct(mean()[j], meanLog()[j]);
    Real a =  Algo::findZero(funct, x0, x1);
    // replace with moment estimate if needed
    if (!Arithmetic<Real>::isFinite(a)) { a =  mean()[j]*mean()[j]/variance()[j];}
    shape()[j] = a;
    scale()[j] = mean()[j]/a;
  }
}

/* Write the parameters on the output stream os */
template<class Data_, class WColVector_>
void ModelGamma_aj_bj<Data_, WColVector_>::writeParametersImpl(ostream& os) const
{
  os << _T("shape = ") << shape();
  os << _T("scale = ") << scale();
}


} // namespace STK

#endif /* STK_JOINTGAMMAMODEL_H */
