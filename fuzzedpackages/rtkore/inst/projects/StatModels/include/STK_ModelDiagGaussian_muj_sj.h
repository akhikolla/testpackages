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

/** @file STK_ModelDiagGaussian_muj_sj.h
 *  @brief In this file we define the class ModelDiagGaussian_muj_sj.
 **/

#ifndef STK_GAUSSIAN_MUJ_SJMODEL_H
#define STK_GAUSSIAN_MUJ_SJMODEL_H

#include <cmath>

#include "STK_Model_Util.h"
#include "STK_IMultiStatModel.h"
#include <STatistiK/include/STK_Law_Normal.h>
#include <STatistiK/include/STK_Stat_Functors.h>

namespace STK
{

// forward declaration
template <class Data_, class WColVector = CVectorX> class ModelDiagGaussian_muj_sj;
struct DiagGaussian_muj_sjParameters;

namespace hidden
{
/** @ingroup hidden
 *  Partial specialization of the StatModelTraits for the ModelDiagGaussian_muj_sj
 *  model
 **/
template<class Data_, class WColVector_>
struct StatModelTraits< ModelDiagGaussian_muj_sj<Data_, WColVector_> >
{
  /** Type of the container storing the data */
  typedef DataBridge<Data_> Data;
  /** Type of the array storing the weights of the data */
  typedef WColVector_ WColVector;
  /** Type of the data in the container */
  typedef typename Data::Type Type;
  /** Type of the parameters of the ModelDiagGaussian_muj_sj */
  typedef DiagGaussian_muj_sjParameters Parameters;
};

} // hidden

/** @ingroup StatModels
 *  Structure encapsulating the parameters of a DiagGaussian_muj_sj model.
 */
struct DiagGaussian_muj_sjParameters
{
    /** default constructor */
    DiagGaussian_muj_sjParameters(): mu_(), sigma_() {}
    /** default constructor */
    DiagGaussian_muj_sjParameters(Range const& range): mu_(range, 0.), sigma_(range, 1.) {}
    /** copy constructor. @param param the parameters to copy. */
    DiagGaussian_muj_sjParameters( DiagGaussian_muj_sjParameters const& param)
                                 : mu_(param.mu_)
                                 , sigma_(param.sigma_)
    {}
    /** destructor */
    ~DiagGaussian_muj_sjParameters() {}

    /** @return the means */
    inline CPointX const& mu() const { return mu_;}
    /** @return the mean of the jth law */
    inline CPointX const& sigma() const { return sigma_;}

    /** @return the mean of the jth law */
    inline Real const mu(int j) const { return mu_[j];}
    /** @return the standard deviation of the jth law */
    inline Real const sigma(int j) const { return sigma_[j];}

    /** resize the parameters only if the range is modified, otherwise, stay
     *  with the current values.
     *  @param range the range of the parameters (= range of the variables of the model)
     **/
    inline void resize(Range const& range)
    {
      if (range != range_)
      {
        mu_.resize(range); sigma_.resize(range);
        range_ = range;
      }
    }
    CPointX mu_;
    CPointX sigma_;
    Range range_;
};

/** @ingroup StatModels
 * A DiagGaussian_muj_sj model is a statistical model of the form:
 * following form
 * \f[
 *     f(\mathbf{x}_i|\theta) =
 *     \prod_{j=1}^p
 *     \frac{1}{\sigma_j\sqrt{2\pi}}
 *            \exp\left(-\frac{\left(x_i-\mu_j\right)^2}{2\sigma_j^2} \right)
 *    \quad j=1,\ldots,p, \quad i=1,\ldots,n.
 * \f]
 *
 **/
template <class Data_, class WColVector_>
class ModelDiagGaussian_muj_sj: public IMultiStatModel< ModelDiagGaussian_muj_sj<Data_, WColVector_> >
{
  public:
    /** Type of the container storing the data */
    typedef DataBridge<Data_> Data;
    typedef typename hidden::Traits<Data_>::Row RowVector;
     /** Type of the array storing the weights of the data */
    typedef WColVector_ WColVector;
    /** Type of the data in the container */
    typedef typename Data::Type Type;
    /** Type of the row vector of the container */
    /** Type of the parameters of the ModelDiagGaussian_muj_sj */
    typedef DiagGaussian_muj_sjParameters Parameters;
    /** Base class */
    typedef IMultiStatModel< ModelDiagGaussian_muj_sj<Data_, WColVector_> > Base;
    using Base::p_data;
    using Base::param;

    /** default constructor. */
    ModelDiagGaussian_muj_sj(): Base() {}
    /** Constructor with data set. */
    ModelDiagGaussian_muj_sj(Data const& data): Base(data) {}
    /** Constructor with a ptr on the data set. */
    ModelDiagGaussian_muj_sj(Data const* p_data): Base(p_data) {}
    /** Copy constructor. */
    ModelDiagGaussian_muj_sj(ModelDiagGaussian_muj_sj const& model): Base(model) {}
    /** destructor */
    ~ModelDiagGaussian_muj_sj() {}

    /** @return the vector of the mean of the observations */
    inline CPointX const& mean() const {return param().mu();}
    /** vector of the mean log of the observations */
    inline CPointX const& sigma() const {return param().sigma();}

    /** compute the number of free parameters */
    inline int computeNbFreeParameters() const { return 2*p_data()->dataij().sizeCols();}
    /** compute the log Likelihood of an observation. */
    Real computeLnLikelihood( RowVector const& rowData) const;
    /** compute the parameters */
    void computeParameters();
    /** compute the weighted parameters */
    void computeParameters( WColVector const& weights);
    /** Write the parameters on the output stream os */
    void writeParametersImpl(ostream& os) const;
};

/* compute the log Likelihood of an observation. */
template<class Data_, class WColVector_>
Real ModelDiagGaussian_muj_sj<Data_, WColVector_>::computeLnLikelihood( RowVector const& rowData) const
{
  Real sum =0.;
  for (Integer j= rowData.begin(); j < rowData.end(); ++j)
  { sum += Law::Normal::lpdf(rowData[j], param().mu(j), param().sigma(j));}
  return sum;
}

/* compute the parameters */
template<class Data_, class WColVector_>
void ModelDiagGaussian_muj_sj<Data_, WColVector_>::computeParameters()
{
  for (int j=p_data()->dataij().beginCols(); j < p_data()->dataij().endCols(); ++j)
  {
    param().mu_[j] = p_data()->dataij().col(j).meanSafe();
    param().sigma_[j] = std::sqrt(p_data()->dataij().col(j).varianceSafe(param().mu(j)));
  }
}
/* compute the weighted parameters */
template<class Data_, class WColVector_>
void ModelDiagGaussian_muj_sj<Data_, WColVector_>::computeParameters( WColVector const& weights)
{
  for (int j=p_data()->dataij().beginCols(); j < p_data()->dataij().endCols(); ++j)
  {
    param().mu_[j] = p_data()->dataij().col(j).wmeanSafe(weights);
    param().sigma_[j] = std::sqrt(p_data()->dataij().col(j).wvarianceSafe(param().mu(j), weights));
  }
}

/* Write the parameters on the output stream os */
template<class Data_, class WColVector_>
void ModelDiagGaussian_muj_sj<Data_, WColVector_>::writeParametersImpl(ostream& os) const
{
  os << _T("mean = ") << mean();
  os << _T("sigma = ")<< sigma();
}

} // namespace STK

#endif /* STK_GAUSSIAN_MUJ_SJMODEL_H */
