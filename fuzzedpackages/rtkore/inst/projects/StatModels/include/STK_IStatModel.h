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
 * Purpose: define the class IStatModel.
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_IStatModel.h
 *  @brief In this file we define the class IStatModel.
 **/

#ifndef STK_ISTATMODEL_H
#define STK_ISTATMODEL_H

#include <cmath>

#include "STK_IStatModelBase.h"
#include "STK_Model_Util.h"
#include <Sdk/include/STK_Macros.h>
#include <STatistiK/include/STK_MultiLaw_IMultiLaw.h>

namespace STK
{

/** @ingroup StatModels
 *
 *  @brief Base class for all Statistical Models [Deprecated], have been
 *  replaced by IStatModel.
 *
 *  A Statistical model, \f$ \mathcal{P}\f$, is a collection of
 *  probability distribution functions or probability density functions
 *  (collectively referred to as ''distributions'' for brevity).
 *  A parametric model is a collection of distributions, each of which is
 *  indexed by a unique finite-dimensional parameter:
 *  \f$\mathcal{P}=\{\mathbb{P}_{\theta} : \theta \in \Theta\}\f$, where
 *  \f$\theta\f$ is a parameter and \f$\Theta \subseteq \mathbb{R}^d\f$ is
 *  the feasible region of parameters, which is a subset of d-dimensional
 * Euclidean space.  A statistical model may be used to describe the set of
 * distributions from which one assumes that a particular data set is sampled.
 * For example, if one assumes that data arise from a univariate Gaussian
 * distribution, then one has assumed a Gaussian model:
 * \f$\mathcal{P}=\{\mathbb{P}(x; \mu, \sigma) = \frac{1}{\sqrt{2 \pi} \sigma}
 * \exp\left\{ -\frac{1}{2\sigma^2}(x-\mu)^2\right\} : \mu \in \mathbb{R}, \sigma > 0\}
 * \f$.
 *
 *  From a computational point of view a statistical model is defined with
 *  the help of two elements
 *  - A data set where the number of samples is the number of rows
 *    and the number of variables is the number of column. This data set is
 *    stored in a Container of type @c Data.
 *  - A probability (density/law) which for each row of the data set can compute
 *  a density/probability.
 *
 *  @tparam Data can be any kind of container for the data set. it should at
 *  least derive from ITContainer, @sa ITContainer, and provide an access to
 *  a single row like the CArray class or the Array2D class.
 **/
template <class Derived>
class IStatModel: public IStatModelBase
{
  public:
    typedef typename hidden::ModelTraits<Derived>::Data Data;
    typedef typename hidden::ModelTraits<Derived>::ParamHandler ParamHandler;

    /** Type of the data contained in the container */
    typedef typename Data::Type Type;
    /** Type of the row of the data container (a sample) */
    typedef typename hidden::Traits<Data>::Row Row;

  protected:
    /** Constructor with data set. */
    IStatModel(Data const& data): IStatModelBase(data.sizeRows(), data.sizeCols())
                                , p_dataij_(&data)
    {}
    /** Constructor with a ptr on the data set. */
    IStatModel(Data const* p_data): IStatModelBase(), p_dataij_(p_data)
    { if (p_data) this->initialize(p_data->sizeRows(), p_data->sizeCols()) ;}

  public:
    /** destructor */
    ~IStatModel() {}
    /** @return the pointer on the main data set */
    inline Data const* const p_dataij() const { return p_dataij_;}
    /** Set the data set of the model.
     *  @param data the data set of the model
     **/
    inline void setData( Data const& data)
    { p_dataij_ = &data;
      this->setNbSample(p_dataij_->sizeRows());
      this->setNbVariable(p_dataij_->sizeCols()) ;
    }
    /** Set the data set of the model.
     *  @param p_data the data set of the model
     **/
    inline void setData( Data const* p_data)
    { p_dataij_ = p_data;
      if (p_dataij_)
      {
        this->setNbSample(p_dataij_->sizeRows());
        this->setNbVariable(p_dataij_->sizeCols()) ;
      }
    }

  protected:
    /** A pointer on the original data set */
    Data const* p_dataij_;
};

} // namespace STK

#endif /* STK_ISTATMODEL_H */
