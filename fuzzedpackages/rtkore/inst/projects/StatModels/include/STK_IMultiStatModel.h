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
 * Purpose: define the class IMultiStatModel.
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_IMultiStatModel.h
 *  @brief In this file we define the class IMultiStatModel.
 **/

#ifndef STK_IMULTISTATMODEL_H
#define STK_IMULTISTATMODEL_H

#include <cmath>

#include "STK_IStatModelBase.h"

#include <Sdk/include/STK_IRunner.h>
#include <Sdk/include/STK_Macros.h>
#include <DManager/include/STK_DataBridge.h>
#include <STatistiK/include/STK_MultiLaw_IMultiLaw.h>


namespace STK
{

namespace hidden
{
/** Policy trait class for (Stat) Model classes */
template<class Derived> struct StatModelTraits;

} // namespace hidden

/** Parameters container class. All statistical models implemented has an Id defined
 *  in STK::Model::StatisticalModels enumeration.
 **/
template <int Id> struct ModelParameters;

/** @ingroup StatModels
 *  @brief Interface base class for all Multivariate Statistical Models.
 *
 *  A Statistical model, \f$ \mathcal{P}\f$, is a collection of multivariate
 *  probability distribution functions or probability density functions
 *  (collectively referred to as ''distributions'' for brevity).
 *  A parametric model is a collection of distributions, each of which is
 *  indexed by a unique finite-dimensional parameter:
 *  \f$\mathcal{P}=\{\mathbb{P}_{\theta} : \theta \in \Theta\}\f$, where
 *  \f$\theta\f$ is a parameter and \f$\Theta \subseteq \mathbb{R}^d\f$ is
 *  the feasible region of parameters, which is a subset of d-dimensional
 *  Euclidean space.
 *
 *  A statistical model may be used to describe the set of
 *  distributions from which one assumes that a particular data set is sampled.
 *  For example, if one assumes that data arise from a multivariate Gaussian
 *  distribution, then one has assumed a Gaussian model:
 *  \f$
 *    \mathcal{P}=\{\mathbb{P}(x; \mu, \Sigma) = \frac{1}{\sqrt{2 \pi |\Sigma|} }
 *    \exp\left\{ -\frac{1}{2}(x-\mu)'\Sigma^{-1}(x-\mu)\right\} : \mu \in \mathbb{R}^p, \Sigma > 0\}
 *  \f$.
 *
 *  From a computational point of view a statistical model is defined with
 *  the help of two elements
 *  - A data set where the number of samples is the number of rows and the number
 *  of variable is the number of columns. This data set is stored in a Container
 *  of type @c Data.
 *  - A set of parameters stored in a class of type @c Parameters. These parameters
 *  can be created using the @c createParameters method with the effect to
 *  call the default constructor of the Parameters class, or can be set to this
 *  class using the method @c setParameters.
 *
 *  The design for this class is the Curious Recursive Pattern. Derived
 *  implementations of this interface have to implement the following
 *  psudo virtual methods:
 *  @code
 *    int computeNbFreeParameters() const;
 *    Real computeLnLikelihood( RowVector const& rowData) const;
 *    void computeParameters();
 *    void computeParameters(ColVector const& weights);
 *    void writeParameters(ostream& os);
 *  @endcode
 *
 *  @tparam Dervived any kind of multidimensionnal model.
 *
 *  @sa IMixtureDensity.
 **/
template <class Derived>
class IMultiStatModel: public IStatModelBase, public IRecursiveTemplate<Derived>
{
  public:
    /** Type of the container with the data */
    typedef typename hidden::StatModelTraits<Derived>::Data Data;
    /** Type of the data in the container */
    typedef typename hidden::StatModelTraits<Derived>::Type Type;
    /** Type of the vector with the weights */
    typedef typename hidden::StatModelTraits<Derived>::WColVector WColVector;
    /** Type of the parameters of the Model */
    typedef typename hidden::StatModelTraits<Derived>::Parameters Parameters;

  protected:
    /** default constructor. */
    IMultiStatModel(): IStatModelBase(), p_data_(0), param_() {}
    /** Constructor with data set. */
    IMultiStatModel( Data const& data): IStatModelBase(), p_data_(&data), param_(data.dataij().cols())
    { this->initialize(data.dataij().sizeRows(), data.dataij().sizeCols());}
    /** Constructor with a ptr on the data set. */
    IMultiStatModel( Data const* p_data): IStatModelBase(), p_data_(p_data), param_()
    {
      if (p_data)
      {
        this->initialize(p_data->sizeRows(), p_data->sizeCols());
        param_.resize(p_data_->cols());
      }
    }
    /** Copy constructor.
     *  @param model the model to copy
     **/
    IMultiStatModel( IMultiStatModel const& model)
                   : IStatModelBase(model)
                   , p_data_(model.p_data_)
                   , param_(param_)
    {}
    /** destructor */
    ~IMultiStatModel() {}

  public:
    /** @return the pointer on the parameters */
    inline Data const* const p_data() const { return (p_data_);}
    /** @return a reference on the parameters */
    inline Parameters const& param() const { return (param_);}
    /** @return the last error message */
    inline String const& error() const { return msg_error_;}
    /** Set the data set. If the state of the derived runner change when a new
     *  data set is set the user have to overload the udpate() method.
     *  @param p_data A pointer on the data set to run
     **/
    inline void setData( Data const* p_data)
    {
      p_data_ = p_data;
      update();
    }
    /** Estimate the parameters of the model and update the model */
    bool run()
    {
#ifdef STK_DEBUG
      if (!p_data())
      { this->msg_error_ = STKERROR_NO_ARG(IMultiStatModel::run,data have not be set);
        return false;
      }
#endif
      try
      {
        // compute parameters
        this->asDerived().computeParameters();
        // compute log-likelihood
        this->setLnLikelihood(computeLnLikelihood());
        // set the number of free parameters
        this->setNbFreeParameter(this->asDerived().computeNbFreeParameters());
      }
      catch (Exception const& e)
      { this->msg_error_ = e.error(); return false;}
      return true;
    }
    /** compute the weighted empirical probability of success based on the observed
     *  variables. The NA values are discarded.
     *  @param weights the weights of the observations
     **/
    bool run(WColVector const& weights)
    {
#ifdef STK_DEBUG
      if (!p_data())
      { this->msg_error_ = STKERROR_NO_ARG(IMultiStatModel::run(weights),data have not be set);
        return false;
      }
#endif
      try
      {
        // compute weighted parameters
        this->asDerived().computeParameters(weights);
        // compute log-likelihood
        this->setLnLikelihood(computeLnLikelihood());
        // set the number of free parameters
        this->setNbFreeParameter(this->asDerived().computeNbFreeParameters());
      }
      catch (Exception const& e)
      { this->msg_error_ = e.error(); return false;}
      return true;
    }
    /** @param os the output stream for the parameters */
    void writeParameters(ostream &os) { this->asDerived().writeParametersImpl(os);}
    /** default implementation of the writeParameters method.
     *  @param os the output stream for the parameters */
    void writeParametersImpl(ostream &os) // default implementation
    {
#ifdef STK_DEBUG
        stk_cout << _T("You should implement this method in your derived class\n");
#endif
    }

  protected:
    /** Pointer on the parameters of the model. */
    Data const* p_data_;
    /** Pointer on the parameters of the model. */
    Parameters param_;
    /** String with the last error message. */
    String msg_error_;
    /** @return the reference on the parameters */
    inline Parameters& param() { return (param_);}
    /** compute the log Likelihood of the statistical model. */
    Real computeLnLikelihood() const
    {
      Real sum = 0.0;
      for (int i= p_data()->dataij().beginRows(); i< p_data()->dataij().endRows(); i++)
      { sum += this->asDerived().computeLnLikelihood(p_data()->dataij().row(i));}
      return(sum);
    }
    /** update the model if a new data set is set */
    void update()
    {
      if (p_data_)
      {
        this->initialize(p_data()->dataij().sizeRows(), p_data()->dataij().sizeCols());
        param_.resize(p_data()->dataij().cols());
      }
      else
      { this->initialize(0,0); }
    }
};

} // namespace STK

#endif /* STK_IMULTISTATMODEL_H */
