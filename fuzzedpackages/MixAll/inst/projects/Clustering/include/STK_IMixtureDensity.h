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
 * Project:  stkpp::Clustering
 * created on: 16 oct. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 * Originally created by Parmeet Bhatia <b..._DOT_p..._AT_gmail_Dot_com>
 **/

/** @file STK_IMixtureDensity.h
 *  @brief In this file we define the main interface class for mixture models.
 **/


#ifndef STK_IMIXTUREDENSITY_H
#define STK_IMIXTUREDENSITY_H

#include "STK_Clust_Util.h"

#include <StatModels/include/STK_Model_Util.h>
#include <STatistiK/include/STK_Law_Categorical.h>

#ifdef STK_MIXTURE_DEBUG
#include "Arrays/include/STK_Display.h"
#endif

namespace STK
{
/**@ingroup Clustering
 *  @brief Base class for all Mixture densities.
 *
 *  Let @e X be an arbitrary measurable space and let
 *  \f$ \mathbf{x} =\{{\mathbf x}_1,...,{\mathbf x}_n\}\f$
 *  be @e n independent vectors in @e X. If each \f${\mathbf x}_i\f$
 *  arises from a probability distribution with density
 *  \f[
 *    f({\mathbf x}_i|\theta) = \sum_{k=1}^K p_k h({\mathbf x}_{i}| \lambda_{k},\alpha)
 *  \f]
 *  where the \f$p_k\f$'s are the mixing proportions, \f$h(\cdot| \lambda_{k},\alpha)\f$
 *  denotes a distribution parameterized by \f$\lambda_k\f$ and \f$\alpha\f$, it
 *  is said that we observe a mixture model.
 *
 *  This interface class is the base class for all distributions part of a
 *  mixture model.
 *
 *  @sa IMixture, IMixtureBridge, MixtureComposer
 *
 * At this level the Parameters struct is created.
 *
 * A call to @e setData trigger a call to @e initializeModel which itself
 * trigger a call to @c initializeModelImpl. This last method can be
 * re-implemented by derived class.
 *
 * The pseudo virtual methods (needed by IMixtureBridge) to implement in concrete
 * derived classes are
 * @code
 *   // return an imputed value for (i,j). In case of
 *   template<class Weights>
 *   Type impute(int i, int j, Weights const& pk) const;
 *   // return a random value for (i,j)
 *   Type rand(int i, int j, int k) const;
 *   Real lnComponentProbability(int i, int k) const;
 *   void randomInit( CArrayXX const* const& p_tik, CPointX const* const& p_tk) ;
 *   bool run( CArrayXX const* const& p_tik, CPointX const* const& p_tk) ;
 *   int computeNbFreeParameters() const;
 *   void initializeModelImpl();
 * @endcode
 *
 * The pseudo virtual methods to re-implement if needed in derived classes are
 * @code
 *   // default implementation "do nothing" provided to all these methods
 *   bool initializeStepImpl(); // return true by default
 *   void finalizeStepImpl();
 * @endcode
 **/
template<class Derived>
class IMixtureDensity: public IRecursiveTemplate<Derived>
{
  public:
    typedef typename hidden::MixtureTraits<Derived>::Array Array;
    typedef typename hidden::MixtureTraits<Derived>::Parameters Parameters;
    typedef typename hidden::Traits<Array>::Type Type;

  protected:
    /** Default constructor.
     *  @param nbCluster the number of cluster
     **/
    IMixtureDensity( int nbCluster);
    /** copy constructor.
     *  - The Parameter class is copied using the copy constructor.
     *  - The pointer on the data set is copied as-is. Check if you should not
     *  change it on the copied object.
     *  @param model the model to copy
     **/
    IMixtureDensity( IMixtureDensity const& model);

  public:
    /** destructor */
    ~IMixtureDensity() {}

    /** @return the number of cluster */
    inline int nbCluster() const { return nbCluster_;}
    /** @return the total available observations */
    inline int nbSample() const { return nbSample_;}
    /** @return the Log of the total available observations */
    inline Real lnNbSample() const
    { return (nbSample_ <= 0) ? -Arithmetic<Real>::infinity() : std::log((Real)nbSample_);}
    /** @return a pointer on the current data set */
    inline Array const* const& p_data() const { return p_dataij_;}
    /** @return a constant reference on the current parameter struct */
    inline Parameters const& param() const {return param_;}

    /** @brief Set the data set.
     *  Setting a (new) data set will trigger the initializeModel() method.
     *  @param data the data set to set
     **/
    void setData(Array const& data);
    /** @brief Set the data set and give dimensions.
     *  Setting a (new) data set will trigger the initializeModel() method.
     *  @param data the data set to set
     *  @param nbRow number of rows of the data set
     *  @param nbCol number of columns of the data set
     *  @param byRow is data vectorized by row or by column ?
     *  @note this method is used for matrix valued data sets where data are
     *  vectorized
     **/
    void setData(Array const& data, int nbRow, int nbCol, bool byRow = true);
    /** @brief This function will be called at the beginning of the estimation
     *  process once the model is created and data is set.
     *  @note a stk++ mixture create and initialize all the containers when the
     *  data is set. Thus the default behavior is @c return true.
     */
    inline bool initializeStep() { return this->asDerived().initializeStepImpl();}
    /** set the parameters obtained with the intermediate results and release
     *  the intermediate results. */
    inline void setParametersStep()
    {
      param_.setParametersStep();
      this->asDerived().setParametersImpl();
    }
    /** @brief This function will be called once the model is estimated.
     *  perform specific model finalization stuff */
    inline void finalizeStep() { this->asDerived().finalizeStepImpl();}

    /** @return a simulated value for the jth variable of the ith sample
     *  @param i,j indexes of the data to impute
     *  @param tk the probabilities of each class for the ith individual
     **/
    template<class Weights>
    inline Type sample(int i, int j, Weights const& tk) const
    { return this->asDerived().rand(i, j, Law::Categorical::rand(tk));}

  protected:
    /** parameters of the derived mixture model. Should be an instance
     *  of the STK::ModelParameters struct.
     *  @sa STK::ModelParameters
     **/
    Parameters param_;

    /** @brief Initialize the model before its first use.
     * This function is triggered when data set is set.
     * In this interface, the @c initializeModel() method
     *  - set the number of samples of the mixture model
     *  - call the derived class implemented method
     * @code
     *   initializeModelImpl()
     * @endcode
     * for initialization of the specific model parameters if needed.
     **/
    void initializeModel();
    // default implementation of the pseudo-virtual methods
    /** default implementation of initializeStepImpl (do nothing and return true) */
    inline bool initializeStepImpl() { return true;/* do nothing*/}
    /** default implementation of finalizeStepImpl (do nothing) */
    inline void finalizeStepImpl() {/* do nothing*/}

    /** Set the number of sample of the model (needed by kernel models)
     *  @param nbSample number of sample of the model
     *  @note to use when data set is not available but only the dimensions
     * */
    inline void setNbSample( int nbSample) { nbSample_ = nbSample;}

  private:
    /** number of cluster. */
    int nbCluster_;
    /** total available samples */
    int nbSample_;
    /** pointer on the data set */
    Array const* p_dataij_;
};


/* Default constructor.
 *  @param nbCluster the number of cluster
 **/
template<class Derived>
IMixtureDensity<Derived>::IMixtureDensity( int nbCluster)
                                         : param_(nbCluster)
                                         , nbCluster_(nbCluster)
                                         , nbSample_(0)
                                         , p_dataij_(0)
{}
/* copy constructor.
 *  - The Parameter class is copied using the copy constructor.
 *  - The pointer on the data set is copied as-is. Check if you should not
 *  change it on the copied object.
 *  @param model the model to copy
 **/
template<class Derived>
IMixtureDensity<Derived>::IMixtureDensity( IMixtureDensity const& model)
                                         : param_(model.param_)
                                         , nbCluster_(model.nbCluster_)
                                         , nbSample_(model.nbSample_)
                                         , p_dataij_(model.p_dataij_)
{}

/* @brief Set the data set.
 *  Setting a (new) data set will trigger the initializeModel() method.
 *  @param data the data set to set
 **/
template<class Derived>
void IMixtureDensity<Derived>::setData(Array const& data)
{
  p_dataij_ = &data;
  initializeModel();
}
/* @brief Set the data set.
 *  Setting a (new) data set will trigger the initializeModel() method.
 *  @param data the data set to set
 **/
template<class Derived>
void IMixtureDensity<Derived>::setData(Array const& data, int nbRow, int nbCol, bool byRow)
{
  p_dataij_ = &data;
  initializeModel();
}


/* @brief Initialize the model before its first use.
 * This function is triggered when data set is set.
 * In this interface, the @c initializeModel() method
 *  - set the number of samples of the mixture model
 *  - call the derived class implemented method
 * @code
 *   initializeModelImpl()
 * @endcode
 * for initialization of the specific model parameters if needed.
 **/
template<class Derived>
void IMixtureDensity<Derived>::initializeModel()
{
  // set dimensions
  this->setNbSample(p_dataij_->sizeRows());
  // call specific model initialization stuff
  this->asDerived().initializeModelImpl();
}

} // namespace STK

#endif /* STK_IMIXTUREDENSITY_H */
