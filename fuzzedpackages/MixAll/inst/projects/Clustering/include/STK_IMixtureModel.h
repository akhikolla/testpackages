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

/** @file STK_IMixtureModel.h
 *  @brief In this file we define the main interface class for mixture models.
 **/


#ifndef STK_IMIXTUREMODEL_H
#define STK_IMIXTUREMODEL_H

#include "STK_IMixtureModelBase.h"
#include <StatModels/include/STK_Model_Util.h>
#include <STatistiK/include/STK_Law_Categorical.h>

#ifdef STK_MIXTURE_DEBUG
#include "Arrays/include/STK_Display.h"
#endif

namespace STK
{
namespace Clust
{
/** Main class for the mixtures traits policy.
 *  The traits struct MixtureTraits must be specialized for any
 *  Mixture deriving from the Interface IMixtureModel.
 **/
template <class Mixture> struct MixtureTraits;

} // namespace Clust

/** Parameters class. All statistical models has an unique Id defined
 *  in STK::Model::StatisticalModels enumeration.
 **/
template <int Id> struct ModelParameters;

/**@ingroup Clustering
 * @brief Main interface class for mixture models.
 * At this level we create the Parameters struct. The call to @e setData
 * trigger the call to @e initializeModel which itself trigger a call to
 * @e initializeModelImpl to the derived class.
 *
 * The methods to implement in derived classes are
 * @code
 *   template<class Weights>
 *   Type impute(int i, int j, Weights const& pk) const;
 *   Type rand(int i, int j, int k) const;
 *   Real lnComponentProbability(int i, int k) const;
 *   void randomInit( CArrayXX const* const& p_tik, CPointX const* const& p_tk) ;
 *   bool run( CArrayXX const* const& p_tik, CPointX const* const& p_tk) ;
 * int computeNbFreeParameters() const;
 * @endcode
 *
 * The pseudo virtual methods to implement if needed in derived classes are
 * @code
 *   // default implementation (do nothing) provided to all these methods
 *   void initializeModelImpl();
 *   bool initializeStepImpl(); // return true by default
 *   void finalizeStepImpl();
 * @endcode
 *
 * @sa IMixtureModelBase, IRecursiveTemplate
 **/
template<class Derived>
class IMixtureModel: public IRecursiveTemplate<Derived>, public IMixtureModelBase
{
  public:
    typedef typename Clust::MixtureTraits<Derived>::Array Array;
    typedef typename Clust::MixtureTraits<Derived>::Parameters Parameters;
    typedef typename Array::Type Type;

  protected:
    /** Default constructor.
     *  @param nbCluster the number of cluster
     **/
    inline IMixtureModel( int nbCluster)
                        : IMixtureModelBase(nbCluster)
                        , param_(nbCluster)
                        , p_dataij_(0)
    {}
    /** copy constructor.
     *  - The Parameter class is copied using the copy constructor.
     *  - The pointer on the data set is copied as-is. Check if you should not
     *  change it on the copied object.
     *  @param model the model to copy
     **/
    IMixtureModel( IMixtureModel const& model)
                 : IMixtureModelBase(model)
                 , param_(model.param_)
                 , p_dataij_(model.p_dataij_)
    {}

  public:
    /** destructor */
    inline ~IMixtureModel() {}
    /** create pattern.  */
    inline IMixtureModel* create() const { return new Derived(this->nbCluster());}
    /** @return a pointer on the current data set */
    inline Array const* p_data() const { return p_dataij_;}

    /** @brief Set the data set.
     *  Setting a (new) data set will trigger the initialization process of the model.
     *  @param data the data set to set
     **/
    inline void setData(Array const& data)
    {
      p_dataij_ = &data;
      initializeModel();
    }
    /** @brief This function will be called once the model is created and data is set.
     *  @note a stk++ mixture create and initialize all the containers when the data
     *  is set. Thus the default behavior is @c return true.
     */
    inline bool initializeStep() { return this->asDerived().initializeStepImpl();}
    /** set the parameters stored in stat_proba_ and release stat_proba_. */
    inline void setParametersStep()
    {
      param_.setParametersStep();
      this->asDerived().setParametersImpl();
    }
    /** @brief This function will be called once the model is estimated.
     *  perform specific model finalization stuff */
    inline void finalizeStep() { this->asDerived().finalizeStepImpl();}

    // default implementation of the pseudo-virtual methods
    /** default implementation of initializeModelImpl (do nothing) */
    inline void initializeModelImpl() {}
    /** default implementation of initializeStepImpl (return true) */
    inline bool initializeStepImpl() { return true;}
    /** default implementation of finalizeStepImpl (do nothing) */
    inline void finalizeStepImpl() {}

    /** @return a simulated value for the jth variable of the ith sample
     *  @param i,j indexes of the data to impute
     *  @param tk the probabilities of each class for the ith individual
     **/
    template<class Weights>
    Type sample(int i, int j, Weights const& tk) const
    { return this->asDerived().rand(i, j, Law::Categorical::rand(tk));}

  protected:
    /** @brief Initialize the model before its first use.
     * This function is triggered when data set is set.
     * In this interface, the @c initializeModel() method
     *  - set the number of samples and variables of the mixture model
     *  - call the derived class implemented method
     * @code
     *   initializeModelImpl()
     * @endcode
     * for initialization of the specific model parameters if needed.
     **/
    void initializeModel()
    {
      // set dimensions
      this->setNbSample(p_dataij_->sizeRows());
      this->setNbVariable(p_dataij_->sizeCols());
      // call specific model initialization stuff
      this->asDerived().initializeModelImpl();
    }
    /** parameters of the derived mixture model */
    Parameters param_;

  private:
    /** pointer on the data set */
    Array const* p_dataij_;
};

} // namespace STK

#endif /* STK_IMIXTUREMODEL_H */
