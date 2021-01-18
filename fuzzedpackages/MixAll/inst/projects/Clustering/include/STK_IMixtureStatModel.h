/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2016  Serge Iovleff

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

/** @file STK_IMixtureStatModel.h
 *  @brief In this file we define the abstract base class for mixture statistical models.
 **/

#ifndef STK_IMIXTURESTATMODEL_H
#define STK_IMIXTURESTATMODEL_H

#include <StatModels/include/STK_IStatModelBase.h>

#include "STK_Clust_Util.h"
#include "STK_IMixtureManager.h"
#include <STatistiK/include/STK_Stat_Functors.h>


namespace STK
{

/** @ingroup StatModels
 *  @brief Interface base class for Mixture (composed) model.
 *
 * In statistics, a mixture model is a probabilistic model for representing
 * the presence of sub-populations within an overall population, without
 * requiring that an observed data-set should identify the sub-population to
 * which an individual observation belongs. Formally a mixture model
 * corresponds to the mixture distribution that represents the probability
 * distribution of observations in the overall population. However, while
 * problems associated with "mixture distributions" relate to deriving the
 * properties of the overall population from those of the sub-populations,
 * "mixture models" are used to make statistical inferences about the
 * properties of the sub-populations given only observations on the pooled
 * population, without sub-population-identity information.
 *
 * Some ways of implementing mixture models involve steps that attribute
 * postulated sub-population-identities to individual observations (or weights
 * towards such sub-populations), in which case these can be regarded as types
 * unsupervised learning or clustering procedures. However not all inference
 * procedures involve such steps.
 *
 * Pure virtual functions to implement in derived class are
 * @code
 *   virtual IMixtureComposer* create() const = 0;
 *   virtual IMixtureComposer* clone() const = 0;
 *   virtual bool randomInit() =0;
 *   virtual Real lnComponentProbability(int i, int k) const = 0;
 * @endcode
 *
 * Virtual functions that can be re-implemented in derived class for a
 * specific behavior are:
 * @code
 *   virtual void imputationStep();
 *   virtual void samplingStep();
 *   virtual void setParametersStep();
 *   virtual void storeIntermediateResults(int iteration);
 *   virtual void releaseIntermediateResults()
 *   virtual void writeParameters(std::ostream& os) const;
 * @endcode
 *
 * Template functions allowing to interact with the composer are
 * @code
 *   template<class Array>
 *   void setMixtureParameters( Array const& tik);
 *   template<class Array, class RowVector>
 *   void setMixtureParameters( Array const& tik, RowVector const& pk);
 *   template<class RowVector>
 *   void setProportions( RowVector const& pk);
 *
 *   template<class Manager, class Parameters>
 *   void setParameters(IMixtureManager<Manager> const& manager, String const& idData, Parameters const& param);
 *
 *   template<class Manager>
 *   void createMixture(IMixtureManager<Manager>& manager);
 *   template<class DataHandler>
 *   IMixture* createMixture(IMixtureManager<Manager>& manager, String const& idData);
 *   template<class DataHandler>
 *   void removeMixture(IMixtureManager<Manager>& manager, String const& idData);
 * @endcode
 *
 * Template functions allowing to get results from the composer are
 * @code
 *   template<class Manager, class Parameters>
 *   void getParameters(IMixtureManager<Manager> const& manager, String const& idData, Parameters& param) const;
 *   template<class Manager, class MissingValues>
 *   void getMissingValues(IMixtureManager<Manager> const& manager, String const& idData, MissingValues& missing) const;
 * @endcode
 *
 * @sa IMixtureComposer, IMixtureLearner, IMixtureManager
 *
 * @note the virtual method @c IMixtureStatModel::initializeStep is called in all
 * the initialization method. Don't forget to called it in the randomInit
 * implementation.
 */
class IMixtureStatModel: public IStatModelBase
{
  protected:
    /** Constructor.
     * @param nbCluster,nbSample number of clusters and samples
     **/
    IMixtureStatModel( int nbSample, int nbCluster);
    /** copy constructor.
     *  @param model the model to clone
     **/
    IMixtureStatModel( IMixtureStatModel const& model);

  public:
    typedef std::vector<IMixture*>::const_iterator ConstMixtIterator;
    typedef std::vector<IMixture*>::iterator MixtIterator;

    /** destructor */
    virtual ~IMixtureStatModel();

    /** @return the number of cluster */
    inline int nbCluster() const { return nbCluster_;}
    /** @return the proportions of each mixtures */
    inline CPointX const& pk() const { return pk_;};
    /** @return the tik probabilities */
    inline CArrayXX const& tik() const { return tik_;};
    /** @return the proportions of individuals */
    inline CPointX const& tk() const { return tk_;};
    /** @return the zi class label */
    inline CVectorXi const& zi() const { return zi_;};
    /** @return a constant reference on the vector of mixture */
    inline std::vector<IMixture*> const& v_mixtures() const { return v_mixtures_;}

    /** @return the computed log-likelihood of the i-th sample.
     *  @param i index of the sample
     **/
    Real computeLnLikelihood(int i) const;
    /** @return the computed likelihood of the i-th sample.
     *  @param i index of the sample
     **/
    Real computeLikelihood(int i) const;
    /** @return the computed log-likelihood. */
    Real computeLnLikelihood() const;
    /** @return the computed ICL criteria. */
    Real computeICL() const;

    /** Utility lookup function allowing to find a Mixture from its idData
     *  @param idData the id name of the data
     *  @return a pointer on the mixture, NULL if the mixture is not found
     **/
     IMixture* getMixture( String const& idData) const;
     /** @brief register a mixture to the composer.
      *  When a mixture is registered, the composer:
      *  - assign composer pointer (itself) to the mixture
      *  - add it to v_mixtures_
      *  - update the number of variables
      * @param p_mixture a pointer on the mixture
      **/
     void registerMixture(IMixture* p_mixture);
     /** @brief release a mixture from the composer.
      *  When a mixture is released, the composer remove it from v_mixtures_.
      *  @note the data set associated with the mixture is still in the manager
      *  list of data set if you use a manager in order to register it.
      *  @param idData the Id of the mixture to release.
      **/
     void releaseMixture(String const& idData);
     /** compute the number of free parameters of the model.
      *  This method is used in IMixtureComposer::initializeStep
      *  in order to give a value to IStatModelBase::nbFreeParameter_.
      *  lookup on the mixtures and sum the nbFreeParameter.
      *  @return the number of free parameters
      **/
     int computeNbFreeParameters() const;
     /** @brief compute the missing values of the model.
      *  lookup on the mixtures and sum the nbMissingValues.
      *  @return the number of missing values
      **/
     int computeNbMissingValues() const;

    // pure virtual
    /** create pattern */
    virtual IMixtureStatModel* create() const = 0;
    /** clone pattern */
    virtual IMixtureStatModel* clone() const = 0;
    /** initialize randomly the parameters of the components of the model */
    virtual void randomInit() = 0;
    /** @return the value of the probability of the ith sample in the kth class.
     *  @param i,k indexes of the sample and of the class
     **/
    virtual Real lnComponentProbability(int i, int k) const = 0;

    // virtual with default implementation
    /** @brief Initialize the model before at its first use.
     *  This function can be overloaded in derived class for initialization of
     *  the specific model parameters. It should be called prior to any used of
     *  the class.
     *  @sa IMixture,MixtureBridge,MixtureLearner
     **/
    virtual void initializeStep();

    /** @brief Impute the missing values.
     *  Default behavior is "do nothing".
     **/
    inline virtual void imputationStep() {}
    /** @brief Simulation of all the latent variables and/or missing data
     *  excluding class labels. Default behavior is "do nothing".
     */
    inline virtual void samplingStep() {}
    /** @brief Utility method allowing to signal to a mixture to set its parameters.
     *  It will be called once enough intermediate results have been stored.
     **/
    inline virtual void setParametersStep() {}
    /**@brief This step can be used to signal to the mixtures that they must
     * store results. This is usually called after a burn-in phase. The composer
     * store the current value of the log-Likelihood.
     **/
    inline virtual void storeIntermediateResults(int iteration) {}
    /**@brief This step can be used to signal to the mixtures that they must
     * release the stored results. This is usually called if the estimation
     * process failed.
     **/
    inline virtual void releaseIntermediateResults() {}
    /** @brief Finalize the estimation of the model.
     *  The default behavior is compute current lnLikelihood.
     **/
    inline virtual void finalizeStep() { setLnLikelihood(computeLnLikelihood());}
    /** write the parameters of the model in the stream os. */
    virtual void writeParameters(ostream& os) const {};

    // initialization of the arrays pk_, tik_, tk_ and zi_
    /** set the mixture parameters using an array of posterior probabilities.
     *  Proportions, numbers in each class and class labels are computed
     *  using these posterior probabilities.
     *  @param tik posterior class probabilities
     **/
    template<class Array>
    void setMixtureParameters( Array const& tik);
    /** set the mixture parameters giving the posterior probabilities and
     *  the proportions.
     *  Numbers in each class and class labels are computed using the
     *  posterior probabilities.
     *  @param tik posterior class probabilities
     *  @param pk prior class proportion
     **/
    template<class Array, class RowVector>
    void setMixtureParameters( Array const& tik, RowVector const& pk);
    /** Set proportions of each classes
     *  @param pk prior class proportion
     **/
    template<class RowVector>
    void setProportions( RowVector const& pk);

    // allow to update specific mixture parameters values
    /** Utility method allowing to set the parameters to a specific mixture.
     *  @param manager the manager with the responsibility of the parameters
     *  @param idData Id of the data we want to set the parameters
     *  @param param structure which contains the parameters
     **/
    template<class Manager, class Parameters>
    void setParameters(IMixtureManager<Manager> const& manager, String const& idData, Parameters const& param);

    // utilities methods
    /** Utility method allowing to create all the mixtures registered in the
     *  data handler of a mixture manager and to register them.
     *  @param manager the manager with the responsibility of the creation.
     **/
    template<class Manager>
    void createMixture(IMixtureManager<Manager>& manager);
    /** Utility method allowing to create a mixture with a given data set
     *  and register it. The Mixture Manager will find the associated model
     *  to use with this data set.
     *  @param manager the manager with the responsibility of the creation.
     *  @param idData the id name of the data to modelize.
     **/
    template<class Manager>
    IMixture* createMixture(IMixtureManager<Manager>& manager, String const& idData);
    /** Utility method allowing to release completely a mixture with its data set.
     *  The MixtureManager will find and release the associated data set.
     *  @param manager the manager with the responsibility of the release.
     *  @param idData the id name of the data to modelize.
     **/
    template<class Manager>
    void removeMixture(IMixtureManager<Manager>& manager, String const& idData);

    /** Utility method allowing to get the parameters of a specific mixture.
     *  @param manager the manager with the responsibility of the parameters
     *  @param idData the Id of the data we want the parameters
     *  @param param the structure which will receive the parameters
     **/
    template<class Manager, class Parameters>
    void getParameters(IMixtureManager<Manager> const& manager, String const& idData, Parameters& param) const;

    /** Utility method allowing to get the missing values of a specific mixture.
     *  @param manager the manager with the responsibility of the parameters
     *  @param idData the Id of the data we want the parameters
     *  @param missing the structure which will receive the missing values
     **/
    template<class Manager, class MissingValues>
    void getMissingValues(IMixtureManager<Manager> const& manager, String const& idData, MissingValues& missing) const;

  protected:
    /** set the number of cluster of the model
     *  @param nbCluster number of cluster of the model
     * */
    inline void setNbCluster( int nbCluster) { nbCluster_ = nbCluster;}

    /** number of cluster. */
    int nbCluster_;
    /** The proportions of each mixtures */
    CPointX pk_;
    /** The tik probabilities */
    CArrayXX tik_;
    /** The sum of the columns of tik_ */
    CPointX tk_;
    /** The zi class label */
    CVectorXi zi_;

    /** vector of pointers to the mixtures components */
    std::vector<IMixture*> v_mixtures_;
};


/* Utility method allowing to create all the mixtures handled by a mixture manager.
 *  @param manager the manager with the responsibility of the creation.
 **/
template<class DataHandler>
void IMixtureStatModel::createMixture(IMixtureManager<DataHandler>& manager)
{
  typedef typename DataHandlerBase<DataHandler>::InfoMap InfoMap;
  typedef typename InfoMap::const_iterator const_iterator;
  for ( const_iterator it=manager.p_handler()->info().begin(); it!=manager.p_handler()->info().end(); ++it)
  {
    IMixture* p_mixture = manager.createMixture(it->first, nbCluster());
#ifdef STK_MIXTURE_DEBUG
  if (!p_mixture)
  { stk_cout << _T("In IMixtureStatModel::createMixture(manager) failed.\n");}
#endif
    if (p_mixture) registerMixture(p_mixture);
  }
}

/* Utility method allowing to create a mixture with a given data set
 *  and register it. The Mixture Manager will find the associated model
 *  to use with this data set.
 *  @param manager the manager with the responsibility of the creation.
 *  @param idData the id name of the data to modelize.
 **/
template<class Manager>
IMixture* IMixtureStatModel::createMixture(IMixtureManager<Manager>& manager, String const& idData)
{
  IMixture* p_mixture = manager.createMixture( idData, nbCluster());
#ifdef STK_MIXTURE_DEBUG
  if (!p_mixture)
  { stk_cout << _T("In IMixtureStatModel::createMixture(manager,")<< idData << _T(") failed.\n");}
#endif
  if (p_mixture) registerMixture(p_mixture);
  return p_mixture;
}

/* Utility method allowing to release completely a mixture with its data set.
 *  The MixtureManager will find and release the associated data set.
 *  @param manager the manager with the responsibility of the release.
 *  @param idData the id name of the data to modelize.
 **/
template<class Manager>
void IMixtureStatModel::removeMixture(IMixtureManager<Manager>& manager, String const& idData)
{
  IMixture* p_mixture = getMixture(idData);
#ifdef STK_MIXTURE_DEBUG
  if (!p_mixture)
  { stk_cout << _T("In IMixtureStatModel::removeMixture(manager,")<< idData << _T(") failed.\n");}
#endif
  if (p_mixture)
  {
    releaseMixture(idData);
    manager.releaseDataBridge( idData);
  }
}

/* Utility method allowing to set the parameters to a specific mixture.
 *  @param manager the manager with the responsibility of the parameters
 *  @param idData Id of the data we want to set the parameters
 *  @param param structure which contains the parameters
 **/
template<class Manager, class Parameters>
void IMixtureStatModel::setParameters(IMixtureManager<Manager> const& manager, String const& idData, Parameters const& param)
{ IMixture* p_mixture= getMixture(idData);
  if (p_mixture) manager.setParameters(p_mixture, param);
}
/* set the mixture parameters using the posterior probabilities.
 **/
template<class Array>
void IMixtureStatModel::setMixtureParameters(Array const& tik)
{
  tik_ = tik;
  tk_ = Stat::sumByCol(tik_);
  pk_ = tk_ / nbSample();
  for (int i = tik_.beginCols(); i< tik_.endCols(); ++i)
  {
    int k;
    tik_.row(i).maxElt(k);
    zi_[i] = k;
  }
}

/* set the mixture parameters using the posterior probabilities.
 **/
template<class Array, class RowVector>
void IMixtureStatModel::setMixtureParameters(Array const& tik, RowVector const& pk)
{
#ifdef STK_MIXTURE_DEBUG
  if (pk_.size() != tik_.sizeCols())
  { STKRUNTIME_ERROR_2ARG(IMixtureLearner::setMixtureParameters,pk_.size(),tik_.sizeCols(),Numbers of class in tik and pk differ);}
#endif
  tik_ = tik;
  pk_  = pk;
  tk_  = Stat::sumByCol(tik_);
  for (int i = tik_.beginCols(); i< tik_.endCols(); ++i)
  {
    int k;
    tik_.row(i).maxElt(k);
    zi_[i] = k;
  }
}

/* set proportions */
template<class RowVector>
void IMixtureStatModel::setProportions( RowVector const& pk)
{
#ifdef STK_MIXTURE_DEBUG
  if (pk_.size() != nbCluster())
  { STKRUNTIME_ERROR_2ARG(IMixtureLearner::setProportions,pk_.size(),nbCluster(),Numbers of class in pk differs);}
#endif
  pk_  = pk;
}

/* Utility method allowing to get the parameters of a specific mixture.
 *  @param manager the manager with the responsibility of the parameters
 *  @param idData the Id of the data we want the parameters
 *  @param param the structure which will receive the parameters
 **/
template<class Manager, class Parameters>
void IMixtureStatModel::getParameters(IMixtureManager<Manager> const& manager, String const& idData, Parameters& param) const
{  IMixture* p_mixture= getMixture(idData);
   if (p_mixture) manager.getParameters(p_mixture, param);
}
/* Utility method allowing to get the missing values of a specific mixture.
 *  @param manager the manager with the responsibility of the parameters
 *  @param idData the Id of the data we want the parameters
 *  @param missing the structure which will receive the missing values
 **/
template<class Manager, class MissingValues>
void IMixtureStatModel::getMissingValues(IMixtureManager<Manager> const& manager, String const& idData, MissingValues& missing) const
{
  IMixture* p_mixture= getMixture(idData);
  if (p_mixture) manager.getMissingValues(p_mixture, missing);
}



} // namespace STK

#endif /* STK_IMIXTURESTATMODEL_H */


