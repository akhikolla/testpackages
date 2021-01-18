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
 * created on: 14 nov. 2013
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 * Originally created by Parmeet Bhatia <b..._DOT_p..._AT_gmail_Dot_com>
 **/

/** @file STK_MixtureComposer.h
 *  @brief In this file we define the class MixtureComposer.
 **/

#ifndef STK_MIXTURECOMPOSER_H
#define STK_MIXTURECOMPOSER_H

#include <list>

#include "STK_IMixtureComposer.h"

namespace STK
{
class IMixture;
/** @ingroup Clustering
 *  Main class for handling composed mixture models.
 *  A composed mixture model on some composed space
 *  \f$  \mathbb{X} = \subset \mathcal{X}^1\times \ldots  \times \mathcal{X}^L \f$
 *  is a density of the form
 * \f[
 *     f(\mathbf{x}|\boldsymbol{\theta})
 *     = \sum_{k=1}^K p_k \prod_{l=1}^L f^l(\mathbf{x}^l;\boldsymbol{\lambda}^l_k,\boldsymbol{\alpha}^l)
 *     \quad \mathbf{x} \in \mathbb{X}.
 * \f]
 * The \f$ p_k > 0\f$ with  \f$ \sum_{k=1}^K p_k =1\f$ are the mixing proportions.
 * The density \e f is called the component of the model. The parameters
 * \f$\boldsymbol{\lambda}^l_k, \, k=1,\ldots K \f$ are the cluster specific parameters
 * and the parameters \f$ \boldsymbol{\alpha}^l \f$ are the shared parameters.
 *
 * The MixtureComposer class is a final class implementing the features requested
 * by the interface class IMixtureComposer.
 *
 * It uses injection dependency in order to create/release mixture and
 * get/set parameters of a mixture. The class responsible of the injection is
 * derived from an STK::IMixtureManager.
 * @sa IMixtureComposer, IMixtureManager, PoissonMixtureManager, DiagGaussianMixtureManager, GammaMixtureManager, CategoricalMixtureManager
 **/
class MixtureComposer: public IMixtureComposer
{
  public:
    typedef std::vector<IMixture*>::const_iterator ConstMixtIterator;
    typedef std::vector<IMixture*>::iterator MixtIterator;
    /** Constructor.
     *  @param nbCluster,nbSample number of clusters and samples.
     */
    MixtureComposer( int nbSample, int nbCluster);
    /** copy constructor.
     *  @param composer the composer to copy
     */
    MixtureComposer(MixtureComposer const& composer);
    /** The registered mixtures will be deleted there.*/
    virtual ~MixtureComposer();

    /** Create a composer, but reinitialize the mixtures parameters. */
    virtual MixtureComposer* create() const;
    /** Create a clone of the current model, with mixtures parameters preserved. */
    virtual MixtureComposer* clone() const;
    /** @return the value of the probability of the i-th sample in the k-th component.
     *  @param i,k indexes of the sample and components
     **/
    virtual Real lnComponentProbability(int i, int k) const;
    /** write the parameters of the model in the stream os. */
    virtual void writeParameters(ostream& os) const;

    /** initialize randomly the parameters of the components of the model */
    virtual void randomInit();
    /** Compute the proportions and the model parameters given the current tik
     *  mixture parameters.
     **/
    virtual void paramUpdateStep();
    /** @brief Impute the missing values.
     **/
    virtual void imputationStep();
    /** @brief Simulation of all the latent variables and/or missing data
     *  excluding class labels.
     */
    virtual void samplingStep();
    /**@brief This step can be used to signal to the mixtures that they must
     * store results. This is usually called after a burn-in phase. The composer
     * store the current value of the log-Likelihood.
     **/
    virtual void storeIntermediateResults(int iteration);
    /**@brief This step can be used to signal to the mixtures that they must
     * release the stored results. This is usually called if the estimation
     * process failed.
     **/
    virtual void releaseIntermediateResults();
    /** @brief Utility method allowing to signal to a mixture to set its parameters.
     *  It will be called once enough intermediate results have been stored. */
    virtual void setParametersStep();
    /**@brief This step can be used by developer to finalize any thing. It will
     *  be called only once after we finish running the estimation algorithm.
     **/
    virtual void finalizeStep();

  protected:
    /** @brief Create the composer using existing data handler and mixtures.
     * This method is essentially used by the create() method and can be
     * reused in derived classes.
     * @sa MixtureComposerFixedProp
     **/
    void createComposer( std::vector<IMixture*> const& v_mixtures_);

    /** averaged lnLikelihood values. Will be used by the
     *  storeIntermediateResults method.
     **/
    Real meanlnLikelihood_;
};

/** @brief specialization of the composer for the fixed proportion case.
 **/
class MixtureComposerFixedProp: public MixtureComposer
{
  public:
    /** Constructor.
     * @param nbCluster,nbSample number of clusters and samples
     */
    MixtureComposerFixedProp( int nbSample, int nbCluster);
    /** copy constructor.
     *  @param model the model to copy
     */
    MixtureComposerFixedProp( MixtureComposer const& model);
    /** destructor */
    inline virtual ~MixtureComposerFixedProp() {}
    /** Create a composer, but reinitialize the mixtures parameters. */
    virtual MixtureComposerFixedProp* create() const;
    /** Create a clone of the current model, with mixtures parameters preserved. */
    virtual MixtureComposerFixedProp* clone() const;
    /** overloading of the computePropotions() method.
     * Let them initialized to 1/K. */
    virtual void pStep();
};


} /* namespace STK */

#endif /* STK_MIXTURECOMPOSER_H */
