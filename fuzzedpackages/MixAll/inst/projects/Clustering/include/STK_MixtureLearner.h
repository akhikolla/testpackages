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

/** @file STK_MixtureLearner.h
 *  @brief In this file we define the class MixtureLearner.
 **/

#ifndef STK_MIXTURELEARNER_H
#define STK_MIXTURELEARNER_H

#include <vector>
#include <list>

#include "STK_IMixtureLearner.h"

namespace STK
{
class IMixture;
/** @ingroup Clustering
 *  Main class for learning mixture models. A composed mixture model on some
 *  composed space
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
 * The MixtureLearner class is a final class implementing the features requested
 * by the interface class IMixtureLearner.
 *
 * It uses injection dependency in order to create/release mixture and
 * get/set parameters of a mixture. The class responsible of the injection is
 * derived from an STK::IMixtureManager.
 * @sa IMixtureLearner, IMixtureManager
 **/
class MixtureLearner: public IMixtureLearner
{
  public:
    typedef std::vector<IMixture*>::const_iterator ConstMixtIterator;
    typedef std::vector<IMixture*>::iterator MixtIterator;
    /** Constructor.
     *  @param nbCluster,nbSample number of clusters and samples.
     */
    MixtureLearner( int nbSample, int nbCluster);
    /** copy constructor.
     *  @param composer the composer to copy
     */
    MixtureLearner(MixtureLearner const& composer);
    /** The registered mixtures will be deleted there.*/
    virtual ~MixtureLearner();

    /** Create a composer, but reinitialize the mixtures parameters. */
    virtual MixtureLearner* create() const;
    /** Create a clone of the current model, with mixtures parameters preserved. */
    virtual MixtureLearner* clone() const;
    /** @return the value of the probability of the i-th sample in the k-th component.
     *  @param i,k indexes of the sample and of the class
     **/
    virtual Real lnComponentProbability(int i, int k) const;
    /** write the parameters of the model in the stream os. */
    virtual void writeParameters(ostream& os) const;
    /** initialize randomly the parameters of the components of the model */
    virtual void randomInit();
    /** Compute the model parameters given the current imputed/simulated
     *  missing values. This method updates then the lnLikelihood_ of the model.
     **/
    virtual void paramUpdateStep();
    /** @brief Impute the data missing values. */
    virtual void imputationStep();
    /** @brief Simulation of all the data missing values.
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
     * @sa MixtureLearnerFixedProp
     **/
    void createLearner( std::vector<IMixture*> const& v_mixtures_);

  private:
    /** averaged lnLikelihood values. Will be used by the
     *  storeIntermediateResults method.
     **/
    Real meanlnLikelihood_;
};

} /* namespace STK */

#endif /* STK_MIXTURELEARNER_H */
