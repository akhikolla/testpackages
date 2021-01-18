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

/** @file STK_IMixtureComposer.h
 *  @brief In this file we define the abstract base class for mixture models.
 **/

#ifndef STK_IMIXTURECOMPOSER_H
#define STK_IMIXTURECOMPOSER_H

#include "STK_IMixtureStatModel.h"

namespace STK
{

/** @ingroup Clustering
 *  @brief Base class for Mixture (composed) model.
 *
 * In this interface we assume there is an underline generative model that will
 * be estimated using either an EM, SEM or CEM algorithm.
 * All mixture parameters are created using the method
 * @code
 *   void initializeParameters();
 * @endcode
 * in the constructor. They can be accessed from the mixtures using constant
 * accessors.
 *
 * The pure virtual function to implement in derived class are
 * @code
 *   virtual IMixtureComposer* create() const = 0;
 *   virtual IMixtureComposer* clone() const = 0;
 *   virtual void randomInit() = 0;
 *   virtual void paramUpdateStep() = 0;
 *   virtual int computeNbFreeParameters() const = 0;
 * @endcode
 *
 * The public virtual function that can be re-implemented in derived class for
 * a specific behavior are:
 * @code
 *   virtual void initializeStep();
 *   virtual void randomClassInit();
 *   virtual void randomFuzzyInit();
 *   virtual int cStep();
 *   virtual int sStep();
 *   virtual Real eStep();
 *   virtual void mapStep();
 *   virtual void finalizeStep();
 *   virtual void pStep();
 * @endcode
 *
 * @sa IMixture
 *
 * @note the virtual method @c IMixtureComposer::initializeStep is called in all
 * the initialization method.  Don't forget to called it in the randomInit
 * implementation.
 */
class IMixtureComposer: public IMixtureStatModel
{
  protected:
    /** Constructor.
     * @param nbCluster,nbSample number of clusters and samples
     **/
    IMixtureComposer( int nbSample, int nbCluster);
    /** copy constructor.
     *  @param model the model to clone
     **/
    IMixtureComposer( IMixtureComposer const& model);

  public:
    /** destructor */
    virtual ~IMixtureComposer();

    /** @return the state of the model*/
    inline Clust::modelState state() const { return state_;}

    /** set the state of the model : should be used by any strategy*/
    inline void setState(Clust::modelState state) { state_ = state;}

    // pure virtual
    /** create pattern */
    virtual IMixtureComposer* create() const = 0;
    /** clone pattern */
    virtual IMixtureComposer* clone() const = 0;
    /** initialize randomly the parameters of the components of the model */
    virtual void randomInit() = 0;
    /** Compute the proportions and the model parameters given the current tik
     *  mixture parameters.
     **/
    virtual void paramUpdateStep() = 0;

    // virtual with default implementation
    /** Replace tik by zik
     *  @return the minimal value of individuals in a class
     **/
    virtual int cStep();
    /** Simulate zi accordingly to tik and replace tik by zik by calling cStep().
     *  @return the minimal value of individuals in a class
     **/
    virtual int sStep();
    /** compute the zi, the lnLikelihood of the current estimates
     *  and the next value of the tik.
     *  @return the minimal value of tk
     **/
    virtual Real eStep();
    /** Compute zi using the Map estimate. */
    virtual void mapStep();
    /** @brief Finalize the estimation of the model.
     *  Compute obtained lnLikelihood and set state to finalized.
     **/
    virtual void finalizeStep();
    /** Compute proportions using the ML estimates, default implementation. Set
     *  as virtual in case we impose fixed proportions in derived Composer.
     **/
    virtual void pStep();
    /** @brief Initialize the model before its first use.
     *  Initialize the values of the mixture parameters @c pk_ and @© tik_ using
     *  virtual method @c initializeMixtureParameters() and compute @c tk_ and
     *  @c zk_ using the virtual methods @c mapStep() and @c pStep().
     **/
    virtual void initializeStep();

    // not virtual
    /** Initialize randomly the labels zi of the model.
     *  Initialize the model parameters using initializeStep() if it has not
     *  been already called. Simulate the zi, compute tik using eStep(), update
     *  the parameters using paramUpdateStep() and terminate using eStep().
     **/
    void randomClassInit();
    /** Initialize randomly the posterior probabilities tik of the model, then
     *  compute the zi values with mapStep, compute the initial parameter values using
     *  paramUpdateStep, and compute the tik using eStep().
     **/
    void randomFuzzyInit();

  protected:
    /** Simulate zi accordingly to tik.
     *  @param i index of the the individual
     **/
    void sStep(int i);
    /** Replace tik by zik
     *  @param i index of the the individual
     **/
    void cStep(int i);
    /** compute one zi and the next value of the tik for i fixed
     *  @param i the individual
     *  @return the contribution of the individual i to the log-likelihood
     **/
    Real eStep(int i);
    /** Compute zi using the Map estimate for i fixed */
    void mapStep(int i);

    // virtual To be re-implemented if some labels or probabilities are known
    /** Create the mixture model parameters pk_ and tik_.
     *  Default implementation is to set pk_ and tik_ arrays to 1/K value. */
    virtual void initializeMixtureParameters();
    /** generate random tik_ */
    virtual int randomTik();
    /** generate random zi_ */
    virtual int randomZi();

  private:
    /** state of the model*/
    Clust::modelState state_;
    /** Auxiliary array used in the eStep */
#ifndef _OPENMP
    CPointX lnComp_;
#endif
};

} // namespace STK

#endif /* STK_IMIXTURECOMPOSER_H */

