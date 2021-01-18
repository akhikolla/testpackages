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
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_MixtureSemiLearner.h
 *  @brief In this file we define the MixtureSemiLearner class for
 *  semi-supervised learning models.
 **/

#ifndef STK_MIXTURESEMILEARNER_H
#define STK_MIXTURESEMILEARNER_H

#include "STK_MixtureComposer.h"

namespace STK
{

/** @ingroup Clustering
 *  @brief Base class for Learning a Mixture mixed model when some classes
 *  are known.
 *
 * This class extend the MixtureComposer class and allow to handle the case
 * where some labels are known (semi-supervised learning).
 *
 * @sa MixtureComposer, MixtureLearner
 */
class MixtureSemiLearner: public MixtureComposer
{
  public:
    /** Constructor.
     *  @param nbCluster,nbSample number of clusters and samples
     **/
    MixtureSemiLearner( int nbSample, int nbCluster);
    /** copy constructor.
     *  @param model the model to clone
     **/
    MixtureSemiLearner( MixtureSemiLearner const& model);
    /** destructor */
    virtual ~MixtureSemiLearner();
    /** Create a composer, but reinitialize the mixtures parameters. */
    virtual MixtureSemiLearner* create() const;
    /** Create a clone of the current model, with mixtures parameters preserved. */
    virtual MixtureSemiLearner* clone() const;
    /** @return the vector of the known labels */
    std::vector<int> const& knownLabels() const { return knownLabels_;}
    /** @return the vector of the unknown labels */
    std::vector<int> const& unknownLabels() const { return unknownLabels_;}
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

    /** Set the known labels. This method is expecting a vector of size
     *  nbSample with NA values for unknown labels.
     *  @param zi the vector with the labels
     **/
    template<class RowVector>
    void setLabels( ExprBase<RowVector> const& zi);

  protected:
    /** Create the mixture model parameters tik_ and pk_. */
    virtual void initializeMixtureParameters();
    /** generate random tik_ */
    virtual int randomTik();
    /** generate random zi_ */
    virtual int randomZi();
    /** Copy the ingredient of the semi-learner in the creation process */
    void createSemiLearner(MixtureSemiLearner const* const p_learner);
    /** indexes of the known labels*/
    std::vector<int> knownLabels_;
    /** indexes of the unknown labels*/
    std::vector<int> unknownLabels_;
};


template<class RowVector>
void MixtureSemiLearner::setLabels( ExprBase<RowVector> const& zi)
{
  if (nbSample() != zi.size())
  { STKRUNTIME_ERROR_NO_ARG(MixtureSemiLearner::setKnownLabels,zi must be of size nbSample);}
  int nbUnknown = zi.isNA().count();
#ifdef STK_MIXTURE_DEBUG
  stk_cout << _T("Entering MixtureSemiLearner::setLabels()\n");
  stk_cout << _T("nbUnknown=") << nbUnknown << _T("\n");
#endif
  // clear previous initialization
  unknownLabels_.clear();
  knownLabels_.clear();
  zi_= baseIdx;
  // reserve and set lables
  unknownLabels_.reserve(nbUnknown);
  knownLabels_.reserve(nbSample() - nbUnknown);
  for (int i=zi.begin(); i < zi.end(); ++i)
  {
    if (isNA(zi[i]))
    { unknownLabels_.push_back(i);}
    else
    {
      knownLabels_.push_back(i); zi_[i] = zi[i];
      tik_.row(i) = 0; tik_(i, zi_[i]) = 1.;
    }
  }
#ifdef STK_MIXTURE_DEBUG
  stk_cout << _T("Terminating MixtureSemiLearner::setLabels()\n");
  stk_cout << _T("unknownLabels_.size()=") << unknownLabels_.size() << _T("\n");
  stk_cout << _T("knownLabels_.size()=") << knownLabels_.size() << _T("\n");
#endif
}

/** @brief specialization of the composer for the fixed proportion case.
 **/
class MixtureSemiLearnerFixedProp: public MixtureSemiLearner
{
  public:
    /** Constructor.
     * @param nbCluster,nbSample number of clusters and samples
     */
    MixtureSemiLearnerFixedProp( int nbSample, int nbCluster);
    /** copy constructor.
     *  @param model the model to copy
     */
    MixtureSemiLearnerFixedProp( MixtureSemiLearnerFixedProp const& model);
    /** destructor */
    inline virtual ~MixtureSemiLearnerFixedProp() {}
    /** Create a composer, but reinitialize the mixtures parameters. */
    virtual MixtureSemiLearnerFixedProp* create() const;
    /** Create a clone of the current model, with mixtures parameters preserved. */
    virtual MixtureSemiLearnerFixedProp* clone() const;
    /** overloading of the pStep() method. Let them initialized to 1/K. */
    virtual void pStep();
};


} // namespace STK

#endif /* STK_MIXTURESEMILEARNER_H */

