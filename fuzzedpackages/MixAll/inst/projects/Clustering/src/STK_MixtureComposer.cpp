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

/** @file STK_MixtureComposer.cpp
 *  @brief In this file we implement the class MixtureComposer.
 **/

#include "../include/STK_MixtureComposer.h"
#include <Arrays/include/STK_Display.h>

namespace STK
{

/* Constructor.
 *  @param nbCluster,nbSample, number of clusters and samples.
 */
MixtureComposer::MixtureComposer( int nbSample, int nbCluster)
                                : IMixtureComposer( nbSample, nbCluster)
                                , meanlnLikelihood_(0.)
{ setNbFreeParameter(nbCluster-1);}

/* copy constructor.
 *  @param composer the composer to copy
 */
MixtureComposer::MixtureComposer( MixtureComposer const& composer)
                                : IMixtureComposer(composer)
                                , meanlnLikelihood_(composer.meanlnLikelihood_)
{}

/* destructor */
MixtureComposer::~MixtureComposer() {}

/* clone pattern */
MixtureComposer* MixtureComposer::clone() const
{ return new MixtureComposer(*this);}

MixtureComposer* MixtureComposer::create() const
{
  // set dimensions
  MixtureComposer* p_composer = new MixtureComposer(nbSample(), nbCluster());
  p_composer->createComposer( v_mixtures_);
  return p_composer;
}

Real MixtureComposer::lnComponentProbability(int i, int k) const
{
  Real sum=0.0;
  for (ConstMixtIterator it = v_mixtures_.begin() ; it != v_mixtures_.end(); ++it)
  { sum += (*it)->lnComponentProbability(i,k);}
  return sum;
}

void MixtureComposer::paramUpdateStep()
{
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("Entering MixtureComposer::paramUpdateStep()\n");
#endif
  for (MixtIterator it = v_mixtures_.begin() ; it != v_mixtures_.end(); ++it)
  { (*it)->paramUpdateStep();}
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("MixtureComposer::paramUpdateStep() done\n");
#endif
}

void MixtureComposer::writeParameters(std::ostream& os) const
{
  os << _T("nbSample = ") << nbSample() << std::endl;
  os << _T("nbCluster = ") << nbCluster() << std::endl;
  os << _T("nbFreeParameter = ") << nbFreeParameter() << std::endl;
  os << _T("nbMissingValues = ") << computeNbMissingValues() << std::endl;
  os << _T("lnLikelihood = ") << lnLikelihood() << std::endl;
  os << _T("v_mixtures_.size() = ") << v_mixtures_.size() << std::endl;
  os << _T("proportions = ") << pk();

  for (ConstMixtIterator it = v_mixtures_.begin(); it != v_mixtures_.end(); ++it)
  {
    os << _T("\nParameters of the mixture: ") << (*it)->idData() << _T("\n");
    (*it)->writeParameters(os);
  }
}

void MixtureComposer::randomInit()
{
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("Entering MixtureComposer::RandomInit()\n");
#endif
  // in case the mixture is not initialized
  if (state() < Clust::modelInitialized_) { initializeStep();}
  // create random tik
  if (randomTik()<2) throw(Clust::randomParamInitFail_);
  // compute zi
  mapStep();
  // set random values to the parameters of the model
  for (MixtIterator it = v_mixtures_.begin(); it != v_mixtures_.end(); ++it)
  { (*it)->randomInit();}
  // perform eStep
  eStep();
  // update state of the model
  setState(Clust::modelParamInitialized_);
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("MixtureComposer::RandomInit() done\n");
#endif
}


void MixtureComposer::imputationStep()
{
  for (MixtIterator it = v_mixtures_.begin(); it != v_mixtures_.end(); ++it)
  { (*it)->imputationStep();}
}

/* @brief Simulation of all the latent variables and/or missing data
 *  excluding class labels.
 */
void MixtureComposer::samplingStep()
{
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("Entering MixtureComposer::samplingStep()\n");
#endif

  for (MixtIterator it = v_mixtures_.begin(); it != v_mixtures_.end(); ++it)
  { (*it)->samplingStep();}
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("MixtureComposer::samplingStep() done\n");
#endif
}

/* store the  intermediate results */
void MixtureComposer::storeIntermediateResults(int iteration)
{
  for (MixtIterator it = v_mixtures_.begin(); it != v_mixtures_.end(); ++it)
  { (*it)->storeIntermediateResults(iteration);}
  meanlnLikelihood_ += (lnLikelihood() - meanlnLikelihood_)/iteration;
}

void MixtureComposer::releaseIntermediateResults()
{
  for (MixtIterator it = v_mixtures_.begin(); it != v_mixtures_.end(); ++it)
  { (*it)->releaseIntermediateResults();}
  meanlnLikelihood_ = 0.;
}

/* Utility method allowing to signal to a mixture to set its parameters */
void MixtureComposer::setParametersStep()
{
  for (MixtIterator it = v_mixtures_.begin(); it != v_mixtures_.end(); ++it)
  { (*it)->setParametersStep();}
  setLnLikelihood(meanlnLikelihood_);
  meanlnLikelihood_ = 0.;
}

/* finalize */
void MixtureComposer::finalizeStep()
{
  for (MixtIterator it = v_mixtures_.begin(); it != v_mixtures_.end(); ++it)
  { (*it)->finalizeStep();}
}

/* @brief Create the composer using existing data handler and ingredients.
 * This method is essentially used by the create() method and can be
 * reused in derived classes. */
void MixtureComposer::createComposer( std::vector<IMixture*> const& v_mixtures)
{
  v_mixtures_.resize( v_mixtures.size());
  for (size_t l = 0; l < v_mixtures_.size(); ++l)
  {
    v_mixtures_[l] = v_mixtures[l]->create();
    v_mixtures_[l]->setMixtureModel(this);
  }
  initializeStep();
}

/* Constructor.
 * @param nbCluster,nbSample number of clusters and samples
 */
MixtureComposerFixedProp::MixtureComposerFixedProp( int nbSample, int nbCluster)
                                                  : MixtureComposer( nbSample, nbCluster)
{ setNbFreeParameter(0); /* remove the count of the pk parameters */}

/* copy constructor.
 *  @param model the model to copy
 */
MixtureComposerFixedProp::MixtureComposerFixedProp( MixtureComposer const& model)
                                                  : MixtureComposer(model) {}

/* Create a composer, but reinitialize the ingredients parameters. */
MixtureComposerFixedProp* MixtureComposerFixedProp::create() const
{
  MixtureComposerFixedProp* p_composer = new MixtureComposerFixedProp(nbSample(), nbCluster());
  p_composer->createComposer(v_mixtures_);
  p_composer->setNbFreeParameter(p_composer->nbFreeParameter()-(nbCluster()-1));
  return p_composer;
}
/* Create a clone of the current model, with ingredients parameters preserved. */
MixtureComposerFixedProp* MixtureComposerFixedProp::clone() const
{ return new MixtureComposerFixedProp(*this);}

/* overloading of the computeProportions() method.
 * Let them initialized to 1/K. */
void MixtureComposerFixedProp::pStep() {}

} /* namespace mixt */
