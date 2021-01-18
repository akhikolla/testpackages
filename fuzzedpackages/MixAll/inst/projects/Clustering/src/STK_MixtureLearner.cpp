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

/** @file STK_MixtureLearner.cpp
 *  @brief In this file we implement the class MixtureLearner.
 **/

#include "../include/STK_MixtureLearner.h"
#include <Arrays/include/STK_Display.h>

namespace STK
{

/* Constructor.
 *  @param nbCluster,nbSample, number of clusters and samples.
 */
MixtureLearner::MixtureLearner( int nbSample, int nbCluster)
                              : IMixtureLearner( nbSample, nbCluster)
                              , meanlnLikelihood_(0.)
{ setNbFreeParameter(nbCluster-1);}

/* copy constructor.
 *  @param composer the composer to copy
 */
MixtureLearner::MixtureLearner( MixtureLearner const& learner)
                              : IMixtureLearner(learner)
                              , meanlnLikelihood_(learner.meanlnLikelihood_)
{}

MixtureLearner::~MixtureLearner() {}

/* clone pattern */
MixtureLearner* MixtureLearner::clone() const
{ return new MixtureLearner(*this);}

MixtureLearner* MixtureLearner::create() const
{
  // set dimensions
  MixtureLearner* p_composer = new MixtureLearner(nbSample(), nbCluster());
  p_composer->createLearner( v_mixtures_);
  return p_composer;
}

Real MixtureLearner::lnComponentProbability(int i, int k) const
{
  Real sum=0.0;
  for (ConstMixtIterator it = v_mixtures_.begin() ; it != v_mixtures_.end(); ++it)
  { sum += (*it)->lnComponentProbability(i,k);}
  return sum;
}

void MixtureLearner::paramUpdateStep()
{
  for (MixtIterator it = v_mixtures_.begin() ; it != v_mixtures_.end(); ++it)
  { (*it)->paramUpdateStep();}
  setLnLikelihood(computeLnLikelihood());
}

void MixtureLearner::writeParameters(std::ostream& os) const
{
  os << _T("Learner nbSample = ") << nbSample() << std::endl;
  os << _T("Learner nbCluster = ") << nbCluster() << std::endl;
  os << _T("Learner nbFreeParameter = ") << nbFreeParameter() << std::endl;
  os << _T("Learner lnLikelihood = ") << lnLikelihood() << std::endl;
  os << _T("Learner proportions = ") << pk();

  for (ConstMixtIterator it = v_mixtures_.begin(); it != v_mixtures_.end(); ++it)
  {
    os << _T("\nParameters of the mixture: ") << (*it)->idData() << _T("\n");
    (*it)->writeParameters(os);
  }
}

void MixtureLearner::randomInit()
{
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("Entering MixtureLearner::RandomInit()\n");
#endif
  if (state() < 2) { initializeStep();}
  for (MixtIterator it = v_mixtures_.begin(); it != v_mixtures_.end(); ++it)
  { (*it)->randomInit();}
  setState(Clust::modelParamInitialized_);
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("MixtureLearner::RandomInit() done\n");
#endif
}


void MixtureLearner::imputationStep()
{
  for (MixtIterator it = v_mixtures_.begin(); it != v_mixtures_.end(); ++it)
  { (*it)->imputationStep();}
}

/* @brief Simulation of all the latent variables and/or missing data
 *  excluding class labels.
 */
void MixtureLearner::samplingStep()
{
  for (MixtIterator it = v_mixtures_.begin(); it != v_mixtures_.end(); ++it)
  { (*it)->samplingStep();}
}

/* store the  intermediate results */
void MixtureLearner::storeIntermediateResults(int iteration)
{
  for (MixtIterator it = v_mixtures_.begin(); it != v_mixtures_.end(); ++it)
  { (*it)->storeIntermediateResults(iteration);}
  meanlnLikelihood_ += (lnLikelihood() - meanlnLikelihood_)/iteration;
}

void MixtureLearner::releaseIntermediateResults()
{
  for (MixtIterator it = v_mixtures_.begin(); it != v_mixtures_.end(); ++it)
  { (*it)->releaseIntermediateResults();}
  meanlnLikelihood_ = 0.;
}

/* Utility method allowing to signal to a mixture to set its parameters */
void MixtureLearner::setParametersStep()
{
  for (MixtIterator it = v_mixtures_.begin(); it != v_mixtures_.end(); ++it)
  { (*it)->setParametersStep();}
  setLnLikelihood(meanlnLikelihood_);
  meanlnLikelihood_ = 0.;
}

/* finalize */
void MixtureLearner::finalizeStep()
{
  for (MixtIterator it = v_mixtures_.begin(); it != v_mixtures_.end(); ++it)
  { (*it)->finalizeStep();}
}

/* @brief Create the composer using existing data handler and ingredients.
 * This method is essentially used by the create() method and can be
 * reused in derived classes. */
void MixtureLearner::createLearner( std::vector<IMixture*> const& v_mixtures)
{
  initialize(nbSample(), 0);
  v_mixtures_.resize( v_mixtures.size());
  for (size_t l = 0; l < v_mixtures_.size(); ++l)
  {
    v_mixtures_[l] = v_mixtures[l]->create();
    v_mixtures_[l]->setMixtureModel(this);
    v_mixtures_[l]->initializeStep();
  }
  // compute number of free parameters
  setNbFreeParameter(computeNbFreeParameters());
}

} /* namespace STK */
