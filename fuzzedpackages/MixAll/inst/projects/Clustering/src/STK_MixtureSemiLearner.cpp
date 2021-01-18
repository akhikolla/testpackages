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
 * created on: 02 June 2016
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_MixtureSemiLearner.cpp
 *  @brief In this file we implement the class MixtureSemiLearner.
 **/

#include <cmath>

#ifdef STK_MIXTURE_DEBUG
#include <Arrays/include/STK_Display.h>
#endif

#include "../include/STK_MixtureSemiLearner.h"

#include <STatistiK/include/STK_Law_Categorical.h>
#include <STatistiK/include/STK_Stat_Functors.h>

namespace STK
{
MixtureSemiLearner::MixtureSemiLearner( int nbSample, int nbCluster)
                                      : MixtureComposer(nbSample, nbCluster)
                                      , knownLabels_()
                                      , unknownLabels_(nbSample)
{
  int i = zi_.begin();
  std::vector<int>::iterator it;
  for (it = unknownLabels_.begin(); it != unknownLabels_.end(); ++it, ++i)
  { *it = i;}
}


/* copy constructor */
MixtureSemiLearner::MixtureSemiLearner( MixtureSemiLearner const& model)
                                      : MixtureComposer(model)
                                      , knownLabels_(model.knownLabels_)
                                      , unknownLabels_(model.unknownLabels_)
{}
/* destructor */
MixtureSemiLearner::~MixtureSemiLearner() {}

/* clone pattern */
MixtureSemiLearner* MixtureSemiLearner::clone() const
{ return new MixtureSemiLearner(*this);}

MixtureSemiLearner* MixtureSemiLearner::create() const
{
  // set dimensions
  MixtureSemiLearner* p_composer = new MixtureSemiLearner(nbSample(), nbCluster());
  p_composer->createSemiLearner(this);
  p_composer->createComposer(v_mixtures());
  return p_composer;
}


/* cStep */
int MixtureSemiLearner::cStep()
{
  std::vector<int>::const_iterator it;
  for (it=unknownLabels_.begin(); it != unknownLabels_.end(); it++)
  { MixtureComposer::cStep(*it);}
  // count the number of individuals in each class
  tk_= Stat::sum(tik_);
  return tk_.minElt();
}

/* simulate zi  */
int MixtureSemiLearner::sStep()
{
  // simulate zi
  std::vector<int>::const_iterator it;
  for (it=unknownLabels_.begin(); it != unknownLabels_.end(); it++)
  { MixtureComposer::sStep(*it);}
  return cStep();
}

/* compute tik, default implementation. */
Real MixtureSemiLearner::eStep()
{
#ifdef STK_MIXTURE_DEBUG
  stk_cout << _T("Entering MixtureSemiLearner::eStep()\n");
#endif
  Real sum = 0.; tk_ =0.;
//#ifdef _OPENMP
//#pragma omp parallel for reduction (+:sum)
//#endif
  std::vector<int>::const_iterator it;
  for(it=unknownLabels_.begin(); it != unknownLabels_.end(); it++)
  { sum += MixtureComposer::eStep(*it);}
  for(it=knownLabels_.begin(); it != knownLabels_.end(); it++)
  { sum += lnComponentProbability(*it, zi_[*it]) + std::log(pk_[zi_[*it]]);}
  // update log-likelihood
  setLnLikelihood(sum);
  // compute proportions
  tk_ = Stat::sumByCol(tik_);
#ifdef STK_MIXTURE_DEBUG
  stk_cout << _T("MixtureSemiLearner::eStep() done\n");
  stk_cout << _T("lnLikelihood =") << sum << _T("\n");
#endif
  return tk_.minElt();
}

/* Compute Zi using the Map estimate, default implementation. */
void MixtureSemiLearner::mapStep()
{
  std::vector<int>::const_iterator it;
  for (it=unknownLabels_.begin(); it != unknownLabels_.end(); it++)
  { MixtureComposer::mapStep(*it);  }
}

/* Create the mixture model parameters. Default implementation is to
 *  set pk_ and tik_ arrays to 1/K value and zi_ to first label value. */
void MixtureSemiLearner::initializeMixtureParameters()
{
  std::vector<int>::const_iterator it;
  for (it=unknownLabels_.begin(); it != unknownLabels_.end(); it++)
  { tik_.row(*it) = 1./nbCluster();}
  for (it=knownLabels_.begin(); it != knownLabels_.end(); it++)
  { tik_.row(*it) = 0.; tik_(*it, zi_[*it]) = 1.;}
}

/* generate random tik_ */
int MixtureSemiLearner::randomTik()
{
  tk_ = 0.;
  std::vector<int>::const_iterator it;
  for (it = unknownLabels_.begin(); it != unknownLabels_.end(); ++it)
  {
    // create a reference on the i-th row
    CPointX tikRowi(tik_.row(*it), true);
    tikRowi.randUnif();
    tikRowi = tikRowi * pk_;
    tikRowi /= tikRowi.sum();
    tk_ += tikRowi;
  }
  return tk_.minElt();
}

/* generate random zi */
int MixtureSemiLearner::randomZi()
{
  Law::Categorical law(pk_);
  std::vector<int>::const_iterator it;
  for (it = unknownLabels_.begin(); it != unknownLabels_.end(); ++it)
  { zi_[*it] = law.rand();}
  // create tik and compute nk
  return cStep();
}

/* Copy the ingredient of the semi-learner in the creation process */
void MixtureSemiLearner::createSemiLearner( MixtureSemiLearner const* const p_learner)
{
  knownLabels_   = p_learner->knownLabels();
  unknownLabels_ = p_learner->unknownLabels();
}

/* Constructor.
 * @param nbCluster,nbSample number of clusters and samples
 */
MixtureSemiLearnerFixedProp::MixtureSemiLearnerFixedProp( int nbSample, int nbCluster)
                                                        : MixtureSemiLearner( nbSample, nbCluster)
{ setNbFreeParameter(0); /* remove the count of the pk parameters */}

/* copy constructor.
 *  @param model the model to copy
 */
MixtureSemiLearnerFixedProp::MixtureSemiLearnerFixedProp( MixtureSemiLearnerFixedProp const& model)
                                                        : MixtureSemiLearner(model) {}

/* Create a composer, but reinitialize the ingredients parameters. */
MixtureSemiLearnerFixedProp* MixtureSemiLearnerFixedProp::create() const
{
  MixtureSemiLearnerFixedProp* p_composer = new MixtureSemiLearnerFixedProp(nbSample(), nbCluster());
  p_composer->createSemiLearner(this);
  p_composer->createComposer(v_mixtures());
  /* remove the count of the pk parameters */
  p_composer->setNbFreeParameter(p_composer->nbFreeParameter()-(nbCluster()-1));
  return p_composer;
}
/* Create a clone of the current model, with ingredients parameters preserved. */
MixtureSemiLearnerFixedProp* MixtureSemiLearnerFixedProp::clone() const
{ return new MixtureSemiLearnerFixedProp(*this);}

/* overloading of the pStep() method. Let them initialized to 1/K. */
void MixtureSemiLearnerFixedProp::pStep() {}



} // namespace STK

