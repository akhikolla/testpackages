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

/** @file STK_IMixtureStatModel.cpp
 *  @brief In this file we implement the abstract base class for mixture statistical models.
 **/

#include "../include/STK_IMixtureStatModel.h"

namespace STK
{


/* Constructor.
 * @param nbCluster,nbSample number of clusters and samples
 **/
IMixtureStatModel::IMixtureStatModel( int nbSample, int nbCluster)
                                    : IStatModelBase(nbSample)
                                    , nbCluster_(nbCluster)
                                    , pk_(nbCluster, 1./nbCluster)
                                    , tik_(nbSample, nbCluster, 1./nbCluster)
                                    , tk_(nbCluster, Real(nbSample)/nbCluster)
                                    , zi_(nbSample, baseIdx)
                                    , v_mixtures_()
{}

/* copy constructor */
IMixtureStatModel::IMixtureStatModel( IMixtureStatModel const& model)
                                    : IStatModelBase(model)
                                    , nbCluster_(model.nbCluster_)
                                    , pk_(model.pk_), tik_(model.tik_)
                                    , tk_(model.tk_), zi_(model.zi_)
                                    , v_mixtures_(model.v_mixtures_.size())
{
  // clone mixtures
  for (size_t l = 0; l < v_mixtures_.size(); ++l)
  {
    v_mixtures_[l] = model.v_mixtures_[l]->clone();
    v_mixtures_[l]->setMixtureModel(this);
  }
}
/* destructor */
IMixtureStatModel::~IMixtureStatModel()
{
  for (MixtIterator it = v_mixtures_.begin() ; it != v_mixtures_.end(); ++it)
  { delete (*it);}
}

Real IMixtureStatModel::computeLikelihood(int i) const
{ return std::exp(computeLnLikelihood(i));}

/* @return the computed likelihood of the i-th sample.
 *  @param i index of the sample
 **/
Real IMixtureStatModel::computeLnLikelihood(int i) const
{
  // get maximal value
  CPointX lnComp(pk_.size());
  for (int k = pk_.begin(); k< pk_.end(); ++k)
  { lnComp[k] = std::log(pk_[k]) + lnComponentProbability(i, k);}
  // compute result
  Real lnCompMax = lnComp.maxElt();
  return std::log((lnComp-lnCompMax).exp().sum())+lnCompMax;
}

/* @return the computed log-likelihood. */
Real IMixtureStatModel::computeLnLikelihood() const
{
  Real res = 0.0;
  for (int i = tik().beginRows(); i< tik().endRows(); ++i)
  { res += computeLnLikelihood(i);}
  return res;
}

/* @return the computed ICL criteria. */
Real IMixtureStatModel::computeICL() const
{
  Real res = 0.0;
  for (int j = tik().beginCols(); j< tik().endCols(); ++j)
  { res += (tik_.col(j) * (tik_.col(j)+1e-15).log()).sum();}
  // compute result
  return (- 2. * lnLikelihood() + nbFreeParameter() * lnNbSample() - 2. * res);
}

/* Utility lookup function allowing to find a Mixture from its idData
 *  @param idData the id name of the mixture we want to get
 *  @return a pointer on the mixture
 **/
IMixture* IMixtureStatModel::getMixture( String const& idData) const
{
  for (ConstMixtIterator it = v_mixtures_.begin(); it != v_mixtures_.end(); ++it)
  { if ((*it)->idData() == idData) return (*it);}
  return 0;
}
/* register the mixture in the composer*/
void IMixtureStatModel::registerMixture(IMixture* p_mixture)
{
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("In IMixtureStatModel::registerMixture, registering mixture: ")
           << p_mixture->idData() << _T("\n");
#endif
  p_mixture->setMixtureModel(this);
  v_mixtures_.push_back(p_mixture);
  // update NbFreeParameters
  setNbFreeParameter(nbFreeParameter()+p_mixture->nbFreeParameter());
}

/* Utility lookup function allowing to find a Mixture from its idData
 *  @param idData the id name of the mixture we want to get
 *  @return a pointer on the mixture
 **/
void IMixtureStatModel::releaseMixture( String const& idData)
{
  for (MixtIterator it = v_mixtures_.begin(); it != v_mixtures_.end(); ++it)
  {
    if ((*it)->idData() == idData)
    {
       setNbFreeParameter(nbFreeParameter()-(*it)->nbFreeParameter());
       // remove mixture
       delete (*it);
       v_mixtures_.erase(it);
       // update log-likelihood
       if (v_mixtures_.size() == 0)
       { setLnLikelihood(-Arithmetic<Real>::infinity());}
       else
       { setLnLikelihood(computeLnLikelihood());}
       // and break
       break;
    }
  }
}

// implement computeNbFreeParameters
int IMixtureStatModel::computeNbFreeParameters() const
{
  int sum = nbCluster_-1; // proportions
  for (ConstMixtIterator it = v_mixtures_.begin(); it != v_mixtures_.end(); ++it)
  { sum+= (*it)->nbFreeParameter();}
  return sum;
}

/** @brief compute the missing values of the model.
 *  lookup on the mixtures and sum the nbMissingValues.
 *  @return the number of missing values
 **/
int IMixtureStatModel::computeNbMissingValues() const
{
  int sum = nbCluster_-1; // proportions
  for (ConstMixtIterator it = v_mixtures_.begin(); it != v_mixtures_.end(); ++it)
  { sum+= (*it)->nbMissingValues();}
  return sum;
}

/* @brief Initialize the model before its first use.
 *  This function can be overloaded in derived class for initialization of
 *  the specific model parameters. It should be called prior to any used of
 *  the class.
 *  @sa IMixture,MixtureBridge,MixtureLearner
 **/
void IMixtureStatModel::initializeStep()
{
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("Entering IMixtureStatModel::initializeStep\n");
#endif
  if (v_mixtures_.size() == 0)
    STKRUNTIME_ERROR_NO_ARG(IMixtureStatModel::initializeStep,no mixture registered);
  setLnLikelihood(-Arithmetic<Real>::infinity());
  // initialize registered mixtures
  for (MixtIterator it = v_mixtures_.begin(); it != v_mixtures_.end(); ++it)
  { (*it)->initializeStep();}
}

} // namespace STK

