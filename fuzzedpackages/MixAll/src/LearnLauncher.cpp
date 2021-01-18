/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2014  Serge Iovleff, University Lille 1, Inria

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU General Public
    License along with this program; if not, write to the
    Free Software Foundation, Inc.,
    59 Temple Place,
    Suite 330,
    Boston, MA 02111-1307
    USA

    Contact : S..._DOT_I..._AT_stkpp.org (see copyright for ...)
*/

/*
 * Project:  MixAll
 * created on: 15 may 2016
 * Author:   iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_LearnLauncher.cpp
 *  @brief In this file we implement the LearnLauncher which
 *  construct properly a mixture model.
 **/


#include "../inst/projects/MixAll/LearnLauncher.h"

using namespace Rcpp;

namespace STK
{
/* facade design pattern.
 * The LearnLauncher allow to create the strategy for estimating a mixture model
 * with less effort
 **/
LearnLauncher::LearnLauncher( Rcpp::S4 model, Rcpp::CharacterVector models, Rcpp::S4 algo )
                            : ILauncher(model, models)
                            , s4_algo_(algo)
                            , criterion_(Rcpp::as<String>(s4_model_.slot("criterionName")))
                            , p_algo_(0)
                            , p_criterion_(0)
                            , p_learner_(0)
                            , isMixedData_(false)

{}
/* facade design pattern.
 * The LearnLauncher allow to create the strategy for estimating a mixture model
 * with less effort
 **/
LearnLauncher::LearnLauncher( Rcpp::S4 model, Rcpp::S4 algo )
                            : ILauncher(model)
                            , s4_algo_(algo)
                            , criterion_(Rcpp::as<String>(s4_model_.slot("criterionName")))
                            , p_algo_(0)
                            , p_criterion_(0)
                            , p_learner_(0)
                            , isMixedData_(true)

{}
/* destructor. */
LearnLauncher::~LearnLauncher()
{
  if (p_algo_)      delete p_algo_;
  if (p_criterion_) delete p_criterion_;
  if (p_learner_)   delete p_learner_;
}

/* run the estimation */
bool LearnLauncher::run()
{
  p_criterion_ = Clust::createCriterion(criterion_);
  if (!p_criterion_)
  { msg_error_ = STKERROR_1ARG(LearnLauncher::run,criterion_,Wrong criterion name);
    return false;
  }
  // create algo runner
  String algoName = s4_algo_.slot("algo");
  Real epsilon    = s4_algo_.slot("epsilon");
  int nbIter      = s4_algo_.slot("nbIteration");
  if (toUpperString(algoName) == "SIMUL")  { p_algo_ = new SimulAlgo();}
  else if (toUpperString(algoName) == "IMPUTE") { p_algo_ = new ImputeAlgo();}
    else
    { msg_error_ = STKERROR_1ARG(LearnLauncher::run,algoName,Wrong algo name);
      return false;
    }
  p_algo_->setEpsilon(epsilon);
  p_algo_->setNbIterMax(nbIter);
  // launch computations
  Real initCriter = s4_model_.slot("criterion");
  Real criter = (isMixedData_) ? selectBestMixedModel()
                               : selectBestSingleModel();

  // release criterion and algo,
  delete p_criterion_; p_criterion_ = 0;
  delete p_algo_; p_algo_ = 0;
  if (!Arithmetic<Real>::isFinite(criter)) return false;

  // get result common part of the estimated model
  s4_model_.slot("criterion")       = criter;
  s4_model_.slot("lnLikelihood")    = p_learner_->lnLikelihood();
  s4_model_.slot("nbFreeParameter") = p_learner_->nbFreeParameter();
  // get zi with predictions
  RVector<int> ziFit = (SEXP)s4_model_.slot("ziFit");
  ziFit = p_learner_->ziPred();
   // compute lnFi
  NumericVector fi = s4_model_.slot("lnFi");
  for (int i=0; i< fi.length(); ++i)
  { fi[i] = p_learner_->computeLnLikelihood(i);}
  // results
  if (criter == initCriter || !Arithmetic<Real>::isFinite(criter)) return false;
  return true;
}

/* get the parameters */
Real LearnLauncher::selectBestSingleModel()
{
  String idDataBestModel;
  // component
  Rcpp::S4 s4_component = s4_model_.slot("component");

  // Get dimensions
  double critValue = s4_model_.slot("criterion");
  int nbSample     = s4_model_.slot("nbSample");
  int K            = s4_model_.slot("nbCluster");
  Rcpp::NumericMatrix r_tik = s4_model_.slot("tik");
  Rcpp::NumericVector r_pk  = s4_model_.slot("pk");
  RMatrix<Real> tik(r_tik);
  RVector<Real> pk(r_pk);
  // interface classes
  IMixtureLearner*   p_current   =0;
  // create criterion
  // loop over all the models and fill
  for (int l=0; l <v_models_.size(); ++l)
  {
    // create idData
    String idData  = "model" + typeToString<int>(l);
    String idModel = Rcpp::as<String>(v_models_[l]);
    // transform R model names to STK++ model names
    // check have been done on the R side so.... Let's go
    bool freeProp;
    Clust::Mixture model           = Clust::stringToMixture(idModel, freeProp);
    Clust::MixtureClass classModel = Clust::mixtureToMixtureClass(model);
    // add Data set with the new model name, m_data is just a pointer on a SEXP
    // structure thus there is no difficulties in doing so
    if ((classModel == Clust::Categorical_)||(classModel == Clust::Poisson_))
    { createDiscreteDataSets(idData, s4_component, model);}
    else
    { createContinuousDataSets(idData, s4_component, model);}
  }

  // start computation
  try
  {
    // start the estimation process, should end with the best model according to
    // the criteria
    // loop over all the models
    for (int l=0; l <v_models_.size(); ++l)
    {
      String idData = "model" + typeToString<int>(l);
      String idModel = Rcpp::as<String>(v_models_[l]);
      // create learner and mixtures
      p_current = new MixtureLearner(nbSample, K);
      p_current->setMixtureParameters(tik, pk);
      createMixtures(static_cast<MixtureLearner*>(p_current));
      // start algo
      p_algo_->setModel(p_current);
      if (p_algo_->run())
      {
        // compute criterion and update model if necessary
        p_criterion_->setModel(p_current);
        p_criterion_->run();
        if (critValue > p_criterion_->value())
        {
          s4_component.slot("modelName") = idModel;
          idDataBestModel = idData;
          critValue = p_criterion_->value();
          // save current value
          std::swap(p_current, p_learner_);
        }
        // release current learner
        if (p_current) { delete p_current; p_current = 0;}
      }
    }
    // get specific parameters
    setParametersToComponent(p_learner_, idDataBestModel, s4_component);
    return critValue;
  }
  catch (Exception const& e)
  {
    if (p_current) delete p_current;
    ::Rf_error(e.error().c_str()) ;
  }
  // failed
  return Arithmetic<Real>::max();
}

/* select best mixed data model */
Real LearnLauncher::selectBestMixedModel()
{
  // list of the component
  Rcpp::List s4_list =s4_model_.slot("lcomponent");
  Real criter  = s4_model_.slot("criterion");
  int nbSample = s4_model_.slot("nbSample");
  int K        = s4_model_.slot("nbCluster");
  Rcpp::NumericMatrix r_tik = s4_model_.slot("tik");
  Rcpp::NumericVector r_pk  = s4_model_.slot("pk");
  RMatrix<Real> tik(r_tik);
  RVector<Real> pk(r_pk);

  // main pointer
  IMixtureLearner* p_current =0;
  try
  {
    bool sameProp = true;
    // loop over the list of component and fil handler_
    for (int l=0; l <s4_list.size(); ++l)
    {
      // component
      Rcpp::S4 s4_component = s4_list[l];
      String idData  = "model" + typeToString<int>(l);
      String idModel = s4_component.slot("modelName");
      // register
      bool freeMixture;
      Clust::Mixture model           = Clust::stringToMixture(idModel, freeMixture);
      Clust::MixtureClass classModel = Clust::mixtureToMixtureClass(model);
      // if one of the model is free proportion, then we use free proportion
      sameProp &= (!freeMixture);
      // add Data set with the new model name, m_data is just a pointer on a SEXP
      // structure thus there is no difficulties in doing so
      if ((classModel == Clust::Categorical_)||(classModel == Clust::Poisson_))
      {
        NumericMatrix m_data = s4_component.slot("data");
        createDataSets(m_data, idData, model);
//        createDiscreteDataSets(idData, s4_component, model);
      }
      else
      {
        IntegerMatrix m_data = s4_component.slot("data");
        createDataSets(m_data, idData, model);
//        createContinuousDataSets(idData, s4_component, model);
      }
    }
    // create learner
    p_current = new MixtureLearner(nbSample, K);
    p_current->setMixtureParameters(tik, pk);
    // create all mixtures
    createMixtures(static_cast<MixtureLearner*>(p_current));
    // start algo
    p_algo_->setModel(p_current);
    if (p_algo_->run())
    {
      // compute criterion and update model if necessary
      p_criterion_->setModel(p_current);
      p_criterion_->run();
      if (criter > p_criterion_->value())
      {
        std::swap(p_current, p_learner_);
        criter = p_criterion_->value();
      }
    }
    // release current learner if necessary (p_crrent should be 0)
    if (p_current) { delete p_current; p_current = 0;}
    // get parameters
    for (int l=0; l <s4_list.size(); ++l)
    {
      // component
      Rcpp::S4 s4_component = s4_list[l];
      // id of the data set and of the model
      String idData  = "model" + typeToString<int>(l);
      setParametersToComponent(p_learner_, idData, s4_component);
    }
    //
    return criter;
  }
  catch (Exception const& e)
  {
    if (p_current) delete p_current;
    ::Rf_error(e.error().c_str()) ;
  }
  // failed
  return Arithmetic<Real>::max();
}


} // namespace STK

