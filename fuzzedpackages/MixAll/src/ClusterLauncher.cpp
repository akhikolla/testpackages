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
 * created on: 4 sept. 2013
 * Author:   iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_ClusterLauncher.cpp
 *  @brief In this file we implement the ClusterLauncher which
 *  construct properly a mixture model.
 **/


#include "../inst/projects/MixAll/ClusterLauncher.h"
#include "../inst/projects/MixAll/ClusterFacade.h"

using namespace Rcpp;

namespace STK
{
/* facade design pattern.
 * The ClusterLauncher allow to create the strategy for estimating a mixture model
 * with less effort
 **/
ClusterLauncher::ClusterLauncher( Rcpp::S4 model, Rcpp::IntegerVector nbCluster, Rcpp::CharacterVector models)
                                : ILauncher(model, models)
                                , s4_strategy_(s4_model_.slot("strategy"))
                                , v_nbCluster_(nbCluster)
                                , criterion_(Rcpp::as<String>(s4_model_.slot("criterionName")))
                                , p_composer_(0)
                                , isMixedData_(false)
{}
/* facade design pattern.
 * The ClusterLauncher allow to create the strategy for estimating a mixture model
 * with less effort
 **/
ClusterLauncher::ClusterLauncher( Rcpp::S4 model, Rcpp::IntegerVector nbCluster)
                                : ILauncher(model)
                                , s4_strategy_(s4_model_.slot("strategy"))
                                , v_nbCluster_(nbCluster)
                                , criterion_(Rcpp::as<String>(s4_model_.slot("criterionName")))
                                , p_composer_(0)
                                , isMixedData_(true)

{}
/* destructor. */
ClusterLauncher::~ClusterLauncher()
{ if (p_composer_) delete p_composer_; p_composer_ = 0;}

/* run the estimation */
bool ClusterLauncher::run()
{
  // compute the best model
  Real initCriter = s4_model_.slot("criterion");
  Real criter = (isMixedData_) ? selectBestMixedModel() : selectBestSingleModel();
  if (!Arithmetic<Real>::isFinite(criter)) return false;

  // get result common part of the estimated model
  s4_model_.slot("criterion")      = criter;
  s4_model_.slot("nbCluster")      = p_composer_->nbCluster();
  s4_model_.slot("lnLikelihood")   = p_composer_->lnLikelihood();
  s4_model_.slot("nbFreeParameter")= p_composer_->nbFreeParameter();
  s4_model_.slot("pk")             = Rcpp::wrap(p_composer_->pk());
  s4_model_.slot("tik")            = Rcpp::wrap(p_composer_->tik());
  s4_model_.slot("zi")             = Rcpp::wrap(p_composer_->zi());
  s4_model_.slot("ziFit")          = Rcpp::wrap(p_composer_->zi());
  NumericVector fi = s4_model_.slot("lnFi");
  IntegerVector zi = s4_model_.slot("zi");
  for (int i=0; i< fi.length(); ++i)
  {
    fi[i] = p_composer_->computeLnLikelihood(i);
    zi[i] += (1 - baseIdx);  // set base 1 for the class labels
  }
  if (criter == initCriter || !Arithmetic<Real>::isFinite(criter)) return false;
  return true;
}

/* get the parameters */
Real ClusterLauncher::selectBestSingleModel()
{
  // component
  Rcpp::S4 s4_component = s4_model_.slot("component");
  double critValue      = s4_model_.slot("criterion");
  int nbSample          = s4_model_.slot("nbSample");
  String idDataBestModel;
  IMixtureComposer*  p_current   =0;
  IMixtureCriterion* p_criterion =0;

  bool freeProp;

  // loop over all the models
  for (int l=0; l <v_models_.size(); ++l)
  {
    // create idData
    String idData  = "model" + typeToString<int>(l);
    String idModel = Rcpp::as<String>(v_models_[l]);
    // transform R model names to STK++ model names
    // check have been done on the R side so.... Let's go
    Clust::Mixture model           = Clust::stringToMixture(idModel, freeProp);
    Clust::MixtureClass classModel = Clust::mixtureToMixtureClass(model);
    // add Data set with the new model name, m_data_* is just a pointer on a SEXP
    // structure thus there is no difficulties in doing so
    if ((classModel == Clust::Categorical_)||(classModel == Clust::Poisson_))
    {
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("In ClusterLauncher::selectBestSingleModel. Creating IntegerMatrix\n");
#endif
      Rcpp::IntegerMatrix m_data_int=s4_component.slot("data");
      createDataSets(m_data_int, idData, model);
    }
    else
    {
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("In ClusterLauncher::selectBestSingleModel. Creating NumericMatrix\n");
#endif
      Rcpp::NumericMatrix m_data_num=s4_component.slot("data");
      createDataSets(m_data_num, idData, model);
    }
  }
  // start computation
  try
  {
    // create criterion
    p_criterion = Clust::createCriterion(criterion_);

    // start the estimation process, should end with the best model according to
    // the criteria
    ClusterFacade facade(p_current);
    facade.createFullStrategy(s4_strategy_);
    // loop over all the models
    for (int l=0; l <v_models_.size(); ++l)
    {
      String idData  = "model" + typeToString<int>(l);
      String idModel = Rcpp::as<String>(v_models_[l]);
      bool freeProp;
      Clust::stringToMixture(idModel, freeProp);
      // for the current model,
      for (int k=0; k <v_nbCluster_.length(); ++k)
      {
        int K = v_nbCluster_[k];
        // create composer
        if (freeProp) { p_current = new MixtureComposer(nbSample, K);}
        else          { p_current = new MixtureComposerFixedProp(nbSample, K);}

        // create current mixture and register it
        createMixtures(static_cast<MixtureComposer*>(p_current));

        // run estimation and get results if possible
        if (!facade.run()) { msg_error_ += facade.error();}
        // compute criterion and update model if necessary
        p_criterion->setModel(p_current);
        p_criterion->run();
        if (critValue > p_criterion->value())
        {
          if (p_composer_) { std::swap(p_current, p_composer_);}
          else             { p_composer_ = p_current; p_current = 0;}
          s4_component.slot("modelName") = idModel;
          idDataBestModel                = idData;
          critValue = p_criterion->value();
        }
        // release current composer
        if (p_current) { delete p_current; p_current = 0;}
      }
    }
    // release
    if (p_criterion) { delete p_criterion; p_criterion = 0;}
    // get specific parameters
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("In ClusterLauncher::selectBestSingleModel. Getting Parameters from best model\n");
#endif
    setParametersToComponent(p_composer_, idDataBestModel, s4_component);
    return critValue;
  }
  catch (Exception const& e)
  {
    if (p_current)   { delete p_current; p_current = 0;}
    if (p_criterion) { delete p_criterion; p_criterion = 0;}
    ::Rf_error(e.error().c_str()) ;
  }
  // failed
  return Arithmetic<Real>::max();
}

/* select best mixed data model */
Real ClusterLauncher::selectBestMixedModel()
{
  // list of the component
  Rcpp::List s4_list =s4_model_.slot("lcomponent");
  Real criter        =s4_model_.slot("criterion");
  int nbSample       =s4_model_.slot("nbSample");
  // main pointers
  IMixtureComposer* p_current =0;
  IMixtureCriterion* p_criterion =0;
  try
  {
    // create facade and strategy
    ClusterFacade facade(p_current);
    facade.createFullStrategy(s4_strategy_);
    bool sameProp = true;
    // loop over the list of component and fill handler_
    for (int l=0; l <s4_list.size(); ++l)
    {
      // component
      Rcpp::S4 s4_component = s4_list[l];
      String idData  = "model" + typeToString<int>(l);
      String idModel = s4_component.slot("modelName");
      // register
      bool freeMixture;
      Clust::Mixture model = Clust::stringToMixture(idModel, freeMixture);
      Clust::MixtureClass classModel = Clust::mixtureToMixtureClass(model);
      // if one of the model is free proportion, then we use free proportion
      sameProp &= (!freeMixture);
      // add Data set with the new model name
      if ((classModel == Clust::Categorical_)||(classModel == Clust::Poisson_))
      {
        Rcpp::IntegerMatrix m_data = s4_component.slot("data");
        createDataSets(m_data, idData, model);
      }
      else
      {
        Rcpp::NumericMatrix m_data = s4_component.slot("data");
        createDataSets(m_data, idData, model);
      }
    }
    // create criterion
    p_criterion = Clust::createCriterion(criterion_);

    for (int k=0; k <v_nbCluster_.length(); ++k)
    {
      int K = v_nbCluster_[k];
      // create composer
      if (!sameProp) { p_current = new MixtureComposer(nbSample, K);}
      else           { p_current = new MixtureComposerFixedProp(nbSample, K);}
      // create all mixtures
      createMixtures(static_cast<MixtureComposer*>(p_current));
      // run estimation and get results if possible
      if (!facade.run()) { msg_error_ += facade.error();}
      // compute criterion and update model if necessary
      p_criterion->setModel(p_current);
      p_criterion->run();
      if (criter > p_criterion->value())
      {
        criter = p_criterion->value();
        if (p_composer_) { std::swap(p_current, p_composer_);}
        else             { p_composer_ = p_current; p_current = 0;}
      }
      // release current composer
      if (p_current) { delete p_current; p_current = 0;}
    }
    // release
    delete p_criterion;
    // get parameters
    for (int l=0; l <s4_list.size(); ++l)
    {
      // get component
      Rcpp::S4 s4_component = s4_list[l];
      // id of the data set and of the model
      String idData  = "model" + typeToString<int>(l);
      setParametersToComponent(p_composer_, idData, s4_component);
    }
    //
    return criter;
  }
  catch (Exception const& e)
  {
    if (p_current) delete p_current;
    if (p_criterion) delete p_criterion;
    ::Rf_error(e.error().c_str()) ;
  }
  // failed
  return Arithmetic<Real>::max();
}

} // namespace STK

