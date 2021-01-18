/*--------------------------------------------------------------------*/
/*     Copyright (C) 2011-2012  Parmeet Singh Bhatia

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

    Author : Parmeet Bhatia
 Contact : bhatia.parmeet@gmail.com , serge.iovleff@stkpp.org

*/

/*
 * Project:  blockcluster::
 * created on: Apr 11, 2012
 * Purpose:  .
 * Author:   modal-laptop, parmeet.bhatia@inria.fr
 *
 **/

/** @file ICoClust.cpp
 *  @brief In this file .
 **/

#include "IDataExchange.h"
//CoClust::Algorithms
#include "coclust/src/Algorithms/EMAlgo.h"
#include "coclust/src/Algorithms/CEMAlgo.h"
#include "coclust/src/Algorithms/SEMAlgo.h"
#include "coclust/src/Algorithms/GibbsAlgo.h"
//CoClust::Initialization
#include "coclust/src/Initialization/CEMInit.h"
#include "coclust/src/Initialization/EMInit.h"
#include "coclust/src/Initialization/RandomInit.h"
//CoClust::strategy
#include "coclust/src/Strategy/XStrategyAlgo.h"
#include "coclust/src/Strategy/XStrategyforSEMAlgo.h"

IDataExchange::~IDataExchange()
{}

void IDataExchange::setInput(Rcpp::S4 & obj)
{
  //Get Strategy
  Rcpp::S4 strategy(obj.slot("strategy"));
  strategy_.Algo_          = S_Algorithm[Rcpp::as<std::string>(strategy.slot("algo"))];
  strategy_.stopcriteria_  = S_StopCriteria[Rcpp::as<std::string>(strategy.slot("stopcriteria"))];
  strategy_.Init_          = S_Init[Rcpp::as<std::string>(strategy.slot("initmethod"))];
  strategy_.DataType_      = S_DataType[Rcpp::as<std::string>(obj.slot("datatype"))];
  strategy_.Model_         = S_Model[Rcpp::as<std::string>(obj.slot("model"))];
  strategy_.SemiSupervised = Rcpp::as<bool>(obj.slot("semisupervised"));

  //Get strategy parameters
  strategyParam_.nbinitmax_  = Rcpp::as<int>(strategy.slot("nbinitmax"));
  strategyParam_.nbiter_xem_ = Rcpp::as<int>(strategy.slot("nbiterationsxem"));
  strategyParam_.nbiter_XEM_ = Rcpp::as<int>(strategy.slot("nbiterationsXEM"));
  strategyParam_.nbtry_      = Rcpp::as<int>(strategy.slot("nbtry"));
  strategyParam_.nbxem_      = Rcpp::as<int>(strategy.slot("nbxem"));

  //get row/column labels if semisupervised
  if(strategy_.SemiSupervised)
  {
    v_rowlabels_ = RVectorInt((SEXP(obj.slot("rowlabels"))));
    v_collabels_ = RVectorInt((SEXP(obj.slot("collabels"))));
  }

  // Set stopping-criteria
  switch (strategy_.stopcriteria_)
  {
    case parameter_:
      strategyParam_.Stop_Criteria = &ICoClustModel::parameterStopCriteria;
      break;
    case likelihood_:
      strategyParam_.Stop_Criteria = &ICoClustModel::likelihoodStopCriteria;
      break;
    default:
      strategyParam_.Stop_Criteria = &ICoClustModel::parameterStopCriteria;
      break;
  }
  //Set various  model parameters

  //get various iterations
  Mparam_.nbinititerations_ = Rcpp::as<int>(strategy.slot("nbinititerations"));
  Mparam_.nbiterations_int_ = Rcpp::as<int>(strategy.slot("nbiterations_int"));

  //get threshold values
  Mparam_.eps_xem_ = Rcpp::as<double>(strategy.slot("epsilonxem"));
  Mparam_.eps_XEM_ = Rcpp::as<double>(strategy.slot("epsilonXEM"));
  Mparam_.initepsilon_ = Rcpp::as<double>(strategy.slot("initepsilon"));
  Mparam_.epsilon_int_ = Rcpp::as<double>(strategy.slot("epsilon_int"));

 // initialize epsilon_
  Mparam_.epsilon_ = Mparam_.eps_xem_;

  //get row and column clusters
  Rcpp::IntegerVector nbcocluster(SEXP(obj.slot("nbcocluster")));
  Mparam_.nbrowclust_ = nbcocluster[0];
  Mparam_.nbcolclust_ = nbcocluster[1];
}

void IDataExchange::initializeParamEnum()
{
  //Datatype
  S_DataType["binary"]      = binary_;
  S_DataType["contingency"] = contingency_;
  S_DataType["continuous"]  = continuous_;
  S_DataType["continuous"]  = categorical_;

  //Algorithm
  S_Algorithm["BEM"]   = bem_;
  S_Algorithm["BCEM"]   = bcem_;
  S_Algorithm["BSEM"]   = bsem_;
  S_Algorithm["BGibbs"] = bgibbs_;

  //StopCriteria
  S_StopCriteria["Parameter"] = parameter_;
  S_StopCriteria["Likelihood"] = likelihood_;

  //Initialization
  S_Init["cemInitStep"] = e_CEMInit_;
  S_Init["emInitStep"] = e_EMInit_;
  S_Init["randomInit"] = e_RandomInit_;

  //Models
  S_Model["pi_rho_epsilon"] = pi_rho_epsilon_;
  S_Model["pik_rhol_epsilon"] = pik_rhol_epsilon_;
  S_Model["pi_rho_epsilonkl"] = pi_rho_epsilonkl_;
  S_Model["pik_rhol_epsilonkl"] = pik_rhol_epsilonkl_;
  S_Model["pi_rho_unknown"] = pi_rho_unknown_;
  S_Model["pik_rhol_unknown"] = pik_rhol_unknown_;
  S_Model["pi_rho_known"] = pi_rho_known_;
  S_Model["pik_rhol_known"] = pik_rhol_known_;
  S_Model["pi_rho_sigma2"] = pi_rho_sigma2_;
  S_Model["pik_rhol_sigma2"] = pik_rhol_sigma2_;
  S_Model["pi_rho_sigma2kl"] = pi_rho_sigma2kl_;
  S_Model["pik_rhol_sigma2kl"] = pik_rhol_sigma2kl_;
  S_Model["pi_rho_multi"] = pi_rho_multi_;
  S_Model["pik_rhol_multi"] = pik_rhol_multi_;
}

void IDataExchange::instantiateAlgo(IAlgo*& algo,IStrategy*& strategy)
{
  switch (strategy_.Algo_)
  {
    case bem_:
      algo = new EMAlgo();
      strategy = new XStrategyAlgo(strategyParam_);
      break;
    case bcem_:
      algo = new CEMAlgo();
      strategy = new XStrategyAlgo(strategyParam_);
      break;
    case bsem_:
      algo = new SEMAlgo();
      strategy = new XStrategyforSEMAlgo(strategyParam_);
      break;
    case bgibbs_:
      algo = new GibbsAlgo();
      strategy = new XStrategyforSEMAlgo(strategyParam_);
      break;
    default:
      algo = new EMAlgo();
      strategy = new XStrategyAlgo(strategyParam_);
      break;
  }
}

void IDataExchange::instantiateInit(IInit*& init)
{
  switch (strategy_.Init_)
  {
    case e_CEMInit_:
      init = new CEMInit();
      break;
    case e_EMInit_:
      init = new emInit();
      break;
    case e_RandomInit_:
      init = new RandomInit();
      break;
    default:
      init = new CEMInit();
      break;
  }
}
