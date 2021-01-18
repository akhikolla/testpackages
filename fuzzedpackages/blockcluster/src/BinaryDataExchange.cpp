/*--------------------------------------------------------------------*/
/*     Copyright (C) 2011-2012  Parmeet Singh Bhatia

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as
 published by the Free Software Foundation; either version 2 of the
 License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public
 License along with this program; if not, write to the
 Free Software Foundation, Inc.,
 59 Temple Place,
 Suite 330,
 Boston, MA 02111-1307
 USA

 Contact : parmeet.bhatia@inria.fr , bhatia.parmeet@gmail.com
 */

/*
 * Project:  RCocluster
 * created on: Feb 22, 2012
 * Author: Parmeet Singh Bhatia
 *
 **/

/** @file BinaryCoClust.cpp
 *  @brief
 **/

#include "BinaryDataExchange.h"
#include "coclust/src/Models/BinaryLBModel.h"
#include "coclust/src/Models/BinaryLBModelequalepsilon.h"

void BinaryDataExchange::dataOutput(Rcpp::S4& obj,ICoClustModel* model,bool successful)
{
  if(!successful)
  {
    obj.slot("successful") = false;
    obj.slot("message") = std::string("Co-Clustering Failed! ") + model->errorMsg();
  }
  else
  {
    obj.slot("successful") = true;
    obj.slot("message") = "Co-Clustering successfully terminated!";
    BinaryLBModel* ptrLBM;
    BinaryLBModelequalepsilon* ptrLBMeq;
    switch (strategy_.Model_)
    {
      case pik_rhol_epsilonkl_:
        ptrLBM = dynamic_cast<BinaryLBModel*>(model);
        obj.slot("classmean")       = STK::wrap(ptrLBM->mean());
        obj.slot("classdispersion") = STK::wrap(ptrLBM->dispersion());
        obj.slot("coclusterdata")   = STK::wrap(ptrLBM->arrangedDataClusters());
        break;
      case pik_rhol_epsilon_:
        ptrLBMeq = dynamic_cast<BinaryLBModelequalepsilon*>(model);
        obj.slot("classmean")       = STK::wrap(ptrLBMeq->mean());
        obj.slot("classdispersion") = STK::wrap(ptrLBMeq->dispersion()*STK::Const::Array<STK::Real>(Mparam_.nbrowclust_,Mparam_.nbcolclust_));
        obj.slot("coclusterdata")   = STK::wrap(ptrLBMeq->arrangedDataClusters());
        break;
      case pi_rho_epsilonkl_:
        ptrLBM = dynamic_cast<BinaryLBModel*>(model);
        obj.slot("classmean")       = STK::wrap(ptrLBM->mean());
        obj.slot("classdispersion") = STK::wrap(ptrLBM->dispersion());
        obj.slot("coclusterdata")   = STK::wrap(ptrLBM->arrangedDataClusters());
        break;
      case pi_rho_epsilon_:
        ptrLBMeq = dynamic_cast<BinaryLBModelequalepsilon*>(model);
        obj.slot("classmean")       = STK::wrap(ptrLBMeq->mean());
        obj.slot("classdispersion") = STK::wrap(ptrLBMeq->dispersion()*STK::Const::Array<STK::Real>(Mparam_.nbrowclust_,Mparam_.nbcolclust_));
        obj.slot("coclusterdata")   = STK::wrap(ptrLBMeq->arrangedDataClusters());
        break;
      default:
        Rcpp::stop("Wrong Model in BinaryDataExchange. Please report Bug.");
        break;
    }
    obj.slot("rowclass")          = STK::wrap(model->rowClassificationVector());
    obj.slot("colclass")          = STK::wrap(model->columnClassificationVector());
    obj.slot("rowproportions")    = STK::wrap(model->rowProportions());
    obj.slot("columnproportions") = STK::wrap(model->colProportions());
    obj.slot("rowposteriorprob")  = STK::wrap(model->rowPosteriorProb());
    obj.slot("colposteriorprob")  = STK::wrap(model->colPosteriorProb());
    obj.slot("likelihood")        = model->likelihood();
    obj.slot("ICLvalue")          = model->iclCriteriaValue();
  }
}

void BinaryDataExchange::dataInput(Rcpp::S4 & obj)
{
  STK::RMatrix<STK::Real> data(SEXP(obj.slot("data")));
  m_Dataij_ = data.cast<bool>();
  Mparam_.nbRow_ = m_Dataij_.sizeRows();
  Mparam_.nbCol_ = m_Dataij_.sizeCols();
  //Get Strategy
  //Rcpp::S4 strategy(obj.slot("strategy"));
  //get hyper-parameters
  Rcpp::NumericVector hyperparam(obj.slot("hyperparam"));
  a_ = hyperparam(0);
  b_ = hyperparam(1);
}

void BinaryDataExchange::instantiateModel(ICoClustModel*& model)
{
  if(!strategy_.SemiSupervised)
  {
    switch (strategy_.Model_)
    {
      case pik_rhol_epsilonkl_:
        Mparam_.fixedproportions_ = false;
        model = new BinaryLBModel(m_Dataij_,Mparam_,a_,b_);
        break;
      case pik_rhol_epsilon_:
        Mparam_.fixedproportions_ = false;
        model = new BinaryLBModelequalepsilon(m_Dataij_,Mparam_,a_,b_);
        break;
      case pi_rho_epsilonkl_:
        Mparam_.fixedproportions_ = true;
        model = new BinaryLBModel(m_Dataij_,Mparam_,a_,b_);
        break;
      case pi_rho_epsilon_:
        Mparam_.fixedproportions_ = true;
        model = new BinaryLBModelequalepsilon(m_Dataij_,Mparam_,a_,b_);
        break;
      default:
        Rcpp::stop("Wrong Model in BinaryDataExchange. Please report Bug.");
        break;
    }
  }
  else
  {
    switch (strategy_.Model_)
    {
      case pik_rhol_epsilonkl_:
        Mparam_.fixedproportions_ = false;
        model = new BinaryLBModel(m_Dataij_,v_rowlabels_, v_collabels_,Mparam_,a_,b_);
        break;
      case pik_rhol_epsilon_:
        Mparam_.fixedproportions_ = false;
        model = new BinaryLBModelequalepsilon(m_Dataij_,v_rowlabels_, v_collabels_,Mparam_,a_,b_);
        break;
      case pi_rho_epsilonkl_:
        Mparam_.fixedproportions_ = true;
        model = new BinaryLBModel(m_Dataij_,v_rowlabels_, v_collabels_,Mparam_,a_,b_);
        break;
      case pi_rho_epsilon_:
        Mparam_.fixedproportions_ =true;
        model = new BinaryLBModelequalepsilon(m_Dataij_,v_rowlabels_, v_collabels_,Mparam_,a_,b_);
        break;
      default:
        Rcpp::stop("Wrong Model in BinaryDataExchange. Please report Bug.");
        break;
    }
  }
}
