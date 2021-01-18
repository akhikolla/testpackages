/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2018  Serge Iovleff, Université Lille 1, Inria

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

/*    Project: stkpp::MixAll
 * created on: Mar 17, 2018
 *     Author: iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file ClusterPredictor.cpp
 *  @brief In this file we implement the ClusterPredictor class
 **/


#include "../inst/projects/MixAll/ClusterPredictor.h"

namespace STK
{
ClusterPredictor::ClusterPredictor( Rcpp::S4 s4_model,  Rcpp::S4 s4_clusterPredict)
                                  : IClusterPredictor(s4_model, s4_clusterPredict)
                                  , s4_component_(s4_model_.slot("component"))
{}
ClusterPredictor::~ClusterPredictor() {}

/* run the estimation */
bool ClusterPredictor::run()
{
  String idModel = s4_component_.slot("modelName");
  Clust::Mixture model           = Clust::stringToMixture(idModel);
  Clust::MixtureClass classModel = Clust::mixtureToMixtureClass(model);
  String idData                  = Clust::mixtureToString(model);

  // put data set to data handler
  int nbSample;
  if ((classModel == Clust::Categorical_)||(classModel == Clust::Poisson_))
  {
    Rcpp::IntegerMatrix r_data_int = s4_clusterPredict_.slot("data");
    nbSample  = r_data_int.rows();
    handler_.addData(r_data_int, idData, idModel);
  }
  else
  {
    Rcpp::NumericMatrix r_data_num = s4_clusterPredict_.slot("data");
    nbSample  = r_data_num.rows();
    handler_.addData(r_data_num, idData, idModel);
  }
  // create composer and mixtures
  int nbCluster = s4_model_.slot("nbCluster");
  p_composer_ = new MixtureComposer(nbSample, nbCluster);
  createMixtures(p_composer_);

#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("In ClusterPredictor::run. Setting model parameters\n");
#endif
  // set proportions parameters
  // get parameters from component and set them to composer
  RVector<double> pk((SEXP)s4_model_.slot("pk"));
  p_composer_->setProportions(pk);
//  ArrayXX params;
//  params.move(getParameters(idData, s4_component_));
//  if (!setParameters(p_composer_, idData, params)) { return false;};
  if (!setParameters(p_composer_, idData, getParameters(idData, s4_component_))) { return false;};

#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("In ClusterPredictor::run. running Algo\n");
#endif
  // run prediction algorithm
  p_algo_->setModel(p_composer_);
  bool flag = p_algo_->run();

#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("In ClusterPredictor::run. Getting parameters\n");
#endif
  // get results
  s4_clusterPredict_.slot("pk")  = Rcpp::wrap(p_composer_->pk());
  s4_clusterPredict_.slot("tik") = Rcpp::wrap(p_composer_->tik());
  s4_clusterPredict_.slot("zi")  = Rcpp::wrap(p_composer_->zi());
  Rcpp::NumericVector fi = s4_clusterPredict_.slot("lnFi");
  Rcpp::IntegerVector zi = s4_clusterPredict_.slot("zi");
  for (int i=0; i< fi.length(); ++i)
  {
    fi[i] = p_composer_->computeLnLikelihood(i);
    zi[i] += (1 - baseIdx);  // set base 1 for the class labels
  }
  // get missing values
  getMissingValues(classModel, idData);
  //
  return flag;
}


} // namespace SK
