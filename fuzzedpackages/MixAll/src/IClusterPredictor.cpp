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

/** @file IClusterPredictor.cpp
 *  @brief In this file we implement the IClusterPredictor class
 **/


#include "../inst/projects/MixAll/IClusterPredictor.h"

namespace STK
{
IClusterPredictor::IClusterPredictor( Rcpp::S4 s4_model,  Rcpp::S4 s4_clusterPredict)
                                  : ILauncher(s4_model)
                                  , s4_clusterPredict_(s4_clusterPredict)
                                  , s4_algo_(s4_clusterPredict_.slot("algo"))
                                  , p_algo_(createAlgo())
                                  , p_composer_(0)
{}
IClusterPredictor::~IClusterPredictor()
{ if (p_algo_) delete p_algo_;
  if (p_composer_) delete p_composer_;
}
/* utility function creating STK algorithm from R algorithm */
IMixtureAlgoPredict* IClusterPredictor::createAlgo()
{
   String algoName = s4_algo_.slot("algo");
   int nbIterBurn = s4_algo_.slot("nbIterBurn");
   int nbIterLong = s4_algo_.slot("nbIterLong");
   double epsilon = s4_algo_.slot("epsilon");
   Clust::algoPredictType algoType = Clust::stringToPredictAlgo(algoName);
   return Clust::createPredictAlgo(algoType,nbIterBurn,nbIterLong,epsilon);
}

/* get missing values */
void IClusterPredictor::getMissingValues(Clust::MixtureClass const& classModel, String const& idData)
{
  switch (classModel)
  {
    case Clust::Categorical_:
    {
  #ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("In IClusterPredictor::getMissingValues. Getting Categorical missing values\n");
  #endif
      RMatrix<Integer> r_data_int = (SEXP)s4_clusterPredict_.slot("data");
      setCategoricalMissingValuesToMatrix(p_composer_,  idData, r_data_int);
    }
    break;
    case Clust::Poisson_:
    {
  #ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("In IClusterPredictor::getMissingValues. Getting Poisson missing values\n");
  #endif
      RMatrix<Integer> r_data_int = (SEXP)s4_clusterPredict_.slot("data");
      setPoissonMissingValuesToMatrix(p_composer_, idData, r_data_int);
    }
    break;
    case Clust::DiagGaussian_:
    {
  #ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("IClusterPredictor::getMissingValues. Getting Diagonal Gaussian missing values\n");
  #endif
      RMatrix<Real> r_data_num = (SEXP)s4_clusterPredict_.slot("data");
      setDiagGaussianMissingValuesToMatrix(p_composer_, idData, r_data_num);
    }
    break;
    case Clust::Gamma_:
    {
  #ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("In IClusterPredictor::getMissingValues. Getting Gamma missing values\n");
  #endif
      RMatrix<Real> r_data_num = (SEXP)s4_clusterPredict_.slot("data");
      setGammaMissingValuesToMatrix(p_composer_, idData, r_data_num);
    }
    break;
    default:
      break;
  }
}

/* get missing values */
void IClusterPredictor::getMissingValues(Clust::MixtureClass const& classModel, String const& idData, int l)
{
  Rcpp::List ldata_ =s4_clusterPredict_.slot("ldata");
  switch (classModel)
  {
    case Clust::Categorical_:
    {
  #ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("In IClusterPredictor::getMissingValues. Getting Categorical missing values\n");
  #endif
      RMatrix<Integer> r_data_int = (SEXP)ldata_[l];
      setCategoricalMissingValuesToMatrix(p_composer_,  idData, r_data_int);
    }
    break;
    case Clust::Poisson_:
    {
  #ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("In IClusterPredictor::getMissingValues. Getting Poisson missing values\n");
  #endif
      RMatrix<Integer> r_data_int = (SEXP)ldata_[l];
      setPoissonMissingValuesToMatrix(p_composer_, idData, r_data_int);
    }
    break;
    case Clust::DiagGaussian_:
    {
  #ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("IClusterPredictor::getMissingValues. Getting Diagonal Gaussian missing values\n");
  #endif
      RMatrix<Real> r_data_num = (SEXP)ldata_[l];
      setDiagGaussianMissingValuesToMatrix(p_composer_, idData, r_data_num);
    }
    break;
    case Clust::Gamma_:
    {
  #ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("In IClusterPredictor::getMissingValues. Getting Gamma missing values\n");
  #endif
      RMatrix<Real> r_data_num = (SEXP)ldata_[l];
      setGammaMissingValuesToMatrix(p_composer_, idData, r_data_num);
    }
    break;
    default:
      break;
  }
}


} // namespace SK
