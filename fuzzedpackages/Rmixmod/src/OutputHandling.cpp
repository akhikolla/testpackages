//
//  OutputHandling.cpp
//  Rmixmod
//
//  Created by Rémi Lebret on 06/02/12.
//  Copyright (c) 2012 CNRS. All rights reserved.
//

#include "Conversion.h"
#include "OutputHandling.h"

#include "mixmod/Kernel/IO/ModelOutput.h"

#include "mixmod/Kernel/IO/LabelDescription.h"
#include "mixmod/Kernel/IO/Label.h"

#include "mixmod/Kernel/IO/ProbaDescription.h"
#include "mixmod/Kernel/IO/Proba.h"
#include "mixmod/Kernel/IO/ParameterDescription.h"

#include "mixmod/Matrix/Matrix.h"

#include <algorithm>

// constructor
OutputHandling::OutputHandling( XEM::ModelOutput* MOutput, Rcpp::S4& xem, const XEM::DataType dataType )
                              : MOutput_(MOutput)
                              , xem_(xem)
                              , nbCluster_((int)MOutput_->getNbCluster())       
{
  
  // add nbCluster value
  xem_.slot("nbCluster") = nbCluster_;
  // add model selected
  xem_.slot("model") = XEM::ModelNameToString(MOutput->getModelType().getModelName());
  // add the error
  xem_.slot("error") = (MOutput->getStrategyRunError()).what();
  
  // fill other slot only if no error
  if ( dynamic_cast<XEM::Exception&>(MOutput->getStrategyRunError()) == XEM::NOERROR ){
    // add likelihood value
    xem_.slot("likelihood") = MOutput_->getLikelihood();
    
    // add parameters
    switch (dataType) {
      case XEM::QuantitativeData:
        setGaussianParameter();
        break;
      case XEM::QualitativeData:
        setMultinomialParameter();
        break;
      case XEM::HeterogeneousData:
        setCompositeParameter();
        break;
      default:
        break;
    }
  }
}
OutputHandling::OutputHandling( XEM::Parameter* par, Rcpp::S4& xem, const XEM::DataType dataType, int nbCluster )
                              : MOutput_(nullptr)
                              , xem_(xem)
                              , nbCluster_(nbCluster)       
{
  
  // add nbCluster value
  //xem_.slot("nbCluster") = nbCluster_;
    
  // add parameters
  switch (dataType) {
  case XEM::QuantitativeData:
    setGaussianParameter(dynamic_cast<XEM::GaussianEDDAParameter const *>(par));
    break;
  case XEM::QualitativeData:
    setMultinomialParameter(dynamic_cast<XEM::BinaryParameter const *>(par));
    break;
  case XEM::HeterogeneousData:
    setCompositeParameter(dynamic_cast<XEM::CompositeParameter *>(par));
    break;
  default:
    break;
  }
  
}

// destructor
OutputHandling::~OutputHandling()
{ }


// set gaussian paramaters 
void OutputHandling::setGaussianParameter(XEM::GaussianEDDAParameter const * parArg)
{
  // get pointer to gaussian parameters
  const XEM::GaussianEDDAParameter * gParam = parArg!=nullptr ? parArg : dynamic_cast<XEM::GaussianEDDAParameter const *>(MOutput_->getParameterDescription()->getParameter());
  // get the number of variables
  nbVariable_ = gParam->getPbDimension();
  
  // create parameter object
  Rcpp::S4 param(xem_.slot("parameters"));
  
  // add proportions values
  param.slot("proportions") = Conversion::CVectorToRcppVector(nbCluster_,gParam->getTabProportion());
  
  // add proportions values
  param.slot("mean") = Conversion::CMatrixToRcppMatrix(nbCluster_,nbVariable_,gParam->getTabMean());
  
  // add variances
  // define list of variance matrix
  Rcpp::GenericVector varianceMatrices(nbCluster_);
  /** Array with matrix variance for each class*/
  XEM::Matrix** matrixVariance = gParam->getTabSigma();
  // loop over clusters
  for (int i=0; i<nbCluster_; i++)
  {
    //varianceMatrices[i] = Conversion::CMatrixToRcppMatrix(nbVariable_,nbVariable_,matrixVariance[i]->storeToArray());
    std::unique_ptr<double*[],XEM::TabDeleter<double>> stoa(matrixVariance[i]->storeToArray(),
                                                            XEM::TabDeleter<double>(matrixVariance[i]->getPbDimension()));
    varianceMatrices[i] = Conversion::CMatrixToRcppMatrix(nbVariable_,nbVariable_,stoa.get());    
  }
  // add variances to S4 object
  param.slot("variance") = varianceMatrices;
  param.slot("nbFreeParam") = gParam->getFreeParameter();
  // add parameters to the output list
  xem_.slot("parameters") = param;
}

// set multinomial parameters 
void OutputHandling::setMultinomialParameter(XEM::BinaryParameter const * parArg)
{
  // get pointer to multinomial parameters
  const XEM::BinaryParameter * bParam = parArg!=nullptr ? parArg : dynamic_cast<XEM::BinaryParameter const *>(MOutput_->getParameterDescription()->getParameter());
  
  // get the number of variables
  nbVariable_ = bParam->getPbDimension();
  // create parameter object
  Rcpp::S4 param(xem_.slot("parameters"));
  
  // add proportions values
  param.slot("proportions") = Conversion::CVectorToRcppVector(nbCluster_,bParam->getTabProportion());
  
  // add means values
  param.slot("center") = Conversion::CMatrixToRcppMatrixForInt(nbCluster_,nbVariable_,bParam->getTabCenter());
  
  //add factor
  param.slot("factor") = Conversion::CVectorToRcppVectorForInt(nbVariable_,bParam->getTabNbModality());
  param.slot("nbFreeParam") = bParam->getFreeParameter();
  
  //-------------------
  // get scatter
  //-------------------
  // get pointer to scatter
  double *** scatter = bParam->scatterToArray();
  // get tab of modalities
  int64_t* tabNbModality = bParam->getTabNbModality();
  // get maximum number of modality
  int64_t max = *max_element(tabNbModality,tabNbModality+nbVariable_);
  
  //VectorMatrix of NumericMatrix
  Rcpp::GenericVector vectorOutput(nbCluster_);
  
  // loop over clusters
  for (int k=0; k<nbCluster_; k++){
    //NumericMatrix for matrix
    Rcpp::NumericMatrix matrixOutput(nbVariable_,max);
    // loop over variables
    for(int j=0; j<nbVariable_; j++){
      // loop over modalities
      for (int h=0; h<tabNbModality[j]; h++) {
        matrixOutput(j,h) = scatter[k][j][h];
      }
    }
    vectorOutput(k) = matrixOutput;
  }  
  // add scatters
  param.slot("scatter") = vectorOutput;
  
  // add parameters to the output list
  xem_.slot("parameters") = param;
}

// set composite parameters
void OutputHandling::setCompositeParameter(XEM::CompositeParameter  *parArg)
{
  XEM::CompositeParameter * cParam =  parArg!=nullptr ? parArg : dynamic_cast<XEM::CompositeParameter  *>(MOutput_->getParameterDescription()->getParameter());
  // get pointer to gaussian parameters
  const XEM::GaussianEDDAParameter * gParam = dynamic_cast<XEM::GaussianEDDAParameter const *>(cParam->getGaussianParameter());
  // get pointer to multinomial parameters
  const XEM::BinaryParameter * bParam = static_cast<XEM::BinaryParameter const *>(cParam->getBinaryParameter());

  // get the number of variables
  int64_t nbVariable_Gaussian = gParam->getPbDimension();
  int64_t nbVariable_Binary = bParam->getPbDimension();

  // create parameter object
  Rcpp::S4 param(xem_.slot("parameters"));
  //create gaussian parameter object
  Rcpp::S4 g_param(param.slot("g_parameter"));
  //create binary parameter object
  Rcpp::S4 b_param(param.slot("m_parameter"));

  //Fill Gaussian parameters slot
  // add proportions values
  g_param.slot("proportions") = Conversion::CVectorToRcppVector(nbCluster_,gParam->getTabProportion());

  // add proportions values
  g_param.slot("mean") = Conversion::CMatrixToRcppMatrix(nbCluster_,nbVariable_Gaussian,gParam->getTabMean());

  // add variances
  // define list of variance matrix
  Rcpp::GenericVector varianceMatrices(nbCluster_);
  /** Array with matrix variance for each class*/
  XEM::Matrix** matrixVariance = gParam->getTabSigma();
  // loop over clusters
  for (int i=0; i<nbCluster_; i++)
  {
    varianceMatrices[i] = Conversion::CMatrixToRcppMatrix(nbVariable_Gaussian,nbVariable_Gaussian,matrixVariance[i]->storeToArray());
  }
  // add variances to S4 object
  g_param.slot("variance") = varianceMatrices;
  g_param.slot("nbFreeParam") = gParam->getFreeParameter();
  //Fill Binary Parameter slots
  // add proportions values
  b_param.slot("proportions") = Conversion::CVectorToRcppVector(nbCluster_,bParam->getTabProportion());

  // add means values
  b_param.slot("center") = Conversion::CMatrixToRcppMatrixForInt(nbCluster_,nbVariable_Binary,bParam->getTabCenter());

  //add factor
  b_param.slot("factor") = Conversion::CVectorToRcppVectorForInt(nbVariable_Binary,bParam->getTabNbModality());


  //-------------------
  // get scatter
  //-------------------
  // get pointer to scatter
  double *** scatter = bParam->scatterToArray();
  // get tab of modalities
  int64_t* tabNbModality = bParam->getTabNbModality();
  // get maximum number of modality
  int64_t max = *max_element(tabNbModality,tabNbModality+nbVariable_Binary);

  //VectorMatrix of NumericMatrix
  Rcpp::GenericVector vectorOutput(nbCluster_);

  // loop over clusters
  for (int k=0; k<nbCluster_; k++){
    //NumericMatrix for matrix
    Rcpp::NumericMatrix matrixOutput(nbVariable_Binary,max);
    // loop over variables
    for(int j=0; j<nbVariable_Binary; j++){
      // loop over modalities
      for (int h=0; h<tabNbModality[j]; h++) {
        matrixOutput(j,h) = scatter[k][j][h];
      }
    }
    vectorOutput(k) = matrixOutput;
  }
  // add scatters
  b_param.slot("scatter") = vectorOutput;
  b_param.slot("nbFreeParam") = bParam->getFreeParameter();

  param.slot("proportions") = g_param.slot("proportions");
  param.slot("nbFreeParam") = cParam->getFreeParameter();
  param.slot("g_parameter") = g_param;
  param.slot("m_parameter") = b_param;
  xem_.slot("parameters") = param;
}
