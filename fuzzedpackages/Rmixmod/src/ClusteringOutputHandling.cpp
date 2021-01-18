//
//  ClusteringOutputHandling.cpp
//  Rmixmod
//
//  Created by Rémi Lebret on 06/02/12.
//  Copyright (c) 2012 CNRS. All rights reserved.
//

#include "Conversion.h"
#include "ClusteringOutputHandling.h"

#include "mixmod/Clustering/ClusteringModelOutput.h"
#include "mixmod/Kernel/IO/LabelDescription.h"
#include "mixmod/Kernel/IO/Label.h"
#include "mixmod/Kernel/IO/ProbaDescription.h"
#include "mixmod/Kernel/IO/Proba.h"

// constructor

ClusteringOutputHandling::ClusteringOutputHandling( XEM::ClusteringModelOutput* cMOutput
                                                  , Rcpp::S4& xem
                                                  , const XEM::DataType dataType
                                                  , Rcpp::CharacterVector const & Rcriterion
                                                  )
                                                  : OutputHandling(cMOutput, xem, dataType)
{  
  // get criterion
  std::vector<std::string> criterionName = Rcpp::as< std::vector<std::string> >(Rcriterion);
  // add criterion name
  xem_.slot("criterion") = criterionName;
  
  // fill other slot only if no error
  if ( dynamic_cast<XEM::Exception&>(cMOutput->getStrategyRunError()) == XEM::NOERROR ){
    
    // declare a vector of criterion
    std::vector<double> criterionValue;
    // loop over criterion names
    for (unsigned int i=0; i<criterionName.size(); i++){
      // BIC criterion
      if (criterionName[i] == "BIC"){ 
        criterionValue.push_back(cMOutput->getCriterionOutput(XEM::BIC).getValue());
      }
      // ICL Criterion
      else if(criterionName[i] == "ICL"){ 
        criterionValue.push_back(cMOutput->getCriterionOutput(XEM::ICL).getValue());
      }
        // NEC Criterion
      else if(criterionName[i] == "NEC"){ 
        criterionValue.push_back(cMOutput->getCriterionOutput(XEM::NEC).getValue());
      }
    }
  
    // add criterion value
    xem_.slot("criterionValue") = criterionValue;
    
    // add labels
    xem_.slot("partition") = Conversion::VectorToRcppVectorForInt(MOutput_->getLabelDescription()->getLabel()->getLabel());
    
    //add proba
    xem_.slot("proba") = Conversion::XEMMatrixToRcppMatrix(MOutput_->getProbaDescription()->getProba()->getProba());
  }
}

ClusteringOutputHandling::ClusteringOutputHandling( XEM::ClusteringModelOutput* cMOutput
                                                    , Rcpp::S4& xem
                                                    , const XEM::DataType dataType
                                                    , std::vector <XEM::CriterionName> const & iCriterion
                                                    )
  : OutputHandling(cMOutput, xem, dataType)
{  
  // get criterion
  std::vector<std::string> criterionName; //= Rcpp::as< std::vector<std::string> >(Rcriterion);
  
  // fill other slot only if no error
  if ( dynamic_cast<XEM::Exception&>(cMOutput->getStrategyRunError()) == XEM::NOERROR ){
    
    // declare a vector of criterion
    std::vector<double> criterionValue;
    // loop over criterion names
    for (unsigned int i=0; i<iCriterion.size(); i++){
      criterionValue.push_back(cMOutput->getCriterionOutput(iCriterion[i]).getValue());
      criterionName.push_back(CriterionNameToString(iCriterion[i]));
    }
  // add criterion name
  xem_.slot("criterion") = criterionName;
  
    // add criterion value
    xem_.slot("criterionValue") = criterionValue;
    
    // add labels
    xem_.slot("partition") = Conversion::VectorToRcppVectorForInt(MOutput_->getLabelDescription()->getLabel()->getLabel());
    
    //add proba
    xem_.slot("proba") = Conversion::XEMMatrixToRcppMatrix(MOutput_->getProbaDescription()->getProba()->getProba());
  }
}

// destructor
ClusteringOutputHandling::~ClusteringOutputHandling()
{ }
