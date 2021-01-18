//
//  LearnOutputHandling.cpp
//  Rmixmod
//
//  Created by Rémi Lebret on 06/02/12.
//  Copyright (c) 2012 CNRS. All rights reserved.
//

#include "Conversion.h"
#include "LearnOutputHandling.h"

#include "mixmod/DiscriminantAnalysis/Learn/LearnModelOutput.h"
#include "mixmod/Kernel/IO/LabelDescription.h"
#include "mixmod/Kernel/IO/Label.h"
#include "mixmod/Kernel/IO/ProbaDescription.h"
#include "mixmod/Kernel/IO/Proba.h"
#include "mixmod/Utilities/Util.h"

// constructor
LearnOutputHandling::LearnOutputHandling( XEM::LearnModelOutput* lMOutput
                                        , Rcpp::S4& xem
                                        , const XEM::DataType dataType
                                        , Rcpp::CharacterVector const & Rcriterion
                                        , std::vector<int64_t> labels
                                        )
                                        : OutputHandling(lMOutput, xem, dataType)
{
  // get criterion
  std::vector<std::string> criterionName = Rcpp::as< std::vector<std::string> >(Rcriterion);
  // add criterion name
  xem_.slot("criterion") = criterionName;
  
  // fill other slot only if no error
  if ( dynamic_cast<XEM::Exception&>(lMOutput->getStrategyRunError()) == XEM::NOERROR ){
	
    // declare a vector of criterion
    std::vector<double> criterionValue;
    // loop over criterion names
    for (unsigned int i=0; i<criterionName.size(); i++){
      // BIC criterion
      if (criterionName[i] == "BIC"){ 
        criterionValue.push_back(lMOutput->getCriterionOutput(XEM::BIC).getValue());
      }
      // CV Criterion
      else if(criterionName[i] == "CV"){ 
        criterionValue.push_back(1-lMOutput->getCriterionOutput(XEM::CV).getValue());
        xem_.slot("CVLabel") = Conversion::VectorToRcppVectorForInt(lMOutput->getCVLabel()->getLabel()->getLabel());
        // define a classification tab
        int64_t** tab = lMOutput->getCVLabel()->getLabel()->getClassificationTab(labels, nbCluster_);
        xem_.slot("CVClassification") = Conversion::CMatrixToRcppMatrixForInt(nbCluster_, nbCluster_, tab);
        // release memory
        for ( int i=0; i<nbCluster_; i++ ) delete [] tab[i];
        delete [] tab;
      }
    }
    
    // add criterion value
    xem_.slot("criterionValue") = criterionValue;    
    
    // add partition
    // get labels size
    const int size = labels.size();
    Rcpp::NumericVector Rlabels(size);
    for ( int i=0; i<size; i++ ) Rlabels[i]=labels[i];
    xem_.slot("partition") = Rlabels;
    
    // add MAP values
    // define a classification tab
    int64_t** tab = lMOutput->getLabelDescription()->getLabel()->getClassificationTab(labels, nbCluster_);
    xem_.slot("MAPClassification") = Conversion::CMatrixToRcppMatrixForInt(nbCluster_, nbCluster_, tab);
    // release memory
    for ( int i=0; i<nbCluster_; i++ ) delete [] tab[i];
    delete [] tab;
    // add MAP error rate
    xem_.slot("MAPErrorRate") = lMOutput->getLabelDescription()->getLabel()->getErrorRate(labels);
  }
}
// constructor

LearnOutputHandling::LearnOutputHandling( XEM::LearnModelOutput* lMOutput
                                        , Rcpp::S4& xem
                                        , const XEM::DataType dataType
                                        , std::vector <XEM::CriterionName> const & iCriterion
                                        , std::vector<int64_t> labels
                                        )
                                        : OutputHandling(lMOutput, xem, dataType)
{
  // get criterion
  std::vector<std::string> criterionName; //= Rcpp::as< std::vector<std::string> >(Rcriterion);
  
  // fill other slot only if no error
  if ( dynamic_cast<XEM::Exception&>(lMOutput->getStrategyRunError()) == XEM::NOERROR ){
	
    // declare a vector of criterion
    std::vector<double> criterionValue;
    // loop over criterion names
    for (unsigned int i=0; i<iCriterion.size(); i++){
      criterionValue.push_back(lMOutput->getCriterionOutput(iCriterion[i]).getValue());
      criterionName.push_back(CriterionNameToString(iCriterion[i]));
    }    
    // add criterion name
    xem_.slot("criterion") = criterionName;
    // add criterion value
    xem_.slot("criterionValue") = criterionValue;    
    
    // add partition
    // get labels size
    const int size = labels.size();
    Rcpp::NumericVector Rlabels(size);
    for ( int i=0; i<size; i++ ) Rlabels[i]=labels[i];
    xem_.slot("partition") = Rlabels;
    
    // add MAP values
    // define a classification tab
    int64_t** tab = lMOutput->getLabelDescription()->getLabel()->getClassificationTab(labels, nbCluster_);
    xem_.slot("MAPClassification") = Conversion::CMatrixToRcppMatrixForInt(nbCluster_, nbCluster_, tab);
    // release memory
    for ( int i=0; i<nbCluster_; i++ ) delete [] tab[i];
    delete [] tab;
    // add MAP error rate
    xem_.slot("MAPErrorRate") = lMOutput->getLabelDescription()->getLabel()->getErrorRate(labels);
  }
}


// destructor
LearnOutputHandling::~LearnOutputHandling()
{ }
