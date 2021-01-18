//
//  ClusteringInputHandling.cpp
//  Rmixmod
//
//  Created by Rémi Lebret on 06/02/12.
//  Copyright (c) 2012 CNRS. All rights reserved.
//

#include "ClusteringInputHandling.h"
#include <vector>

#include "mixmod/Clustering/ClusteringInput.h"
#include "mixmod/Clustering/ClusteringStrategy.h"
#include "mixmod/Clustering/ClusteringStrategyInit.h"
#include "mixmod/Kernel/Model/ModelType.h"
#include "mixmod/Kernel/Algo/Algo.h"

#include "mixmod/Utilities/Util.h"
#include "Conversion.h"
//default constructor
ClusteringInputHandling::ClusteringInputHandling( XEM::ClusteringInput* cInput
                                                , Rcpp::S4& algoOptions
                                                ) 
                                                : InputHandling(cInput)
                                                , cStrategy_(cInput->getStrategy())
                                                , algoOptions_(algoOptions)
{ run();}

//destructor
ClusteringInputHandling::~ClusteringInputHandling()
{ }


///setAlgo
void ClusteringInputHandling::setAlgo()
{
  // check if it is a NULL value
  if (Rf_isNull( algoOptions_.slot("algo") )) return;
  // check if it is a NULL value
  if (Rf_isNull( algoOptions_.slot("epsilonInAlgo") )) return;
  // check if it is a NULL value
  if (Rf_isNull( algoOptions_.slot("nbIterationInAlgo") )) return;
  
  // get algo
  std::vector<std::string> algo = Rcpp::as< std::vector<std::string> >(algoOptions_.slot("algo"));
  // get int value
  std::vector<int> nbIterationInAlgo = Rcpp::as< std::vector<int> >(algoOptions_.slot("nbIterationInAlgo"));
  // get double value
  std::vector<double> epsilonInAlgo = Rcpp::as< std::vector<double> >(algoOptions_.slot("epsilonInAlgo"));
  
  // remove default algo
  cStrategy_->removeAlgo(0);

  // loop over algorithm
  for (unsigned int i=0; i<algo.size(); i++){
    // EM Algo
    if (algo[i] == "EM")
    { 
      cStrategy_->addAlgo(XEM::EM);
      // set epsilon value
      if ( (epsilonInAlgo[i]>0) & (epsilonInAlgo[i]<1) )
      { cStrategy_->setAlgoEpsilon(i, epsilonInAlgo[i]); }
    }
    // SEM Algo
    else if (algo[i] == "SEM")
    { cStrategy_->addAlgo(XEM::SEM); }
    // CEM Algo
    else if (algo[i] == "CEM")
    { 
      cStrategy_->addAlgo(XEM::CEM);
      // set epsilon value
      if ( (epsilonInAlgo[i]>0) & (epsilonInAlgo[i]<1) )
      { cStrategy_->setAlgoEpsilon(i, epsilonInAlgo[i]); }
    }
    else{
      //exception because wrong criterionName
      Rcpp::stop("In InputHandling::setAlgo Invalid algo");
    }
    
    // set nbIteration value
    if (nbIterationInAlgo[i] > 0)
    { cStrategy_->setAlgoIteration(i, nbIterationInAlgo[i]); }
  }
}

/// setInitAlgo
void ClusteringInputHandling::setInitAlgo()
{
  // check if it is a NULL value
  if (Rf_isNull( algoOptions_.slot("initMethod") )) return;
  // get initialization algorithm name
  std::string init = Rcpp::as<std::string>(algoOptions_.slot("initMethod"));
  // random initialization
  if (init == "random")
  {
    setNbTry();
    setNbTryInInit();
    cStrategy_->setStrategyInitName(XEM::RANDOM);    
    return;
  }
  else if (init == "smallEM")
  {
    setNbTry();
    setNbTryInInit();
    setNbIterationInInit();
    setEpsilonInInit();
    cStrategy_->setStrategyInitName(XEM::SMALL_EM);
    return;
  }
  else if (init == "CEM")
  {
    setNbTry();
    setNbTryInInit();
    cStrategy_->setStrategyInitName(XEM::CEM_INIT);
    return;
  }
  else if (init == "SEMMax")
  {
    setNbTry();
    setNbIterationInInit();
    cStrategy_->setStrategyInitName(XEM::SEM_MAX);
    return;
  }
  else if (init == "parameter")
  {
    //throw(std::runtime_error("Error in InputHandling::run"
    //                         " parameter initialization is not implemented"));
    //std::cout<<" parameter initialization is not implemented" <<std::endl;
    cStrategy_->setStrategyInitName(XEM::USER);
    Rcpp::S4 userParam(algoOptions_.slot("parameter"));
    XEM::Parameter** paramArray1 = nullptr;
    int nbCluster = 0;
    int pbDim = 0;
    Rcpp::NumericVector rProportions(userParam.slot("proportions"));
    std::unique_ptr<double[]> proportions(new double[rProportions.size()]);
    //Rcpp::NumericVector::iterator it;
    Rcpp::NumericVector::iterator it;
    int i; 
    for(it = rProportions.begin(), i=0; it != rProportions.end(); ++it, i++)
      { proportions[i] = (double)(*it);}

    if(cInput_->getDataType() == XEM::QuantitativeData) {
      Rcpp::NumericMatrix rMean(SEXP(userParam.slot("mean")));
      Rcpp::List rVariance(userParam.slot("variance"));
      nbCluster = rMean.nrow();
      pbDim = rMean.ncol();
      std::unique_ptr<double*[], XEM::TabDeleter<double>> cMean(Conversion::RcppMatrixToC2DArray(rMean),XEM::TabDeleter<double>(nbCluster));      
      std::unique_ptr<double**[], XEM::VectTabDeleter<double>> cVariance(Conversion::RcppListOfMatrixToC3DArray(rVariance), XEM::VectTabDeleter<double>(nbCluster, pbDim)); 
      paramArray1 = new XEM::Parameter*[1];
      paramArray1[0] = new XEM::GaussianGeneralParameter(nbCluster, pbDim, cInput_->getModelType()[0], proportions.get(), cMean.get(), cVariance.get() );
      cStrategy_->setTabInitParameter (paramArray1, 1);
    } else if(cInput_->getDataType() == XEM::QualitativeData) {
      Rcpp::IntegerVector rFactor(userParam.slot("factor"));
      std::unique_ptr<int64_t[]> cFactor(new int64_t[rFactor.size()]);
      //Rcpp::NumericVector::iterator it;
      Rcpp::IntegerVector::iterator it2;
      int k; 
      for(it2 = rFactor.begin(), k=0; it2 != rFactor.end(); ++it2, k++)
        { cFactor[k] = (int64_t)(*it2);}
        
      Rcpp::NumericMatrix rCenter(SEXP(userParam.slot("center")));
      Rcpp::List rScatter(userParam.slot("scatter"));
      nbCluster = rCenter.nrow();
      pbDim = rCenter.ncol();
      std::unique_ptr<double*[], XEM::TabDeleter<double>> cCenter(Conversion::RcppMatrixToC2DArray(rCenter),XEM::TabDeleter<double>(nbCluster));      
      std::unique_ptr<double**[], XEM::VectTabDeleter<double>> cScatter(Conversion::RcppListOfMatrixToC3DArray(rScatter), XEM::VectTabDeleter<double>(nbCluster, pbDim)); 
      paramArray1 = new XEM::Parameter*[1];
      paramArray1[0] = new XEM::BinaryEkjhParameter(nbCluster, pbDim, cInput_->getModelType()[0], cFactor.get(), proportions.get(), cCenter.get(), cScatter.get() );
      cStrategy_->setTabInitParameter (paramArray1, 1);
    } else {
      Rcpp::stop("Error in InputHandling::run"
                             " parameter initialization for heterogeneous data  his not implemented");      
    }
  }
  else if (init == "partition")
  {
      Rcpp::IntegerVector rLabels(algoOptions_.slot("labels"));
      int64_t nbSample = cInput_->getNbSample();
      if(rLabels.size() != nbSample){
        Rcpp::stop("Error in InputHandling::run"
                                 "The size of labels vector doesn't match "
                                 "the number of samples ");      
      }
      std::vector<int64_t> cLabelArray;
      Rcpp::IntegerVector::iterator it;
      int k; 
      for(it = rLabels.begin(), k=0; it != rLabels.end(); ++it, k++)
        { cLabelArray.push_back((int64_t)(*it));}
      cStrategy_->setStrategyInitName(XEM::USER_PARTITION);
      std::unique_ptr<XEM::Label> cLabel(new XEM::Label(nbSample));
      cLabel->setLabel(cLabelArray, nbSample);
      XEM::Partition* cPartition = new XEM::Partition(cLabel.get(),cInput_->getNbCluster(0));
      cStrategy_->setInitPartition(cPartition, 0);
        
  }
  //
  //  XEMStrategyInitName initName = cInput_->getStrategy()->getStrategyInit()->getStrategyInitName();
  //  cInput_->getStrategy()->setStrategyInitName(initName);
}

/// setNbTry
void ClusteringInputHandling::setNbTry()
{
  // check if it is a NULL value
  if (Rf_isNull( algoOptions_.slot("nbTry") )) return;
  // get int value
  int64_t nbTry = (int64_t)Rcpp::as<int>(algoOptions_.slot("nbTry"));
  // set value
  if (nbTry > 0)
  { cStrategy_->setNbTry(nbTry); }
}

/// setInitParameter
void ClusteringInputHandling::setInitParameter(XEM::Parameter** parameter)
{
  //get initName
  XEM::StrategyInitName initName = cStrategy_->getStrategyInit()->getStrategyInitName();
  //get nbInitParameter value
  int64_t nbInitParameter = cStrategy_->getStrategyInit()->getNbInitParameter();
  //get parameter
  //XEMParameter** parameter = cInput_->getStrategy()->getStrategyInit()->getTabInitParameter();
  //get nbCluster value
  int64_t nbCluster = cInput_->getNbCluster(0);
  if (initName == XEM::USER)
  {
    if (nbCluster == 1)
    {
      cStrategy_->setTabInitParameter(parameter, nbInitParameter);
    }
  }
}

/// setInitPartition
void ClusteringInputHandling::setInitPartition( XEM::Partition ** tabPartition
                                     , int64_t nbPartition
                                     )
{ cStrategy_->setTabPartition(tabPartition, nbPartition);}

/// setNbTryInInit
void ClusteringInputHandling::setNbTryInInit()
{
  // check if it is a NULL value
  if (Rf_isNull( algoOptions_.slot("nbTryInInit"))) return;
  // get int value
  int64_t nbTryInInit = (int64_t)Rcpp::as<int>(algoOptions_.slot("nbTryInInit"));
  // set value
  if (nbTryInInit > 0)
  { cStrategy_->setNbTryInInit(nbTryInInit); }
}

///setNbIterationInInit
void ClusteringInputHandling::setNbIterationInInit()
{
  // check if it is a NULL value
  if (Rf_isNull( algoOptions_.slot("nbIterationInInit") )) return;
  // get int value
  int64_t nbIterationInInit = (int64_t)Rcpp::as<int>(algoOptions_.slot("nbIterationInInit"));
  // set value
  if (nbIterationInInit > 0)
  { cStrategy_->setNbIterationInInit(nbIterationInInit); }
}

///setEpsilonInInit
void ClusteringInputHandling::setEpsilonInInit()
{
  // check if it is a NULL value
  if (Rf_isNull( algoOptions_.slot("epsilonInInit") )) return;
  // get double value
  double epsilonInInit = Rcpp::as<double>(algoOptions_.slot("epsilonInInit"));
  // set value
  if (epsilonInInit > 0)
  { cStrategy_->setEpsilonInInit(epsilonInInit); }
}

// run method
void ClusteringInputHandling::run()
{
  try
  {
    // set Algo name
    setAlgo();
    // set init Algo name
    setInitAlgo();
  }
  catch(...)
  {
    Rcpp::stop("InputHandling::run get a mixmod error");
  }
}

