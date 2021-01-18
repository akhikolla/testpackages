/**
 * Project:  Rmixmod
 * created on: 4 avr. 2011
 * Purpose:  Create the main for the mixmod call.
 * Author:   iovleff, serge.iovleff@stkpp.org
 *
 **/

/***************************************************************************
                             ClusteringMain.cpp  description
                             ------------------------
    copyright            : (C) MIXMOD Team - 2001-2003-2004-2005-2006-2007-2008-2009
    email                : mixmod@univ-fcomte.fr
 ***************************************************************************/

/***************************************************************************
    This file is part of MIXMOD

    MIXMOD is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MIXMOD is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MIXMOD.  If not, see <http://www.gnu.org/licenses/>.

 ***************************************************************************/
/** @file ClusteringMain.cpp
 *  @brief In this file we implement the wrapper .
 **/

/*
 *
 ***** File "rcpp_hello_world.h"
 *
 *
 * #ifndef _Rmixmod_RCPP_HELLO_WORLD_H
 * #define _Rmixmod_RCPP_HELLO_WORLD_H
 *
 * #include <Rcpp.h>
 *
 * note : RcppExport is an alias to `extern "C"` defined by Rcpp.
 *
 * It gives C calling convention to the rcpp_hello_world function so that
 * it can be called from .Call in R. Otherwise, the C++ compiler mangles the
 * name of the function and .Call can't find it.
 *
 * It is only useful to use RcppExport when the function is intended to be called
 * by .Call. See the thread http://thread.gmane.org/gmane.comp.lang.r.rcpp/649/focus=672
 * on Rcpp-devel for a misuse of RcppExport
 *
 * RcppExport SEXP rcpp_hello_world() ;
 *
 * #endif
 *
 *
 ***** File "rcpp_hello_world.cpp"
 *
 * #include "rcpp_hello_world.h"
 *
 * SEXP rcpp_hello_world(){
 *   using namespace Rcpp ;
 *
 *   CharacterVector x = CharacterVector::create( "foo", "bar" )  ;
 *   NumericVector y   = NumericVector::create( 0.0, 1.0 ) ;
 *   List z            = List::create( x, y ) ;
 *
 *   return z ;
 * }
 *
 */

#include <vector>
#include <string>

#include "mixmod/Utilities/Util.h"
#include "mixmod/Utilities/exceptions/OtherException.h"
#include "mixmod/Clustering/ClusteringMain.h"
#include "mixmod/Clustering/ClusteringInput.h"
#include "mixmod/Clustering/ClusteringOutput.h"
#include "mixmod/Clustering/ClusteringModelOutput.h"
#include "mixmod/Clustering/ClusteringStrategy.h"
#include "mixmod/Clustering/ClusteringStrategyInit.h"
#include "mixmod/Kernel/IO/ParameterDescription.h"
#include "mixmod/Kernel/Parameter/GaussianParameter.h"
#include "mixmod/Kernel/IO/GaussianData.h"
#include "mixmod/Kernel/IO/BinaryData.h"
#include "mixmod/Kernel/IO/CompositeData.h"
#include "mixmod/Kernel/Algo/Algo.h"
#include "mixmod/Matrix/Matrix.h"
#ifdef RMIXMOD_XML
#include "mixmod_iostream/IOStreamUtil.h"
//#include "mixmod_iostream/DomClusteringProject.h"
#endif
#include <Rcpp.h>

#include "Conversion.h"
#include "ClusteringInputHandling.h"
#include "ClusteringOutputHandling.h"


/** This is the main method doing the interface between R and Mixmod for Clustering.
 *  The method will create a matrix in the format of mixmod and copy the data
 *  inside
 *
 * @param xem R S4 object containing data and options for clustering
 */
RcppExport SEXP clusteringMain( SEXP xem )
{
  BEGIN_RCPP
  std::unique_ptr<XEM::ClusteringMain> cMain; //(cInput.get());
 std::unique_ptr<XEM::ClusteringInput> cInput; //(new XEM::ClusteringInput(clusterArray, *dataDescription)); 
 // wrap S4 object
  Rcpp::S4 mixmodClustering(xem);
  // define a new MixmodResults object
  Rcpp::S4 iResults(mixmodClustering.slot("bestResult"));
  //Rcpp::StringVector XmlIn(mixmodClustering.slot("xmlIn"));
  Rcpp::S4 XmlInObj(mixmodClustering.slot("xmlIn"));
  bool conversionOnly = false; // actually conversionOnly could be true only when an xmlIn in given
  XEM::IoMode iomode = XEM::IoMode::NUMERIC;  
  std::string xmlInput = "";
  //Rcpp::S4 sfile = XmlInObj.slot("file");
  std::vector<std::string>  xmlv =  Rcpp::as< std::vector<std::string> >(XmlInObj.slot("file"));
  if(xmlv.size()>0){
    xmlInput = xmlv[0];
  }
  //if(XmlIn.size()>0)  xmlInput = Rcpp::as< std::vector<std::string> >(XmlIn)[0];
  Rcpp::StringVector XmlOut(mixmodClustering.slot("xmlOut"));
  std::string xmlOutput = "";
  if(XmlOut.size()>0)  xmlOutput = Rcpp::as< std::vector<std::string> >(XmlOut)[0];
  //std::string xmlInput = Rcpp::as<std::string>(XmlIn);
  //std::cout << "xmlIn:::::: " << XmlIn << std::endl;
  int seed = Rcpp::as<int>(mixmodClustering.slot("seed"));
  int trace = Rcpp::as<int>(mixmodClustering.slot("trace"));
  int massiccc = Rcpp::as<int>(mixmodClustering.slot("massiccc"));    
  if(xmlInput==""){
    // wrap data in Rcpp matrix
    Rcpp::NumericMatrix RData(SEXP(mixmodClustering.slot("data"))); // creates Rcpp matrix from SEXP
    // wrap clusters in Rcpp
    Rcpp::NumericVector RnbCluster(mixmodClustering.slot("nbCluster")); // keep a copy (as the classic version does)
    // wrap partition matrix in Rcpp matrix
    Rcpp::NumericVector RPartition(mixmodClustering.slot("knownLabels"));
    // wrap list algoOptions
    Rcpp::S4 RalgoOptions(mixmodClustering.slot("strategy"));
    // wrap criterion
    Rcpp::CharacterVector Rcriterion(mixmodClustering.slot("criterion"));
    // wrap models
    Rcpp::S4 Rmodel(mixmodClustering.slot("models"));
    // wrap type of data
    Rcpp::StringVector Rtype(mixmodClustering.slot("dataType"));
    // wrap weight in Rcpp
    Rcpp::NumericVector Rweight(mixmodClustering.slot("weight")); 
    // wrap factor in Rcpp
    Rcpp::NumericVector Rfactor(mixmodClustering.slot("factor"));

    // gaussian, binary or heterogeneous ?
    XEM::DataType dataType;
    if(Rcpp::as<std::string>(Rtype) == std::string("quantitative"))
      dataType = XEM::QuantitativeData;
    else if(Rcpp::as<std::string>(Rtype) == std::string("qualitative"))
      dataType = XEM::QualitativeData;
    else
      dataType = XEM::HeterogeneousData;
    // data descritor
    XEM::DataDescription* dataDescription;
    XEM::GaussianData* gData(0);
    XEM::BinaryData* bData(0);
    XEM::CompositeData* cData(0);
    switch (dataType) {
    case XEM::QuantitativeData:
      /*===============================================*/
      /*  Create XEM gaussian data set and description */
      gData = Conversion::DataToXemGaussianData(RData);
      dataDescription = new XEM::DataDescription(gData);
      break;
    case XEM::QualitativeData:
      /*===============================================*/
      /* Create XEM binary data set and description */
      bData = Conversion::DataToXemBinaryData(RData, Rfactor);
      dataDescription = new XEM::DataDescription(bData);
      break;
    case XEM::HeterogeneousData:
      /*===============================================*/
      /* Create XEM binary data set and description */
      cData = Conversion::DataToXemCompositeData(RData, Rfactor);
      dataDescription = new XEM::DataDescription(cData);
      break;
    default:
      break;
    }

    /*===============================================*/
    /*  Create XEM cluster array */
    std::vector<int64_t> clusterArray;
    Rcpp::NumericVector::iterator it;
    for( it = RnbCluster.begin(); it != RnbCluster.end(); ++it)
      { clusterArray.push_back((int64_t)(*it));}
  
    /*===============================================*/
    /* Initialize input parameters in MIXMOD         */
    // create XEMClusteringInput
    //XEM::ClusteringInput * cInput =  new XEM::ClusteringInput(clusterArray, *dataDescription);
    //std::unique_ptr<XEM::ClusteringInput> cInput(new XEM::ClusteringInput(clusterArray, *dataDescription));
    cInput.reset(new XEM::ClusteringInput(clusterArray, *dataDescription));
    // initialize the parameters using user defined values (see Algo.R)
    ClusteringInputHandling initInput(cInput.get(), RalgoOptions);

    // initialize criterion name (BIC, AIC, ...)
    initInput.setCriterionName(Rcriterion);

    // initialize the model (gaussian_...)
    initInput.setModel(Rmodel);

    // set weight
    initInput.setWeight(Rweight);
  
    // set knownPartition
    initInput.setKnownPartition(RPartition);
  
    /* XEM will check if the option are corrects (should be the case) */
    cInput->finalize();
    cMain.reset(new XEM::ClusteringMain(cInput.get()));
    // check if it is a NULL value
    if (!Rf_isNull( RalgoOptions.slot("seed") )){
      //seed = -1;
      //}else{
      int strategy_seed = (int)Rcpp::as<int>(RalgoOptions.slot("seed"));
      if(seed!=-1 && strategy_seed!=-1 && seed!=strategy_seed){
        THROW(XEM::OtherException, XEM::internalMixmodError);
      }
      seed = strategy_seed!=-1?strategy_seed:seed;
    }
    // release memory
    //if (dataDescription) delete dataDescription;
    //if (gData) delete gData;
    //if (bData) delete bData;
    //if (cData) delete cData;

  } else { // i.e. with xml
#ifdef RMIXMOD_XML
    //std::cout<< "works with xmlIn"<< std::endl;
    std::string numFormat = Rcpp::as<std::string>(XmlInObj.slot("numFormat"));
    if(numFormat == "hexBinary") iomode = XEM::IoMode::BINARY;
    conversionOnly = Rcpp::as<bool>(XmlInObj.slot("conversionOnly"));
    cMain.reset(XEM::IStream(xmlInput, XEM::IOStreamFormat::XML, false, iomode));
    //cMain.reset(XEM::IStream_XML(xmlInput, false));
    //cMain.reset(XEM::IStream(xmlInput, XEM::IOStreamFormat::XML, false, XEM::IoMode::BINARY));
    cInput.reset(dynamic_cast<XEM::ClusteringInput*>(cMain->getInput()));
    Rcpp::S4 params(iResults.slot("parameters"));
    switch(cInput->getDataType()){
    case XEM::QuantitativeData:
      iResults.slot("parameters") = params.slot("g_parameter");
      mixmodClustering.slot("dataType") = Rcpp::StringVector("quantitative");
      mixmodClustering.slot("models") = XmlInObj.slot("gaussianModel");
      break;
    case XEM::QualitativeData:
      iResults.slot("parameters") = params.slot("m_parameter");
      mixmodClustering.slot("dataType") = Rcpp::StringVector("qualitative");
      mixmodClustering.slot("models") = XmlInObj.slot("multinomialModel");
      break;
    case XEM::HeterogeneousData:
      mixmodClustering.slot("dataType") = Rcpp::StringVector("composite");
      mixmodClustering.slot("models") = XmlInObj.slot("compositeModel");
      break;
    }
    Rcpp::NumericVector vectorNbCluster(cInput->getNbClusterSize());
    for(int i=0;i<cInput->getNbClusterSize();i++) vectorNbCluster(i) = cInput->getNbCluster(i);
    mixmodClustering.slot("nbCluster") = vectorNbCluster;
    Rcpp::StringVector vectorCrit(cInput->getNbCriterion());
    for(int i=0;i<cInput->getNbCriterion();i++) vectorCrit(i) = CriterionNameToString(cInput->getCriterionName(i));
    mixmodClustering.slot("criterion") =  vectorCrit;
    mixmodClustering.slot("nbSample") = Rcpp::NumericVector(1,cInput->getNbSample());
    mixmodClustering.slot("nbVariable") = Rcpp::NumericVector(1,cInput->getPbDimension());
    Rcpp::S4 strategy(mixmodClustering.slot("strategy"));
    Rcpp::StringVector vectorAlgo(cInput->getStrategy()->getNbAlgo());
    for(int i=0;i<cInput->getStrategy()->getNbAlgo();i++) vectorAlgo(i) = AlgoNameToString(cInput->getStrategy()->getAlgo(i)->getAlgoName());
    strategy.slot("algo") = vectorAlgo;
    strategy.slot("nbTry") = Rcpp::NumericVector(1,cInput->getStrategy()->getNbTry());
    strategy.slot("initMethod") = Rcpp::StringVector(StrategyInitNameToStringApp(cInput->getStrategy()->getStrategyInit()->getStrategyInitName()));
    strategy.slot("nbTryInInit") = Rcpp::NumericVector(1,cInput->getStrategy()->getNbTryInInit());
    strategy.slot("nbIterationInInit") = Rcpp::NumericVector(1,cInput->getStrategy()->getNbIterationInInit());
    strategy.slot("epsilonInInit") = Rcpp::NumericVector(1,cInput->getStrategy()->getEpsilonInInit());
    strategy.slot("nbIterationInAlgo") = Rcpp::NumericVector(1,cInput->getStrategy()->getAlgo(0)->getNbIteration());
    strategy.slot("epsilonInAlgo") = Rcpp::NumericVector(1,cInput->getStrategy()->getAlgo(0)->getEpsilon());
    if(seed>=0){
      strategy.slot("seed") = mixmodClustering.slot("seed");
    }
    mixmodClustering.slot("strategy") = strategy;
    // labels ...
    if (cInput->getKnownLabelDescription()) {
      if (cInput->getKnownPartition()) {
        // knowPartition and knownLabelDescription can't be both not NULL !
        THROW(XEM::OtherException, XEM::internalMixmodError);
      }
      const std::vector<int64_t> & labels = cInput->getKnownLabelDescription()->getLabel()->getLabel();
      Rcpp::NumericVector labelsNV(labels.size());
      for(int64_t i=0;i<labels.size();++i){
        labelsNV(i) = labels[i];
      }

      mixmodClustering.slot("knownLabels") = labelsNV;
    }
    // ... or partition
    if (cInput->getKnownPartition()) {
      //std::vector<int64_t> labels(cInput->getNbSample());
      Rcpp::NumericVector labelsNV(cInput->getNbSample());
      int64_t** tv = cInput->getKnownPartition()->getTabValue();
      for(int64_t i=0;i<cInput->getNbSample();i++){
        int64_t label_i = 0;
        for(int64_t j=0;j<cInput->getNbCluster()[0];j++){
          label_i += tv[i][j]*(j+1);
        }
        labelsNV(i) = label_i;
      }
      mixmodClustering.slot("knownLabels") = labelsNV;
    }
    // weight
    if(cInput->getData()->getWeight()){
      Rcpp::NumericVector weightNV(cInput->getNbSample());
      for(int64_t i=0;i<weightNV.size();++i){
        weightNV(i) = cInput->getData()->getWeightI(i);
      }
      mixmodClustering.slot("weight") = weightNV;
    }

    // models
    std::vector<XEM::ModelType*> modelTypeV = cInput->getModelType();
    
    Rcpp::StringVector modelsV(modelTypeV.size());
    for(int64_t i=0;i<modelTypeV.size();++i){
      modelsV(i) = XEM::ModelNameToString(modelTypeV[i]->getModelName());
    }
    Rcpp::S4 slotModels(mixmodClustering.slot("models"));
    slotModels.slot("listModels") = modelsV;

    /*
#c@strategy     c@data         c@factor       c@knownLabels  c@nbVariable   c@criterion    c@error        c@xmlIn
    #c@bestResult   c@dataType     c@nbCluster    c@weight       c@nbSample     c@models       c@results      c@xmlOut
    > c@strategy@
c@strategy@algo               c@strategy@epsilonInInit
c@strategy@nbTry              c@strategy@epsilonInAlgo
c@strategy@initMethod         c@strategy@seed
c@strategy@nbTryInInit        c@strategy@parameter
c@strategy@nbIterationInInit  c@strategy@labels
c@strategy@nbIterationInAlgo  
    */
    //http://www.hep.by/gnu/r-patched/r-exts/R-exts_11.html
    //return mixmodClustering;
#else
    THROW(XEM::OtherException, XEM::xmlFeaturesNotAvailable);
#endif    
  } //end else
//  cInput->edit(std::cout);
//  std::cout << "Input finished" << std::endl;
  /*===============================================*/
  /* Computation of the estimates */
  // XEMClusteringMain
  //XEM::ClusteringMain cMain(cInput.get());
  
  //cout<<"seed is..... "<<seed<<endl;
  //cout<<"trace is..... "<<trace<<endl;  
  // xmain run
  if(!conversionOnly){
    try{
      cMain->run(seed, XEM::IoMode::NUMERIC, trace, massiccc);
    }
    catch(XEM::OtherException & e){
      // The following "if" (based on what()) is FRAGILE! SO:
      // TO DO ASAP (when mixmodLib were available):
      // * activate (uncomment) if(e.getErrorType()...) instead of:
      /// if(std::string(e.what())=...)
      if(std::string(e.what())=="All models got errors"){
      //if(e.getErrorType()==XEM::AllModelsGotErros){
        mixmodClustering.slot("error") = true;
        return mixmodClustering;
      }
      throw;
    }
  } //else {
    //cout<<"********************************** CONVERSION ONLY ********************************"<<endl;  
  //}
  if(xmlOutput!=""){
#ifdef RMIXMOD_XML
    XEM::OStream(xmlOutput, XEM::IOStreamFormat::XML, cMain.get(), XEM::IoMode::NUMERIC);
#else
    THROW(XEM::OtherException, XEM::xmlFeaturesNotAvailable);
#endif    
  }
  //std::cout << "Run finished" << std::endl;
  /*===============================================*/
  // get XEMClusteringOutput object from cMain
  XEM::ClusteringOutput * cOutput = cMain->getOutput();
  // Treatment : sort with the first criterion (there is only one criterion)  
  cOutput->sort(cMain->getInput()->getCriterionName(0));
  //std::cout << "Sort finished" << std::endl;
  
  
  /*===============================================*/
  /* get output parameters from MIXMOD             */
  
  if ( cOutput->atLeastOneEstimationNoError() )
  {
    // create a list of results
    Rcpp::List resList(cOutput->getNbClusteringModelOutput());
    // loop over all estimation
    for ( int i=0; i<cOutput->getNbClusteringModelOutput(); i++ )
      {
        // create the output object for R
        //ClusteringOutputHandling(cOutput->getClusteringModelOutput(i), iResults,cMain->getInput()->getDataType(), Rcriterion);
        ClusteringOutputHandling(cOutput->getClusteringModelOutput(i), iResults,cMain->getInput()->getDataType(), cMain->getInput()->getCriterionName());      
        // add those results to the list
        resList[i] = Rcpp::clone(iResults);
      } 
    // add all results
    mixmodClustering.slot("results") = resList;
    // TODO: if no criterion...
    
    // add the best results
    //ClusteringOutputHandling(cOutput->getClusteringModelOutput().front(), iResults,cMain->getInput()->getDataType() , Rcriterion);
    ClusteringOutputHandling(cOutput->getClusteringModelOutput().front(), iResults,cMain->getInput()->getDataType() ,cMain->getInput()->getCriterionName());
    mixmodClustering.slot("bestResult") = Rcpp::clone(iResults);
  }
  
  // add error
  mixmodClustering.slot("error") = !cOutput->atLeastOneEstimationNoError();
  
  //std::cout << "Output finished" << std::endl;

  // return final output
  return mixmodClustering;
  END_RCPP  

}
