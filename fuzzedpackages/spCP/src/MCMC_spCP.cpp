#include <RcppArmadillo.h>
#include "MCMC_spCP.h"

//This function is being exported to R for use in this package exclusively...
//not for use by users.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List spCP_Rcpp(Rcpp::List DatObj_List,  Rcpp::List HyPara_List,
                     Rcpp::List MetrObj_List, Rcpp::List Para_List,
                     Rcpp::List DatAug_List,  Rcpp::List McmcObj_List,
                     arma::mat RawSamples, bool Interactive) {

  //Convet Rcpp::Lists to C++ structs
  datobj DatObj = ConvertDatObj(DatObj_List);
  hypara HyPara = ConvertHyPara(HyPara_List);
  metrobj MetrObj = ConvertMetrObj(MetrObj_List);
  para Para = ConvertPara(Para_List);
  dataug DatAug = ConvertDatAug(DatAug_List);
  mcmcobj McmcObj = ConvertMcmcObj(McmcObj_List);

  //Set objects to be used in MCMC sampler
  int FamilyInd = DatObj.FamilyInd;
  int NTotal = McmcObj.NTotal;
  int NBurn = McmcObj.NBurn;
  int NTrunc = DatAug.NBelow + DatAug.NAbove;
  arma::vec WhichPilotAdapt = McmcObj.WhichPilotAdapt;
  arma::vec WhichKeep = McmcObj.WhichKeep;
  arma::vec WhichBurnInProgress = McmcObj.WhichBurnInProgress;
  arma::vec WhichBurnInProgressInt = McmcObj.WhichBurnInProgressInt;
  arma::vec WhichSamplerProgress = McmcObj.WhichSamplerProgress;
  std::pair<para, metrobj> Update;

  //User output
  BeginBurnInProgress(McmcObj, Interactive);

  //Begin MCMC Sampler
  for (int s = 1; s < NTotal + 1; s++) {
  // for (int s = 1; s < 3; s++) {

    //Check for user interrupt every 500 iterations
    if (s % 500 == 0) Rcpp::checkUserInterrupt();

    // Data Augmentation Step
    if ((FamilyInd != 0) & (NTrunc > 0)) DatObj = SampleY(DatObj, Para, DatAug);

    //Gibbs step for Delta
    Para = SampleDelta(DatObj, Para, HyPara);

    //Gibbs step for Beta (i.e. Beta0 and Beta1)
    Para = SampleBeta(DatObj, Para);

    //Metropolis step for Lambda0
    Update = SampleLambda0(DatObj, Para, MetrObj);
    Para = Update.first;
    MetrObj = Update.second;

    //Metropolis step for Lambda1
    Update = SampleLambda1(DatObj, Para, MetrObj);
    Para = Update.first;
    MetrObj = Update.second;

    //Metropolis step for Eta
    Update = SampleEta(DatObj, Para, MetrObj);
    Para = Update.first;
    MetrObj = Update.second;

    //Gibbs sampler step for Sigma
    Para = SampleSigma(DatObj, Para, HyPara);

    //Metropolis step for Alpha
    Update = SampleAlpha(DatObj, Para, HyPara, MetrObj);
    Para = Update.first;
    MetrObj = Update.second;

    //Pilot adaptation
    if (std::find(WhichPilotAdapt.begin(), WhichPilotAdapt.end(), s) != WhichPilotAdapt.end())
      MetrObj = PilotAdaptation(DatObj, MetrObj, McmcObj);

    //Store raw samples
    if (std::find(WhichKeep.begin(), WhichKeep.end(), s) != WhichKeep.end())
      RawSamples.cols(find(s == WhichKeep)) = StoreSamples(DatObj, Para);

    //Update burn-in progress bar
    if (Interactive) if (std::find(WhichBurnInProgress.begin(), WhichBurnInProgress.end(), s) != WhichBurnInProgress.end())
      UpdateBurnInBar(s, McmcObj);
    if (!Interactive) if (std::find(WhichBurnInProgressInt.begin(), WhichBurnInProgressInt.end(), s) != WhichBurnInProgressInt.end())
      UpdateBurnInBarInt(s, McmcObj);

    //Post burn-in progress
    if (s == NBurn) Rcpp::Rcout << std::fixed << "\nSampler progress:  0%..  ";
    if (std::find(WhichSamplerProgress.begin(), WhichSamplerProgress.end(), s) != WhichSamplerProgress.end())
       SamplerProgress(s, McmcObj);

  //End MCMC Sampler
  }

  //Output Metropolis object for summary
  Rcpp::List Metropolis = OutputMetrObj(MetrObj);

  //Return raw samples
  return Rcpp::List::create(Rcpp::Named("rawsamples") = RawSamples,
                            Rcpp::Named("metropolis") = Metropolis);

//End MCMC sampler function
}
