#include <RcppArmadillo.h>
#include "MCMC_bfa_sp.h"

//This function is being exported to R for use in this package exclusively...
//not for use by users.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List bfa_sp_Rcpp(Rcpp::List DatObj_List,  Rcpp::List HyPara_List,
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
  arma::Col<int> FamilyInd = DatObj.FamilyInd;
  int SpCorInd = DatObj.SpCorInd;
  int NTotal = McmcObj.NTotal;
  int NBurn = McmcObj.NBurn;
  int LInf = DatObj.LInf;
  int BNP = DatObj.CL;
  arma::vec WhichPilotAdapt = McmcObj.WhichPilotAdapt;
  arma::vec WhichKeep = McmcObj.WhichKeep;
  arma::vec WhichBurnInProgress = McmcObj.WhichBurnInProgress;
  arma::vec WhichBurnInProgressInt = McmcObj.WhichBurnInProgressInt;
  arma::vec WhichSamplerProgress = McmcObj.WhichSamplerProgress;
  std::pair<datobj, para> DAUpdate;
  std::pair<para, metrobj> Update;
  
  //User output
  BeginBurnInProgress(McmcObj, Interactive);

  //Begin MCMC Sampler
  for (int s = 1; s < NTotal + 1; s++) {

    // Rcpp::Rcout << std::fixed << s << std::endl;
    
    //Check for user interrupt every 50 iterations
    if (s % 50 == 0) Rcpp::checkUserInterrupt();
    
    // Rcpp::Rcout << std::fixed << "Sampling U:" << std::endl;
    
    //Gibbs step for Latent U
    if (LInf == 1) Para = SampleU(DatObj, Para);
    
    // Rcpp::Rcout << std::fixed << "Sampling Y:" << std::endl;
    
    // Data Augmentation Step
    if (any(FamilyInd != 0)) {
      DAUpdate = SampleY(DatObj, Para, DatAug);
      DatObj = DAUpdate.first;
      Para = DAUpdate.second;
    }
    
    // Rcpp::Rcout << std::fixed << "Sampling Theta:" << std::endl;
    
    //Gibbs step for Theta
    if (BNP == 1) Para = SampleTheta(DatObj, Para);

    // Rcpp::Rcout << std::fixed << "Sampling Xi:" << std::endl;
    
    //Gibbs step for Xi
    if (BNP == 1) Para = SampleXi(DatObj, Para);
    
    // Rcpp::Rcout << std::fixed << "Sampling Z:" << std::endl;
    
    //Gibbs step for Z
    if (BNP == 1) Para = SampleZ(DatObj, Para);

    // Rcpp::Rcout << std::fixed << "Sampling Alpha:" << std::endl;
    
    //Gibbs step for Alpha
    Para = SampleAlpha(DatObj, Para);
    
    // Rcpp::Rcout << std::fixed << "Sampling Kappa:" << std::endl;
    
    //Gibbs step for Kappa
    Para = SampleKappa(DatObj, Para, HyPara);
    
    // Rcpp::Rcout << std::fixed << "Sampling Rho:" << std::endl;
    
    //Metropolis step for Rho
    if ((SpCorInd == 0) & (DatObj.IS == 1)) {
      Update = SampleRho(DatObj, Para, HyPara, MetrObj);
      Para = Update.first;
      MetrObj = Update.second;
    }
  
    // Rcpp::Rcout << std::fixed << "Sampling Delta:" << std::endl;
    
    //Gibbs step for Delta
    if (BNP == 1) Para = SampleDelta(DatObj, Para, HyPara);
    
    // Rcpp::Rcout << std::fixed << "Sampling Eta:" << std::endl;
    
    //Gibbs step for Eta
    Para = SampleEta(DatObj, Para, HyPara);
    
    // Rcpp::Rcout << std::fixed << "Sampling Upsilon:" << std::endl;
    
    //Gibbs step for Upsilon
    Para = SampleUpsilon(DatObj, Para, HyPara);
    
    // Rcpp::Rcout << std::fixed << "Sampling Psi:" << std::endl;
    
    // Metropolis step for Psi
    Update = SamplePsi(DatObj, Para, HyPara, MetrObj);
    Para = Update.first;
    MetrObj = Update.second;
    
    // Rcpp::Rcout << std::fixed << "Sampling Sigma2:" << std::endl;
    
    //Gibbs step for Sigma2
    if (any(FamilyInd != 3)) {
      Para = SampleSigma2(DatObj, Para, HyPara);
    }
    
    // Rcpp::Rcout << std::fixed << "Sampling Beta:" << std::endl;
    
    //Gibbs step for Beta
    if (DatObj.P > 0) {
      Para = SampleBeta(DatObj, Para, HyPara);
    }
    
    //Pilot adaptation
    if (std::find(WhichPilotAdapt.begin(), WhichPilotAdapt.end(), s) != WhichPilotAdapt.end())
      MetrObj = PilotAdaptation(MetrObj, McmcObj, DatObj);

    //Store raw samples
    if (std::find(WhichKeep.begin(), WhichKeep.end(), s) != WhichKeep.end())
      RawSamples.cols(find(s == WhichKeep)) = StoreSamples(DatObj, Para);
    
    //Update burn-in progress bar
    if (Interactive) if (std::find(WhichBurnInProgress.begin(), WhichBurnInProgress.end(), s) != WhichBurnInProgress.end())
      UpdateBurnInBar(s, McmcObj);
    if (!Interactive) if (std::find(WhichBurnInProgressInt.begin(), WhichBurnInProgressInt.end(), s) != WhichBurnInProgressInt.end())
      UpdateBurnInBarInt(s, McmcObj);

    //Post burn-in progress
    if (s == NBurn) Rcpp::Rcout << std::fixed << "\nSampler progress:  0%.. ";
    if (std::find(WhichSamplerProgress.begin(), WhichSamplerProgress.end(), s) != WhichSamplerProgress.end())
       SamplerProgress(s, McmcObj);

  //End MCMC Sampler
  }

  //Output Metropolis object for summary
  Rcpp::List Metropolis = OutputMetrObj(MetrObj, DatObj);

  //Return raw samples
  return Rcpp::List::create(Rcpp::Named("rawsamples") = RawSamples,
                            Rcpp::Named("metropolis") = Metropolis);

//End MCMC sampler function
}
