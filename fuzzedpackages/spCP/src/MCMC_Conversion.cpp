#include <RcppArmadillo.h>
#include "MCMC_spCP.h"

//Function to convert Rcpp::List DatObj to a custom C++ struct datobj--------------------------------------------------
datobj ConvertDatObj(Rcpp::List DatObj_List) {

  //Set objects from List
  double Rho = DatObj_List["Rho"];
  double ScaleY = DatObj_List["ScaleY"];
  double ScaleDM = DatObj_List["ScaleDM"];
  double tNu = DatObj_List["tNu"];
  double t1 = DatObj_List["t1"];
  int N = DatObj_List["N"];
  int M = DatObj_List["M"];
  int Nu = DatObj_List["Nu"];
  int FamilyInd = DatObj_List["FamilyInd"];
  int WeightsInd = DatObj_List["WeightsInd"];
  arma::vec YStar = DatObj_List["YStar"];
  arma::mat YStarWide = DatObj_List["YStarWide"];
  arma::vec DM = DatObj_List["DM"];
  arma::Mat<int> W = DatObj_List["W"];
  arma::vec Time = DatObj_List["Time"];
  arma::vec TimeVec = DatObj_List["TimeVec"];
  arma::vec OneM = DatObj_List["OneM"];
  arma::vec OneNu = DatObj_List["OneNu"];
  arma::vec OneN = DatObj_List["OneN"];
  arma::mat EyeM = DatObj_List["EyeM"];
  arma::mat EyeNu = DatObj_List["EyeNu"];
  arma::mat EyeN = DatObj_List["EyeN"];
  arma::mat Eye5 = DatObj_List["Eye5"];
  arma::mat Eye5M = DatObj_List["Eye5M"];
  arma::mat ZDelta = DatObj_List["ZDelta"];
  arma::vec DMLong = DatObj_List["DMLong"];
  arma::umat AdjacentEdgesBoolean = DatObj_List["AdjacentEdgesBoolean"];
  arma::umat PhiIndeces = DatObj_List["PhiIndeces"];
  arma::uvec XThetaInd = DatObj_List["XThetaInd"];

  //Convert to C++ struct
  datobj DatObj;
  DatObj.Rho = Rho;
  DatObj.ScaleY = ScaleY;
  DatObj.ScaleDM = ScaleDM;
  DatObj.tNu = tNu;
  DatObj.t1 = t1;
  DatObj.N = N;
  DatObj.M = M;
  DatObj.Nu = Nu;
  DatObj.FamilyInd = FamilyInd;
  DatObj.WeightsInd = WeightsInd;
  DatObj.YStar = YStar;
  DatObj.YStarWide = YStarWide;
  DatObj.DM = DM;
  DatObj.W = W;
  DatObj.Time = Time;
  DatObj.TimeVec = TimeVec;
  DatObj.OneM = OneM;
  DatObj.OneNu = OneNu;
  DatObj.OneN = OneN;
  DatObj.EyeM = EyeM;
  DatObj.EyeNu = EyeNu;
  DatObj.EyeN = EyeN;
  DatObj.Eye5 = Eye5;
  DatObj.Eye5M = Eye5M;
  DatObj.ZDelta = ZDelta;
  DatObj.DMLong = DMLong;
  DatObj.AdjacentEdgesBoolean = AdjacentEdgesBoolean;
  DatObj.PhiIndeces = PhiIndeces;
  DatObj.XThetaInd = XThetaInd;
  return DatObj;

}



//Function to convert Rcpp::List HyPara to a custom C++ struct hypara--------------------------------------------------
hypara ConvertHyPara(Rcpp::List HyPara_List) {

  //Set objects from List
  double Kappa2 = HyPara_List["Kappa2"];
  double Xi = HyPara_List["Xi"];
  arma::mat Psi = HyPara_List["Psi"];
  double AAlpha = HyPara_List["AAlpha"];
  double BAlpha = HyPara_List["BAlpha"];

  //Convert to C++ struct
  hypara HyPara;
  HyPara.Kappa2 = Kappa2;
  HyPara.Xi = Xi;
  HyPara.Psi = Psi;
  HyPara.AAlpha = AAlpha;
  HyPara.BAlpha = BAlpha;
  return HyPara;

}



//Function to convert Rcpp::List MetrObj to a custom C++ struct metrobj-----------------------------------------------
metrobj ConvertMetrObj(Rcpp::List MetrObj_List) {

  //Set objects from List
  arma::vec MetropLambda0Vec = MetrObj_List["MetropLambda0Vec"];
  arma::vec AcceptanceLambda0Vec = MetrObj_List["AcceptanceLambda0Vec"];
  arma::vec MetropLambda1Vec = MetrObj_List["MetropLambda1Vec"];
  arma::vec AcceptanceLambda1Vec = MetrObj_List["AcceptanceLambda1Vec"];
  arma::vec MetropEtaVec = MetrObj_List["MetropEtaVec"];
  arma::vec AcceptanceEtaVec = MetrObj_List["AcceptanceEtaVec"];
  double MetropAlpha = MetrObj_List["MetropAlpha"];
  double AcceptanceAlpha = MetrObj_List["AcceptanceAlpha"];

  //Convert to C++ struct
  metrobj MetrObj;
  MetrObj.MetropLambda0Vec = MetropLambda0Vec;
  MetrObj.AcceptanceLambda0Vec = AcceptanceLambda0Vec;
  MetrObj.MetropLambda1Vec = MetropLambda1Vec;
  MetrObj.AcceptanceLambda1Vec = AcceptanceLambda1Vec;
  MetrObj.MetropEtaVec = MetropEtaVec;
  MetrObj.AcceptanceEtaVec = AcceptanceEtaVec;
  MetrObj.MetropAlpha = MetropAlpha;
  MetrObj.AcceptanceAlpha = AcceptanceAlpha;
  return MetrObj;

}



//Function to convert Rcpp::List Para to a custom C++ struct para-----------------------------------------------------
para ConvertPara(Rcpp::List Para_List) {

  //Set objects from List
  arma::vec Beta = Para_List["Beta"];
  arma::vec Lambda = Para_List["Lambda"];
  arma::vec Eta = Para_List["Eta"];
  arma::vec Delta = Para_List["Delta"];
  double Alpha = Para_List["Alpha"];
  arma::mat Sigma = Para_List["Sigma"];
  arma::vec Sigma2 = Para_List["Sigma2"];
  arma::mat Omega = Para_List["Omega"];
  arma::mat OmegaInv = Para_List["OmegaInv"];
  arma::mat WAlpha = Para_List["WAlpha"];
  arma::mat QInv = Para_List["QInv"];
  arma::mat Q = Para_List["Q"];
  arma::mat SigmaInv = Para_List["SigmaInv"];
  arma::vec Theta = Para_List["Theta"];
  arma::mat XTheta = Para_List["XTheta"];
  arma::vec Mu = Para_List["Mu"];
  arma::vec Phi = Para_List["Phi"];
  arma::mat PhiPrec = Para_List["PhiPrec"];
  arma::mat PhiCov = Para_List["PhiCov"];
  arma::vec PhiMean = Para_List["PhiMean"];

  //Convert to C++ struct
  para Para;
  Para.Beta = Beta;
  Para.Lambda = Lambda;
  Para.Eta = Eta;
  Para.Delta = Delta;
  Para.Alpha = Alpha;
  Para.Sigma = Sigma;
  Para.Sigma2 = Sigma2;
  Para.Omega = Omega;
  Para.OmegaInv = OmegaInv;
  Para.WAlpha = WAlpha;
  Para.QInv = QInv;
  Para.Q = Q;
  Para.SigmaInv = SigmaInv;
  Para.Theta = Theta;
  Para.XTheta = XTheta;
  Para.Mu = Mu;
  Para.Phi = Phi;
  Para.PhiPrec = PhiPrec;
  Para.PhiCov = PhiCov;
  Para.PhiMean = PhiMean;
  return Para;
}



//Function to convert Rcpp::List DatAug to a custom C++ struct dataug-----------------------------------------------------
dataug ConvertDatAug(Rcpp::List DatAug_List) {

  //Set objects from List
  int NBelow = DatAug_List["NBelow"];
  int NAbove = DatAug_List["NAbove"];
  arma::uvec WhichBelow = DatAug_List["WhichBelow"];
  arma::uvec WhichAbove = DatAug_List["WhichAbove"];

  //Convert to C++ struct
  dataug DatAug;
  DatAug.NBelow = NBelow;
  DatAug.NAbove = NAbove;
  DatAug.WhichBelow = WhichBelow;
  DatAug.WhichAbove = WhichAbove;
  return DatAug;
}



//Function to convert Rcpp::List McmcObj to a custom C++ struct mcmcmobj-----------------------------------------------------
mcmcobj ConvertMcmcObj(Rcpp::List McmcObj_List) {

  //Set objects from List
  int NBurn = McmcObj_List["NBurn"];
  int NSims = McmcObj_List["NSims"];
  int NThin = McmcObj_List["NThin"];
  int NPilot = McmcObj_List["NPilot"];
  int NTotal = McmcObj_List["NTotal"];
  int NKeep = McmcObj_List["NKeep"];
  arma::vec WhichKeep = McmcObj_List["WhichKeep"];
  arma::vec WhichPilotAdapt = McmcObj_List["WhichPilotAdapt"];
  arma::vec WhichBurnInProgress = McmcObj_List["WhichBurnInProgress"];
  arma::vec WhichBurnInProgressInt = McmcObj_List["WhichBurnInProgressInt"];
  arma::vec WhichSamplerProgress = McmcObj_List["WhichSamplerProgress"];
  arma::vec BurnInProgress = McmcObj_List["BurnInProgress"];
  int BarLength = McmcObj_List["BarLength"];
  int PilotAdaptDenominator = McmcObj_List["PilotAdaptDenominator"];

  //Convert to C++ struct
  mcmcobj McmcObj;
  McmcObj.NBurn = NBurn;
  McmcObj.NSims = NSims;
  McmcObj.NThin = NThin;
  McmcObj.NPilot = NPilot;
  McmcObj.NTotal = NTotal;
  McmcObj.NKeep = NKeep;
  McmcObj.WhichKeep = WhichKeep;
  McmcObj.WhichPilotAdapt = WhichPilotAdapt;
  McmcObj.WhichBurnInProgress = WhichBurnInProgress;
  McmcObj.WhichBurnInProgressInt = WhichBurnInProgressInt;
  McmcObj.WhichSamplerProgress = WhichSamplerProgress;
  McmcObj.BurnInProgress = BurnInProgress;
  McmcObj.PilotAdaptDenominator = PilotAdaptDenominator;
  McmcObj.BarLength = BarLength;
  return McmcObj;
}
