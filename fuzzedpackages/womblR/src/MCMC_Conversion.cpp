#include <RcppArmadillo.h>
#include "MCMC_STBDwDM.h"

//Function to convert Rcpp::List DatObj to a custom C++ struct datobj--------------------------------------------------
datobj ConvertDatObj(Rcpp::List DatObj_List) {

  //Set objects from List
  double ScaleY = DatObj_List["ScaleY"];
  double ScaleDM = DatObj_List["ScaleDM"];
  arma::mat YStarWide = DatObj_List["YStarWide"];
  arma::vec DM = DatObj_List["DM"];
  arma::mat W = DatObj_List["W"];
  arma::mat UpperW = arma::trimatu(W);
  arma::umat AdjacentEdgesBoolean = find(UpperW == 1);
  arma::mat TimeDist = DatObj_List["TimeDist"];
  double Rho = DatObj_List["Rho"];
  int N = DatObj_List["N"];
  int M = DatObj_List["M"];
  int Nu = DatObj_List["Nu"];
  int NTheta = DatObj_List["NTheta"];
  arma::vec OneM = DatObj_List["OneM"];
  arma::vec OneNu = DatObj_List["OneNu"];
  arma::mat EyeM = DatObj_List["EyeM"];
  arma::mat EyeNu = DatObj_List["EyeNu"];
  arma::mat EyeNTheta = DatObj_List["EyeNTheta"];
  arma::mat Eye3 = DatObj_List["Eye3"];
  arma::mat ZDelta = DatObj_List["ZDelta"];
  arma::vec Z = DatObj_List["Z"];
  int TempCorInd = DatObj_List["TempCorInd"];
  int FamilyInd = DatObj_List["FamilyInd"];
  int WeightsInd = DatObj_List["WeightsInd"];

  //Convert to C++ struct
  datobj DatObj;
  DatObj.ScaleY = ScaleY;
  DatObj.ScaleDM = ScaleDM;
  DatObj.YStarWide = YStarWide;
  DatObj.DM = DM;
  DatObj.W = W;
  DatObj.AdjacentEdgesBoolean = AdjacentEdgesBoolean;
  DatObj.TimeDist = TimeDist;
  DatObj.Rho = Rho;
  DatObj.N = N;
  DatObj.M = M;
  DatObj.Nu = Nu;
  DatObj.NTheta = NTheta;
  DatObj.OneM = OneM;
  DatObj.OneNu = OneNu;
  DatObj.EyeM = EyeM;
  DatObj.EyeNu = EyeNu;
  DatObj.EyeNTheta = EyeNTheta;
  DatObj.Eye3 = Eye3;
  DatObj.ZDelta = ZDelta;
  DatObj.Z = Z;
  DatObj.TempCorInd = TempCorInd;
  DatObj.FamilyInd = FamilyInd;
  DatObj.WeightsInd = WeightsInd;
  return DatObj;

}



//Function to convert Rcpp::List HyPara to a custom C++ struct hypara--------------------------------------------------
hypara ConvertHyPara(Rcpp::List HyPara_List) {

  //Set objects from List
  arma::vec OmegaDeltaInvMuDelta = HyPara_List["OmegaDeltaInvMuDelta"];
  arma::mat OmegaDeltaInv = HyPara_List["OmegaDeltaInv"];
  double Xi = HyPara_List["Xi"];
  arma::mat Psi = HyPara_List["Psi"];
  double APhi = HyPara_List["APhi"];
  double BPhi = HyPara_List["BPhi"];

  //Convert to C++ struct
  hypara HyPara;
  HyPara.OmegaDeltaInvMuDelta = OmegaDeltaInvMuDelta;
  HyPara.OmegaDeltaInv = OmegaDeltaInv;
  HyPara.Xi = Xi;
  HyPara.Psi = Psi;
  HyPara.APhi = APhi;
  HyPara.BPhi = BPhi;
  return HyPara;

}



//Function to convert Rcpp::List MetrObj to a custom C++ struct metrobj-----------------------------------------------
metrobj ConvertMetrObj(Rcpp::List MetrObj_List) {

  //Set objects from List
  arma::vec MetropTheta2 = MetrObj_List["MetropTheta2"];
  arma::vec AcceptanceTheta2 = MetrObj_List["AcceptanceTheta2"];
  arma::vec MetropTheta3 = MetrObj_List["MetropTheta3"];
  arma::vec AcceptanceTheta3 = MetrObj_List["AcceptanceTheta3"];
  double MetropPhi = MetrObj_List["MetropPhi"];
  double AcceptancePhi = MetrObj_List["AcceptancePhi"];

  //Convert to C++ struct
  metrobj MetrObj;
  MetrObj.MetropTheta2 = MetropTheta2;
  MetrObj.AcceptanceTheta2 = AcceptanceTheta2;
  MetrObj.MetropTheta3 = MetropTheta3;
  MetrObj.AcceptanceTheta3 = AcceptanceTheta3;
  MetrObj.MetropPhi = MetropPhi;
  MetrObj.AcceptancePhi = AcceptancePhi;
  return MetrObj;

}



//Function to convert Rcpp::List Para to a custom C++ struct para-----------------------------------------------------
para ConvertPara(Rcpp::List Para_List) {

  //Set objects from List
  arma::vec Mu = Para_List["Mu"];
  arma::vec Tau2 = Para_List["Tau2"];
  arma::vec Alpha = Para_List["Alpha"];
  arma::cube WAlphas = Para_List["WAlphas"];
  arma::cube JointCovariances = Para_List["JointCovariances"];
  arma::cube RootiLikelihoods = Para_List["RootiLikelihoods"];
  arma::vec VecTheta = Para_List["VecTheta"];
  arma::mat Theta = Para_List["Theta"];
  arma::vec Delta = Para_List["Delta"];
  arma::vec MeanTheta = Para_List["MeanTheta"];
  arma::mat T = Para_List["T"];
  arma::mat TInv = Para_List["TInv"];
  double Phi = Para_List["Phi"];
  arma::mat SIGMAPhi = Para_List["SIGMAPhi"];
  arma::mat SIGMAPhiInv = Para_List["SIGMAPhiInv"];
  arma::mat CovThetaInv = Para_List["CovThetaInv"];
  arma::mat RootiTheta = Para_List["RootiTheta"];
  arma::mat MMat = Para_List["MMat"];

  //Convert to C++ struct
  para Para;
  Para.Mu = Mu;
  Para.Tau2 = Tau2;
  Para.Alpha = Alpha;
  Para.WAlphas = WAlphas;
  Para.JointCovariances = JointCovariances;
  Para.RootiLikelihoods = RootiLikelihoods;
  Para.VecTheta = VecTheta;
  Para.Theta = Theta;
  Para.Delta = Delta;
  Para.MeanTheta = MeanTheta;
  Para.T = T;
  Para.TInv = TInv;
  Para.Phi = Phi;
  Para.SIGMAPhi = SIGMAPhi;
  Para.SIGMAPhiInv = SIGMAPhiInv;
  Para.CovThetaInv = CovThetaInv;
  Para.RootiTheta = RootiTheta;
  Para.MMat = MMat;
  return Para;
}



//Function to convert Rcpp::List DatAug to a custom C++ struct dataug-----------------------------------------------------
dataug ConvertDatAug(Rcpp::List DatAug_List) {

  //Set objects from List
  int NBelow = DatAug_List["NBelow"];
  int NAbove = DatAug_List["NAbove"];
  arma::mat TobitIndeces = DatAug_List["TobitIndeces"];
  arma::mat ProbitIndeces = DatAug_List["ProbitIndeces"];

  //Convert to C++ struct
  dataug DatAug;
  DatAug.NBelow = NBelow;
  DatAug.NAbove = NAbove;
  DatAug.TobitIndeces = TobitIndeces;
  DatAug.ProbitIndeces = ProbitIndeces;
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



