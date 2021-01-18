#include <RcppArmadillo.h>
#include "MCMC_bfa_sp.h"

//Function to convert Rcpp::List DatObj to a custom C++ struct datobj--------------------------------------------------
datobj ConvertDatObj(Rcpp::List DatObj_List) {

  //Set objects from List
  arma::mat YStarWide = DatObj_List["YStarWide"];
  arma::colvec YStar = DatObj_List["YStar"];
  arma::cube YObserved = DatObj_List["YObserved"];
  arma::mat SpDist = DatObj_List["SpDist"];
  arma::mat TimeDist = DatObj_List["TimeDist"];
  arma::colvec Time = DatObj_List["Time"];
  int N = DatObj_List["N"];
  int K = DatObj_List["K"];
  int L = DatObj_List["L"];
  int M = DatObj_List["M"];
  int Nu = DatObj_List["Nu"];
  int O = DatObj_List["O"];
  int C = DatObj_List["C"];
  int TempCorInd = DatObj_List["TempCorInd"];
  int SpCorInd = DatObj_List["SpCorInd"];
  int GS = DatObj_List["GS"];
  int IS = DatObj_List["IS"];
  int CL = DatObj_List["CL"];
  arma::Col<int> FamilyInd = DatObj_List["FamilyInd"];
  int LInf = DatObj_List["LInf"];
  arma::mat EyeNu = DatObj_List["EyeNu"];
  arma::Col<int> SeqL = DatObj_List["SeqL"];
  arma::mat EyeM = DatObj_List["EyeM"];
  arma::mat EyeO = DatObj_List["EyeO"];
  arma::mat EyeK = DatObj_List["EyeK"];
  arma::mat EyeOM = DatObj_List["EyeOM"];
  arma::mat EyeKbyNu = DatObj_List["EyeKbyNu"];
  arma::colvec ZeroKbyNu = DatObj_List["ZeroKbyNu"];
  arma::colvec ZeroM = DatObj_List["ZeroM"];
  arma::colvec ZeroOM = DatObj_List["ZeroOM"];
  arma::colvec OneNu = DatObj_List["OneNu"];
  arma::colvec OneO = DatObj_List["OneO"];
  arma::cube Trials = DatObj_List["Trials"];
  arma::cube Chi = DatObj_List["Chi"];
  arma::mat X = DatObj_List["X"];
  int P = DatObj_List["P"];
  arma::Col<int> Indeces = DatObj_List["Indeces"];

  //Convert to C++ struct
  datobj DatObj;
  DatObj.YStarWide = YStarWide;
  DatObj.YStar = YStar;
  DatObj.YObserved = YObserved;
  DatObj.SpDist = SpDist;
  DatObj.TimeDist = TimeDist;
  DatObj.Time = Time;
  DatObj.N = N;
  DatObj.K = K;
  DatObj.L = L;
  DatObj.M = M;
  DatObj.Nu = Nu;
  DatObj.O = O;
  DatObj.C = C;
  DatObj.GS = GS;
  DatObj.IS = IS;
  DatObj.CL = CL;
  DatObj.TempCorInd = TempCorInd;
  DatObj.SpCorInd  = SpCorInd;
  DatObj.FamilyInd = FamilyInd;
  DatObj.LInf = LInf;
  DatObj.EyeNu = EyeNu;
  DatObj.EyeO = EyeO;
  DatObj.SeqL = SeqL;
  DatObj.EyeM = EyeM;
  DatObj.EyeK = EyeK;
  DatObj.EyeOM = EyeOM;
  DatObj.EyeKbyNu = EyeKbyNu;
  DatObj.ZeroOM = ZeroOM;
  DatObj.ZeroKbyNu = ZeroKbyNu;
  DatObj.ZeroM = ZeroM;
  DatObj.OneNu = OneNu;
  DatObj.OneO = OneO;
  DatObj.Trials = Trials;
  DatObj.Chi = Chi;
  DatObj.P = P;
  DatObj.X = X;
  DatObj.Indeces = Indeces;
  return DatObj;

}



//Function to convert Rcpp::List HyPara to a custom C++ struct hypara--------------------------------------------------
hypara ConvertHyPara(Rcpp::List HyPara_List) {

  //Set objects from List
  double A = HyPara_List["A"];
  double B = HyPara_List["B"];
  double SmallUpsilon = HyPara_List["SmallUpsilon"];
  double A1 = HyPara_List["A1"];
  double A2 = HyPara_List["A2"];
  double APsi = HyPara_List["APsi"];
  double BPsi = HyPara_List["BPsi"];
  double ARho = HyPara_List["ARho"];
  double BRho = HyPara_List["BRho"];
  double Gamma = HyPara_List["Gamma"];
  double Beta = HyPara_List["Beta"];
  double Zeta = HyPara_List["Zeta"];
  arma::mat Omega = HyPara_List["Omega"];
  arma::mat BigTheta = HyPara_List["BigTheta"];
  arma::colvec SigmaBetaInvMuBeta = HyPara_List["SigmaBetaInvMuBeta"];
  arma::mat SigmaBetaInv = HyPara_List["SigmaBetaInv"];
  
  //Convert to C++ struct
  hypara HyPara;
  HyPara.A = A;
  HyPara.B = B;
  HyPara.SmallUpsilon = SmallUpsilon;
  HyPara.A1 = A1;
  HyPara.A2 = A2;
  HyPara.APsi = APsi;
  HyPara.BPsi = BPsi;
  HyPara.ARho = ARho;
  HyPara.BRho = BRho;
  HyPara.Gamma = Gamma;
  HyPara.Beta = Beta;
  HyPara.Zeta = Zeta;
  HyPara.Omega = Omega;
  HyPara.BigTheta = BigTheta;
  HyPara.SigmaBetaInvMuBeta = SigmaBetaInvMuBeta;
  HyPara.SigmaBetaInv = SigmaBetaInv;
  return HyPara;

}



//Function to convert Rcpp::List MetrObj to a custom C++ struct metrobj-----------------------------------------------
metrobj ConvertMetrObj(Rcpp::List MetrObj_List) {

  //Set objects from List
  double MetropPsi = MetrObj_List["MetropPsi"];
  int AcceptancePsi = MetrObj_List["AcceptancePsi"];
  double MetropRho = MetrObj_List["MetropRho"];
  int AcceptanceRho = MetrObj_List["AcceptanceRho"];
  arma::vec OriginalTuners = MetrObj_List["OriginalTuners"];
  
  //Convert to C++ struct
  metrobj MetrObj;
  MetrObj.MetropPsi = MetropPsi;
  MetrObj.AcceptancePsi = AcceptancePsi;
  MetrObj.MetropRho = MetropRho;
  MetrObj.AcceptanceRho = AcceptanceRho;
  MetrObj.OriginalTuners = OriginalTuners;
  return MetrObj;

}



//Function to convert Rcpp::List Para to a custom C++ struct para-----------------------------------------------------
para ConvertPara(Rcpp::List Para_List) {

  //Set objects from List
  arma::mat Sigma2 = Para_List["Sigma2"];
  arma::mat Kappa = Para_List["Kappa"];
  double Rho = Para_List["Rho"];
  arma::colvec Delta = Para_List["Delta"];
  arma::colvec Beta = Para_List["Beta"];
  double Psi = Para_List["Psi"];
  arma::mat Upsilon = Para_List["Upsilon"];
  arma::mat UpsilonInv = Para_List["UpsilonInv"];
  arma::umat Xi = Para_List["Xi"];
  arma::mat Theta = Para_List["Theta"];
  arma::mat Lambda = Para_List["Lambda"];
  arma::colvec Tau = Para_List["Tau"];
  arma::mat BigPhi = Para_List["BigPhi"];
  arma::colvec Eta = Para_List["Eta"];
  arma::cube Alpha = Para_List["Alpha"];
  arma::cube Z = Para_List["Z"];
  arma::mat HPsi = Para_List["HPsi"];
  arma::mat CholHPsi = Para_List["CholHPsi"];
  arma::mat HPsiInv = Para_List["HPsiInv"];
  arma::colvec Mean = Para_List["Mean"];
  arma::cube Weights = Para_List["Weights"];
  arma::cube logWeights = Para_List["logWeights"];
  arma::mat U = Para_List["U"];
  arma::colvec LStarJ = Para_List["LStarJ"];
  arma::mat SpCov = Para_List["SpCov"];
  arma::mat SpCovInv = Para_List["SpCovInv"];
  arma::mat CholSpCov = Para_List["CholSpCov"];
  arma::mat CholKappa = Para_List["CholKappa"];
  arma::mat KappaInv = Para_List["KappaInv"];
  arma::cube Cov = Para_List["Cov"];
  arma::colvec XBeta = Para_List["XBeta"];

  //Convert to C++ struct
  para Para;
  Para.Sigma2 = Sigma2;
  Para.Kappa = Kappa;
  Para.Rho = Rho;
  Para.Delta = Delta;
  Para.Beta = Beta;
  Para.Psi = Psi;
  Para.Upsilon = Upsilon;
  Para.UpsilonInv = UpsilonInv;
  Para.Xi = Xi;
  Para.Theta = Theta;
  Para.Lambda = Lambda;
  Para.Tau = Tau;
  Para.BigPhi = BigPhi;
  Para.Eta = Eta;
  Para.Alpha = Alpha;
  Para.Z = Z;
  Para.HPsi = HPsi;
  Para.CholHPsi = CholHPsi;
  Para.HPsiInv = HPsiInv;
  Para.Mean = Mean;
  Para.Weights = Weights;
  Para.logWeights = logWeights;
  Para.U = U;
  Para.LStarJ = LStarJ;
  Para.SpCov = SpCov;
  Para.SpCovInv = SpCovInv;
  Para.CholSpCov = CholSpCov;
  Para.CholKappa = CholKappa;
  Para.KappaInv = KappaInv;
  Para.Cov = Cov;
  Para.XBeta = XBeta;
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



