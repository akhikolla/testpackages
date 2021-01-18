//
//  functions used in STBDwDM package
//

#ifndef __STBDwDM__
#define __STBDwDM__

//MCMC Sampler
Rcpp::List STBDwDM_Rcpp(Rcpp::List DatObj_List,  Rcpp::List HyPara_List,
                        Rcpp::List MetrObj_List, Rcpp::List Para_List,
                        Rcpp::List DatAug_List,  Rcpp::List McmcObj_List,
                        arma::mat RawSamples, bool Interactive);

//STRUCT DEFINITIONS
struct datobj {
  double Rho;
  double ScaleY;
  double ScaleDM;
  int N;
  int M;
  int Nu;
  int NTheta;
  int TempCorInd;
  int FamilyInd;
  int WeightsInd;
  arma::mat YStarWide;
  arma::vec DM;
  arma::mat W;
  arma::mat TimeDist;
  arma::vec OneM;
  arma::vec OneNu;
  arma::mat EyeM;
  arma::mat EyeNu;
  arma::mat EyeNTheta;
  arma::mat Eye3;
  arma::mat ZDelta;
  arma::vec Z;
  arma::umat AdjacentEdgesBoolean;
};
struct hypara {
  arma::vec OmegaDeltaInvMuDelta;
  arma::mat OmegaDeltaInv;
  double Xi;
  arma::mat Psi;
  double APhi;
  double BPhi;
};
struct metrobj {
  arma::vec MetropTheta2;
  arma::vec MetropTheta3;
  arma::vec AcceptanceTheta2;
  arma::vec AcceptanceTheta3;
  double MetropPhi;
  double AcceptancePhi;
};
struct para {
  arma::vec Mu;
  arma::vec Tau2;
  arma::vec Alpha;
  arma::cube WAlphas;
  arma::cube JointCovariances;
  arma::cube RootiLikelihoods;
  arma::vec VecTheta;
  arma::mat Theta;
  arma::vec Delta;
  arma::vec MeanTheta;
  arma::mat T;
  arma::mat TInv;
  double Phi;
  arma::mat SIGMAPhi;
  arma::mat SIGMAPhiInv;
  arma::mat CovThetaInv;
  arma::mat RootiTheta;
  arma::mat MMat;
};
struct dataug {
  int NBelow;
  int NAbove;
  arma::mat TobitIndeces;
  arma::mat ProbitIndeces;
};
struct mcmcobj {
  int NBurn;
  int NSims;
  int NThin;
  int NPilot;
  int NTotal;
  int NKeep;
  arma::vec WhichKeep;
  arma::vec WhichPilotAdapt;
  arma::vec WhichBurnInProgress;
  arma::vec WhichBurnInProgressInt;
  arma::vec WhichSamplerProgress;
  arma::vec BurnInProgress;
  int BarLength;
  int PilotAdaptDenominator;
};

//COVARIANCE FUNCTIONS
arma::mat GetRooti(arma::mat const& Cov, arma::mat const& Eye);
arma::cube JointCovarianceCube(arma::cube const& WAlphas, arma::vec const& Tau2, arma::mat const& EyeM, double Rho, int M, int Nu);
arma::mat JointCovarianceMatrix(arma::mat const& WAlpha, double tau2, arma::mat const& EyeM, double Rho, int M);
arma::cube JointPrecisionCube(arma::cube const& WAlphas, arma::vec const& Tau2, arma::mat const& EyeM, double Rho, int M, int Nu);
arma::mat SIGMA(double Phi, int TempCorInd, arma::mat const& TimeDist, int Nu);
arma::cube WAlphaCube(arma::vec const& Alpha, arma::colvec const& Z, arma::mat const& W, int M, int Nu, int WeightsInd);
arma::mat WAlphaMatrix(double alpha, arma::colvec const& Z, arma::umat const& AdjacentEdgesBoolean, arma::mat const& W, int M, int WeightsInd);

//DISTRIBUTION FUNCTIONS
double lndMvn(arma::vec const& Y, arma::vec const& Mu, arma::mat const& Rooti);
double randuRcpp();
arma::mat rmvnormRcpp(int n, arma::vec const& mean, arma::mat const& sigma);
arma::vec rnormRcpp(int n, double mean, double sd);
double rtnormRcpp(double mean, double sd, bool Above);
double rtnormRcppMSM(double mean, double sd, double lower, double upper);
arma::mat rwishRcpp(double n, arma::mat const& V);

//MCMC CONVERSION FUNCTIONS
datobj ConvertDatObj(Rcpp::List DatObj_List);
hypara ConvertHyPara(Rcpp::List HyPara_List);
metrobj ConvertMetrObj(Rcpp::List MetrObj_List);
para ConvertPara(Rcpp::List Para_List);
mcmcobj ConvertMcmcObj(Rcpp::List McmcObj_List);
dataug ConvertDatAug(Rcpp::List DatAug_List);

//MCMC SAMPLER FUNCTIONS
para SampleDelta(datobj DatObj, para Para, hypara HyPara);
std::pair<para, metrobj> SamplePhi(datobj DatObj, para Para, hypara HyPara, metrobj MetrObj);
para SampleT(datobj DatObj, para Para, hypara HyPara);
para SampleTheta1(datobj DatObj, para Para);
std::pair<para, metrobj> SampleTheta2(datobj DatObj, para Para, metrobj MetrObj);
std::pair<para, metrobj> SampleTheta3(datobj DatObj, para Para, metrobj MetrObj);
arma::mat SampleProbit(datobj DatObj, para Para, dataug DatAug);
arma::mat SampleTobit(datobj DatObj, para Para, dataug DatAug);
datobj SampleY(datobj DatObj, para Para, dataug DatAug);

//MCMC UTILITY FUNCTIONS
void BeginBurnInProgress(mcmcobj McmcObj, bool Interactive);
Rcpp::List OutputMetrObj(metrobj MetrObj);
metrobj PilotAdaptation(datobj DatObj, metrobj MetrObj, mcmcobj McmcObj);
void SamplerProgress(int s, mcmcobj McmcObj);
arma::colvec StoreSamples(datobj DatObj, para Para);
void UpdateBurnInBar(int s, mcmcobj McmcObj);
void UpdateBurnInBarInt(int s, mcmcobj McmcObj);

//UTILITY FUNCTIONS
arma::mat CholInv(arma::mat const& Cov);
arma::mat Inv2(arma::mat const& A);
arma::mat Inv3(arma::mat const& A);
bool rows_equal(arma::mat const& lhs, arma::mat const& rhs, double tol);

#endif // __STBDwDM__
