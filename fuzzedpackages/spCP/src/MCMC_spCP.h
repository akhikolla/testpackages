//
//  functions used in spCP package
//

#ifndef __spCP__
#define __spCP__

//MCMC Sampler
Rcpp::List spCP_Rcpp(Rcpp::List DatObj_List,  Rcpp::List HyPara_List,
                        Rcpp::List MetrObj_List, Rcpp::List Para_List,
                        Rcpp::List DatAug_List,  Rcpp::List McmcObj_List,
                        arma::mat RawSamples, bool Interactive);

//STRUCT DEFINITIONS
struct datobj {
  double Rho;
  double ScaleY;
  double ScaleDM;
  double tNu;
  double t1;
  int N;
  int M;
  int Nu;
  int FamilyInd;
  int WeightsInd;
  arma::vec YStar;
  arma::mat YStarWide;
  arma::vec DM;
  arma::Mat<int> W;
  arma::vec Time;
  arma::vec TimeVec;
  arma::vec OneM;
  arma::vec OneNu;
  arma::vec OneN;
  arma::mat EyeM;
  arma::mat EyeNu;
  arma::mat EyeN;
  arma::mat Eye5;
  arma::mat Eye5M;
  arma::mat ZDelta;
  arma::vec DMLong;
  arma::umat AdjacentEdgesBoolean;
  arma::umat PhiIndeces;
  arma::uvec XThetaInd;
};
struct hypara {
  double Kappa2;
  double Xi;
  arma::mat Psi;
  double AAlpha;
  double BAlpha;
};
struct metrobj {
  arma::vec MetropLambda0Vec;
  arma::vec AcceptanceLambda0Vec;
  arma::vec MetropLambda1Vec;
  arma::vec AcceptanceLambda1Vec;
  arma::vec MetropEtaVec;
  arma::vec AcceptanceEtaVec;
  double MetropAlpha;
  double AcceptanceAlpha;
};
struct para {
  arma::vec Beta;
  arma::vec Lambda;
  arma::vec Eta;
  arma::vec Delta;
  double Alpha;
  arma::mat Sigma;
  arma::vec Sigma2;
  arma::mat Omega;
  arma::mat OmegaInv;
  arma::mat WAlpha;
  arma::mat QInv;
  arma::mat Q;
  arma::mat SigmaInv;
  arma::vec Theta;
  arma::mat XTheta;
  arma::vec Mu;
  arma::vec Phi;
  arma::mat PhiPrec;
  arma::mat PhiCov;
  arma::vec PhiMean;
};
struct dataug {
  int NBelow;
  int NAbove;
  arma::uvec WhichAbove;
  arma::uvec WhichBelow;
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
arma::mat WAlphaFnc(double Alpha, arma::colvec const& DMLong, arma::umat const& AdjacentEdgesBoolean, arma::Mat<int> const& W, int M, int WeightsInd);
arma::mat QFnc(arma::mat const& WAlpha, arma::mat const& EyeM, double Rho, int M);
arma::mat QInvFnc(arma::mat const& WAlpha, arma::mat const& EyeM, double Rho, int M);

//DISTRIBUTION FUNCTIONS
double lndMvn(arma::vec const& Y, arma::vec const& Mu, arma::mat const& Rooti);
double randuRcpp();
arma::mat rmvnormRcpp(int n, arma::vec const& mean, arma::mat const& sigma);
arma::vec rnormRcpp(int n, double mean, double sd);
double rtnormRcpp(double mean, double sd, bool Above);
double rtnormRcppMSM(double mean, double sd, double lower, double upper);
// arma::vec rtnormRcppMSM(int N, arma::vec const& mean, arma::vec const& sd, double lower, double upper);
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
std::pair<para, metrobj> SampleAlpha(datobj DatObj, para Para, hypara HyPara, metrobj MetrObj);
para SampleSigma(datobj DatObj, para Para, hypara HyPara);
para SampleBeta(datobj DatObj, para Para);
std::pair<para, metrobj> SampleLambda0(datobj DatObj, para Para, metrobj MetrObj);
std::pair<para, metrobj> SampleLambda1(datobj DatObj, para Para, metrobj MetrObj);
std::pair<para, metrobj> SampleEta(datobj DatObj, para Para, metrobj MetrObj);
arma::vec SampleProbit(datobj DatObj, para Para, dataug DatAug);
arma::vec SampleTobit(datobj DatObj, para Para, dataug DatAug);
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
arma::mat GetXThetaLoc(double ThetaLoc, arma::vec const& Time, arma::vec const& OneNu, int Nu);
arma::vec CreatePhi(arma::vec const& Beta, arma::vec const& Lambda, arma::vec const& Eta, int M);
arma::mat CholInv(arma::mat const& Cov);
arma::mat Inv2(arma::mat const& A);
arma::mat Inv3(arma::mat const& A);
bool rows_equal(arma::mat const& lhs, arma::mat const& rhs, double tol);

#endif // __spCP__
