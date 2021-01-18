//
//  functions used in spBFA package
//

#ifndef __spBFA__
#define __spBFA__

//MCMC Sampler
Rcpp::List bfa_sp_Rcpp(Rcpp::List DatObj_List,  Rcpp::List HyPara_List,
                       Rcpp::List MetrObj_List, Rcpp::List Para_List,
                       Rcpp::List DatAug_List,  Rcpp::List McmcObj_List,
                       arma::mat RawSamples, bool Interactive);

//STRUCT DEFINITIONS
struct datobj {
  int N;
  int M;
  int Nu;
  int K;
  int L;
  int O;
  int C;
  int P;
  int GS;
  int IS;
  int CL;
  arma::Col<int> FamilyInd;
  arma::Col<int> Indeces;
  int TempCorInd;
  int SpCorInd;
  int LInf;
  arma::colvec YStar;
  arma::cube YObserved;
  arma::mat YStarWide;
  arma::mat SpDist;
  arma::colvec Time;
  arma::mat TimeDist;
  arma::mat EyeNu;
  arma::Col<int> SeqL;
  arma::mat EyeM;
  arma::mat EyeK;
  arma::mat EyeO;
  arma::mat EyeOM;
  arma::mat EyeKbyNu;
  arma::mat X;
  arma::colvec ZeroKbyNu;
  arma::colvec ZeroM;
  arma::colvec ZeroOM;
  arma::colvec OneNu;
  arma::colvec OneO;
  arma::cube Trials;
  arma::cube Chi;
};
struct hypara {
  double A;
  double B;
  double SmallUpsilon;
  double A1;
  double A2;
  double APsi;
  double BPsi;
  double ARho;
  double BRho;
  double Gamma;
  double Beta;
  double Zeta;
  arma::mat Omega;
  arma::mat BigTheta;
  arma::colvec SigmaBetaInvMuBeta;
  arma::mat SigmaBetaInv;
};
struct metrobj {
  double MetropPsi;
  double MetropRho;
  int AcceptanceRho;
  int AcceptancePsi;
  arma::vec OriginalTuners;
};
struct para {
  arma::mat Sigma2;
  arma::colvec Delta;
  arma::mat Kappa;
  double Psi;
  double Rho;
  arma::colvec Beta;
  arma::mat Upsilon;
  arma::mat UpsilonInv;
  arma::umat Xi;
  arma::mat Theta;
  arma::mat Lambda;
  arma::colvec Tau;
  arma::mat BigPhi;
  arma::colvec Eta;
  arma::cube Alpha;
  arma::cube Z;
  arma::mat HPsi;
  arma::mat CholHPsi;
  arma::mat HPsiInv;
  arma::colvec Mean;
  arma::cube Weights;
  arma::cube logWeights;
  arma::mat U;
  arma::colvec LStarJ;
  arma::mat SpCov;
  arma::mat SpCovInv;
  arma::mat CholSpCov;
  arma::mat CholKappa;
  arma::mat KappaInv;
  arma::cube Cov;
  arma::colvec XBeta;
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
arma::mat H(double Psi, int TempCorInd, arma::mat const& TimeDist, int Nu);
arma::mat SpEXP(double Rho, arma::mat const& SpDist, int M);

//DISTRIBUTION FUNCTIONS
arma::vec rnormRcpp(int n, double mean, double sd);
arma::vec sampleRcpp(arma::Col<int> const& x, int size, bool replace, arma::vec const& prob);
double rtnormRcppMSM(double mean, double sd, double lower, double upper);
arma::mat rmvnormRcpp(int n, arma::vec const& mean, arma::mat const& sigma);
double pnormRcpp(double q);
double lpnormRcpp(double q);
double UpperpnormRcpp(double q);
double lUpperpnormRcpp(double q);
double rigammaRcpp(double Alpha, double Theta);
double rgammaRcpp(double Alpha, double Theta);
arma::mat rwishRcpp(double n, arma::mat const& V);
double lndMvn(arma::vec const& Y, arma::vec const& Mu, arma::mat const& Rooti);
double randuRcpp();
double rtnormRcpp(double mean, double sd, bool Above);
arma::vec rtnormRcppMSM(int N, arma::vec const& mean, arma::vec const& sd, double lower, double upper);
arma::vec pgRcpp(arma::vec const& b, arma::vec const& c);

//MCMC CONVERSION FUNCTIONS
datobj ConvertDatObj(Rcpp::List DatObj_List);
hypara ConvertHyPara(Rcpp::List HyPara_List);
metrobj ConvertMetrObj(Rcpp::List MetrObj_List);
para ConvertPara(Rcpp::List Para_List);
mcmcobj ConvertMcmcObj(Rcpp::List McmcObj_List);
dataug ConvertDatAug(Rcpp::List DatAug_List);

//MCMC SAMPLER FUNCTIONS
para SampleTheta(datobj DatObj, para Para);
para SampleXi(datobj DatObj, para Para);
para SampleZ(datobj DatObj, para Para);
para SampleAlpha(datobj DatObj, para Para);
para SampleKappa(datobj DatObj, para Para, hypara HyPara);
para SampleDelta(datobj DatObj, para Para, hypara HyPara);
para SampleEta(datobj DatObj, para Para, hypara HyPara);
para SampleUpsilon(datobj DatObj, para Para, hypara HyPara);
std::pair<para, metrobj> SamplePsi(datobj DatObj, para Para, hypara HyPara, metrobj MetrObj);
para SampleSigma2(datobj DatObj, para Para, hypara HyPara);
std::pair<datobj, para> SampleY(datobj DatObj, para Para, dataug DatAug);
para SampleU(datobj DatObj, para Para);
std::pair<para, metrobj> SampleRho(datobj DatObj, para Para, hypara HyPara, metrobj MetrObj);
para SampleBeta(datobj DatObj, para Para, hypara HyPara);
  
//MCMC UTILITY FUNCTIONS
void BeginBurnInProgress(mcmcobj McmcObj, bool Interactive);
Rcpp::List OutputMetrObj(metrobj MetrObj, datobj DatObj);
metrobj PilotAdaptation(metrobj MetrObj, mcmcobj McmcObj, datobj DatObj);
void SamplerProgress(int s, mcmcobj McmcObj);
arma::colvec StoreSamples(datobj DatObj, para Para);
void UpdateBurnInBar(int s, mcmcobj McmcObj);
void UpdateBurnInBarInt(int s, mcmcobj McmcObj);

//UTILITY FUNCTIONS
arma::mat GetRooti(arma::mat const& Cov, arma::mat const& Eye);
arma::mat CholInv(arma::mat const& Cov);
arma::mat GetLambda(arma::mat const& Theta, arma::umat const& Xi, int K, int M, int O);
arma::cube GetlogWeights(arma::cube const& Alpha, int K, int M, int L, int O);
arma::cube GetWeights(arma::cube const& Alpha, int K, int M, int L, int O);
bool rows_equal(arma::mat const& lhs, arma::mat const& rhs, double tol);
arma::colvec GetLStarJ(arma::mat const& U, arma::cube const& Weights, int K, int M, int O);

#endif // __spBFA__
