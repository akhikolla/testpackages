//
//  functions used in spCP package for dianostics
//

#ifndef __diagnostics__
#define __diagnostics__

//STRUCT DEFINITIONS
struct datobjDIAG {
  double Rho;
  double ScaleY;
  double ScaleDM;
  double tNu;
  double t1;
  int M;
  int Nu;
  int N;
  int FamilyInd;
  arma::vec YObserved;
  arma::mat YStarWide;
  arma::vec OneM;
  arma::mat EyeM;
  arma::uvec XThetaInd;
  arma::vec OneN;
  arma::mat EyeN;
  arma::vec TimeVec;
  arma::vec OneNu;
  arma::mat EyeNu;
};
struct paraDIAG {
  arma::mat Beta0;
  arma::mat Beta1;
  arma::mat Lambda0;
  arma::mat Lambda1;
  arma::mat Eta;
  arma::vec MuMean;
  arma::vec Sigma2Mean;
};

//COVARIANCE FUNCTIONS
arma::mat GetXTheta(arma::vec const& Theta, arma::uvec const& XThetaInd, arma::vec const& TimeVec, arma::vec const& OneNu, arma::vec const& OneN, double tNu, int N, int M);
arma::mat GetXTheta_lmc(arma::vec const& Theta, arma::uvec const& XThetaInd, arma::vec const& TimeVec, arma::vec const& OneNu, arma::vec const& OneN, double tNu, int N, int M);

//DIAGNOSTIC FUNCTIONS
arma::colvec GetLogLik(Rcpp::List DatObj_List, Rcpp::List Para_List, int NKeep);
double GetLogLikMean(Rcpp::List DatObj_List, Rcpp::List Para_List);
arma::colvec NormalLogLik(datobjDIAG DatObj, paraDIAG Para, int NKeep);
double NormalLogLikMean(datobjDIAG DatObj, paraDIAG Para);
arma::colvec TobitLogLik(datobjDIAG DatObj, paraDIAG Para, int NKeep);
double TobitLogLikMean(datobjDIAG DatObj, paraDIAG Para);
arma::colvec ProbitLogLik(datobjDIAG DatObj, paraDIAG Para, int NKeep);
double ProbitLogLikMean(datobjDIAG DatObj, paraDIAG Para);

//DISTRIBUTION FUNCTIONS
double dlnorm(double x, double mu, double sigma2);
double lndMvn(arma::vec const& Y, arma::vec const& Mu, arma::mat const& Rooti);
double pnormRcpp(double q);
arma::mat rmvnormRcpp(int n, arma::vec const& mean, arma::mat const& sigma);

//DIAGNOSTIC CONVERSION FUNCTIONS
datobjDIAG ConvertDatObjDIAG(Rcpp::List DatObj_List);
paraDIAG ConvertParaDIAG(Rcpp::List Para_List);

//UTILITY FUNCTIONS
arma::mat CholInv(arma::mat const& Cov);

#endif // __diagnostics__
