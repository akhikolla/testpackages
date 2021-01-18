//
//  functions used in spBDwDM package for predictions
//

#ifndef __predictions__
#define __predictions__

//STRUCT DEFINITIONS
struct datobjPRED {
  int M;
  int K;
  int Nu;
  int O;
  int C;
  int P;
  arma::Col<int> FamilyInd;
  int TempCorInd;
  int NNewVisits;
  arma::mat EyeK;
  arma::mat TimeDist;
  arma::uvec NewVisits;
  arma::uvec OriginalVisits;
  arma::cube Trials;
  arma::mat NewX;
};
struct paraPRED {
  arma::mat Psi;
  arma::mat Upsilon;
  arma::mat Lambda;
  arma::mat Sigma2;
  arma::mat Eta;
  arma::mat Beta;
};

//COVARIANCE FUNCTIONS
arma::mat H(double Psi, int TempCorInd, arma::mat const& TimeDist, int Nu);
arma::mat SpEXP(double Rho, arma::mat const& SpDist, int M);

//PREDICTION FUNCTIONS
arma::mat EtaKrigging(Rcpp::List DatObj_List, Rcpp::List Para_List, int NKeep, bool Verbose);
arma::cube YKrigging(Rcpp::List DatObj_List, Rcpp::List Para_List, arma::mat EtaKrig, int NKeep, bool Verbose);
  
//DISTRIBUTION FUNCTIONS
arma::mat rmvnormRcpp(int n, arma::vec const& mean, arma::mat const& sigma);
double rbinomRcpp(double n, double p);
arma::vec rnormVecRcpp(arma::vec const& mean, arma::vec const& sd);

//PREDICTION CONVERSION FUNCTIONS
datobjPRED ConvertDatObjPRED(Rcpp::List DatObj_List);
paraPRED ConvertParaPRED(Rcpp::List Para_List);

//UTILITY FUNCTIONS
arma::mat CholInv(arma::mat const& Cov);

#endif // __predictions__
