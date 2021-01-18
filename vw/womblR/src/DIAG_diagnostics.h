//
//  functions used in spBDwDM package for dianostics
//

#ifndef __diagnostics__
#define __diagnostics__

//STRUCT DEFINITIONS
struct datobjDIAG {
  double Rho;
  double ScaleY;
  double ScaleDM;
  int M;
  int Nu;
  int FamilyInd;
  int WeightsInd;
  arma::vec YObserved;
  arma::mat YStarWide;
  arma::mat W;
  arma::umat AdjacentEdgesBoolean;
  arma::vec OneM;
  arma::mat EyeM;
  arma::vec Z;
};
struct paraDIAG {
  arma::mat Mu;
  arma::mat Tau2;
  arma::mat Alpha;
  arma::vec MuMean;
  arma::cube CovMean;
};
struct dataugDIAG {
  int NBelow;
  arma::vec NBelowCount;
  arma::field<arma::uvec> NBelowBoolean;
  arma::field<arma::uvec> NAboveBoolean;
  arma::field<arma::vec> YStarNonZero;
};

//COVARIANCE FUNCTIONS
arma::mat GetRooti(arma::mat const& Cov, arma::mat const& Eye);
arma::cube JointCovarianceCube(arma::cube const& WAlphas, arma::vec const& Tau2, arma::mat const& EyeM, double Rho, int M, int Nu);
arma::cube WAlphaCube(arma::vec const& Alpha, arma::colvec const& Z, arma::mat const& W, int M, int Nu, int WeightsInd);

//DIAGNOSTIC FUNCTIONS
arma::colvec GetLogLik(Rcpp::List DatObj_List, Rcpp::List Para_List, Rcpp::List DatAug_List, int NKeep);
double GetLogLikMean(Rcpp::List DatObj_List, Rcpp::List Para_List, Rcpp::List DatAug_List);
arma::colvec NormalLogLik(datobjDIAG DatObj, paraDIAG Para, int NKeep);
double NormalLogLikMean(datobjDIAG DatObj, paraDIAG Para);
arma::colvec TobitLogLik(datobjDIAG DatObj, paraDIAG Para, dataugDIAG DatAug, int NKeep);
double TobitLogLikMean(datobjDIAG DatObj, paraDIAG Para, dataugDIAG DatAug);
arma::mat SamplePPD(Rcpp::List DatObj_List, Rcpp::List Para_List, int NKeep);

//DISTRIBUTION FUNCTIONS
double lndMvn(arma::vec const& Y, arma::vec const& Mu, arma::mat const& Rooti);
double pmvnormRcpp(int NBelowVisit, arma::vec const& CondMean, arma::mat const& CondCovRound);
arma::mat rmvnormRcpp(int n, arma::vec const& mean, arma::mat const& sigma);

//DIAGNOSTIC CONVERSION FUNCTIONS
datobjDIAG ConvertDatObjDIAG(Rcpp::List DatObj_List);
paraDIAG ConvertParaDIAG(Rcpp::List Para_List);
dataugDIAG ConvertDatAugDIAG(Rcpp::List DatAug_List, datobjDIAG DatObj);

//UTILITY FUNCTIONS
arma::mat CholInv(arma::mat const& Cov);

#endif // __diagnostics__
