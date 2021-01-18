//
//  functions used in spBDwDM package for predictions
//

#ifndef __predictions__
#define __predictions__

//STRUCT DEFINITIONS
struct datobjPRED {
  double Rho;
  double ScaleY;
  double ScaleDM;
  int M;
  int Nu;
  int FamilyInd;
  int TempCorInd;
  int NNewVisits;
  int WeightsInd;
  arma::vec YObserved;
  arma::mat YStarWide;
  arma::mat W;
  arma::umat AdjacentEdgesBoolean;
  arma::vec OneNu;
  arma::vec OneM;
  arma::mat EyeM;
  arma::mat TimeDist;
  arma::uvec NewVisits;
  arma::uvec OriginalVisits;
  arma::vec Z;
};
struct paraPRED {
  arma::mat Mu;
  arma::mat Tau2;
  arma::mat Alpha;
  arma::mat Delta;
  arma::mat T;
  arma::mat Phi;
};

//COVARIANCE FUNCTIONS
arma::cube JointCovarianceCube(arma::cube const& WAlphas, arma::vec const& Tau2, arma::mat const& EyeM, double Rho, int M, int Nu);
arma::cube WAlphaCube(arma::vec const& Alpha, arma::colvec const& Z, arma::mat const& W, int M, int Nu, int WeightsInd);
arma::mat SIGMA(double Phi, int TempCorInd, arma::mat const& DayDist, int Nu);

//PREDICTION FUNCTIONS
arma::mat ThetaKrigging(Rcpp::List DatObj_List, Rcpp::List Para_List, int NKeep);
arma::mat YKrigging(Rcpp::List DatObj_List, arma::mat ThetaKrig, int NKeep);

//DISTRIBUTION FUNCTIONS
arma::mat rmvnormRcpp(int n, arma::vec const& mean, arma::mat const& sigma);

//PREDICTION CONVERSION FUNCTIONS
datobjPRED ConvertDatObjPRED(Rcpp::List DatObj_List);
paraPRED ConvertParaPRED(Rcpp::List Para_List);

//UTILITY FUNCTIONS
arma::mat CholInv(arma::mat const& Cov);

#endif // __predictions__
