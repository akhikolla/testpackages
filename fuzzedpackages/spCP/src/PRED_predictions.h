//
//  functions used in spCP package for predictions
//

#ifndef __predictions__
#define __predictions__

//STRUCT DEFINITIONS
struct datobjPRED {
  double Rho;
  double ScaleY;
  double ScaleDM;
  double tNu;
  double t1;
  int M;
  int Nu;
  int FamilyInd;
  int NNewTimes;
  arma::vec YObserved;
  arma::mat YStarWide;
  arma::mat W;
  arma::uvec AdjacentEdgesBoolean;
  arma::vec OneNu;
  arma::vec OneM;
  arma::mat EyeM;
  arma::vec NewTimes;
  arma::vec TimeVec;
  arma::uvec XThetaInd;
};
struct paraPRED {
  arma::mat Beta0;
  arma::mat Beta1;
  arma::mat Lambda0;
  arma::mat Lambda1;
  arma::mat Eta;
};

//PREDICTION FUNCTIONS
arma::cube PredictFuture(Rcpp::List DatObj_List, Rcpp::List Para_List, int NKeep);

//DISTRIBUTION FUNCTIONS
arma::vec rnormRcpp(int n, double mean, double sd);

//PREDICTION CONVERSION FUNCTIONS
datobjPRED ConvertDatObjPRED(Rcpp::List DatObj_List);
paraPRED ConvertParaPRED(Rcpp::List Para_List);

//UTILITY FUNCTIONS
arma::mat CholInv(arma::mat const& Cov);

#endif // __predictions__
