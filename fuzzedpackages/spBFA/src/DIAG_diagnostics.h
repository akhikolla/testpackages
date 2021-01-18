//
//  functions used in spCP package for dianostics
//

#ifndef __diagnostics__
#define __diagnostics__

//STRUCT DEFINITIONS
struct datobjDIAG {
  int M;
  int O;
  int Nu;
  int N;
  int K;
  int P;
  arma::vec YObserved;
  arma::mat X;
  arma::mat EyeNu;
  arma::cube Trials;
  arma::Col<int> FamilyInd;
};
struct paraDIAG {
  arma::mat Lambda;
  arma::mat Eta;
  arma::mat Beta;
  arma::mat Sigma2;
  arma::cube MuMean;
  arma::cube CovMean;
};

//DIAGNOSTIC FUNCTIONS
arma::colvec GetLogLik(Rcpp::List DatObj_List, Rcpp::List Para_List, int NKeep, bool Verbose);
double GetLogLikMean(Rcpp::List DatObj_List, Rcpp::List Para_List);

//DISTRIBUTION FUNCTIONS
double dlnorm(double x, double mu, double sigma2);
arma::colvec dlnormRcpp(arma::vec const& x, arma::vec const& mu, arma::vec const& sigma2);
double pnormRcpp(double q);
double dlbinom(int x, int n, double pi);
arma::vec rnormVecRcpp(arma::vec const& mean, arma::vec const& sd);
double rbinomRcpp(double n, double p);

//DIAGNOSTIC CONVERSION FUNCTIONS
datobjDIAG ConvertDatObjDIAG(Rcpp::List DatObj_List);
paraDIAG ConvertParaDIAG(Rcpp::List Para_List);

//UTILITY FUNCTIONS
arma::mat CholInv(arma::mat const& Cov);

#endif // __diagnostics__
