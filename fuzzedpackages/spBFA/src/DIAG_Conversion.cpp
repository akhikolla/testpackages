#include <RcppArmadillo.h>
#include "DIAG_diagnostics.h"

//Function to convert Rcpp::List DatObj to a custom C++ struct datobj--------------------------------------------------
datobjDIAG ConvertDatObjDIAG(Rcpp::List DatObj_List) {

  //Set objects from List
  int M = DatObj_List["M"];
  int O = DatObj_List["O"];
  int Nu = DatObj_List["Nu"];
  int N = DatObj_List["N"];
  int K = DatObj_List["K"];
  int P = DatObj_List["P"];
  arma::vec YObserved  = DatObj_List["YObserved"];
  arma::mat X = DatObj_List["X"];
  arma::mat EyeNu = DatObj_List["EyeNu"];
  arma::cube Trials = DatObj_List["Trials"];
  arma::Col<int> FamilyInd = DatObj_List["FamilyInd"];
  
  //Convert to C++ struct
  datobjDIAG DatObj;
  DatObj.M = M;
  DatObj.Nu = Nu;
  DatObj.FamilyInd = FamilyInd;
  DatObj.YObserved = YObserved;
  DatObj.N = N;
  DatObj.EyeNu = EyeNu;
  DatObj.O = O;
  DatObj.K = K;
  DatObj.P = P;
  DatObj.X = X;
  DatObj.Trials = Trials;
  return DatObj;

}



//Function to convert Rcpp::List Para to a custom C++ struct para-----------------------------------------------------
paraDIAG ConvertParaDIAG(Rcpp::List Para_List) {

  //Set objects from List
  arma::mat Eta = Para_List["Eta"];
  arma::mat Beta = Para_List["Beta"];
  arma::mat Lambda = Para_List["Lambda"];
  arma::mat Sigma2 = Para_List["Sigma2"];
  arma::cube MuMean = Para_List["MuMean"];
  arma::cube CovMean = Para_List["CovMean"];

  //Convert to C++ struct
  paraDIAG Para;
  Para.Beta = Beta;
  Para.Lambda = Lambda;
  Para.Eta = Eta;
  Para.MuMean = MuMean;
  Para.CovMean = CovMean;
  Para.Sigma2 = Sigma2;
  return Para;
}
