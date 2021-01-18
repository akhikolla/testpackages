#include <RcppArmadillo.h>
#include "PRED_predictions.h"

//Function to convert Rcpp::List DatObj to a custom C++ struct datobj--------------------------------------------------
datobjPRED ConvertDatObjPRED(Rcpp::List DatObj_List) {

  //Set objects from List
  int M = DatObj_List["M"];
  int K = DatObj_List["K"];
  int Nu = DatObj_List["Nu"];
  int O = DatObj_List["O"];
  int C = DatObj_List["C"];
  int P = DatObj_List["P"];
  arma::Col<int> FamilyInd = DatObj_List["FamilyInd"];
  int NNewVisits = DatObj_List["NNewVisits"];
  int TempCorInd = DatObj_List["TempCorInd"];
  arma::mat EyeK = DatObj_List["EyeK"];
  arma::mat TimeDist = DatObj_List["TimeDist"];
  arma::uvec NewVisits = DatObj_List["NewVisits"];
  arma::uvec OriginalVisits = DatObj_List["OriginalVisits"];
  arma::cube Trials = DatObj_List["Trials"];
  arma::mat NewX = DatObj_List["NewX"];

  //Convert to C++ struct
  datobjPRED DatObj;
  DatObj.M = M;
  DatObj.K = K;
  DatObj.Nu = Nu;
  DatObj.P = P;
  DatObj.FamilyInd = FamilyInd;
  DatObj.O = O;
  DatObj.C = C;
  DatObj.EyeK = EyeK;
  DatObj.TimeDist = TimeDist;
  DatObj.TempCorInd = TempCorInd;
  DatObj.NNewVisits = NNewVisits;
  DatObj.NewVisits = NewVisits;
  DatObj.OriginalVisits = OriginalVisits;
  DatObj.Trials = Trials;
  DatObj.NewX = NewX;
  return DatObj;

}



//Function to convert Rcpp::List Para to a custom C++ struct para-----------------------------------------------------
paraPRED ConvertParaPRED(Rcpp::List Para_List) {

  //Set objects from List
  arma::mat Psi = Para_List["Psi"];
  arma::mat Upsilon = Para_List["Upsilon"];
  arma::mat Lambda = Para_List["Lambda"];
  arma::mat Sigma2 = Para_List["Sigma2"];
  arma::mat Eta = Para_List["Eta"];
  arma::mat Beta = Para_List["Beta"];

  //Convert to C++ struct
  paraPRED Para;
  Para.Psi = Psi;
  Para.Upsilon = Upsilon;
  Para.Lambda = Lambda;
  Para.Sigma2 = Sigma2;
  Para.Eta = Eta;
  Para.Beta = Beta;
  return Para;
}
