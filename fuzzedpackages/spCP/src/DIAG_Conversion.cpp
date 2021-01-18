#include <RcppArmadillo.h>
#include "DIAG_diagnostics.h"

//Function to convert Rcpp::List DatObj to a custom C++ struct datobj--------------------------------------------------
datobjDIAG ConvertDatObjDIAG(Rcpp::List DatObj_List) {

  //Set objects from List
  double tNu = DatObj_List["tNu"];
  double t1 = DatObj_List["t1"];
  double ScaleY = DatObj_List["ScaleY"];
  double ScaleDM = DatObj_List["ScaleDM"];
  int M = DatObj_List["M"];
  int N = DatObj_List["N"];
  int Nu = DatObj_List["Nu"];
  int FamilyInd = DatObj_List["FamilyInd"];
  arma::vec YObserved = DatObj_List["YObserved"];
  arma::mat YStarWide = DatObj_List["YStarWide"];
  arma::vec OneM = DatObj_List["OneM"];
  arma::mat EyeM = DatObj_List["EyeM"];
  arma::uvec XThetaInd = DatObj_List["XThetaInd"];
  arma::mat OneN = DatObj_List["OneN"];
  arma::mat EyeN = DatObj_List["EyeN"];
  arma::mat EyeNu = DatObj_List["EyeNu"];
  arma::vec TimeVec = DatObj_List["TimeVec"];
  arma::vec OneNu = DatObj_List["OneNu"];

  //Convert to C++ struct
  datobjDIAG DatObj;
  DatObj.OneN = OneN;
  DatObj.EyeN = EyeN;
  DatObj.TimeVec = TimeVec;
  DatObj.OneNu = OneNu;
  DatObj.ScaleY = ScaleY;
  DatObj.ScaleDM = ScaleDM;
  DatObj.M = M;
  DatObj.Nu = Nu;
  DatObj.FamilyInd = FamilyInd;
  DatObj.YObserved = YObserved;
  DatObj.YStarWide = YStarWide;
  DatObj.OneM = OneM;
  DatObj.EyeM = EyeM;
  DatObj.tNu = tNu;
  DatObj.t1 = t1;
  DatObj.XThetaInd = XThetaInd;
  DatObj.N = N;
  DatObj.EyeNu = EyeNu;
  return DatObj;

}



//Function to convert Rcpp::List Para to a custom C++ struct para-----------------------------------------------------
paraDIAG ConvertParaDIAG(Rcpp::List Para_List) {

  //Set objects from List
  arma::mat Beta0 = Para_List["Beta0"];
  arma::mat Beta1 = Para_List["Beta1"];
  arma::mat Lambda0 = Para_List["Lambda0"];
  arma::mat Lambda1 = Para_List["Lambda1"];
  arma::mat Eta = Para_List["Eta"];
  arma::mat MuMean = Para_List["MuMean"];
  arma::mat Sigma2Mean = Para_List["Sigma2Mean"];

  //Convert to C++ struct
  paraDIAG Para;
  Para.Beta0 = Beta0;
  Para.Beta1 = Beta1;
  Para.Lambda0 = Lambda0;
  Para.Lambda1 = Lambda1;
  Para.Eta = Eta;
  Para.MuMean = MuMean;
  Para.Sigma2Mean = Sigma2Mean;
  return Para;
}
