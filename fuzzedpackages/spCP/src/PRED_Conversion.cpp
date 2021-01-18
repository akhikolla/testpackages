#include <RcppArmadillo.h>
#include "PRED_predictions.h"

//Function to convert Rcpp::List DatObj to a custom C++ struct datobj--------------------------------------------------
datobjPRED ConvertDatObjPRED(Rcpp::List DatObj_List) {

  //Set objects from List
  double Rho = DatObj_List["Rho"];
  double ScaleY = DatObj_List["ScaleY"];
  double ScaleDM = DatObj_List["ScaleDM"];
  double t1 = DatObj_List["t1"];
  double tNu = DatObj_List["tNu"];
  int M = DatObj_List["M"];
  int Nu = DatObj_List["Nu"];
  int FamilyInd = DatObj_List["FamilyInd"];
  int NNewTimes = DatObj_List["NNewTimes"];
  arma::vec YObserved = DatObj_List["YObserved"];
  arma::mat YStarWide = DatObj_List["YStarWide"];
  arma::mat W = DatObj_List["W"];
  arma::mat UpperW = arma::trimatu(W);
  arma::uvec AdjacentEdgesBoolean = find(UpperW == 1);
  arma::vec OneM = DatObj_List["OneM"];
  arma::vec OneNu = DatObj_List["OneNu"];
  arma::mat EyeM = DatObj_List["EyeM"];
  arma::vec NewTimes = DatObj_List["NewTimes"];
  arma::vec TimeVec = DatObj_List["TimeVec"];
  arma::uvec XThetaInd = DatObj_List["XThetaInd"];

  //Convert to C++ struct
  datobjPRED DatObj;
  DatObj.Rho = Rho;
  DatObj.ScaleY = ScaleY;
  DatObj.ScaleDM = ScaleDM;
  DatObj.M = M;
  DatObj.Nu = Nu;
  DatObj.FamilyInd = FamilyInd;
  DatObj.YObserved = YObserved;
  DatObj.YStarWide = YStarWide;
  DatObj.W = W;
  DatObj.AdjacentEdgesBoolean = AdjacentEdgesBoolean;
  DatObj.OneM = OneM;
  DatObj.OneNu = OneNu;
  DatObj.EyeM = EyeM;
  DatObj.NNewTimes = NNewTimes;
  DatObj.NewTimes = NewTimes;
  DatObj.t1 = t1;
  DatObj.tNu = tNu;
  DatObj.TimeVec = TimeVec;
  DatObj.XThetaInd = XThetaInd;
  return DatObj;

}



//Function to convert Rcpp::List Para to a custom C++ struct para-----------------------------------------------------
paraPRED ConvertParaPRED(Rcpp::List Para_List) {

  //Set objects from List
  arma::mat Beta0 = Para_List["Beta0"];
  arma::mat Beta1 = Para_List["Beta1"];
  arma::mat Lambda0 = Para_List["Lambda0"];
  arma::mat Lambda1 = Para_List["Lambda1"];
  arma::mat Eta = Para_List["Eta"];

  //Convert to C++ struct
  paraPRED Para;
  Para.Beta0 = Beta0;
  Para.Beta1 = Beta1;
  Para.Lambda0 = Lambda0;
  Para.Lambda1 = Lambda1;
  Para.Eta = Eta;
  return Para;

}
