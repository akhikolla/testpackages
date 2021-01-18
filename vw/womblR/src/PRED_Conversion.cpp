#include <RcppArmadillo.h>
#include "PRED_predictions.h"

//Function to convert Rcpp::List DatObj to a custom C++ struct datobj--------------------------------------------------
datobjPRED ConvertDatObjPRED(Rcpp::List DatObj_List) {

  //Set objects from List
  double Rho = DatObj_List["Rho"];
  double ScaleY = DatObj_List["ScaleY"];
  double ScaleDM = DatObj_List["ScaleDM"];
  int M = DatObj_List["M"];
  int Nu = DatObj_List["Nu"];
  int FamilyInd = DatObj_List["FamilyInd"];
  int NNewVisits = DatObj_List["NNewVisits"];
  int TempCorInd = DatObj_List["TempCorInd"];
  int WeightsInd = DatObj_List["WeightsInd"];
  arma::vec YObserved = DatObj_List["YObserved"];
  arma::mat YStarWide = DatObj_List["YStarWide"];
  arma::mat W = DatObj_List["W"];
  arma::mat UpperW = arma::trimatu(W);
  arma::umat AdjacentEdgesBoolean = find(UpperW == 1);
  arma::vec OneM = DatObj_List["OneM"];
  arma::vec OneNu = DatObj_List["OneNu"];
  arma::mat EyeM = DatObj_List["EyeM"];
  arma::vec Z = DatObj_List["Z"];
  arma::mat TimeDist = DatObj_List["TimeDist"];
  arma::uvec NewVisits = DatObj_List["NewVisits"];
  arma::uvec OriginalVisits = DatObj_List["OriginalVisits"];

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
  DatObj.Z = Z;
  DatObj.TimeDist = TimeDist;
  DatObj.TempCorInd = TempCorInd;
  DatObj.WeightsInd = WeightsInd;
  DatObj.NNewVisits = NNewVisits;
  DatObj.NewVisits = NewVisits;
  DatObj.OriginalVisits = OriginalVisits;
  return DatObj;

}



//Function to convert Rcpp::List Para to a custom C++ struct para-----------------------------------------------------
paraPRED ConvertParaPRED(Rcpp::List Para_List) {

  //Set objects from List
  arma::mat Mu = Para_List["Mu"];
  arma::mat Tau2 = Para_List["Tau2"];
  arma::mat Alpha = Para_List["Alpha"];
  arma::mat Delta = Para_List["Delta"];
  arma::mat T = Para_List["T"];
  arma::mat Phi = Para_List["Phi"];

  //Convert to C++ struct
  paraPRED Para;
  Para.Mu = Mu;
  Para.Tau2 = Tau2;
  Para.Alpha = Alpha;
  Para.Delta = Delta;
  Para.T = T;
  Para.Phi = Phi;
  return Para;
}
