#include <RcppArmadillo.h>
#include "DIAG_diagnostics.h"

//Function to convert Rcpp::List DatObj to a custom C++ struct datobj--------------------------------------------------
datobjDIAG ConvertDatObjDIAG(Rcpp::List DatObj_List) {

  //Set objects from List
  double Rho = DatObj_List["Rho"];
  double ScaleY = DatObj_List["ScaleY"];
  double ScaleDM = DatObj_List["ScaleDM"];
  int M = DatObj_List["M"];
  int Nu = DatObj_List["Nu"];
  int FamilyInd = DatObj_List["FamilyInd"];
  int WeightsInd = DatObj_List["WeightsInd"];
  arma::vec YObserved = DatObj_List["YObserved"];
  arma::mat YStarWide = DatObj_List["YStarWide"];
  arma::mat W = DatObj_List["W"];
  arma::mat UpperW = arma::trimatu(W);
  arma::umat AdjacentEdgesBoolean = find(UpperW == 1);
  arma::vec OneM = DatObj_List["OneM"];
  arma::mat EyeM = DatObj_List["EyeM"];
  arma::vec Z = DatObj_List["Z"];

  //Convert to C++ struct
  datobjDIAG DatObj;
  DatObj.Rho = Rho;
  DatObj.ScaleY = ScaleY;
  DatObj.ScaleDM = ScaleDM;
  DatObj.M = M;
  DatObj.Nu = Nu;
  DatObj.FamilyInd = FamilyInd;
  DatObj.WeightsInd = WeightsInd;
  DatObj.YObserved = YObserved;
  DatObj.YStarWide = YStarWide;
  DatObj.W = W;
  DatObj.AdjacentEdgesBoolean = AdjacentEdgesBoolean;
  DatObj.OneM = OneM;
  DatObj.EyeM = EyeM;
  DatObj.Z = Z;
  return DatObj;

}



//Function to convert Rcpp::List Para to a custom C++ struct para-----------------------------------------------------
paraDIAG ConvertParaDIAG(Rcpp::List Para_List) {

  //Set objects from List
  arma::mat Mu = Para_List["Mu"];
  arma::mat Tau2 = Para_List["Tau2"];
  arma::mat Alpha = Para_List["Alpha"];
  arma::vec MuMean = Para_List["MuMean"];
  arma::cube CovMean = Para_List["CovMean"];

  //Convert to C++ struct
  paraDIAG Para;
  Para.Mu = Mu;
  Para.Tau2 = Tau2;
  Para.Alpha = Alpha;
  Para.MuMean = MuMean;
  Para.CovMean = CovMean;
  return Para;
}



//Function to convert Rcpp::List DatAug to a custom C++ struct dataug-----------------------------------------------------
dataugDIAG ConvertDatAugDIAG(Rcpp::List DatAug_List, datobjDIAG DatObj) {

  //Set data objects
  arma::mat YStarWide = DatObj.YStarWide;
  int Nu = DatObj.Nu;

  //Set objects from List
  int NBelow = DatAug_List["NBelow"];
  arma::vec NBelowCount = DatAug_List["NBelowCount"];
  Rcpp::List YStarNonZeroList = DatAug_List["YStarNonZero"];
  arma::field<arma::vec> YStarNonZero(Nu, 1);
  arma::field<arma::uvec> NBelowBoolean(Nu, 1), NAboveBoolean(Nu, 1);
  for (int i = 0; i < Nu; i++) {
    YStarNonZero(i, 0) = Rcpp::as<arma::vec>(YStarNonZeroList[i]);
    NBelowBoolean(i, 0) = find(YStarWide.col(i) <= 0);
    NAboveBoolean(i, 0) = find(YStarWide.col(i) > 0);
  }

  //Convert to C++ struct
  dataugDIAG DatAug;
  DatAug.NBelow = NBelow;
  DatAug.NBelowCount = NBelowCount;
  DatAug.NBelowBoolean = NBelowBoolean;
  DatAug.NAboveBoolean = NAboveBoolean;
  DatAug.YStarNonZero = YStarNonZero;
  return DatAug;
}
