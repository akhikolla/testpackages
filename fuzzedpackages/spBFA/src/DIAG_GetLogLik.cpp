#include <RcppArmadillo.h>
#include "DIAG_diagnostics.h"
#include <math.h>

//Function that computes the log-likelihood for spBDwDM model----------------------------------------------------------------
//This function is being exported to R for use in this package exclusively...
//not for use by users.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::colvec GetLogLik(Rcpp::List DatObj_List, Rcpp::List Para_List, int NKeep, bool Verbose) {

  //Convert Rcpp::Lists to C++ structs
  datobjDIAG DatObj = ConvertDatObjDIAG(DatObj_List);
  paraDIAG Para = ConvertParaDIAG(Para_List);

  //Set data object
  arma::Col<int> FamilyInd = DatObj.FamilyInd;
  int M = DatObj.M;
  int O = DatObj.O;
  int Nu = DatObj.Nu;
  int N = DatObj.N;
  int K = DatObj.K;
  int P = DatObj.P;
  arma::vec YObserved = DatObj.YObserved;
  arma::mat X = DatObj.X;
  arma::mat EyeNu = DatObj.EyeNu;
  arma::cube Trials = DatObj.Trials;
  
  //Set parameters
  arma::mat BetaMat = Para.Beta;
  arma::mat LambdaMat = Para.Lambda;
  arma::mat EtaMat = Para.Eta;
  arma::mat Sigma2Mat = Para.Sigma2;

  //Declarations
  arma::cube MeanOut(N, 1, 1), YCube(N, 1, 1);
  arma::mat Lambda(M * O, K);
  arma::colvec Beta(P), Eta(Nu * K), Mean(N), MeanVec(Nu * M), YVec(Nu * M), LogLik(NKeep), PiVec(Nu * M), TrialsVec(Nu * M);
  int count = 0;
  double YS, SDS, MeanS, LogLiks, Pi;
  int YSInt, TrialInt;
  YCube(arma::span::all, arma::span(0, 0), arma::span(0, 0)) = YObserved;
  arma::cube Y = arma::reshape(YCube, M, O, Nu);
  arma::vec PNormS(1);
  
  //Verbose output
  arma::vec VerboseSeq;
  if (Verbose) {
    VerboseSeq << 0.25 << 0.50 << 0.75;
    VerboseSeq *= NKeep;
    Rcpp::Rcout << std::fixed << "Calculating Log-Lik: 0%.. ";
  }
    
  //Loop over scans
  for (int s = 0; s < NKeep; s++) {
    
    //Compute moments and get log-likelihood
    Eta = EtaMat.row(s).t();
    Lambda = arma::reshape(LambdaMat.row(s), K, M * O).t();
    Beta = BetaMat.row(s).t();
    Mean = arma::kron(EyeNu, Lambda) * Eta + X * Beta;
    arma::cube MeanOut(N, 1, 1);
    MeanOut(arma::span::all, arma::span(0, 0), arma::span(0, 0)) = Mean;
    MeanOut = arma::reshape(MeanOut, M, O, Nu);
    LogLiks = 0; count = 0;
    
    //Loop over spatial observation types
    for (arma::uword f = 0; f < O; f++) {
      
      //Set family type
      int FamilyType = FamilyInd[f];
      
      //Normal
      if (FamilyType == 0) {
        arma::mat Sigma2 = arma::reshape(Sigma2Mat.row(s), O, M);
        arma::colvec VarVec = Sigma2.row(count).t();
        arma::colvec Var = arma::repmat(VarVec, Nu, 1);
        count++;
        MeanVec = arma::vectorise(MeanOut(arma::span::all, arma::span(f, f), arma::span::all));
        YVec = arma::vectorise(Y(arma::span::all, arma::span(f, f), arma::span::all));
        LogLiks += arma::sum(dlnormRcpp(YVec, MeanVec, Var));
      }
      
      //Probit
      if (FamilyType == 1) {
        arma::mat Sigma2 = arma::reshape(Sigma2Mat.row(s), O, M);
        arma::colvec SDVec = arma::sqrt(Sigma2.row(count).t());
        arma::colvec SD = arma::repmat(SDVec, Nu, 1);
        count++;
        MeanVec = arma::vectorise(MeanOut(arma::span::all, arma::span(f, f), arma::span::all));
        YVec = arma::vectorise(Y(arma::span::all, arma::span(f, f), arma::span::all));
        for (int i = 0; i < (M * Nu); i++) {
          YS = YVec(i);
          SDS = SD(i);
          MeanS = MeanVec(i);
          PNormS(0) = pnormRcpp(-MeanS / SDS);
          if (YS == 0) LogLiks += arma::as_scalar(arma::log(PNormS));
          if (YS == 1) LogLiks += arma::as_scalar(arma::log(1 - PNormS));
        }
      }
      
      //Tobit
      if (FamilyType == 2) {
        arma::mat Sigma2 = arma::reshape(Sigma2Mat.row(s), O, M);
        arma::colvec SDVec = arma::sqrt(Sigma2.row(count).t());
        arma::colvec SD = arma::repmat(SDVec, Nu, 1);
        count++;
        MeanVec = arma::vectorise(MeanOut(arma::span::all, arma::span(f, f), arma::span::all));
        YVec = arma::vectorise(Y(arma::span::all, arma::span(f, f), arma::span::all));
        for (int i = 0; i < (M * Nu); i++) {
          YS = YVec(i);
          SDS = SD(i);
          MeanS = MeanVec(i);
          PNormS(0) = pnormRcpp(-MeanS / SDS);
          if (YS == 0) LogLiks += arma::as_scalar(arma::log(PNormS));
          if (YS > 0) LogLiks += dlnorm(YS, MeanS, SDS * SDS);
        }
      }
      
      //Categorical
      if (FamilyType == 3) {
        MeanVec = arma::vectorise(MeanOut(arma::span::all, arma::span(f, f), arma::span::all));
        PiVec = arma::exp(MeanVec) / (1 + arma::exp(MeanVec));
        TrialsVec = arma::vectorise(Trials(arma::span::all, arma::span(f, f), arma::span::all));
        YVec = arma::vectorise(Y(arma::span::all, arma::span(f, f), arma::span::all));
        for (int i = 0; i < (M * Nu); i++) {
          YSInt = YVec(i);
          TrialInt = TrialsVec(i);
          Pi = PiVec(i);
          LogLiks += dlbinom(YSInt, TrialInt, Pi);
        }
      }
      
    //End loop over spatial observation types 
    }

    //Save log-likelihood
    LogLik(s) = LogLiks;

    //Add a new percentage
    if (Verbose) {
      Rcpp::Rcout.precision(0);
      if (std::find(VerboseSeq.begin(), VerboseSeq.end(), s) != VerboseSeq.end())
        Rcpp::Rcout << std::fixed << 100 * (s) / NKeep << "%.. ";
    }
    
  //End loop over scans
  }
  
  //Output final percentage
  if (Verbose) Rcpp::Rcout << std::fixed << "100%.. Done!" << std::endl;
  
  //Return log-likelihood
  return LogLik;

}



//Function that computes the log-likelihood for STBDwDM model for the mean parameters------------------------------------------
//This function is being exported to R for use in this package exclusively...
//not for use by users.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double GetLogLikMean(Rcpp::List DatObj_List, Rcpp::List Para_List) {

  //Convert Rcpp::Lists to C++ structs
  datobjDIAG DatObj = ConvertDatObjDIAG(DatObj_List);
  paraDIAG Para = ConvertParaDIAG(Para_List);

  //Set data object
  arma::Col<int> FamilyInd = DatObj.FamilyInd;
  arma::colvec YObserved = DatObj.YObserved;
  int O = DatObj.O;
  int M = DatObj.M;
  int Nu = DatObj.Nu;
  int N = DatObj.N;
  arma::cube Trials = DatObj.Trials;
  
  //Set parameter objects
  arma::cube MuMean = Para.MuMean;
  arma::cube CovMean = Para.CovMean;
  
  //Declarations
  int YSInt, TrialInt;
  double LogLikMean = 0, YS, MuS, SDS, Pi;
  arma::mat MuVec(M, Nu), CovVec(M, Nu), YVec(M, Nu), PiVec(M, Nu);
  arma::vec PNormS(1), SDVec(M * Nu), TrialsVec(Nu * M);
  arma::cube YCube(N, 1, 1);
  YCube(arma::span::all, arma::span(0, 0), arma::span(0, 0)) = YObserved;
  arma::cube Y = arma::reshape(YCube, M, O, Nu);
  
  //Loop over spatial observation types
  for (arma::uword f = 0; f < O; f++) {
    
    //Set family type and extract objects
    int FamilyType = FamilyInd[f];
    MuVec = arma::vectorise(MuMean(arma::span::all, arma::span(f, f), arma::span::all));
    CovVec = arma::vectorise(CovMean(arma::span::all, arma::span(f, f), arma::span::all));
    YVec = arma::vectorise(Y(arma::span::all, arma::span(f, f), arma::span::all));
    
    //Normal
    if (FamilyType == 0) LogLikMean += arma::sum(dlnormRcpp(YVec, MuVec, CovVec));
    if (FamilyType == 1) { //Probit
      for (int i = 0; i < (M * Nu); i++) {
        SDVec = arma::sqrt(arma::vectorise(CovMean(arma::span::all, arma::span(f, f), arma::span::all)));
        YS = YVec(i);
        MuS = MuVec(i);
        SDS = SDVec(i);
        PNormS(0) = pnormRcpp(-MuS / SDS);
        if (YS == 0) LogLikMean += arma::as_scalar(arma::log(PNormS));
        if (YS == 1) LogLikMean += arma::as_scalar(arma::log(1 - PNormS));
      }
    }
    if (FamilyType == 2) { //Tobit
      for (int i = 0; i < (M * Nu); i++) {
        SDVec = arma::sqrt(arma::vectorise(CovMean(arma::span::all, arma::span(f, f), arma::span::all)));
        YS = YVec(i);
        MuS = MuVec(i);
        SDS = SDVec(i);
        PNormS(0) = pnormRcpp(-MuS / SDS);
        if (YS == 0) LogLikMean += arma::as_scalar(arma::log(PNormS));
        if (YS > 1) LogLikMean += dlnorm(YS, MuS, SDS * SDS);
      }
    }
    if (FamilyType == 3) { //Categorical
      PiVec = arma::exp(MuVec) / (1 + arma::exp(MuVec));
      TrialsVec = arma::vectorise(Trials(arma::span::all, arma::span(f, f), arma::span::all));
      for (int i = 0; i < (M * Nu); i++) {
        YSInt = YVec(i);
        TrialInt = TrialsVec(i);
        Pi = PiVec(i);
        LogLikMean += dlbinom(YSInt, TrialInt, Pi);
      }
    }
    
  //End loop over spatial types
  }

  //Return log-likelihood
  return LogLikMean;

}
