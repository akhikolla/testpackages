#include <RcppArmadillo.h>
#include "DIAG_diagnostics.h"

//Function that samples from the posterior predictive distribution for the spBDwDM model------------------
//This function is being exported to R for use in this package exclusively...
//not for use by users.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat SamplePPD(Rcpp::List DatObj_List, Rcpp::List Para_List, int NKeep, bool Verbose) {

  //Convet Rcpp::Lists to C++ structs
  datobjDIAG DatObj = ConvertDatObjDIAG(DatObj_List);
  paraDIAG Para = ConvertParaDIAG(Para_List);

  //Set data objects
  arma::vec YObserved = DatObj.YObserved;
  arma::Col<int> FamilyInd = DatObj.FamilyInd;
  arma::mat X = DatObj.X;
  int K = DatObj.K;
  int M = DatObj.M;
  int O = DatObj.O;
  int N = DatObj.N;
  int Nu = DatObj.Nu;
  int P = DatObj.P;
  arma::mat EyeNu = DatObj.EyeNu;
  arma::cube Trials = DatObj.Trials;
  
  //Set parameters
  arma::mat BetaMat = Para.Beta;
  arma::mat EtaMat = Para.Eta;
  arma::mat LambdaMat = Para.Lambda;
  arma::mat Sigma2Mat = Para.Sigma2;

  //Verbose output
  arma::vec VerboseSeq;
  if (Verbose) {
    VerboseSeq << 0.25 << 0.50 << 0.75;
    VerboseSeq *= NKeep;
    Rcpp::Rcout << std::fixed << "Calculating PPD: 0%.. ";
  }
  
  //Initialize objects
  arma::colvec Mean(N), MeanVec(M * Nu), YStar(M * Nu), YMax(M * Nu, arma::fill::zeros), Eta(K * Nu), Beta(P);
  arma::cube YStarCube(M, O, Nu);
  arma::mat Lambda(M * O, Nu), MeanMat(M, Nu), TrialsMat(M, Nu), PiMat(M, Nu), YStarMat(M, Nu), PPD(N, NKeep);
  int count;
  arma::umat ProbitOnes;
  
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
    count = 0;
    
    //Loop over spatial observation types
    for (arma::uword f = 0; f < O; f++) {
      
      //Set family type
      int FamilyType = FamilyInd[f];
      
      //Normal based
      if ((FamilyType == 0) | (FamilyType == 1) | (FamilyType == 2)) {
        arma::mat Sigma2 = arma::reshape(Sigma2Mat.row(s), O, M);
        arma::colvec SDVec = arma::sqrt(Sigma2.row(count).t());
        arma::colvec SD = arma::repmat(SDVec, Nu, 1);
        count++;
        MeanVec = arma::vectorise(MeanOut(arma::span::all, arma::span(f, f), arma::span::all));
        YStar = rnormVecRcpp(MeanVec, SD);
        if (FamilyType == 1) { // Probit
          YStar = arma::max(YStar, YMax);
          ProbitOnes = find(YStar > 0);
          YStar(ProbitOnes) = arma::ones<arma::vec>(ProbitOnes.size());
        }
        if (FamilyType == 2) YStar = arma::max(YStar, YMax); // Tobit
        YStarCube(arma::span::all, arma::span(f, f), arma::span::all) = arma::reshape(YStar, M, Nu);
      }
      
      //Categorical
      if (FamilyType == 3) {
        MeanMat = MeanOut(arma::span::all, arma::span(f, f), arma::span::all);
        PiMat = arma::exp(MeanMat) / (1 + arma::exp(MeanMat));
        TrialsMat = Trials(arma::span::all, arma::span(f, f), arma::span::all);
        for (arma::uword t = 0; t < Nu; t++) {
          for (arma::uword i = 0; i < M; i++) {
            YStarMat(i, t) = rbinomRcpp(TrialsMat(i, t), PiMat(i, t));
          }
        }
        YStarCube(arma::span::all, arma::span(f, f), arma::span::all) = YStarMat;
      }
      
    //End loop over spatial observation types 
    }
    

    //Reshape and save ppd
    PPD.col(s) = arma::vectorise(YStarCube);

    //Add a new percentage
    if (Verbose) {
      Rcpp::Rcout.precision(0);
      if (std::find(VerboseSeq.begin(), VerboseSeq.end(), s) != VerboseSeq.end())
        Rcpp::Rcout << std::fixed << 100 * (s) / NKeep << "%.. ";
    }
    
  //End loop
  }

  //Output final percentage
  if (Verbose) Rcpp::Rcout << std::fixed << "100%.. Done!" << std::endl;

  //Return PPD
  return PPD;

}
