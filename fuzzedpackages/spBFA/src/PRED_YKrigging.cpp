#include <RcppArmadillo.h>
#include "PRED_predictions.h"

//Function that samples from the posterior predictive distribution for the spBDwDM model------------------
//This function is being exported to R for use in this package exclusively...
//not for use by users.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube YKrigging(Rcpp::List DatObj_List, Rcpp::List Para_List, arma::mat EtaKrig, int NKeep, bool Verbose) {

  //Convet Rcpp::Lists to C++ structs
  datobjPRED DatObj = ConvertDatObjPRED(DatObj_List);
  paraPRED Para = ConvertParaPRED(Para_List);
  
  //Set data object
  arma::Col<int> FamilyInd = DatObj.FamilyInd;
  int M = DatObj.M;
  int O = DatObj.O;
  int K = DatObj.K;
  int P = DatObj.P;
  int NNewVisits = DatObj.NNewVisits;
  arma::cube Trials = DatObj.Trials;
  arma::mat NewX = DatObj.NewX;

  //Set parameters
  arma::mat LambdaMat = Para.Lambda;
  arma::mat Sigma2Mat = Para.Sigma2;
  arma::mat BetaMat = Para.Beta;
  
  //Verbose output
  arma::vec VerboseSeq;
  if (Verbose) {
    VerboseSeq << 0.25 << 0.50 << 0.75;
    VerboseSeq *= NKeep;
    Rcpp::Rcout << std::fixed << "Krigging Y: 0%.. ";
  }
  
  //Initialize objects
  arma::cube Out(M * O, NNewVisits, NKeep);
  arma::mat YMax(M, O, arma::fill::zeros), Lambda(M * O, K), OutTemp(M * O, NNewVisits);
  arma::mat BigPhi(K, NNewVisits), Mean(M * O, NNewVisits), MeanMat(M, O), PredMat(M, O);
  arma::colvec MeanVec(M), Pi(M), TrialsVec(M), Beta(P);
  arma::umat ProbitOnes;
  int FamilyType, count;
  
  //Loop over scans
  for (arma::uword s = 0; s < NKeep; s++) {

    //Extract level 1 parameters
    BigPhi = arma::reshape(EtaKrig.col(s), K, NNewVisits);
    Lambda = arma::reshape(LambdaMat.row(s), K, M * O).t();
    Beta = BetaMat.row(s).t();
    
    //Get joint moments
    Mean = Lambda * BigPhi + arma::reshape(NewX * Beta, M * O, NNewVisits);
    
    //Loop over predictions
    for (arma::uword n = 0; n < NNewVisits; n++) {
      
      //Mean for particular prediction
      MeanMat = arma::reshape(Mean.col(n), M, O);
      arma::mat TrialsMat = Trials.slice(n);
      
      //Loop over family dimension
      count = 0;
      for (arma::uword f = 0; f < O; f++) {
      
        //Family type
        FamilyType = FamilyInd(f);
        
        //Non-categorical
        if (FamilyType != 3) {
          arma::mat Sigma2 = arma::reshape(Sigma2Mat.row(s), O, M);
          arma::vec SD = arma::sqrt(Sigma2.row(count).t());
          MeanVec = MeanMat.col(f);
          PredMat.col(f) = rnormVecRcpp(MeanVec, SD);
          count++;
          if (FamilyType == 1) { //Probit
            PredMat = arma::max(PredMat, YMax);
            ProbitOnes = find(PredMat > 0);
            PredMat(ProbitOnes) = arma::ones<arma::vec>(ProbitOnes.size());
          }
          if (FamilyType == 2) { //Tobit
            PredMat = arma::max(PredMat, YMax);
          }
        }
        
        //Categorical
        if (FamilyType == 3) {
          MeanVec = MeanMat.col(f);
          TrialsVec = TrialsMat.col(f - count);
          Pi = arma::exp(MeanVec) / (1 + arma::exp(MeanVec));
          for (arma::uword i = 0; i < M; i++) PredMat(i, f) = rbinomRcpp(TrialsVec(i), Pi(i));
        }
        
      //End loop over family type
      }
      
      //Save output for a prediction
      OutTemp.col(n) = arma::vectorise(PredMat);
      
    //End loop over predictions 
    }
    
    //Save output
    Out.slice(s) = OutTemp;

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
  return Out;

}
