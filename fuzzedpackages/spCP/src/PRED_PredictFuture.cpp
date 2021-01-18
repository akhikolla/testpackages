#include <RcppArmadillo.h>
#include "PRED_predictions.h"
#include <math.h>

//Function that samples from the posterior distribution of a future time point--------------------------------------------------
//This function is being exported to R for use in this package exclusively...
//not for use by users.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube PredictFuture(Rcpp::List DatObj_List, Rcpp::List Para_List, int NKeep) {

  //Convet Rcpp::Lists to C++ structs
  datobjPRED DatObj = ConvertDatObjPRED(DatObj_List);
  paraPRED Para = ConvertParaPRED(Para_List);

  //Set data object
  int M = DatObj.M;
  int NNewTimes = DatObj.NNewTimes;
  double t1 = DatObj.t1;
  double tNu = DatObj.tNu;
  arma::vec NewTimes = DatObj.NewTimes;
  arma::vec OneM = DatObj.OneM;

  //Set parameters
  arma::mat Beta0Mat = arma::trans(Para.Beta0);
  arma::mat Beta1Mat = arma::trans(Para.Beta1);
  arma::mat Lambda0Mat = arma::trans(Para.Lambda0);
  arma::mat Lambda1Mat = arma::trans(Para.Lambda1);
  arma::mat EtaMat = arma::trans(Para.Eta);

  //Verbose output
  arma::vec VerboseSeq;
  VerboseSeq << 0.25 << 0.50 << 0.75;
  VerboseSeq *= NKeep;
  Rcpp::Rcout << std::fixed << "Future Predictions: 0%.. ";

  //Initialize objects
  arma::colvec Beta0(M), Beta1(M), Lambda0(M), Lambda1(M), Eta(M), Theta(M);
  arma::cube YPred(NKeep, M, NNewTimes);
  double NewTime, CP, Mu, Sigma;

  //Loop over scans
  for (arma::uword i = 0; i < NKeep; i++) {

    //Specify parameters at scan s
    Beta0 = Beta0Mat.col(i);
    Beta1 = Beta1Mat.col(i);
    Lambda0 = Lambda0Mat.col(i);
    Lambda1 = Lambda1Mat.col(i);
    Eta = EtaMat.col(i);
    Theta = arma::max(arma::min(tNu * OneM, Eta), t1 * OneM);

    //Loop over locations and new times
    for (arma::uword t = 0; t < NNewTimes; t++) {
      for (arma::uword s = 0; s < M; s++) {
        NewTime = NewTimes(t);
        CP = Theta(s);
        if (CP <= NewTime) Mu = Beta0(s);
        if (CP > NewTime) Mu = Beta0(s) + Beta1(s) * (NewTime - CP);
        if (CP <= NewTime) Sigma = exp(Lambda0(s));
        if (CP > NewTime) Sigma = exp(Lambda0(s) + Lambda1(s) * (NewTime - CP));
        YPred.subcube(i, s, t, i, s, t) = fmax(0, arma::as_scalar(rnormRcpp(1, Mu, Sigma)));
      }
    }

    //Add a new percentage
    Rcpp::Rcout.precision(0);
    if (std::find(VerboseSeq.begin(), VerboseSeq.end(), i) != VerboseSeq.end())
      Rcpp::Rcout << std::fixed << 100 * (i) / NKeep << "%.. ";

  //End loop
  }

  //Output final percentage
  Rcpp::Rcout << std::fixed << "100%.. Done!" << std::endl;

  //Return PPD
  return YPred;

}
