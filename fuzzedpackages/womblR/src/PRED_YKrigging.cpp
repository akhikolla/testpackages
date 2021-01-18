#include <RcppArmadillo.h>
#include "PRED_predictions.h"

//Function that samples from the posterior predictive distribution for the spBDwDM model------------------
//This function is being exported to R for use in this package exclusively...
//not for use by users.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat YKrigging(Rcpp::List DatObj_List, arma::mat ThetaKrig, int NKeep) {

  //Convet Rcpp::Lists to C++ structs
  datobjPRED DatObj = ConvertDatObjPRED(DatObj_List);

  //Set data object
  int FamilyInd = DatObj.FamilyInd;
  int M = DatObj.M;
  int NNewVisits = DatObj.NNewVisits;
  int WeightsInd = DatObj.WeightsInd;
  double ScaleY = DatObj.ScaleY;
  double Rho = DatObj.Rho;
  arma::vec Z = DatObj.Z;
  arma::vec OneM = DatObj.OneM;
  arma::mat EyeM = DatObj.EyeM;
  arma::mat W = DatObj.W;

  //Verbose output
  arma::vec VerboseSeq;
  VerboseSeq << 0.25 << 0.50 << 0.75;
  VerboseSeq *= NKeep;
  Rcpp::Rcout << std::fixed << "Krigging Y: 0%.. ";

  //Initialize objects
  arma::mat ThetaT(NNewVisits, 3);
  arma::rowvec ThetaRow(3 * NNewVisits);
  arma::mat PPD(M * NNewVisits, NKeep), YTemp(M, NNewVisits), YMax(M, NNewVisits, arma::fill::zeros);
  arma::colvec Mu(NNewVisits), Tau2(NNewVisits), Alpha(NNewVisits);
  arma::cube WAlphas(M, M, NNewVisits), JointCovariances(M, M, NNewVisits);
  arma::umat ProbitOnes;

  //Loop over scans
  for (arma::uword s = 0; s < NKeep; s++) {

    //Extract level 1 parameters
    ThetaRow = ThetaKrig.row(s);
    ThetaT = arma::reshape(ThetaRow, 3, NNewVisits).t();
    Mu = ThetaT.col(0);
    Tau2 = square(exp(ThetaT.col(1)));
    Alpha = exp(ThetaT.col(2));

    //Get joint moments
    WAlphas = WAlphaCube(Alpha, Z, W, M, NNewVisits, WeightsInd);
    JointCovariances = JointCovarianceCube(WAlphas, Tau2, EyeM, Rho, M, NNewVisits);

    //Sample from posterior predictive distribution
    for (arma::uword i = 0; i < NNewVisits; i++) YTemp.col(i) = rmvnormRcpp(1, Mu(i) * OneM, JointCovariances.slice(i));

    //Adjust for probit or tobit censoring
    if (FamilyInd == 1) {
      YTemp = arma::max(YTemp, YMax);
      ProbitOnes = find(YTemp > 0);
      YTemp(ProbitOnes) = arma::ones<arma::vec>(ProbitOnes.size());
    }
    if (FamilyInd == 2) {
      YTemp = arma::max(YTemp, YMax);
    }

    //Scale the samples
    YTemp *= ScaleY;

    //Reshape and save ppd
    PPD.col(s) = arma::reshape(YTemp, M * NNewVisits, 1);

    //Add a new percentage
    Rcpp::Rcout.precision(0);
    if (std::find(VerboseSeq.begin(), VerboseSeq.end(), s) != VerboseSeq.end())
      Rcpp::Rcout << std::fixed << 100 * (s) / NKeep << "%.. ";

    //End loop
  }

  //Output final percentage
  Rcpp::Rcout << std::fixed << "100%.. Done!" << std::endl;

  //Return PPD
  return PPD;

}
