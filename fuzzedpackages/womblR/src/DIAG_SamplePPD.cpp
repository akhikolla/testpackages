#include <RcppArmadillo.h>
#include "DIAG_diagnostics.h"

//Function that samples from the posterior predictive distribution for the spBDwDM model------------------
//This function is being exported to R for use in this package exclusively...
//not for use by users.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat SamplePPD(Rcpp::List DatObj_List, Rcpp::List Para_List, int NKeep) {

  //Convet Rcpp::Lists to C++ structs
  datobjDIAG DatObj = ConvertDatObjDIAG(DatObj_List);
  paraDIAG Para = ConvertParaDIAG(Para_List);

  //Set data object
  int FamilyInd = DatObj.FamilyInd;
  int Nu = DatObj.Nu;
  int M = DatObj.M;
  int WeightsInd = DatObj.WeightsInd;
  double ScaleY = DatObj.ScaleY;
  double Rho = DatObj.Rho;
  arma::vec Z = DatObj.Z;
  arma::vec OneM = DatObj.OneM;
  arma::mat EyeM = DatObj.EyeM;
  arma::mat W = DatObj.W;

  //Set parameters
  arma::mat MuMat = Para.Mu.t();
  arma::mat AlphaMat = Para.Alpha.t();
  arma::mat Tau2Mat = Para.Tau2.t();

  //Verbose output
  arma::vec VerboseSeq;
  VerboseSeq << 0.25 << 0.50 << 0.75;
  VerboseSeq *= NKeep;
  Rcpp::Rcout << std::fixed << "Calculating PPD: 0%.. ";

  //Initialize objects
  arma::mat PPD(M * Nu, NKeep), YTemp(M, Nu), YMax(M, Nu, arma::fill::zeros);
  arma::colvec Mu(Nu), Tau2(Nu), Alpha(Nu);
  arma::cube WAlphas(M, M, Nu), JointCovariances(M, M, Nu);
  arma::umat ProbitOnes;

  //Loop over scans
  for (int s = 0; s < NKeep; s++) {

    //Get joint moments
    Mu = MuMat.col(s);
    Tau2 = Tau2Mat.col(s);
    Alpha = AlphaMat.col(s);
    WAlphas = WAlphaCube(Alpha, Z, W, M, Nu, WeightsInd);
    JointCovariances = JointCovarianceCube(WAlphas, Tau2, EyeM, Rho, M, Nu);

    //Sample from posterior predictive distribution
    for (int i = 0; i < Nu; i++) YTemp.col(i) = rmvnormRcpp(1, Mu(i) * OneM, JointCovariances.slice(i));

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
    PPD.col(s) = arma::reshape(YTemp, M * Nu, 1);

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
