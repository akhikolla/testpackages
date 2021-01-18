#include <RcppArmadillo.h>
#include "DIAG_diagnostics.h"

//Function to compute log likelihood for tobit model--------------------------------------------------------------------------
arma::colvec NormalLogLik(datobjDIAG DatObj, paraDIAG Para, int NKeep) {

  //Set base round function
  //Rcpp::Environment base("package:base"); //This is equivalent to library(base)
  Rcpp::Environment base = Rcpp::Environment::namespace_env("base"); //This is equivalent to PACKAGE::FUNCTION()
  Rcpp::Function round = base["round"];

  //Set data objects
  int Nu = DatObj.Nu;
  int M = DatObj.M;
  int WeightsInd = DatObj.WeightsInd;
  arma::mat EyeM = DatObj.EyeM;
  arma::vec OneM = DatObj.OneM;
  arma::vec Z = DatObj.Z;
  arma::mat W = DatObj.W;
  double Rho = DatObj.Rho;
  arma::mat YStarWide = DatObj.YStarWide;

  //Set parameters
  arma::mat MuMat = Para.Mu;
  arma::mat Tau2Mat = Para.Tau2;
  arma::mat AlphaMat = Para.Alpha;

  //Verbose output
  arma::vec VerboseSeq;
  VerboseSeq << 0.25 << 0.50 << 0.75;
  VerboseSeq *= NKeep;

  //Initialize objects
  arma::colvec LogLik(NKeep);
  arma::colvec Mu(Nu), Tau2(Nu), Alpha(Nu);
  arma::cube WAlphas(M, M, Nu), JointCovariances(M, M, Nu);
  Rcpp::Rcout << std::fixed << "Calculating Log-Lik: 0%.. ";

  //Initialize objects
  double loglik;
  arma::mat Rooti(M, M);
  arma::colvec Mean(M), Y(M);

  //Loop over scans
  for (int s = 0; s < NKeep; s++) {

    //Compute moments and get log-likelihood
    Mu = MuMat.row(s).t();
    Tau2 = Tau2Mat.row(s).t();
    Alpha = AlphaMat.row(s).t();
    WAlphas = WAlphaCube(Alpha, Z, W, M, Nu, WeightsInd);
    JointCovariances = JointCovarianceCube(WAlphas, Tau2, EyeM, Rho, M, Nu);

    //Loop over visits and return log-likelihood
    loglik = 0;
    for (int i = 0; i < Nu; i++) {
      Y = YStarWide.col(i);
      Mean = Mu(i) * OneM;
      Rooti = GetRooti(JointCovariances.slice(i), EyeM);
      loglik += lndMvn(Y, Mean, Rooti);
    }
    LogLik(s) = loglik;

    //Add a new percentage
    Rcpp::Rcout.precision(0);
    if (std::find(VerboseSeq.begin(), VerboseSeq.end(), s) != VerboseSeq.end())
      Rcpp::Rcout << std::fixed << 100 * (s) / NKeep << "%.. ";

    //End loop over scans
  }

  //Output final percentage
  Rcpp::Rcout << std::fixed << "100%.. Done!" << std::endl;

  //Return log-likelihood
  return LogLik;

  //End function
}



//Function to compute log likelihood for tobit model--------------------------------------------------------------------------
double NormalLogLikMean(datobjDIAG DatObj, paraDIAG Para) {

  //Set base round function
  //Rcpp::Environment base("package:base"); //This is equivalent to library(base)
  Rcpp::Environment base = Rcpp::Environment::namespace_env("base"); //This is equivalent to PACKAGE::FUNCTION()
  Rcpp::Function round = base["round"];

  //Set data objects
  int Nu = DatObj.Nu;
  int M = DatObj.M;
  arma::mat EyeM = DatObj.EyeM;
  arma::vec OneM = DatObj.OneM;
  arma::mat YStarWide = DatObj.YStarWide;

  //Set parameters
  arma::vec MuMean = Para.MuMean;
  arma::cube CovMean = Para.CovMean;

  //Initialize objects
  double loglik = 0;
  arma::mat Rooti(M, M);
  arma::colvec Mean(M), Y(M);

  //Loop over visits and return log-likelihood
  for (int i = 0; i < Nu; i++) {
    Y = YStarWide.col(i);
    Mean = arma::as_scalar(MuMean(i)) * OneM;
    Rooti = GetRooti(CovMean.slice(i), EyeM);
    loglik += lndMvn(Y, Mean, Rooti);
  }

  //Return log-likelihood
  return loglik;

  //End function
}
