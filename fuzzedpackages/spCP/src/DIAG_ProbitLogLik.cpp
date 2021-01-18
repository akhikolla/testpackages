#include <RcppArmadillo.h>
#include "DIAG_diagnostics.h"

//Function to compute log likelihood for tobit model--------------------------------------------------------------------------
arma::colvec ProbitLogLik(datobjDIAG DatObj, paraDIAG Para, int NKeep) {

  //Set data objects
  int M = DatObj.M;
  int N = DatObj.N;
  double t1 = DatObj.t1;
  double tNu = DatObj.tNu;
  arma::uvec XThetaInd = DatObj.XThetaInd;
  arma::mat EyeM = DatObj.EyeM;
  arma::vec OneM = DatObj.OneM;
  arma::mat OneN = DatObj.OneN;
  arma::mat EyeN = DatObj.EyeN;
  arma::vec TimeVec = DatObj.TimeVec;
  arma::vec OneNu = DatObj.OneNu;
  arma::vec YObserved = DatObj.YObserved;

  //Set parameters
  arma::mat Beta0Mat = Para.Beta0;
  arma::mat Beta1Mat = Para.Beta1;
  arma::mat Lambda0Mat = Para.Lambda0;
  arma::mat Lambda1Mat = Para.Lambda1;
  arma::mat EtaMat = arma::trans(Para.Eta);

  //Verbose output
  arma::vec VerboseSeq;
  VerboseSeq << 0.25 << 0.50 << 0.75;
  VerboseSeq *= NKeep;

  //Initialize objects
  arma::colvec LogLik(NKeep), Theta(M), Mu(N), Sigma2(N);
  arma::rowvec Beta0(M), Beta1(M), Lambda0(M), Lambda1(M);
  arma::colvec Eta(M), Beta(2 * M), Lambda(2 * M);
  arma::mat XTheta(N, 2 * M), BetaMat(2, M), LambdaMat(2, M), OmegaChol(N, N);
  double loglik, Y;
  Rcpp::Rcout << std::fixed << "Calculating Log-Lik: 0%.. ";

  //Loop over scans
  for (int s = 0; s < NKeep; s++) {

    //Compute moments and get log-likelihood
    Beta0 = Beta0Mat.row(s);
    Beta1 = Beta1Mat.row(s);
    BetaMat.row(0) = Beta0;
    BetaMat.row(1) = Beta1;
    Beta = arma::vectorise(BetaMat);
    Lambda0 = Lambda0Mat.row(s);
    Lambda1 = Lambda1Mat.row(s);
    LambdaMat.row(0) = Lambda0;
    LambdaMat.row(1) = Lambda1;
    Lambda = arma::vectorise(LambdaMat);
    Eta = EtaMat.col(s);
    Theta = arma::max(arma::min(tNu * OneM, Eta), t1 * OneM);
    XTheta = GetXTheta(Theta, XThetaInd, TimeVec, OneNu, OneN, tNu, N, M);
    Mu = XTheta * Beta;
    Sigma2 = arma::exp(2 * (XTheta * Lambda));

    //Loop over observations
    loglik = 0;
    for (int i = 0; i < N; i++) {
      Y = YObserved(i);
      if (Y == 0) loglik += log(pnormRcpp(-Mu(i) / sqrt(Sigma2(i))));
      if (Y == 1) loglik += log(1 - pnormRcpp(-Mu(i) / sqrt(Sigma2(i))));
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



//Function to compute log likelihood for probit model--------------------------------------------------------------------------
double ProbitLogLikMean(datobjDIAG DatObj, paraDIAG Para) {

  //Set data objects
  arma::mat EyeN = DatObj.EyeN;
  arma::vec YObserved = DatObj.YObserved;
  int N = DatObj.N;

  //Set parameters
  arma::mat MuMean = Para.MuMean;
  arma::mat Sigma2Mean = Para.Sigma2Mean;

  //Initialize objects
  double LogLik = 0, Y;

  //Loop over observations
  LogLik = 0;
  for (int i = 0; i < N; i++) {
    Y = YObserved(i);
    if (Y == 1) LogLik += log(1 - pnormRcpp(-MuMean(i) / sqrt(Sigma2Mean(i))));
    if (Y == 0) LogLik += log(pnormRcpp(-MuMean(i) / sqrt(Sigma2Mean(i))));
  }

  //Return log-likelihood
  return LogLik;

  //End function
}
