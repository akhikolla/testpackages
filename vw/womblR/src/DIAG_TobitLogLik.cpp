#include <RcppArmadillo.h>
#include "DIAG_diagnostics.h"

//Function to compute log likelihood for tobit model--------------------------------------------------------------------------
arma::colvec TobitLogLik(datobjDIAG DatObj, paraDIAG Para, dataugDIAG DatAug, int NKeep) {

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

  //Set parameters
  arma::mat MuMat = Para.Mu;
  arma::mat Tau2Mat = Para.Tau2;
  arma::mat AlphaMat = Para.Alpha;

  //Verbose output
  arma::vec VerboseSeq;
  VerboseSeq << 0.25 << 0.50 << 0.75;
  VerboseSeq *= NKeep;

  //Set data augmentation objects
  int NBelow = DatAug.NBelow;
  arma::vec NBelowCount = DatAug.NBelowCount;
  arma::field<arma::uvec> NBelowBoolean = DatAug.NBelowBoolean;
  arma::field<arma::uvec> NAboveBoolean = DatAug.NAboveBoolean;
  arma::field<arma::vec> YStarNonZero = DatAug.YStarNonZero;

  //Initialize objects
  arma::colvec LogLik(NKeep);
  arma::colvec Mu(Nu), Tau2(Nu), Alpha(Nu);
  arma::cube WAlphas(M, M, Nu), JointCovariances(M, M, Nu);
  Rcpp::Rcout << std::fixed << "Calculating Log-Lik: 0%.. ";

  //Loop over scans
  for (int s = 0; s < NKeep; s++) {

    //Compute moments and get log-likelihood
    Mu = MuMat.row(s).t();
    Tau2 = Tau2Mat.row(s).t();
    Alpha = AlphaMat.row(s).t();
    WAlphas = WAlphaCube(Alpha, Z, W, M, Nu, WeightsInd);
    JointCovariances = JointCovarianceCube(WAlphas, Tau2, EyeM, Rho, M, Nu);

    //No zero truncation
    if (NBelow == 0) {

      //Initialize objects
      double loglik = 0;
      arma::mat Rooti(M, M);
      arma::colvec Mean(M), Y(M);

      //Loop over visits and return log-likelihood
      for (int i = 0; i < Nu; i++) {
        Y = YStarNonZero(i, 0);
        Mean = Mu(i) * OneM;
        Rooti = GetRooti(JointCovariances.slice(i), EyeM);
        loglik += lndMvn(Y, Mean, Rooti);
      }
      LogLik(s) = loglik;

    }

    //Any zero truncation
    else {

      //Get log-likelihood at each visit independently
      double logLik = 0;
      for (int i = 0; i < Nu; i++) {

        //Set parameters
        arma::colvec JointMean = Mu(i) * OneM;
        arma::mat JointCovariance = JointCovariances.slice(i);
        arma::uvec BelowBoolean = NBelowBoolean(i, 0);
        arma::uvec AboveBoolean = NAboveBoolean(i, 0);
        arma::colvec YStarVisit = YStarNonZero(i, 0);
        int NBelowVisit = NBelowCount(i);
        int NAboveVisit = M - NBelowVisit;

        //If particular time point has no zeros
        if (NBelowVisit == 0) logLik += lndMvn(YStarVisit, JointMean, GetRooti(JointCovariance, EyeM));

        //Otherwise proceed...
        else {

          //Initialize objects
          double logLik1, logLik2;
          arma::mat Cov11(NBelowVisit, NBelowVisit), Cov12(NBelowVisit, NAboveVisit);
          arma::mat Cov22(NAboveVisit, NAboveVisit), Cov22Inv(NAboveVisit, NAboveVisit);
          arma::mat CondCov(NBelowVisit, NBelowVisit), CondCovRound(NBelowVisit, NBelowVisit);
          arma::mat EyeNAboveVisit(NAboveVisit, NAboveVisit, arma::fill::eye) ;
          arma::mat CovStar(NBelowVisit, NAboveVisit);
          arma::colvec Mean1(NBelowVisit), Mean2(NAboveVisit), CondMean(NBelowVisit);
          arma::uvec AnyFalse;

          //Subset to block moments
          Cov11 = JointCovariance(BelowBoolean, BelowBoolean);
          Cov12 = JointCovariance(BelowBoolean, AboveBoolean);
          Cov22 = JointCovariance(AboveBoolean, AboveBoolean);
          Mean1 = JointMean(BelowBoolean);
          Mean2 = JointMean(AboveBoolean);
          Cov22Inv = CholInv(Cov22);

          //Compute conditional moments
          CovStar = Cov12 * Cov22Inv;
          CondCov = Cov11 - CovStar * arma::trans(Cov12);
          CondMean = Mean1 + CovStar * (YStarVisit - Mean2);

          //Resolve numerical issue
          int digit = 32;
          SEXP roundSEXP = round(Rcpp::Named("x", CondCov), Rcpp::Named("digits", digit));
          CondCovRound = Rcpp::as<arma::mat>(roundSEXP);
          AnyFalse = find(arma::trans(CondCovRound) != CondCovRound);
          while (AnyFalse.n_elem > 0) {
            digit--;
            SEXP roundSEXP = round(Rcpp::Named("x", CondCov), Rcpp::Named("digits", digit));
            CondCovRound = Rcpp::as<arma::mat>(roundSEXP);
            AnyFalse = find(arma::trans(CondCovRound) != CondCovRound);
          }

          //Compute deviance
          logLik1 = lndMvn(YStarVisit, Mean2, GetRooti(Cov22, EyeNAboveVisit));
          logLik2 = log(pmvnormRcpp(NBelowVisit, CondMean, CondCovRound));
          logLik += (logLik1 + logLik2);

        //End otherwise proceed
        }

      //End loop over visits
      }

      //Save log-likelhood
      LogLik(s) = logLik;

    //End if statetment of zero observations
    }

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
double TobitLogLikMean(datobjDIAG DatObj, paraDIAG Para, dataugDIAG DatAug) {

  //Set base round function
  //Rcpp::Environment base("package:base"); //This is equivalent to library(base)
  Rcpp::Environment base = Rcpp::Environment::namespace_env("base"); //This is equivalent to PACKAGE::FUNCTION()
  Rcpp::Function round = base["round"];

  //Set data objects
  int Nu = DatObj.Nu;
  int M = DatObj.M;
  arma::mat EyeM = DatObj.EyeM;
  arma::vec OneM = DatObj.OneM;

  //Set parameters
  arma::vec MuMean = Para.MuMean;
  arma::cube CovMean = Para.CovMean;

  //Set data augmentation objects
  int NBelow = DatAug.NBelow;
  arma::vec NBelowCount = DatAug.NBelowCount;
  arma::field<arma::uvec> NBelowBoolean = DatAug.NBelowBoolean;
  arma::field<arma::uvec> NAboveBoolean = DatAug.NAboveBoolean;
  arma::field<arma::vec> YStarNonZero = DatAug.YStarNonZero;

  //Initialize output object
  double LogLik = 0;

  //No zero truncation
  if (NBelow == 0) {

    //Initialize objects
    double loglik = 0;
    arma::mat Rooti(M, M);
    arma::colvec Mean(M), Y(M);

    //Loop over visits and return log-likelihood
    for (int i = 0; i < Nu; i++) {
      Y = YStarNonZero(i, 0);
      Mean = MuMean(i) * OneM;
      Rooti = GetRooti(CovMean.slice(i), EyeM);
      loglik += lndMvn(Y, Mean, Rooti);
    }
    LogLik = loglik;

  }

  //Any zero truncation
  else {

    //Get log-likelihood at each visit independently
    double logLik = 0;
    for (int i = 0; i < Nu; i++) {

      //Set parameters
      arma::colvec JointMean = MuMean(i) * OneM;
      arma::mat JointCovariance = CovMean.slice(i);
      arma::uvec BelowBoolean = NBelowBoolean(i, 0);
      arma::uvec AboveBoolean = NAboveBoolean(i, 0);
      arma::colvec YStarVisit = YStarNonZero(i, 0);
      int NBelowVisit = NBelowCount(i);
      int NAboveVisit = M - NBelowVisit;

      //If particular time point has no zeros
      if (NBelowVisit == 0) logLik += lndMvn(YStarVisit, JointMean, GetRooti(JointCovariance, EyeM));

      //Otherwise proceed...
      else {

        //Initialize objects
        double logLik1, logLik2;
        arma::mat Cov11(NBelowVisit, NBelowVisit), Cov12(NBelowVisit, NAboveVisit);
        arma::mat Cov22(NAboveVisit, NAboveVisit), Cov22Inv(NAboveVisit, NAboveVisit);
        arma::mat CondCov(NBelowVisit, NBelowVisit), CondCovRound(NBelowVisit, NBelowVisit);
        arma::mat EyeNAboveVisit(NAboveVisit, NAboveVisit, arma::fill::eye) ;
        arma::mat CovStar(NBelowVisit, NAboveVisit);
        arma::colvec Mean1(NBelowVisit), Mean2(NAboveVisit), CondMean(NBelowVisit);
        arma::uvec AnyFalse;

        //Subset to block moments
        Cov11 = JointCovariance(BelowBoolean, BelowBoolean);
        Cov12 = JointCovariance(BelowBoolean, AboveBoolean);
        Cov22 = JointCovariance(AboveBoolean, AboveBoolean);
        Mean1 = JointMean(BelowBoolean);
        Mean2 = JointMean(AboveBoolean);
        Cov22Inv = CholInv(Cov22);

        //Compute conditional moments
        CovStar = Cov12 * Cov22Inv;
        CondCov = Cov11 - CovStar * arma::trans(Cov12);
        CondMean = Mean1 + CovStar * (YStarVisit - Mean2);

        //Resolve numerical issue
        int digit = 32;
        SEXP roundSEXP = round(Rcpp::Named("x", CondCov), Rcpp::Named("digits", digit));
        CondCovRound = Rcpp::as<arma::mat>(roundSEXP);
        AnyFalse = find(arma::trans(CondCovRound) != CondCovRound);
        while (AnyFalse.n_elem > 0) {
          digit--;
          SEXP roundSEXP = round(Rcpp::Named("x", CondCov), Rcpp::Named("digits", digit));
          CondCovRound = Rcpp::as<arma::mat>(roundSEXP);
          AnyFalse = find(arma::trans(CondCovRound) != CondCovRound);
        }

        //Compute deviance
        logLik1 = lndMvn(YStarVisit, Mean2, GetRooti(Cov22, EyeNAboveVisit));
        logLik2 = log(pmvnormRcpp(NBelowVisit, CondMean, CondCovRound));
        logLik += (logLik1 + logLik2);

      //End otherwise proceed
      }

    //End loop over visits
    }

    //Save log-likelhood
    LogLik = logLik;

  //End if statetment of zero observations
  }

  //Return log-likelihood
  return LogLik;

//End function
}
