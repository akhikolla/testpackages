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

  //Set data objects
  int M = DatObj.M;
  int N = DatObj.N;
  double t1 = DatObj.t1;
  double tNu = DatObj.tNu;
  double ScaleY = DatObj.ScaleY;
  arma::uvec XThetaInd = DatObj.XThetaInd;
  arma::mat EyeM = DatObj.EyeM;
  arma::vec OneM = DatObj.OneM;
  arma::mat OneN = DatObj.OneN;
  arma::mat EyeN = DatObj.EyeN;
  arma::vec TimeVec = DatObj.TimeVec;
  arma::vec OneNu = DatObj.OneNu;
  arma::vec YObserved = DatObj.YObserved;
  int FamilyInd = DatObj.FamilyInd;

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
  Rcpp::Rcout << std::fixed << "Calculating PPD: 0%.. ";

  //Initialize objects
  arma::colvec YStar(N), Y(N), Theta(M), Mu(N), Sigma2(N);
  arma::rowvec Beta0(M), Beta1(M), Lambda0(M), Lambda1(M);
  arma::colvec Eta(M), Beta(2 * M), Lambda(2 * M), YMax(N, arma::fill::zeros);;
  arma::mat XTheta(N, 2 * M), BetaMat(2, M), LambdaMat(2, M), OmegaChol(N, N);
  arma::mat PPD(N, NKeep);
  arma::umat ProbitOnes;

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

    //Sample from posterior predictive distribution
    YStar = rmvnormRcpp(1, Mu, arma::diagmat(Sigma2));

    //Adjust for probit or tobit censoring
    if (FamilyInd == 1) {
      Y = arma::max(YStar, YMax);
      ProbitOnes = find(YStar > 0);
      Y(ProbitOnes) = arma::ones<arma::vec>(ProbitOnes.size());
    }
    if (FamilyInd == 2) {
      Y = arma::max(YStar, YMax);
    }

    //Scale the samples
    Y *= ScaleY;

    //Reshape and save ppd
    PPD.col(s) = Y;

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
