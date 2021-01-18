#include <RcppArmadillo.h>
#include "PRED_predictions.h"

//Function that samples from the krigged distribution of theta--------------------------------------------------
//This function is being exported to R for use in this package exclusively...
//not for use by users.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat EtaKrigging(Rcpp::List DatObj_List, Rcpp::List Para_List, int NKeep, bool Verbose) {

  //Convet Rcpp::Lists to C++ structs
  datobjPRED DatObj = ConvertDatObjPRED(DatObj_List);
  paraPRED Para = ConvertParaPRED(Para_List);

  //Set data object
  int Nu = DatObj.Nu;
  int K = DatObj.K;
  int NNewVisits = DatObj.NNewVisits;
  int TempCorInd = DatObj.TempCorInd;
  arma::mat TimeDist = DatObj.TimeDist;
  arma::uvec NewVisits = DatObj.NewVisits;
  arma::uvec OriginalVisits = DatObj.OriginalVisits;
  arma::mat EyeK = DatObj.EyeK;
  
  //Set parameters objects
  arma::mat PsiMat = Para.Psi;
  arma::mat UpsilonMat = Para.Upsilon;
  arma::mat EtaMat = Para.Eta;
  
  //Verbose output
  arma::vec VerboseSeq;
  if (Verbose) {
    VerboseSeq << 0.25 << 0.50 << 0.75;
    VerboseSeq *= NKeep;
    Rcpp::Rcout << std::fixed << "Krigging Eta: 0%.. ";
  }
  
  //Initialize objects
  arma::mat Upsilon(K, K), Sigma(K * NNewVisits, K * NNewVisits);
  arma::rowvec UpsilonRow(K);
  arma::colvec Eta(K * Nu), Mean(K * NNewVisits);
  arma::mat HPsiFull(Nu + NNewVisits, Nu + NNewVisits), KrigOut(K * NNewVisits, NKeep);
  arma::mat CovNewOrig(NNewVisits, Nu), CovOrigOrig(Nu, Nu), CovPlus(NNewVisits, Nu);
  arma::mat CovNewNew(NNewVisits, NNewVisits), CovStar(NNewVisits, NNewVisits);
  double Psi;
  int counter;

  //Loop over scans
  for (arma::uword s = 0; s < NKeep; s++) {

    //Recover parameters (includes reshaping Upsilon)
    Eta = EtaMat.row(s).t();
    arma::rowvec UpsilonRow = UpsilonMat.row(s);
    counter = 0;
    for (arma::uword i = 0; i < K; i++) {
      for (arma::uword j = 0; j <= i; j++) {
        Upsilon(i, j) = UpsilonRow(counter);
        counter++;
      }
    }
    Upsilon = arma::symmatl(Upsilon);
    Psi = arma::as_scalar(PsiMat.row(s));

    //Calculate covariance components
    HPsiFull = H(Psi, TempCorInd, TimeDist, Nu + NNewVisits);
    CovNewOrig = HPsiFull(NewVisits, OriginalVisits);
    CovOrigOrig = HPsiFull(OriginalVisits, OriginalVisits);
    CovPlus = CovNewOrig * CholInv(CovOrigOrig);
    CovNewNew = HPsiFull(NewVisits, NewVisits);
    CovStar = CovNewNew - CovPlus * arma::trans(CovNewOrig);

    //Compute moments
    Sigma = arma::kron(CovStar, Upsilon);
    Mean = arma::kron(CovPlus, EyeK) * Eta;

    //Sample from posterior predictive distribution
    KrigOut.col(s) = rmvnormRcpp(1, Mean, Sigma);

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
  return KrigOut;

}
