#include <RcppArmadillo.h>
#include "PRED_predictions.h"

//Function that samples from the krigged distribution of theta--------------------------------------------------
//This function is being exported to R for use in this package exclusively...
//not for use by users.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat ThetaKrigging(Rcpp::List DatObj_List, Rcpp::List Para_List, int NKeep) {

  //Convet Rcpp::Lists to C++ structs
  datobjPRED DatObj = ConvertDatObjPRED(DatObj_List);
  paraPRED Para = ConvertParaPRED(Para_List);

  //Set data object
  int Nu = DatObj.Nu;
  int NNewVisits = DatObj.NNewVisits;
  int TempCorInd = DatObj.TempCorInd;
  arma::vec OneNu = DatObj.OneNu;
  arma::mat TimeDist = DatObj.TimeDist;
  arma::uvec NewVisits = DatObj.NewVisits;
  arma::uvec OriginalVisits = DatObj.OriginalVisits;

  //Set parameters objects
  arma::mat Theta1 = Para.Mu;
  arma::mat Tau2Mat = Para.Tau2;
  arma::mat Theta2 = log(sqrt(Tau2Mat));
  arma::mat AlphaMat = Para.Alpha;
  arma::mat Theta3 = log(AlphaMat);
  arma::mat DeltaMat= Para.Delta.t();
  arma::mat TMat = Para.T;
  arma::mat PhiMat = Para.Phi;

  //Verbose output
  arma::vec VerboseSeq;
  VerboseSeq << 0.25 << 0.50 << 0.75;
  VerboseSeq *= NKeep;
  Rcpp::Rcout << std::fixed << "Krigging Thetas: 0%.. ";

  //Initialize objects
  arma::colvec Delta(3), ThetaDiff(3 * Nu), Kron1(3 * NNewVisits);
  arma::mat T(3, 3), Theta(3, Nu);
  arma::mat SIGMAPhiFull(Nu + NNewVisits, Nu + NNewVisits);
  arma::mat CovNewOrig(NNewVisits, Nu), CovOrigOrig(Nu, Nu), CovPlus(NNewVisits, Nu);
  arma::mat CovNewNew(NNewVisits, NNewVisits), CovStar(NNewVisits, NNewVisits);
  arma::mat Sigma(NNewVisits * 3, NNewVisits * 3), Kron2(3 * NNewVisits, 3 * Nu);
  arma::colvec Mean(3 * NNewVisits), OneNNewVisits(NNewVisits, arma::fill::ones);
  arma::rowvec TMat2(6);
  double Phi;
  int counter;
  arma::mat KrigOut(3 * NNewVisits, NKeep);

  //Loop over scans
  for (arma::uword s = 0; s < NKeep; s++) {

    //Recover parameters (includes reshaping T)
    Theta.row(0) = Theta1.row(s);
    Theta.row(1) = Theta2.row(s);
    Theta.row(2) = Theta3.row(s);
    Delta = DeltaMat.col(s);
    TMat2 = TMat.row(s);
    counter = 0;
    for (arma::uword i = 0; i < 3; i++) {
      for (arma::uword j = 0; j <= i; j++) {
        T(i, j) = TMat2(counter);
        counter++;
      }
    }
    T = arma::symmatl(T);
    Phi = arma::as_scalar(PhiMat.row(s));

    //Calculate covariance components
    SIGMAPhiFull = SIGMA(Phi, TempCorInd, TimeDist, Nu + NNewVisits);
    CovNewOrig = SIGMAPhiFull(NewVisits, OriginalVisits);
    CovOrigOrig = SIGMAPhiFull(OriginalVisits, OriginalVisits);
    CovPlus = CovNewOrig * CholInv(CovOrigOrig);
    CovNewNew = SIGMAPhiFull(NewVisits, NewVisits);
    CovStar = CovNewNew - CovPlus * arma::trans(CovNewOrig);

    //Compute moments
    Sigma = arma::kron(CovStar, T);
    ThetaDiff = arma::vectorise(Theta) - arma::kron(OneNu, Delta);
    Kron1 = arma::kron(OneNNewVisits, Delta);
    Kron2 = arma::kron(CovPlus, T);
    Mean = Kron1 + Kron2 * ThetaDiff;

    //Sample from posterior predictive distribution
    KrigOut.col(s) = rmvnormRcpp(1, Mean, Sigma);

    //Add a new percentage
    Rcpp::Rcout.precision(0);
    if (std::find(VerboseSeq.begin(), VerboseSeq.end(), s) != VerboseSeq.end())
      Rcpp::Rcout << std::fixed << 100 * (s) / NKeep << "%.. ";

  //End loop
  }

  //Output final percentage
  Rcpp::Rcout << std::fixed << "100%.. Done!" << std::endl;

  //Return PPD
  return arma::trans(KrigOut);

}
