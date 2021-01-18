#include <RcppArmadillo.h>
#include "MCMC_STBDwDM.h"

//Initiate burn-in progress bar--------------------------------------------------------------------------
void BeginBurnInProgress(mcmcobj McmcObj, bool Interactive) {

  //Set MCMC object
  int BarLength = McmcObj.BarLength;

  //Initialize burn-in bar
  if (Interactive) {
    Rcpp::Rcout << std::fixed << "Burn-in progress:  |";
    for (int i = 0; i < BarLength - 1; i++) Rcpp::Rcout << std::fixed << " ";
    Rcpp::Rcout << std::fixed <<  "|" << std::fixed;
  }
  if (!Interactive) {
    Rcpp::Rcout << std::fixed << "Burn-in progress:  0%..  ";
  }

}



//Function to pilot adapt tuning parameter--------------------------------------------------------------
double PilotAdaptFunc(double TuningParameter, double AcceptancePct) {

  //Adjust tuning parameter using scaling based on size of acceptance rate
  if (AcceptancePct >= 0.90) TuningParameter *= 1.3;
  if ( (AcceptancePct >= 0.75 ) & (AcceptancePct < 0.90 ) ) TuningParameter *= 1.2;
  if ( (AcceptancePct >= 0.45 ) & (AcceptancePct < 0.75 ) ) TuningParameter *= 1.1;
  if ( (AcceptancePct <= 0.25 ) & (AcceptancePct > 0.15 ) ) TuningParameter *= 0.9;
  if ( (AcceptancePct <= 0.15 ) & (AcceptancePct > 0.10 ) ) TuningParameter *= 0.8;
  if (AcceptancePct <= 0.10) TuningParameter *= 0.7;
  return TuningParameter;

}



//Function for implementing pilot adaptation in MCMC sampler--------------------------------------------
metrobj PilotAdaptation(datobj DatObj, metrobj MetrObj, mcmcobj McmcObj) {

  //Set data objects
  int Nu = DatObj.Nu;

  //Set Metropolis objects
  arma::vec MetropTheta2 = MetrObj.MetropTheta2;
  arma::vec AcceptanceTheta2 = MetrObj.AcceptanceTheta2;
  arma::vec MetropTheta3 = MetrObj.MetropTheta3;
  arma::vec AcceptanceTheta3 = MetrObj.AcceptanceTheta3;
  double MetropPhi = MetrObj.MetropPhi;
  double AcceptancePhi = MetrObj.AcceptancePhi;

  //Set MCMC objects
  int PilotAdaptDenominator = McmcObj.PilotAdaptDenominator;

  //Get acceptance percentages
  arma::vec PctTheta2 = AcceptanceTheta2 / double(PilotAdaptDenominator);
  arma::vec PctTheta3 = AcceptanceTheta3 / double(PilotAdaptDenominator);
  double PctPhi = AcceptancePhi / double(PilotAdaptDenominator);

  //Update Tuning Parameter
  for (int i = 0; i < Nu; i++) MetropTheta2(i) = PilotAdaptFunc(MetropTheta2(i), PctTheta2(i));
  for (int i = 0; i < Nu; i++) MetropTheta3(i) = PilotAdaptFunc(MetropTheta3(i), PctTheta3(i));
  MetropPhi = PilotAdaptFunc(MetropPhi, PctPhi);
  MetrObj.MetropTheta2 = MetropTheta2;
  MetrObj.MetropTheta3 = MetropTheta3;
  MetrObj.MetropPhi = MetropPhi;

  //Zero the acceptance counters
  AcceptanceTheta2.zeros();
  AcceptanceTheta3.zeros();
  AcceptancePhi = 0;
  MetrObj.AcceptanceTheta2 = AcceptanceTheta2;
  MetrObj.AcceptanceTheta3 = AcceptanceTheta3;
  MetrObj.AcceptancePhi = AcceptancePhi;
  return MetrObj;

}



//Output Metropolis object for summary-------------------------------------------------------------------
Rcpp::List OutputMetrObj(metrobj MetrObj) {

  Rcpp::List Out = Rcpp::List::create(Rcpp::Named("AcceptanceTheta2") = MetrObj.AcceptanceTheta2,
                                      Rcpp::Named("MetropTheta2") = MetrObj.MetropTheta2,
                                      Rcpp::Named("AcceptanceTheta3") = MetrObj.AcceptanceTheta3,
                                      Rcpp::Named("MetropTheta3") = MetrObj.MetropTheta3,
                                      Rcpp::Named("AcceptancePhi") = MetrObj.AcceptancePhi,
                                      Rcpp::Named("MetropPhi") = MetrObj.MetropPhi);
  return Out;

}



//Initiate burn-in progress bar-------------------------------------------------------------------------------------
void SamplerProgress(int s, mcmcobj McmcObj) {

  //Set MCMC object
  int NSims = McmcObj.NSims;
  int NBurn = McmcObj.NBurn;

  //Add a new percentage
  Rcpp::Rcout.precision(0);
  Rcpp::Rcout << std::fixed << 100 * (s - NBurn) / NSims << "%..  ";

}



//Function for storing raw MCMC samples to to an object in memory
arma::colvec StoreSamples(datobj DatObj, para Para) {

  //Set data object
  int Nu = DatObj.Nu;

  //Set parameter objects
  arma::vec Mu = Para.Mu;
  arma::vec Tau2 = Para.Tau2;
  arma::vec Alpha = Para.Alpha;
  arma::vec Delta = Para.Delta;
  arma::mat T = Para.T;
  double Phi = Para.Phi;

  //Save raw samples
  int counter = 0;
  arma::colvec col(3 * Nu + 3 + 6 + 1);
  for (int i = 0; i < Nu; i++) col(i) = Mu(i);
  for (int i = 0; i < Nu; i++) col(i + Nu) = Tau2(i);
  for (int i = 0; i < Nu; i++) col(i + 2 * Nu) = Alpha(i);
  for (int i = 0; i < 3; i++) col(i +  3 * Nu) = Delta(i);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j <= i; j++) {
      col(3 * Nu + 3 + counter) = T(i, j);
      counter++;
    }
  }
  col(3 * Nu + 3 + 6) = Phi;
  return col;
}



// //Update burn-in progress bar----------------------------------------------------------------------------
// void UpdateBurnInBar(int s, mcmcobj McmcObj) {
//
//   //Set MCMC object
//   arma::vec WhichBurnInProgress = McmcObj.WhichBurnInProgress;
//   int BarLength = McmcObj.BarLength;
//
//   //Add a new star
//   arma::uvec NewStarBoolean = find(s == WhichBurnInProgress);
//   arma::vec NewStarBooleanVec = arma::conv_to<arma::vec>::from(NewStarBoolean);
//   int NewStar = NewStarBooleanVec(0);
//   for (int i = 0; i < (BarLength + 1 - NewStar); i++) Rcpp::Rcout << std::fixed << "\b";
//   Rcpp::Rcout << std::fixed << "*";
//   for (int i = 0; i < (BarLength - 1 - NewStar); i++) Rcpp::Rcout << std::fixed << " ";
//   Rcpp::Rcout << std::fixed << "|";
//
// }




//Update burn-in progress bar----------------------------------------------------------------------------
void UpdateBurnInBarInt(int s, mcmcobj McmcObj) {

  //Set MCMC object
  arma::vec WhichBurnInProgressInt = McmcObj.WhichBurnInProgressInt;
  arma::uvec NewStarBoolean = find(s == WhichBurnInProgressInt);
  arma::vec NewStarBooleanVec = arma::conv_to<arma::vec>::from(NewStarBoolean);
  int NewStar = NewStarBooleanVec(0);

  //Add percentage to submited job mode
  Rcpp::Rcout.precision(0);
  if (NewStar == 0) Rcpp::Rcout << std::fixed << "10%..  ";
  if (NewStar == 1) Rcpp::Rcout << std::fixed << "20%..  ";
  if (NewStar == 2) Rcpp::Rcout << std::fixed << "30%..  ";
  if (NewStar == 3) Rcpp::Rcout << std::fixed << "40%..  ";
  if (NewStar == 4) Rcpp::Rcout << std::fixed << "50%..  ";
  if (NewStar == 5) Rcpp::Rcout << std::fixed << "60%..  ";
  if (NewStar == 6) Rcpp::Rcout << std::fixed << "70%..  ";
  if (NewStar == 7) Rcpp::Rcout << std::fixed << "80%..  ";
  if (NewStar == 8) Rcpp::Rcout << std::fixed << "90%..  ";
  if (NewStar == 9) Rcpp::Rcout << std::fixed << "100%!  ";

}



// //Update burn-in progress bar----------------------------------------------------------------------------
// void UpdateBurnInBar(int s, mcmcobj McmcObj) {
//
//   //Set MCMC object
//   arma::vec WhichBurnInProgress = McmcObj.WhichBurnInProgress;
//   int BarLength = McmcObj.BarLength;
//
//   //Number of new star in interactive mode
//   arma::uvec NewStarBoolean = find(s == WhichBurnInProgress);
//   arma::vec NewStarBooleanVec = arma::conv_to<arma::vec>::from(NewStarBoolean);
//   int NewStar = NewStarBooleanVec(0);
//   for (int i = 0; i < (BarLength + 1 - NewStar); i++) Rcpp::Rcout << std::fixed << "\b";
//   Rcpp::Rcout << std::fixed << "*";
//   for (int i = 0; i < (BarLength - 1 - NewStar); i++) Rcpp::Rcout << std::fixed << " ";
//   Rcpp::Rcout << std::fixed << "|";
//
// }



//Update burn-in progress bar----------------------------------------------------------------------------
void UpdateBurnInBar(int s, mcmcobj McmcObj) {

  //Set MCMC object
  arma::vec WhichBurnInProgress = McmcObj.WhichBurnInProgress;
  int BarLength = McmcObj.BarLength;

  //Add a new star
  arma::uvec NewStarBoolean = find(s == WhichBurnInProgress);
  arma::vec NewStarBooleanVec = arma::conv_to<arma::vec>::from(NewStarBoolean);
  int NewStar = NewStarBooleanVec(0);
  Rcpp::Rcout << std::fixed << "\rBurn-in progress:  |";
  for (int i = 0; i < NewStar; i++) Rcpp::Rcout << std::fixed << "*";
  for (int i = 0; i < (BarLength - 1 - NewStar); i++) Rcpp::Rcout << std::fixed << " ";
  Rcpp::Rcout << std::fixed << "|";

}

