#include <RcppArmadillo.h>
#include "MCMC_bfa_sp.h"

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
metrobj PilotAdaptation(metrobj MetrObj, mcmcobj McmcObj, datobj DatObj) {

  //Set Data objects
  int SpCorInd = DatObj.SpCorInd;
  
  //Set Metropolis objects
  double MetropPsi = MetrObj.MetropPsi;
  double AcceptancePsi = MetrObj.AcceptancePsi;

  //Set MCMC objects
  int PilotAdaptDenominator = McmcObj.PilotAdaptDenominator;

  //Get acceptance percentages
  double PctPsi = AcceptancePsi / double(PilotAdaptDenominator);
  
  //Update Tuning Parameter
  MetropPsi = PilotAdaptFunc(MetropPsi, PctPsi);
  MetrObj.MetropPsi = MetropPsi;
  
  //Zero the acceptance counters
  AcceptancePsi = 0;
  MetrObj.AcceptancePsi = AcceptancePsi;
  
  //Update Rho
  if ((SpCorInd == 0) & (DatObj.IS == 1)) {
    double MetropRho = MetrObj.MetropRho;
    double AcceptanceRho = MetrObj.AcceptanceRho;
    double PctRho = AcceptanceRho / double(PilotAdaptDenominator);
    MetropRho = PilotAdaptFunc(MetropRho, PctRho);
    MetrObj.MetropRho = MetropRho;
    AcceptanceRho = 0;
    MetrObj.AcceptanceRho = AcceptanceRho;
  }  
  return MetrObj;

}



//Output Metropolis object for summary-------------------------------------------------------------------
Rcpp::List OutputMetrObj(metrobj MetrObj, datobj DatObj) {

  Rcpp::List Out;
  if (DatObj.SpCorInd == 1) Out = Rcpp::List::create(Rcpp::Named("AcceptancePsi") = MetrObj.AcceptancePsi,
                                                     Rcpp::Named("MetropPsi") = MetrObj.MetropPsi);
  if (DatObj.SpCorInd == 0) Out = Rcpp::List::create(Rcpp::Named("AcceptancePsi") = MetrObj.AcceptancePsi,
                                                     Rcpp::Named("MetropPsi") = MetrObj.MetropPsi,
                                                     Rcpp::Named("AcceptanceRho") = MetrObj.AcceptanceRho,
                                                     Rcpp::Named("MetropRho") = MetrObj.MetropRho);
  return Out;

}



//Initiate burn-in progress bar-------------------------------------------------------------------------------------
void SamplerProgress(int s, mcmcobj McmcObj) {

  //Set MCMC object
  int NSims = McmcObj.NSims;
  int NBurn = McmcObj.NBurn;

  //Add a new percentage
  Rcpp::Rcout.precision(0);
  if (s < NSims + NBurn)Rcpp::Rcout << std::fixed << 100 * (s - NBurn) / NSims << "%.. ";
  if (s == NSims + NBurn) Rcpp::Rcout << std::fixed << 100 * (s - NBurn) / NSims << "%!";

}



//Function for storing raw MCMC samples to to an object in memory
arma::colvec StoreSamples(datobj DatObj, para Para) {

  //Set data object
  int M = DatObj.M;
  int K = DatObj.K;
  int Nu = DatObj.Nu;
  int O = DatObj.O;
  int C = DatObj.C;
  int P = DatObj.P;

  //Set parameter objects
  arma::mat Lambda = Para.Lambda;
  arma::mat BigPhi = Para.BigPhi;
  arma::mat Sigma2 = Para.Sigma2;
  arma::mat Kappa = Para.Kappa;
  arma::colvec Delta = Para.Delta;
  arma::mat Upsilon = Para.Upsilon;
  double Psi = Para.Psi;
  arma::umat Xi = Para.Xi;
  double Rho = Para.Rho;
  arma::colvec Beta = Para.Beta;
  
  //Save raw samples
  int counter = 0;
  arma::uword Index;
  arma::colvec col(O * M * K + K * Nu + M * (O - C) + ((O + 1) * O) / 2 + K + ((K + 1) * K) / 2 + 1 + O * M * K + 1 + P);
  for (arma::uword o = 0; o < O; o++) {
    for (arma::uword i = 0; i < M; i++) {
      Index = i + M * o;
      for (arma::uword j = 0; j < K; j++) {
        col(counter) = Lambda(Index, j);
        counter++;
      }
    }
  }
  for (arma::uword t = 0; t < Nu; t++) {
    for (arma::uword j = 0; j < K; j++) {
      col(counter) = BigPhi(j, t);
      counter++;
    }
  }
  for (arma::uword c = 0; c < (O - C); c++) {
    for (arma::uword i = 0; i < M; i++) {
      col(counter) = Sigma2(i, c);
      counter++;
    }
  }
  for (arma::uword i = 0; i < O; i++) {
    for (arma::uword j = 0; j <= i; j++) {
      col(counter) = Kappa(i, j);
      counter++;
    }
  }
  for (arma::uword j = 0; j < K; j++) {
    col(counter) = Delta(j);
    counter++;
  }
  for (arma::uword i = 0; i < K; i++) {
    for (arma::uword j = 0; j <= i; j++) {
      col(counter) = Upsilon(i, j);
      counter++;
    }
  }
  col(counter) = Psi;
  counter++;
  for (arma::uword o = 0; o < O; o++) {
    for (arma::uword i = 0; i < M; i++) {
      Index = i + M * o;
      for (arma::uword j = 0; j < K; j++) {
        col(counter) = Xi(Index, j);
        counter++;
      }
    }
  }
  col(counter) = Rho;
  counter++;
  for (arma::uword p = 0; p < P; p++) {
    col(counter) = Beta(p);
    counter++;
  }
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
  if (NewStar == 0) Rcpp::Rcout << std::fixed << "10%.. ";
  if (NewStar == 1) Rcpp::Rcout << std::fixed << "20%.. ";
  if (NewStar == 2) Rcpp::Rcout << std::fixed << "30%.. ";
  if (NewStar == 3) Rcpp::Rcout << std::fixed << "40%.. ";
  if (NewStar == 4) Rcpp::Rcout << std::fixed << "50%.. ";
  if (NewStar == 5) Rcpp::Rcout << std::fixed << "60%.. ";
  if (NewStar == 6) Rcpp::Rcout << std::fixed << "70%.. ";
  if (NewStar == 7) Rcpp::Rcout << std::fixed << "80%.. ";
  if (NewStar == 8) Rcpp::Rcout << std::fixed << "90%.. ";
  if (NewStar == 9) Rcpp::Rcout << std::fixed << "100%!";

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

