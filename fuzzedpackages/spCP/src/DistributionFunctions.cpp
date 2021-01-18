#include <RcppArmadillo.h>

//Log density of a multivariate normal-----------------------------------------------------
double lndMvn(arma::vec const& Y, arma::vec const& Mu, arma::mat const& Rooti){
  arma::vec Z = arma::vectorise( arma::trans(Rooti) * (Y - Mu) );
  return ( -( Y.size() / 2.0 ) * log( 2 * M_PI ) - 0.5*(arma::trans(Z) * Z) + sum( log( arma::diagvec(Rooti) ) ) )[0];
}



//Log density of a normal distribution
double dlnorm(double x, double mu, double sigma2) {
  return -0.5 * log(2 * M_PI * sigma2) - 0.5 * ((x - mu) * (x - mu)) / sigma2;
}



//Function for sampling random standard uniform variable-----------------------------------
double randuRcpp() {
  return R::runif(0, 1);
}



//Function for sampling random chi square variable-----------------------------------------
double rchisqRcpp(double df) {
  return R::rchisq(df);
}



//Distribution function of a standard normal-----------------------------------------------
double pnormRcpp(double q) {
  return R::pnorm5(q, 0, 1, 1, 0);
}



//Inverse distribution function of a standard normal---------------------------------------
double qnormRcpp(double p) {
  return R::qnorm5(p, 0, 1, 1, 0);
}



//Function for sampling from a standard normal distribution--------------------------------
arma::vec rnormSNRcpp(int n) {
  arma::vec out(n);
  for (int i = 0; i < n; i++) out(i) = R::rnorm(0, 1);
  return out;
}



//Sample from a standard normal distribution-----------------------------------------------
arma::vec rnormRcpp(int n, double mean, double sd) {
  arma::vec muvec(1);
  muvec(0) = mean;
  return arma::repmat(muvec, n, 1) + rnormSNRcpp(n) * sd;
}



//Sample from a multivariate normal distribution-------------------------------------------
arma::mat rmvnormRcpp(int n, arma::vec const& mean, arma::mat const& sigma) {
  int ncols = sigma.n_cols;
  arma::vec Yvec(n * ncols);
  Yvec = rnormSNRcpp(n * ncols);
  arma::mat Y = arma::reshape(Yvec, n, ncols);
  return arma::trans(arma::repmat(mean, 1, n).t() + Y * arma::chol(sigma));
}



//Sample from a normal distribution truncated below or above by zero-----------------------
// This first version uses the inverse sampling from a CDF. It has problems
// when the data get to far from zero due to precision issues with qnorm and pnorm.
// See "Simulation of truncated normal variables" C. Robert (2009). As such I also
// include the truncated normal from the msm R package. This implementation uses an
// accept-reject method that is more efficient than the inverse sampling version
// that I currently use for values far from the truncation value (i.e, zero).
double rtnormRcpp(double mean, double sd, bool Above) {

  //Declare Variables
  double RandU = randuRcpp();
  double ScaledMean = -mean / sd;

  //Truncation Above by Zero
  if (Above) return sd * qnormRcpp(RandU * pnormRcpp(ScaledMean)) + mean;

  //Truncation Below by Zero
  else return sd * qnormRcpp(RandU - pnormRcpp(ScaledMean) * (RandU - 1)) + mean;

}
// arma::vec rtnormRcppMSM(int N, arma::vec const& mean, arma::vec const& sd, double lower, double upper) {
//
//   //Set truncated normal function
//   // Rcpp::Environment msm("package:msm"); //this is equivalent to library(msm)
//   Rcpp::Environment msm = Rcpp::Environment::namespace_env("msm"); //This is equivalent to PACKAGE::FUNCTION()
//   Rcpp::Function rtnorm = msm["rtnorm"];
//
//   //Evaluate pmvnorm
//   SEXP rtnormSEXP;
//   rtnormSEXP = rtnorm(Rcpp::Named("n", N),
//                       Rcpp::Named("mean", mean),
//                       Rcpp::Named("sd", sd),
//                       Rcpp::Named("lower", lower),
//                       Rcpp::Named("upper", upper));
//
//   //Convert output to double
//   return Rcpp::as<arma::vec>(rtnormSEXP);
//
// }
double rtnormRcppMSM(double mean, double sd, double lower, double upper) {

  //Set truncated normal function
  // Rcpp::Environment msm("package:msm"); //this is equivalent to library(msm)
  Rcpp::Environment msm = Rcpp::Environment::namespace_env("msm"); //This is equivalent to PACKAGE::FUNCTION()
  Rcpp::Function rtnorm = msm["rtnorm"];

  //Evaluate pmvnorm
  SEXP rtnormSEXP;
  rtnormSEXP = rtnorm(Rcpp::Named("n", 1),
                      Rcpp::Named("mean", mean),
                      Rcpp::Named("sd", sd),
                      Rcpp::Named("lower", lower),
                      Rcpp::Named("upper", upper));

  //Convert output to double
  return Rcpp::as<double>(rtnormSEXP);

}



//Sample from a Wishart distribution using the Bartlett decomposition----------------------
arma::mat rwishRcpp(double n, arma::mat const& V) {
  int p = V.n_rows;
  arma::mat L = arma::chol(V);
  arma::mat A(p, p, arma::fill::zeros);
  for (int i = 0; i < p; i++) A(i, i) = sqrt(rchisqRcpp(n - i));
  if (p > 1) {
    arma::vec RandSN = rnormSNRcpp(p * (p - 1) / 2);
    int counter = 0;
    for (int i = 0; i < p; i++) {
      for (int j = 0; j < i; j++){
        A(j, i) = RandSN(counter);
        counter++;
      }
    }
  }
  arma::mat AL = A * L;
  return arma::trans(AL) * AL;
}



//Sample from an Inverse-Wishart distribution----------------------------------------------
arma::mat riwishRcpp(double n, arma::mat const& V) {
  return arma::inv_sympd(rwishRcpp(n, arma::inv_sympd(V)));
}



//Multivariate normal CDF------------------------------------------------------------------
double pmvnormRcpp(int NBelowVisit, arma::vec const& CondMean, arma::mat const& CondCovRound) {

    //Set multivariate normal CDF function
    // Rcpp::Environment mvtnorm("package:mvtnorm"); //This is equivalent to library(mvtnorm)
    Rcpp::Environment mvtnorm = Rcpp::Environment::namespace_env("mvtnorm"); //This is equivalent to PACKAGE::FUNCTION()
    Rcpp::Function pmvnorm = mvtnorm["pmvnorm"];

    //Evaluate pmvnorm
    Rcpp::NumericVector Upper(NBelowVisit);
    Rcpp::NumericVector Mean = Rcpp::NumericVector(CondMean.begin(), CondMean.end());
    SEXP pmvnormSEXP = pmvnorm(Rcpp::Named("upper", Upper),
                               Rcpp::Named("mean", Mean),
                               Rcpp::Named("sigma", CondCovRound));

    //Convert output to double
    return Rcpp::as<double>(pmvnormSEXP);

}
