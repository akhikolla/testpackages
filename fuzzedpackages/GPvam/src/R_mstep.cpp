#include <RcppArmadillo.h>
RcppExport SEXP R_mstep_cpp ( SEXP invsqrtW_, SEXP JYp_, SEXP loopsize_,
SEXP patternlength_, SEXP rownumber_,SEXP ybetas_,SEXP etahat_, SEXP tempmatR_,
SEXP JXpi_, SEXP JXpp_, SEXP JXpx_, SEXP JXpdim_,
SEXP JZpi_, SEXP JZpp_,SEXP JZpx_,SEXP JZpdim_){
arma::mat invsqrtW = Rcpp::as<arma::mat>(invsqrtW_);
arma::mat etahat = Rcpp::as<arma::mat>(etahat_);
arma::mat tempmatR = Rcpp::as<arma::mat>(tempmatR_);
arma::mat JY = Rcpp::as<arma::mat>(JYp_);
arma::umat JXpi = Rcpp::as<arma::umat>(JXpi_);
arma::umat JXpp = Rcpp::as<arma::umat>(JXpp_);
arma::mat JXpx = Rcpp::as<arma::mat>(JXpx_);
arma::mat JXpdim = Rcpp::as<arma::mat>(JXpdim_);
arma::sp_mat JX(JXpi,JXpp,JXpx,JXpdim(0,0),JXpdim(1,0));
arma::umat JZpi = Rcpp::as<arma::umat>(JZpi_);
arma::umat JZpp = Rcpp::as<arma::umat>(JZpp_);
arma::mat JZpx = Rcpp::as<arma::mat>(JZpx_);
arma::mat JZpdim = Rcpp::as<arma::mat>(JZpdim_);
arma::sp_mat JZ(JZpi,JZpp,JZpx,JZpdim(0,0),JZpdim(1,0));
int loopsize =Rcpp::as<int>(loopsize_);
int patternlength =Rcpp::as<int>(patternlength_);
arma::umat rownumber = Rcpp::as<arma::umat>(rownumber_);
arma::mat ybetas = Rcpp::as<arma::mat>(ybetas_);
arma::mat JXt, JZt, JYt,tempt,patternsum,W;
arma::umat rownumbert;
patternsum.zeros(patternlength,patternlength);
for (int i=0; i<loopsize; i++) {
 JXt=JX.rows( i * patternlength, (i+1) * patternlength-1 );
 JYt=JY.rows( i * patternlength, (i+1) * patternlength-1 );
 JZt=JZ.rows( i * patternlength, (i+1) * patternlength-1 );
 rownumbert=rownumber.rows( i * patternlength, (i+1) * patternlength-1 );
 W.zeros(patternlength,patternlength);
 W.diag()= invsqrtW.submat(rownumbert(0,0)-1,0,rownumbert(patternlength-1,0)-1,0);
 tempt=JYt-JXt*ybetas;
 patternsum+=W*(tempt*tempt.t()-tempt*(JZt*etahat).t()-JZt*etahat*tempt.t()+JZt*tempmatR*JZt.t())*W;
// Rcpp::Rcout << patternsum << std::endl;
  };
 return(Rcpp::wrap(patternsum));
 }


