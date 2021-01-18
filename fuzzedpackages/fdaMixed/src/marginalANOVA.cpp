#include <RcppArmadillo.h>

// Arguments:
//   Gmat           : covariance matrix for random effect
//   Ymat           : matrix of observations
//   gammaMat       : marginal design matrix of fixed effect
//   Zmat           : marginal design matrix of random effect
//   condRes_       : matrix to contain conditional residuals
//   betaHat_       : matrix to contain estimate of fixed effect
//   uBLUP_         : matrix to contain prediction of random effect
//   uBLUPinvG_     : matrix to contain rotation of uBLUP
//   Cbeta_         : matrix to contain marginal variance of betaHat-beta
//   Cu_            : matrix to contain marginal variance of uBLUP-u
//
// Argument types, from which the dimensions (pp,qq,NN,MM) are extracted:
//   Gmat           : matrix of dimension (qq,qq)
//   Ymat           : matrix of dimension (NN,MM)
//   gammaMat       : matrix of dimension (MM,pp)
//   Zmat           : matrix of dimension (MM,qq)
//   condRes_       : matrix of dimension (NN,MM)
//   betaHat_       : matrix of dimension (NN,pp)
//   uBLUP_         : matrix of dimension (NN,qq)
//   uBLUPinvG_     : matrix of dimension (NN,qq)
//   Cbeta_         : matrix of dimension (pp,pp)
//   Cu_            : matrix of dimension (qq,qq)
//
// Value:
//   1) estimates of fixed effect retured in variable betaHat
//   2) prediction of random effect retured in variable uBLUP
//   3) conditional residuals returned in variable condRes
//   4) rotation of uBLUP returned in variable uBLUPinvG
//   5) variance matrices returned in variables Cbeta and Cu

extern "C" SEXP marginalANOVA(SEXP Gmat,
                              SEXP Ymat,
                              SEXP gammaMat, SEXP Zmat,
                              SEXP condRes_,
                              SEXP betaHat_, SEXP uBLUP_, SEXP uBLUPinvG_,
                              SEXP Cbeta_,   SEXP Cu_
                              ) {
  try{
    // creates Rcpp objects from SEXP's
    Rcpp::NumericMatrix Gr(Gmat);
    Rcpp::NumericMatrix Yr(Ymat);
    Rcpp::NumericMatrix gammaR(gammaMat);
    Rcpp::NumericMatrix Zr(Zmat);
    Rcpp::NumericMatrix condResR(condRes_);
    Rcpp::NumericMatrix betaR(betaHat_);
    Rcpp::NumericMatrix uR(uBLUP_);
    Rcpp::NumericMatrix uGr(uBLUPinvG_);
    Rcpp::NumericMatrix CbetaR(Cbeta_);
    Rcpp::NumericMatrix CuR(Cu_);

    // extract dimensions from defined objects
    int NN = Yr.nrow();
    int MM = Yr.ncol();
    int pp = betaR.ncol();
    int qq = uR.ncol();

    // introduce iterators to handle situation with pp=0 and/or qq=0
    arma::mat dummy_mat = arma::zeros<arma::mat>(1,1);
    arma::mat::iterator Gp     = dummy_mat.begin();
    arma::mat::iterator gammaP = dummy_mat.begin();
    arma::mat::iterator Zp     = dummy_mat.begin();
    arma::mat::iterator betaP  = dummy_mat.begin();
    arma::mat::iterator uP     = dummy_mat.begin();
    arma::mat::iterator uGp    = dummy_mat.begin();
    arma::mat::iterator CbetaP = dummy_mat.begin();
    arma::mat::iterator CuP    = dummy_mat.begin();
    if (pp > 0) {
      gammaP = gammaR.begin();
      betaP  = betaR.begin();
      CbetaP = CbetaR.begin();
    }
    if (qq > 0) {
      Gp  = Gr.begin();
      Zp  = Zr.begin();
      uP  = uR.begin();
      uGp = uGr.begin();
      CuP = CuR.begin();
    }

    // reuses memory and avoids extra copies  
    arma::mat Y(Yr.begin(),NN,MM,false);
    arma::mat condRes(condResR.begin(),NN,MM,false);
    arma::mat gamma(gammaP,MM,pp,false);
    arma::mat betaHat(betaP,NN,pp,false);
    arma::mat Cbeta(CbetaP,pp,pp,false);
    arma::mat G(Gp,qq,qq,false);
    arma::mat Z(Zp,MM,qq,false);
    arma::mat uBLUP(uP,NN,qq,false);
    arma::mat uBLUPinvG(uGp,NN,qq,false);
    arma::mat Cu(CuP,qq,qq,false);

    // Make covariance and projection matrices
    arma::mat Ginv     = arma::zeros<arma::mat>(qq,qq);
    arma::mat uProj    = arma::zeros<arma::mat>(MM,qq);
    arma::mat betaProj = arma::zeros<arma::mat>(MM,pp);
    if (qq!=0) {
      Ginv  = arma::inv(G);
      Cu    = arma::inv(Ginv+arma::trans(Z)*Z);
      uProj = Z*Cu;
      if (pp!=0) {
	Cbeta = arma::inv(arma::trans(gamma)*gamma-arma::trans(gamma)*uProj*arma::trans(Z)*gamma);
	betaProj = (gamma-uProj*arma::trans(Z)*gamma)*Cbeta;
      }
    } else {
      if (pp!=0) {
	Cbeta    = arma::inv(arma::trans(gamma)*gamma);
	betaProj = gamma*Cbeta;
      }
    }

    // Do projections
    condRes = Y;

    if (pp!=0) {
      for (int nn=0; nn<NN; nn++) {
	betaHat.row(nn)  = arma::reshape(condRes.row(nn)*betaProj,1,pp);
    	condRes.row(nn) -= betaHat.row(nn)*arma::trans(gamma);
      }
    }

    if (qq!=0) {
      for (int nn=0; nn<NN; nn++) {
	uBLUP.row(nn)      = condRes.row(nn)*uProj;
	condRes.row(nn)   -= uBLUP.row(nn)*arma::trans(Z);
	uBLUPinvG.row(nn)  = uBLUP.row(nn)*Ginv;
      }
    }

    // Result is saved and returned in the external variables:
    //   betaHat, Cbeta, uBLUP, uBLUPinvG, Cu, condRes

    // Below possible exceptions are taken care off
  } catch( std::exception &ex ) {
      forward_exception_to_r( ex );
  } catch(...) { 
      ::Rf_error( "c++ exception (unknown reason)"); 
  }

  // return to R
  return R_NilValue;
}
