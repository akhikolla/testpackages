#include <RcppArmadillo.h>

// Arguments:
//   left_limit_    : left  limit of sampling interval (including border)
//   right_limit_   : right limit of sampling interval (including border)
//   alpha_2k_      : leading coefficient of L-operator
//   eta_left_real  : vector of real      part of left most eigenvalues
//   eta_left_imag  : vector of imaginary part of left most eigenvalues 
//   eta_right_real : vector of real      part of right most eigenvalues
//   eta_right_imag : vector of imaginary part of right most eigenvalues
//   FleftMat       : matrix of left  boundary conditions
//   FrightMat      : matrix of right boundary conditions
//   GreenVec       : vector to contain diagonal of Greens function
//   acfVec         : vector to contain the averaged auto covariance function
//
// Argument types, from which the dimensions (kk,order,pp,qq, NN,MM,lag) are extracted:
//   left_limit_    : double
//   right_limit_   : double
//   alpha_2k_      : double
//   eta_left_real  : vector of length kk
//   eta_left_imag  : vector of length kk
//   eta_right_real : vector of length kk
//   eta_right_imag : vector of length kk
//   FleftMat       : matrix of dimension (kk,2*kk)
//   FrightMat      : matrix of dimension (kk,2*kk)
//   GreenVec       : vector of length NN
//   acfVec         : vector of length lag
//
// Value:
//   1) 1-Diagonal of Green function returned in vector GreenVec
//   2) Average auto covariance function returned in vector acfVec
//
// Implicit assumptions:
//   1) left_limit_ < right_limit_
//   2) Sampling is equidistant as described in Markussen (2013)
//
// Remarks:
//   1) If eta_left_real <=0 and eta_right_real >= 0, then the algorithm 
//      is numerically stable.

extern "C" SEXP diagGreen(SEXP left_limit_, SEXP right_limit_, 
                          SEXP alpha_2k_, 
                          SEXP eta_left_real,  SEXP eta_left_imag,
                          SEXP eta_right_real, SEXP eta_right_imag,
                          SEXP FleftMat, SEXP FrightMat,
                          SEXP GreenVec, SEXP acfVec
                          ) {
  try{
    // copy data to armadillo structures
    double left_limit  = Rcpp::as<double>(left_limit_);
    double right_limit = Rcpp::as<double>(right_limit_);
    double alpha_2k    = Rcpp::as<double>(alpha_2k_);
    arma::cx_rowvec eta_left = 
      arma::cx_rowvec(Rcpp::as<arma::rowvec>(eta_left_real),
                      Rcpp::as<arma::rowvec>(eta_left_imag));
    arma::cx_rowvec eta_right = 
      arma::cx_rowvec(Rcpp::as<arma::rowvec>(eta_right_real),
                      Rcpp::as<arma::rowvec>(eta_right_imag));
    arma::mat Fleft  = Rcpp::as<arma::mat>(FleftMat);
    arma::mat Fright = Rcpp::as<arma::mat>(FrightMat);

    // creates Rcpp objects from SEXP's
    Rcpp::NumericVector GreenR(GreenVec);
    Rcpp::NumericVector acfR(acfVec);

    // extract dimensions from defined objects
    int kk = Fleft.n_rows;
    int NN = GreenR.size();
    int ll = acfR.size();

    // reuses memory and avoids extra copies
    arma::vec Green(GreenR.begin(),NN,false);
    arma::vec acf(acfR.begin(),ll,false);

    // define variables from Markussen (2013), Proposition 3
    double Delta = (right_limit-left_limit)/NN;
    arma::cx_mat Wleft = arma::ones<arma::cx_mat>(2*kk,kk);
    for (int ii=1; ii<2*kk; ii++) {
      Wleft.row(ii) = eta_left % Wleft.row(ii-1);   
    }
    arma::cx_mat Wright = arma::ones<arma::cx_mat>(2*kk,kk);
    for (int ii=1; ii<2*kk; ii++) {
      Wright.row(ii) = eta_right % Wright.row(ii-1);
    }
    arma::cx_mat W         = arma::join_rows(Wleft,Wright);
    arma::cx_mat invW      = arma::pinv(W);
    arma::rowvec v1        = arma::ones<arma::rowvec>(kk);
    arma::cx_vec vleft     = invW.submat(0, 2*kk-1,  kk-1,2*kk-1);
    arma::cx_vec vright    = invW.submat(kk,2*kk-1,2*kk-1,2*kk-1);
    
    // make some pre-computations
    arma::cx_mat Fleft_W    = arma::pinv(Fleft *Wleft) *Fleft *Wright;
    arma::cx_mat Fright_W   = arma::pinv(Fright*Wright)*Fright*Wleft;
    arma::cx_mat exp_FW_exp = arma::zeros<arma::cx_mat>(kk,kk);

    // compute discounters
    arma::cx_vec discounter_left       = arma::trans(exp(     Delta*eta_left));
    arma::cx_vec discounter_right      = arma::trans(exp(    -Delta*eta_right));
    arma::cx_vec prod_discounter_left  = arma::zeros<arma::cx_vec>(kk);
    arma::cx_vec prod_discounter_right = arma::zeros<arma::cx_vec>(kk);
    arma::cx_mat st_left               = arma::eye<arma::cx_mat>(kk,kk);
    arma::cx_mat st_right              = arma::eye<arma::cx_mat>(kk,kk);

    // backward loop over ll
    // last step ends with the diagonal as also required for output
    for (int ii=ll; ii>0; ii--) {
      // initialize product discounters
      prod_discounter_left  = arma::trans(exp( 0.5*Delta*eta_left));
      prod_discounter_right = arma::trans(exp(-0.5*Delta*eta_right));
      st_left               = arma::diagmat(exp( (ii-1)*Delta*eta_left));
      st_right              = arma::diagmat(exp(-(ii-1)*Delta*eta_right));

      // loop over s=nn*Delta
      for (int nn=0; nn<=NN-ii; nn++) {
        // compute part of phi: t=(nn+ii-1)*Delta
        exp_FW_exp = arma::diagmat(exp((nn+ii-0.5-NN)*Delta*eta_right))*Fright_W*
                     arma::diagmat(exp((NN-nn-ii+0.5)*Delta*eta_left));

        // update diagonal in Greens function: t=s+(ii-1)*Delta
        Green(nn) = arma::as_scalar(arma::real((v1-v1*exp_FW_exp)*
	  arma::pinv(arma::eye<arma::cx_mat>(kk,kk)-
                     st_left *arma::diagmat(prod_discounter_left)*Fleft_W*
                     st_right*arma::diagmat(prod_discounter_right)*exp_FW_exp)*
		     st_left*
  	            (vleft+arma::diagmat(prod_discounter_left)*Fleft_W*
		           arma::diagmat(prod_discounter_right)*vright)))/alpha_2k;

        // update discounters
        prod_discounter_left  = prod_discounter_left  % discounter_left;
        prod_discounter_right = prod_discounter_right % discounter_right;
      }

      // Substract from identity matrix
      Green *= -1;
      if (ii==1) Green += 1;

      // Mean covariance
      acf(ii-1) = arma::mean(Green.rows(0,NN-ii));
    }

    // Result is saved and returned in the external variables:
    //   GreenVec, acfVec

    // Below possible exceptions are taken care off
  } catch( std::exception &ex ) {
      forward_exception_to_r( ex );
  } catch(...) { 
      ::Rf_error( "c++ exception (unknown reason)"); 
  }

  // return to R
  return R_NilValue;
}
