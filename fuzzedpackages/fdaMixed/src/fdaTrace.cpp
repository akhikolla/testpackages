#include <RcppArmadillo.h>

// Arguments:
// left_limit_    : left  limit of sampling interval
// right_limit_   : right limit of sampling interval
// tau            : leading coefficient of L_lambda operator
// eta_left_real  : matrix of      real part of left most eigenvalues
// eta_left_imag  : matrix of imaginary part of left most eigenvalues
// eta_right_real : matrix of      real part of right most eigenvalues
// eta_right_imag : matrix of imaginary part of right most eigenvalues
// FleftMat       : matrix of left  boundary conditions
// FrightMat      : matrix of right boundary conditions
// N_             : sample length (number of rows in Y-matrix)
//
// Argument types, from which the dimensions (kk,order,NN,MM) are extracted:
// left_limit_    : double
// right_limit_   : double
// tau            : double
// eta_left_real  : matrix of dimension (J,k)
// eta_left_imag  : matrix of dimension (J,k)
// eta_right_real : matrix of dimension (J,k)
// eta_right_imag : matrix of dimension (J,k)
// FleftMat       : matrix of dimension (k,2*k)
// FrightMat      : matrix of dimension (k,2*k)
// N_             : double
//
// Value:
// vector (length J) of diagonal integrals of the Greens functions 
//
// Implicit assumptions:
// 1) left_limit_ < right_limit_
// 2) Sampling is equidistant as described in Markussen (2013)
//
// Remarks:
// 1) If eta_left_real <=0 and eta_right_real >= 0, then the algorithm is
//    numerically stable.


extern "C" SEXP fdaTrace(SEXP left_limit_, SEXP right_limit_, 
                         SEXP tau, 
                         SEXP eta_left_real,  SEXP eta_left_imag,
                         SEXP eta_right_real, SEXP eta_right_imag,
                         SEXP FleftMat, SEXP FrightMat,
                         SEXP N_
                         ) {
  try{
    // copy data to armadillo structures
    double  left_limit = Rcpp::as<double>( left_limit_);
    double right_limit = Rcpp::as<double>(right_limit_);
    arma::cx_mat eta_left = 
      arma::cx_mat(Rcpp::as<arma::mat>(eta_left_real),
                   Rcpp::as<arma::mat>(eta_left_imag));
    arma::cx_mat eta_right = 
      arma::cx_mat(Rcpp::as<arma::mat>(eta_right_real),
                   Rcpp::as<arma::mat>(eta_right_imag));
    // remark: eta_left and eta_right defined as row vectors
    arma::mat Fleft    = Rcpp::as<arma::mat>(FleftMat);
    arma::mat Fright   = Rcpp::as<arma::mat>(FrightMat);
    int kk       = eta_left.n_cols;
    int JJ       = eta_left.n_rows;
    int NN       = Rcpp::as<int>(N_);
    double Delta = (right_limit-left_limit)/NN;

    // initialize vector to contain result
    arma::vec traceGreen = arma::zeros<arma::vec>(JJ);

    // initialize variables
    // remark: first rows of Wleft and Wright remain fixed at initial values
    //         vector v1 remain fixed at initial value
    arma::cx_mat    Wleft  = arma::ones<arma::cx_mat>(2*kk,kk);
    arma::cx_mat    Wright = arma::ones<arma::cx_mat>(2*kk,kk);
    arma::cx_rowvec v1     = arma::ones<arma::cx_rowvec>(kk);
    arma::cx_mat    W            = arma::zeros<arma::cx_mat>(2*kk,2*kk);
    arma::cx_mat    invW         = arma::zeros<arma::cx_mat>(2*kk,2*kk);
    arma::cx_vec    vleft        = arma::zeros<arma::cx_vec>(kk);
    arma::cx_vec    vright       = arma::zeros<arma::cx_vec>(kk);
    arma::cx_rowvec exp_left     = arma::zeros<arma::cx_rowvec>(kk);
    arma::cx_rowvec exp_right    = arma::zeros<arma::cx_rowvec>(kk);
    arma::cx_mat    Aleft_left   = arma::zeros<arma::cx_mat>(kk,kk);
    arma::cx_mat    Aright_right = arma::zeros<arma::cx_mat>(kk,kk);
    arma::cx_mat    Aleft_right  = arma::zeros<arma::cx_mat>(kk,kk);
    arma::cx_mat    Aright_left  = arma::zeros<arma::cx_mat>(kk,kk);
    arma::cx_mat    Bmat         = arma::zeros<arma::cx_mat>(kk,kk);
    arma::cx_mat    Fleft_W      = arma::zeros<arma::cx_mat>(kk,kk);
    arma::cx_mat    Fright_W     = arma::zeros<arma::cx_mat>(kk,kk);

    // loop over sets of eigenvalues
    for (int jj=0; jj<JJ; jj++) {
      // define variables from Markussen (2013), Theorem 2
      for (int ii=1; ii<2*kk; ii++) {
        Wleft.row(ii) = eta_left.row(jj) % Wleft.row(ii-1);   
      }
      for (int ii=1; ii<2*kk; ii++) {
        Wright.row(ii) = eta_right.row(jj) % Wright.row(ii-1);
      }
      W      = arma::join_rows(Wleft,Wright);
      invW   = arma::pinv(W);
      vleft  = invW.submat(0, 2*kk-1,  kk-1,2*kk-1);
      vright = invW.submat(kk,2*kk-1,2*kk-1,2*kk-1);

      // make some pre-computations
      exp_left  = exp((right_limit-left_limit)*eta_left.row(jj));
      exp_right = exp((left_limit-right_limit)*eta_right.row(jj));
      Fleft_W   = arma::pinv(Fleft *Wleft) *Fleft *Wright;
      Fright_W  = arma::pinv(Fright*Wright)*Fright*Wleft;

      // compute Aleft_left
      // remark: Bmat is temporarily used before it is computed below
      Aleft_left          = arma::repmat(arma::trans(exp_left),1,kk)-arma::repmat(exp_left,kk,1);
      Bmat                = Delta*(arma::repmat(arma::trans(eta_left.row(jj)),1,kk)-
                                   arma::repmat(eta_left.row(jj),kk,1));
      Aleft_left.diag()   = NN*exp_left;
      Bmat.diag()         = arma::ones<arma::cx_vec>(kk);
      Aleft_left          = Aleft_left / Bmat;

      // compute Aright_right
      // remark: Bmat is temporarily used before it is computed below
      Aright_right        = arma::repmat(arma::trans(exp_right),1,kk)-arma::repmat(exp_right,kk,1);
      Bmat                = Delta*(arma::repmat(eta_right.row(jj),kk,1)-
                                   arma::repmat(arma::trans(eta_right.row(jj)),1,kk));
      Aright_right.diag() = NN*exp_right;
      Bmat.diag()         = arma::ones<arma::cx_vec>(kk);
      Aright_right        = Aright_right / Bmat;

      // compute Aleft_right
      // remark: Bmat is temporarily used before it is computed below
      Bmat        = arma::repmat(arma::trans(eta_left.row(jj)),1,kk)-arma::repmat(eta_right.row(jj),kk,1);
      Aleft_right = (exp((right_limit-left_limit)*Bmat)-arma::ones<arma::cx_mat>(kk,kk)) / (Delta*Bmat);

      // compute Aright_left
      Aright_left = arma::trans(Aleft_right);

      // compute Bmat
      Bmat = Fleft_W*arma::diagmat(exp_right)*Fright_W*
             arma::pinv(arma::eye<arma::cx_mat>(kk,kk)-
	                arma::diagmat(exp_left)*Fleft_W*arma::diagmat(exp_right)*Fright_W);

      // compute trace
      traceGreen(jj)  = arma::as_scalar(arma::real(NN*v1*vleft));
      traceGreen(jj) += arma::as_scalar(arma::real(v1*(Fleft_W %Aleft_right)*vright));
      traceGreen(jj) -= arma::as_scalar(arma::real(v1*(Fright_W%Aright_left)*vleft));
      traceGreen(jj) -= arma::as_scalar(arma::real(v1*((Fright_W*arma::diagmat(exp_left)*
                                                        Fleft_W)%Aright_right)*vright));
      traceGreen(jj) += arma::as_scalar(arma::real(v1*(Bmat%Aleft_left)*vleft));
      traceGreen(jj) += arma::as_scalar(arma::real(v1*((Bmat*arma::diagmat(exp_left)*Fleft_W)%
                                                       Aleft_right)*vright));
      traceGreen(jj) -= arma::as_scalar(arma::real(v1*((Fright_W*arma::diagmat(exp_left)*Bmat)%
                                                       Aright_left)*vleft));
      traceGreen(jj) -= arma::as_scalar(arma::real(v1*((Fright_W*arma::diagmat(exp_left)*Bmat*
                                                       arma::diagmat(exp_left)*Fleft_W)%Aright_right)*vright));
    }

    // Rescale trace
    traceGreen = traceGreen/Rcpp::as<double>(tau);

    // return it to R
    return Rcpp::wrap(traceGreen);

    // Below possible exceptions are taken care off
  } catch( std::exception &ex ) {
      forward_exception_to_r( ex );
  } catch(...) { 
      ::Rf_error( "c++ exception (unknown reason)"); 
  }

  // return to R
  return R_NilValue;
}
