#include <RcppArmadillo.h>

// Arguments:
//   boxcox         : Box-Cox parameter
//   geometricMean  : Geometric mean of untransformed observations
//   Yoriginal      : matrix of untransformed observations
//   Ymat           : matrix to contain transformed observations
//
// Argument types, from which the dimensions (NN,MM) are extracted:
//   boxcox         : double
//   geometricMean  : double
//   Yoriginal      : matrix of dimension (NN,MM)
//   Ymat           : matrix of dimension (NN,MM*(1+pp+qq))
//
// Value:
//   Box-Cox transformed observations are copied to first MM columns of Ymat

extern "C" SEXP boxcoxTransform(SEXP boxcox,
                                SEXP geometricMean,
                                SEXP Yoriginal,
                                SEXP Ymat
				) {
  try{
    // copy data to armadillo structures
    double mu      = Rcpp::as<double>(boxcox);
    double geoMean = Rcpp::as<double>(geometricMean);

    // creates Rcpp objects from SEXP's
    Rcpp::NumericMatrix YorigR(Yoriginal);
    Rcpp::NumericMatrix Yr(Ymat);

    // extract dimensions from defined objects
    int NN = YorigR.nrow();
    int MM = YorigR.ncol();

    // reuses memory and avoids extra copies  
    arma::mat Yorig(YorigR.begin(),NN,MM,false);
    arma::mat Y(Yr.begin(),NN,Yr.ncol(),false);

    // BoxCox transformation
    if (mu==0) {
      Y.cols(0,MM-1) = geoMean*arma::log(Yorig);
    } else {
      Y.cols(0,MM-1) = (arma::pow(Yorig,mu)-1)/(mu*pow(geoMean,mu-1));
    }

    // Result is saved and returned in the external variables:
    //   Ymat

    // Below possible exceptions are taken care off
  } catch( std::exception &ex ) {
      forward_exception_to_r( ex );
  } catch(...) { 
      ::Rf_error( "c++ exception (unknown reason)"); 
  }

  // return to R
  return R_NilValue;
}
