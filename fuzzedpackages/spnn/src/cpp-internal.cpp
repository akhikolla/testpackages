#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.spnn_predict_cpp)]]
arma::mat spnn_predict_cpp(arma::cube& setCube,
                           arma::mat& newData,
                           arma::mat& sigmaInverse){

  arma::mat final(newData.n_rows, setCube.n_slices); // final matrix of probabilities

  for(unsigned int i = 0; i < newData.n_rows; i++){ // iterate through newData rows
    
    arma::rowvec results(final.n_cols); // final probability vector
    
    for(unsigned int j = 0; j < setCube.n_slices; j++){ // iterate through setCube slices
      
      arma::mat slice = setCube.slice(j); // subset slice for each class
      arma::uvec nan_cutoff_vec = arma::find_nonfinite(slice.col(0));
      int nan_cutoff = arma::conv_to<int>::from(nan_cutoff_vec.head(1));
      slice = slice.submat(0, 0, nan_cutoff - 1, slice.n_cols - 1); // cut off row where NaNs begin
      
      double f; // scale invariant estimate
      double F = 0; // density for given class j
      
      for(unsigned int k = 0; k < slice.n_rows; k++){ // iterate through class observations
        
        arma::mat xd = newData.row(i) - slice.row(k);
        f = std::exp(-0.5*std::abs(arma::as_scalar(xd*sigmaInverse*xd.t()))); // scale invariant estimate
        if(ISNAN(f)){ f = 0; } // error handler
        
        F += f; // add to density
      }
      
      results(j) = F / slice.n_rows;
    }
    
    final.row(i) = results / arma::accu(results); // final output layer
  }

  return final;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(.cspnn_predict_cpp)]]
arma::mat cspnn_predict_cpp(arma::mat& xr,
                            arma::mat& newData,
                            arma::mat& sigmaInverse){

  arma::mat final(newData.n_rows, xr.n_rows); // final matrix of probabilities

  for(unsigned int i = 0; i < newData.n_rows; i++){
    
    arma::rowvec results(xr.n_rows);
    
    for(unsigned int j = 0; j < xr.n_rows; j++){
      
      arma::mat xd = newData.row(i)-xr.row(j);
      double f = std::exp(-0.5*std::abs(arma::as_scalar(xd*sigmaInverse*xd.t()))); // scale invariant estimate
      if(ISNAN(f)){ f = 0; } // error handler
      
      results(j) = f;
    }
    
    final.row(i) = results / arma::accu(results); // final output layer
  }

  return final;
}
