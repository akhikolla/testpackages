#include <RcppArmadillo.h>
#include <iostream>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::cube cmosaic_compose(Rcpp::StringVector file, arma::mat xy_pos, int xmax, int ymax, int zmax){
  
  // define size of the final cube using data from the first
  std::string fname(file[0]);
  arma::cube A;
  A.load(fname, arma::arma_binary);
  if(zmax == -1){
    zmax = A.n_slices;
  }
  arma::cube out(A.n_rows * (xmax+1), A.n_cols * (ymax+1), zmax);
  
  // get the fpa size
  int fpa = A.n_rows;
  
  // loading each cube and inserting them in out
  for(int i = 0; i <= file.size() - 1 ; ++i)
  {
    // reading a cube
    std::string fname(file[i]);
    arma::cube A;
    A.load(fname, arma::arma_binary);
    
    // copying it to the cube subview of out.
    // we are using the .tube subview:
    int first_row = xy_pos(i, 0) * fpa;
    int first_col = xy_pos(i, 1) * fpa;
    out.tube(first_row, first_col, first_row + fpa - 1, first_col + fpa - 1) = A.slices(0, zmax-1);
  }
  
  return out;
}
