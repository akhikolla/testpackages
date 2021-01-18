#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat sam_internal(arma::mat& x, arma::mat& em) {
  int ns = x.n_rows;
  int nend = em.n_rows;
  arma::mat out(ns, nend);
  out.zeros();
  arma::mat normEM = sum(pow(em,2),1);
  arma::mat normT  = sum(pow(x,2),1);

  for(int j = 0; j < nend; ++j){
    for(int i = 0; i < ns; ++i){
      out(i,j) = acos(accu(x.row(i) % em.row(j)) / sqrt(normT(i,0) * normEM(j,0)));
    }
  }
  return out;
}

// [[Rcpp::export]]
arma::cube sam_main(arma::cube myCube, arma::mat em){

  arma::cube out(myCube.n_rows, myCube.n_cols, em.n_rows);
  out.zeros();

  int crows = myCube.n_rows;

  for(int i = 0; i < crows; ++i){
    arma::mat cubeSlice(myCube.row(i));
    //out.row(i) = cubeSlice;
    out.row(i) = sam_internal(cubeSlice, em);
  }

  // Transpose the cube, to match uFTIR output
  // is transposing and flipping the rows
  //for (size_t s = 0; s < out.n_slices; ++s)
  //     out.slice(s) = arma::flipud(out.slice(s).t());

  return out;
}

// [[Rcpp::export]]
arma::ucube sam_match(arma::cube myCube)
{
  if(myCube.has_nan())
  {
    myCube.replace(arma::datum::nan, arma::datum::inf);
  }

  arma::ucube out(size(myCube));

  for(arma::uword i = 0; i < myCube.n_rows; ++i)
  {
    for(arma::uword j = 0; j < myCube.n_cols; ++j)
    {
      arma::vec match = myCube.tube(i, j);
      out.tube(i, j) = arma::sort_index(match, "ascend");
    }
  }

  // Since R starts counting from 1
  out += 1;

  return out;
}

// [[Rcpp::export]]
Rcpp::List ctile_sam(arma::cube myCube, arma::mat em)
{
  myCube = sam_main(myCube, em);

  arma::ucube myCubeSort(size(myCube));
  myCubeSort = sam_match(myCube);

  return Rcpp::List::create(Rcpp::Named("raw_sam") = myCube,
                            Rcpp::Named("match_sam") = myCubeSort);

}
