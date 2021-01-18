#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace std;
using namespace arma;

////////////////////// Residuals from conditional logistic regression //////////////////


void conditResClust(const vec theta, int ysum, int clustsize, int min_idx, int inv, const vec & y, const mat & X, int nparams, vec & res, mat & dres){

  vec yi = y.rows(min_idx - 1, min_idx + clustsize - 2);

  if (ysum == 1) {

    mat Xi  = X.rows( min_idx - 1, min_idx + clustsize - 2 );
    vec wi =  exp( Xi * theta );

    double sumwi = sum( wi );
    vec yihat = wi / sumwi;

    mat yiXi = Xi; 
    yiXi.each_col() %= yihat; 

    mat dresi = -Xi ;
    dresi.each_row() += sum(yiXi, 0);
    dresi.each_col() %= yihat;

    res.rows(min_idx - 1, min_idx + clustsize - 2) = yi - yihat;
    dres.rows(min_idx - 1, min_idx + clustsize - 2) = dresi  ;

  } else {

    mat Xi(clustsize, nparams);

    if(inv) {
      Xi  = -X.rows( min_idx - 1, min_idx + clustsize - 2 );
      ysum = clustsize - ysum;
    } else {
      Xi  = X.rows( min_idx - 1, min_idx + clustsize - 2 );
    }

    vec wi =  exp( Xi * theta );

    vec b(ysum + 1, fill::zeros);
    b(0) = 1;
    b(1) = wi(0);

    mat Db(ysum + 1, nparams, fill::zeros);
    Db.row(1) = Xi.row(0) * wi(0);

    for (int j = 1; j < clustsize; j++) {

      // Set element of the first row
      // b(1) = b(1) + wi(j);
      int maxrow = min(j + 1, ysum); 

      mat Db_plus_bXi = Db.rows(0, ysum - 1) + b.rows(0, ysum - 1) * Xi.row(j);
      Db.rows(1, ysum) = Db.rows(1, ysum) + Db_plus_bXi * wi(j);
      b.subvec(1, maxrow) = b.subvec(1, maxrow) + b.subvec(0, maxrow - 1) * wi(j);
    }

    mat Pwi(clustsize, ysum + 1, fill::ones);

    for(int l = 1; l < ysum + 1; l++) {
      Pwi.col(l) = Pwi.col(l - 1) % -wi;
    }

    vec cwi = - Pwi.tail_cols(ysum) * flipud(b.rows(0, ysum - 1));
    vec yihat = cwi / b( ysum );

    mat Dc(clustsize, nparams, fill::zeros); 
    mat Dc_sum_add(clustsize, nparams); 

    for (int l = 0; l < ysum; l++) {
      Dc_sum_add = Xi * (ysum - 1 - l) * b(l);
      Dc_sum_add.each_row() += Db.row(l);  
      Dc_sum_add.each_col() %= Pwi.col(ysum - 1 - l);  
      Dc = Dc + Dc_sum_add;
    }

    mat Xiyihat = Xi;
    Xiyihat.each_col() %= yihat;
    mat Dcwi = Dc;
    Dcwi.each_col() %= wi;
    mat yihatDb = yihat * Db.tail_rows(1);

    if(inv) {
      res.rows(min_idx - 1, min_idx + clustsize - 2) = yi + yihat - 1;
      dres.rows(min_idx - 1, min_idx + clustsize - 2) = Xiyihat - ( yihatDb - Dcwi ) / b(ysum);
    } else {
      res.rows(min_idx - 1, min_idx + clustsize - 2) = yi - yihat;
      dres.rows(min_idx - 1, min_idx + clustsize - 2) = -Xiyihat + ( yihatDb - Dcwi ) / b(ysum);
    }
  }
}

// [[Rcpp::export]]
RcppExport SEXP _conditRes(SEXP thetahat, SEXP ysums, SEXP clustsizes, SEXP minidx, SEXP inv, SEXP yin, SEXP Xcent){
  IntegerVector y_sums(ysums);
  IntegerVector clust_sizes(clustsizes);
  IntegerVector min_idx(minidx);
  IntegerVector invert(inv);
  NumericMatrix X_cent(Xcent);
  NumericVector thetatmp(thetahat);
  NumericVector ytmp(yin);

  int ndisc = y_sums.length(), tot_rows = X_cent.nrow(), nparams =  thetatmp.size();
  vec theta = as<vec>(thetatmp);
  vec y = as<vec>(ytmp);
  vec res(tot_rows, fill::zeros);
  mat dres(tot_rows, nparams, fill::zeros);
  mat X = as<mat>(X_cent);

  for(int j=0; j < ndisc; j++){
    conditResClust(theta, y_sums[j], clust_sizes[j], min_idx[j], invert[j], y, X, nparams, res, dres);
  }

  return List::create( Named("res") = res,
  		       Named("dres") = dres
  		      );
}

// [[Rcpp::export]]
RcppExport SEXP _center(SEXP Uin, SEXP ID) {
  NumericMatrix U_tmp(Uin);
  mat U = as<mat>(U_tmp);

  IntegerVector id_tmp(ID);
  uvec id = as<uvec>(id_tmp);
  uvec uid = unique(id);

  int n_obs = U.n_rows;
  int n_col = U.n_cols;
  // int n_clust = uid.n_elem;
  // int n_clust(nclust);

  // mat Uc_means(size(U));
  mat Uout(U);

  // Initialize the cluster identifier flag
  unsigned int clust_id = id[0];
  // Initialize the cluster index flags
  unsigned int clust_start_idx = 0;
  unsigned int clust_end_idx = 0;

  // double c_sum = 0;
  // Temporary vector to store the sums
  // for each column in the matrix Uin
  // for each cluster
  vec u_sums(n_col);
  u_sums.fill(0.0);
  // int c_size = 0; 

  // Loop over the rows of the matrix Uin
  for(int k = 0; k < n_obs; ++k){

    // If id[k] is a new cluster
    if( id[k] != clust_id ){

      // Update the out matrix
      clust_end_idx = k - 1;

      // Loop over the columns of U for the cluster id[k]
      // to create centered elements in Uout
      for(int l = 0; l < n_col; ++l){
	Uout(span(clust_start_idx, clust_end_idx),l) -= u_sums(l) / (clust_end_idx - clust_start_idx + 1);
	// Uout(clust_start_idx,l) = u_sums(l) / (clust_end_idx - clust_start_idx + 1);
	// Uout(clust_end_idx,l) = u_sums(l) / (clust_end_idx - clust_start_idx + 1);
	// Uout(clust_start_idx,l) = u_sums(l) ;
	// Uout(clust_end_idx,l) = u_sums(l) ;
	// Uout(clust_end_idx,l) = k;
      }

      // Update the cluster_id to the new one
      clust_id = id[k];
      // Update the cluster_start_idx to the new one
      clust_start_idx = k;
      // Update the cluster sums
      u_sums = U.row(k).t();
      // // Reset the cluster sums
      // u_sums.fill(0.0);

    } else {

      // Update the cluster sums
      u_sums += U.row(k).t();

    }

  }

  // Center the last cluster
  clust_end_idx = n_obs - 1;

  for(int l = 0; l < n_col; ++l){
    Uout(span(clust_start_idx, clust_end_idx),l) -= u_sums(l) / (clust_end_idx - clust_start_idx + 1);
    // Uout(clust_start_idx,l) -= u_sums(l) / (clust_end_idx - clust_start_idx + 1);
    // Uout(clust_end_idx,l) -= u_sums(l);

  }

  // Uout.rows(clust_start_idx, clust_end_idx).each_col() -= u_sums / (clust_end_idx - clust_start_idx + 1);

  return wrap(Uout);

}

static const R_CallMethodDef CallEntries[] = {
    {"center", (DL_FUNC) &_center, 2},
    {"conditRes", (DL_FUNC) &_conditRes, 7},
    {NULL, NULL, 0}
};
 
// Register routines and disable symbol search
RcppExport void R_init_drgee(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

