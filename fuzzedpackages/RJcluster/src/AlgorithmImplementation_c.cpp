#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec getIndividualList(const arma::colvec& classification, const arma::colvec& temp_index, int i)
{
  // get vector of desired size
  arma::vec x = temp_index.elem( find(classification == (i + 1)) );

  return(x);
}


// [[Rcpp::export]]
List getCCmatrix_c(const arma::colvec& classification, const arma::colvec& temp_index, int G)
{
  List CC(G);

  for (int i = 0; i < G; i++)
  {
    // get the values to replace
    arma::vec temp = getIndividualList(classification, temp_index, i);

    // create new entry in CC
    CC[i] = temp;
  }

  return(CC);
}

// [[Rcpp::export]]
arma::mat getMatrixMeans_c(const List CC, const arma::mat X, const int d)
{
  // get dimensions
  int p = X.n_cols;

  arma::mat Cmeans(d, p);

  for(int i = 0; i < d; i++)
  {
    // subset matrix X
    arma::uvec indices = CC(i);

    // subtract one from indices
    for (int j = 0; j < indices.n_elem; j++)
    {
      indices[j] = indices[j] - 1;
    }

    // subset X
    arma::mat X_temp = X.rows(indices);
    //
    // // get vector of column means
    arma::rowvec Cmeans_temp = sum(X_temp);

    for (int j = 0; j < Cmeans_temp.n_elem; j++)
    {
      Cmeans_temp[j] = Cmeans_temp[j] / X_temp.n_rows;
    }

    // // add to means matrix
    Cmeans.row(i) = Cmeans_temp;
  }

  return(Cmeans);
}

// [[Rcpp::export]]
arma::vec assignGroups_c(const int N, const int G, const arma::colvec& classification, const List CC)
{
  // Group to be returned
  arma::vec Group = arma::zeros<arma::vec>(N);


  for (int i = 0; i < G; i++)
  {
    arma::uvec Cluster1 = find(classification == (i + 1));

    for (int j = 0; j < Cluster1.n_elem; j++)
    {
      // find elemetns in CC that fit that assignment
      arma::uvec indices = CC[Cluster1[j]];

      // subtract 1 from indices
      for (int k = 0; k < indices.n_elem; k++)
      {
        indices[k] = indices[k] - 1;
      }

      // assign group at those indices
      arma::vec vals(indices.n_elem);
      vals.fill(i + 1);
      Group.elem(indices) = vals;
    }
  }

  return(Group);
}

// [[Rcpp::export]]
arma::mat getFinalMeans_c(const int G, const arma::vec& Group, const arma::mat& X)
{
  int p = X.n_cols;

  arma::mat Lmark(G, p);

  for (int i = 0; i < G; i++)
  {
    // get indices
    arma::uvec indices = find( Group == i + 1);

    // subset X
    arma::mat X_temp = X.rows(indices);

    // get vector of column means
    arma::rowvec means_temp = sum(X_temp);

    for (int j = 0; j < means_temp.n_elem; j++)
    {
      means_temp[j] = means_temp[j] / X_temp.n_rows;
    }

    // add to final matrix
    Lmark.row(i) = means_temp;
  }


  return(Lmark);
}
