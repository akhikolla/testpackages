// matrix_utils.cpp

#include "dmbc.h"

// Matrix column sums
void colsums(double* colsums, const double* A, int nrows, int ncols){
  for(int j = 0; j < ncols; j++){
    for(int i = 0; i < nrows; i++){
      colsums[j] += A[i + nrows*j];
    }
  }
}

// Matrix row sums
void rowsums(double* rowsums, const double* A, int nrows, int ncols){
  for(int i = 0; i < nrows; i++){
    for(int j = 0; j < ncols; j++){
      rowsums[i] += A[i + nrows*j];
    }
  }
}

// Conversion of the dissimilarities from matrix to vector
arma::vec dissM2V(const arma::mat& d){
	int n = d.n_rows;
	int m = n*(n - 1)/2;
	int count = 0;
	arma::vec dvec(m, arma::fill::zeros);
	
	for(int j = 0; j < n; j++){
		for(int i = (j + 1); i < n; i++){
			dvec(count) = d(i, j);
			count++;
		}
	}

	return dvec;
}

// Matrix to vector
arma::vec mat2vec(const arma::mat& A, const int& j){
	int nrows = A.n_rows;
	arma::vec v(nrows, arma::fill::zeros);
	
	for(int i = 0; i < nrows; i++){
		v(i) = A(i, j);
	}

	return v;
}

// Vector to matrix
arma::mat vec2mat(const arma::mat& A, const int& j, const arma::vec& v){
	int nrows = A.n_rows;
	arma::mat B = A;

	for(int i = 0; i < nrows; i++){
		B(i, j) = v(i);
	}

	return B;
}

// Computation of the Mahalanobis distance
arma::vec mahalanobis(const arma::mat& x, const arma::vec& center, const arma::mat& sigma){
  int nrows = x.n_rows;
  int ncols = x.n_cols;
  arma::mat x_centered = x;
  arma::mat sigmainv = arma::inv_sympd(sigma);

  for(int i = 0; i < nrows; i++){
    for(int j = 0; j < ncols; j++){
      x_centered(i, j) -= center(j);
    }
  }

  arma::vec m = arma::sum(x_centered % (x_centered * sigmainv), 1);
  
  return m;
}
