#include <Rcpp.h>
using namespace Rcpp;

//' fast Euclidean distance matrix
//'
//' @param x matrix with sample rows for which the distance matrix is computed (to use with vectors, use \code{as.matrix(x)})
//' @examples
//' #require(microbenchmark)
//' #x = rnorm(100)
//' #microbenchmark(fastdist(as.matrix(x)),as.matrix(dist(x)))
//' @export
// [[Rcpp::export]]
NumericMatrix fastdist (const NumericMatrix & x){
  unsigned int outrows = x.nrow(), i = 0, j = 0;
  double d;
  Rcpp::NumericMatrix out(outrows,outrows);

  for (i = 0; i < outrows - 1; i++){
    Rcpp::NumericVector v1 = x.row(i);
    for (j = i + 1; j < outrows ; j ++){
      d = sqrt(sum(pow(v1-x.row(j), 2.0)));
      out(j,i)=d;
      out(i,j)=d;
    }
  }

  return out;
}

//' double center a symmetric matrix
//'
//' @param x symmetric matrix
//' @param normalize boolean. If \code{TRUE} the matrix will be normalized to mean 1.
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix doubleCenterSymMat(const NumericMatrix & x, bool & normalize) {
  int i, j;
  int N = x.nrow();
  NumericVector colmeans(N);
  NumericMatrix out(N, N);
  double fullmean;
  double tmp;

  for (i=0; i<N; i++) {
    colmeans(i) = sum(x(i,_))/(double)(N);
  }
  fullmean = sum(colmeans)/N;

  if ( (fullmean == 0) | !normalize) {
    // unnormalized centering
    if (fullmean == 0) warning("It seems that one variable is constant. Constants are always independent. \n");
    for (i=0; i<N; i++)
      for (j=i; j<N; j++) {
        tmp = - x(i, j) + colmeans(i) + colmeans(j) - fullmean;
        out(j, i) = tmp;
        out(i, j) = tmp;
      }
  } else { // normalized centering
    for (i=0; i<N; i++)
      for (j=i; j<N; j++) {
        tmp = (- x(i, j) + colmeans(i) + colmeans(j) - fullmean)/fullmean;
        out(j, i) = tmp;
        out(i, j) = tmp;
      }
  }

  return out;
}


//' fast centered Euclidean distance matrix
//'
//' @param x matrix with sample rows for which the distance matrix is computed (to use with vectors, use \code{as.matrix(x)})
//' @param normalize boolean. If \code{TRUE} the matrix will be normalized to mean 1.
//' @export
// [[Rcpp::export]]
NumericMatrix fastEuclideanCdm (const NumericMatrix & x, bool & normalize){
  unsigned const int N = x.nrow();
  unsigned int i = 0, j = 0;
  NumericMatrix out(N,N);
  NumericVector colmeans(N);
  double tmp, m;


  for (i = 0; i < N - 1; i++){ // row
    NumericVector v1 = x.row(i);
    for (j = i + 1; j < N ; j ++){ // column
      tmp = sqrt(sum(pow(v1-x.row(j), 2.0)));
      out(i,j) = tmp;
      colmeans(i) += tmp;
      colmeans(j) += tmp;
    }
  }

  colmeans = colmeans/N;
  m = sum(colmeans)/(double) N;

  if ( (m == 0) | !normalize) {
    // unnormalized centering
    if (m == 0) warning("It seems that one variable is constant. Constants are always independent. \n");
    for (i = 0; i < N-1; i++){ // row
      for (j = i+1 ; j < N ; j++){ // column
        tmp = -out(i,j) + colmeans(i) + colmeans(j) - m;
        out(i,j) = tmp;
        out(j,i) = tmp;

      }
    }


    for (i = 0; i < N ; i++){ // diag
      out(i,i) =  2*colmeans(i) - m;
    }
  } else {
    // normalized centering
    for (i = 0; i < N-1; i++){ // row
      for (j = i+1 ; j < N ; j++){ // column
        tmp = (-out(i,j) + colmeans(i) + colmeans(j) - m)/m;
        out(i,j) = tmp;
        out(j,i) = tmp;
      }
    }


    for (i = 0; i < N ; i++){ // diag
      out(i,i) =  (2*colmeans(i) - m)/m;
    }
  }
  /* */
   return out;
}


//' for the fast detection of the full dependence structure
//'
//' Returns the row indicies of matrix A which match with B
//'
//' @param A matrix
//' @param B matrix whose rows are subset of A
//'
//' @examples
//' # A = t(utils::combn(10,3))
//' # B = A[sort(sample.int(nrow(A),10)),]
//' # match_rows(A,B)
//'
//' @keywords internal
// [[Rcpp::export]]

NumericVector match_rows(NumericMatrix & A,NumericMatrix &B){
  int i = 0, k;
  NumericVector res (B.nrow());

  for (k = 0; k < B.nrow(); k++) {
    while( is_true(any(A.row(i) != B.row(k)))) {
      i++;
    }
    res(k) = i;
  }
  return res+1;
}

//' bias corrected double centering
//' # included for speed comparison
//'
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix doubleCenterBiasCorrected(const NumericMatrix & x) {
  int i, j;
  int N = x.nrow();
  NumericVector colmeans(N);
  NumericMatrix out(N, N);
  double fullmean;
  double tmp;

  for (i=0; i<N; i++) {
    colmeans(i) = sum(x(i,_))/(double)(N-2);
  }
  fullmean = sum(colmeans)/(N-1);

  // unnormalized centering
  for (i=0; i<N; i++)
    for (j=i+1; j<N; j++) {
      tmp = - x(i, j) + colmeans(i) + colmeans(j) - fullmean;
      out(j, i) = tmp;
      out(i, j) = tmp;
    }

  return out;
}

//' bias corrected double centering with normalizing
//' # included for speed comparison
//'
//' @keywords internal
// [[Rcpp::export]]
NumericVector doubleCenterBiasCorrectedUpper(const NumericMatrix & x) {
  // returning a vector
  int i, j, k = 0;
  int N = x.nrow();
  NumericVector colmeans(N);
  NumericVector out(N*(N-1)/2);
  const NumericVector zero(N*(N-1)/2);
  double fullmean;
  double tmp;
  double centeredmean = 0.0;

  for (i=0; i<N; i++) {
    colmeans(i) = sum(x(i,_))/(double)(N-2);
  }
  fullmean = sum(colmeans)/(N-1);

  // unnormalized centering
  for (i=0; i<N; i++)
    for (j=i+1; j<N; j++) {
      tmp = - x(i, j) + colmeans(i) + colmeans(j) - fullmean;
      centeredmean += pow(tmp, 2);
      out(k) = tmp;
      k += 1;
    }

  centeredmean = 2*centeredmean/(N*(N-3));

  if (centeredmean < DBL_EPSILON) {
//    warning("A doubly centered distance matrix seems to have only zero entries. Maybe one variable is constant. Constants are always independent. \n");
    return zero;
  } else {
    return out/pow(centeredmean, 0.5);
  }
}

//' bias corrected double centering with normalizing constants for upper and lower bound
//'
//' @keywords internal
// [[Rcpp::export]]
List doubleCenterBiasCorrectedUpperLower(const NumericMatrix & x, const int & n) {
  // returning a list, with vector instead of matrix
  int i, j, k = 0;
  int N = x.nrow();
  NumericVector colmeans(N);
  NumericVector out(N*(N-1)/2);
  const NumericVector zero(N*(N-1)/2);
  double fullmean;
  double tmp;
  double upper = 0.0, lower = 0.0;
  List ret;

  for (i=0; i<N; i++) {
    colmeans(i) = sum(x(i,_))/(double)(N-2);
  }
  fullmean = sum(colmeans)/(N-1);

  // unnormalized centering
  for (i=0; i<N; i++)
    for (j=i+1; j<N; j++) {
      tmp = - x(i, j) + colmeans(i) + colmeans(j) - fullmean;
      upper += pow(tmp, 2);
      lower += fabs(pow(tmp, n));
      out(k) = tmp;
      k += 1;
    }



  upper = 2*upper/(N*(N-3));

  if (upper < DBL_EPSILON) {
//    warning("A doubly centered distance matrix seems to have only zero entries. Maybe one variable is constant. Constants are always independent. \n");
    out = zero;
    upper = 1;  // to prevent division by 0
    lower = 1;  // to prevent division by 0
  }

  ret["out"] = out;

  ret["upper"] =  pow(upper, 0.5);

  lower = 2*lower/(N*(N-3));
  ret["lower"] =  pow(lower, (double)1/n);

  return ret;
}




/*** R
*/
