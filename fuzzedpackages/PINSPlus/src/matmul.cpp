// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppParallel.h>
#include <RcppArmadillo.h>



using namespace Rcpp;
using namespace RcppParallel;
using namespace arma;


// arma::mat stra(arma::mat A, arma::mat B) {
//   double m1 = (A(0,0) + A(1,1))*(B(0,0)+B(1,1));
//   double m2 = (A(1,0) + A(1,1))*B(0,0);
//   double m3 = A(0,0)*(B(0,1) - B(1,1));
//   double m4 = A(1,1)*(B(1,0) - B(0,0));
//   double m5 = (A(0,0) + A(0,1))*B(1,1);
//   double m6 = (A(1,0) - A(0,0))*(B(0,0) + B(0,1));
//   double m7 = (A(0,1) - A(1,1))*(B(1,0) + B(1,1));
//   
//   armat:mat C(2, 2);
//   C(0,0) = m1 + m4 - m5 + m7;
//   C(0,1) = m3 + m5;
//   C(1,0) = m2 + m4;
//   C(1,1) = m1 - m2 + m3 + m6;
//   return C;
// }
// 
// // [[Rcpp::export]]
// arma::mat fast_mat(arma::mat& A, arma::mat& B) {
//   int n = A.n_rows;
//   int m = B.n_cols;
//   int l = A.n_cols;
//   armat:mat C(n, m, arma::fill::zeros);
//   int i = 0;
//   int j = 0;
//   int k = 0;
// 
//   while(i + 3 < n){
//     while(j + 3 < m) {
//       while(k + 3 < l)
//       {
//         C.submat(i,j,i+3,j+3) += A.submat(i, k, i + 3, k + 3)*B.submat(k, j, k + 3, j + 3); //stra(A.submat(i, k, i + 1, k + 1), B.submat(k, j, k + 1, j + 1));
//         k += 4;
//       }
//       k = 0;
//       j += 4;
//     }
//     j = 0;
//     i += 4;
//   }
//   return C;
// }

struct mat_mul : public Worker {
  
  // input matrix to read from
  arma::mat*  A;
  arma::mat*  B;
  arma::mat*  res;
  
  // output matrix to write to
  
  
  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  mat_mul(arma::mat* A, arma::mat* B, arma::mat* res)
    : A(A), B(B), res(res) {}
  
  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    res->rows(begin, end-1) = A->rows(begin, end-1) * (*B);
    
  }
  
};

// [[Rcpp::export]]
arma::mat mat_mul_para(arma::mat& A, arma::mat& B) {
  arma::mat res(A.n_rows, B.n_cols, arma::fill::zeros);
  
  // create the worker
  mat_mul mat_mul(&A, &B, &res);
  
  // call it with parallelFor
  parallelFor(0, A.n_rows, mat_mul, 20);
  
  return res;
}

// [[Rcpp::export]]
arma::mat c_transpose(arma::mat& A) {

  return A.t();
} 

